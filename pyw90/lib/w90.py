import re, os, copy
import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
from typing import Dict, List
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy import interpolate
from scipy.signal import savgol_filter

from lib.config import Config

import logging
logger = logging.getLogger(__name__)

class W90():
    r'''
    This is main class of `Auto-Wannier90-Fit` and contains necessary `dis_windows` parameters for `Wannier90`. The class also offers methods including `dis_window` suggestion and evaluating the quality of Wannier Functions.
    '''
    def __init__(self, config:Config, path:str='.', nbnds_excl:int=None, ndeg:int=1):
        '''
        Init

        :param path: path of the working folder.
        :param nbnds_excl: Number of bands which is not including below the lowest tight-binding model energy band (It doesn't mean the absolute value of lowest tight-bingding model energy band is higher than the excluded band from VASP ). The concept is similar to `exclued_bands` block in `.win` input for `Wannier90` to generate the overlap and projection of wave functions such as `.mmn` file.
        :param ndeg: Degeneracy of bands. Only used in calcution of number Wannier functions from number of projections, i.e. #WFs = ndeg * #projs
        '''
        self.config = config
        self._win   = config.win
        self._sys   = self._win[:-4]
        self._fname = self._sys + '.eig'
        self._dname = path    # the directory containing the input file
        self.efermi = config.efermi
        self.nbnds_excl = nbnds_excl
        self.nwann  = config.nwann
        self.ndeg   = ndeg      # denegeracy of bands, actually need to set 2 only meets Kramers degeneracy.
        self.eps    = 4.0e-3     # tolerance for generate `dis_windows`
        self.inf    = 1.0e+6     # evaluate self.inf when job meets illegal input `dis_windows`
        self._block = '█'

        self.read_eigenval()

    def read_eigenval(self):
        r'''
        Read band energies from VASP EIGENVAL file or `self._sys.eig` file.
        '''
        if self._fname[-3:] == 'eig':
            self.nspin = 1
            self.nelect = None
            data = np.loadtxt(f'{self._dname}/{self._fname}')
            self.ir_kpath = None
            self.ir_kptwt = None
            self.ir_nkpts = int(data[:, 1].max())
            self.nbnds = int(data[:, 0].max())
            self.ir_ebands = data[:, 2].reshape((self.ir_nkpts, 1, -1)).swapaxes(0, 1)
        else:
            with open(f'{self._dname}/{self._fname}') as inp:
                # read all the data.
                dat = np.array([line.strip() for line in inp if line.strip()])

            # extract the needed info in the header.
            self.nspin = int(dat[0].split()[-1])
            self.nelect, self.ir_nkpts, self.nbnds = map(int, dat[5].split())

            # remove the header
            dat = dat[6:]

            # extract the k-points info
            dump = np.array(
                [xx.split() for xx in dat[::self.nspin * self.nbnds + 1]], dtype=float)
            self.ir_kpath = dump[:self.ir_nkpts, :3]
            self.ir_kptwt = dump[:self.ir_nkpts, -1]

            # extract the ebands info
            ebands_flag = np.ones(dat.size, dtype=bool)
            ebands_flag[::self.nspin * self.nbnds + 1] = 0
            if self.nspin == 1:
                ebands = np.array([xx.split()[1] for xx in dat[ebands_flag]],
                                    dtype=float)
            else:
                ebands = np.array([xx.split()[1:3] for xx in dat[ebands_flag]],
                                    dtype=float)
            ebands.shape = (self.ir_nkpts, self.nspin, self.nbnds)
            self.ir_ebands = ebands.swapaxes(0, 1)

        self.emax = self.ir_ebands.max()
        self.emin = self.ir_ebands.min()
        self.eband_max = np.max(self.ir_ebands, axis=1)[0]
        self.eband_min = np.min(self.ir_ebands, axis=1)[0]

    def read_wannier90_win(self):
        r'''
        Read parameters from `self._win` file.
        '''
        with open(self._win) as f:
            dat = [line.strip() for line in f if line.strip()]

        for l in dat:
            r = re.match('num_wann\s*=\s*(\d+)', l)
            if r:
                self.nwann = eval(r.group(1))
            r = re.match('dis_win_max\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.win_max = eval(r.group(1))
            r = re.match('dis_froz_max\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.froz_max = eval(r.group(1))
            r = re.match('dis_froz_min\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.froz_min = eval(r.group(1))
            r = re.match('dis_win_min\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.win_min = eval(r.group(1))

    def plot_eigenval(self, erange:List[float, float]=None, separate:bool=False, savefig:str='eigenval_dis.png'):
        r'''
        Plot the eigenvalue distribution of each bands. When `separate` is `False`, we will merge the distribution of bands when there is no global gap.

        :param erange: The energy interval you concerned about.
        :param separate: Control whether to merge the bands without the global gap or not. Default is `False`.
        :param savefig: The file name of saved image result.
        '''
        fig, ax = plt.subplots(figsize=(8, 6))

        def label_bar(string, height, rect):
            """Attach a text label on top of bar."""
            ax.annotate(f'{string}',
                        xy=(rect.get_x() + rect.get_width()/2, height),
                        xytext=(0, 0),  # 4 points vertical offset.
                        textcoords='offset points',
                        ha='center', va='bottom')

        # merge bands without global gap
        if not separate:
            idx = self.eband_min[1:] > self.eband_max[:-1]
            idx = np.nonzero(idx)[0]

            if self.nbnds_excl:
                idx = idx[idx > (self.nbnds_excl - 1)]
                idx_min, idx_max = [self.nbnds_excl]+list(idx+1), list(idx)+[self.nbnds-1]
            else:
                idx_min, idx_max = [0]+list(idx+1), list(idx)+[self.nbnds-1]
            
            eplot_min, eplot_max = self.eband_min[idx_min], self.eband_max[idx_max]

            ymin, ymax = erange if erange else (-1e4, 1e4)

            for i, (emin, emax) in enumerate(zip(eplot_min, eplot_max)):
                if emax > ymin and emin < ymax:
                    rect = ax.bar(i, emax-emin, width=0.6, bottom=emin-self.efermi, color='b')
                    label_bar(f'{idx_min[i]-self.nbnds_excl}', emax, rect[0])

        # dont merge bands without global gap
        else:
            idx = np.arange(self.nbnds_excl, self.nbnds-1, self.ndeg) if self.nbnds_excl else np.arange(0, self.nbnds-1, self.ndeg)
            ymin, ymax = erange if erange else (-1e4, 1e4)
            for i in idx:
                emin, emax = self.eband_min[i], self.eband_max[i]
                if emax > ymin and emin < ymax:
                    rect = ax.bar(i, emax-emin, width=0.6, bottom=emin-self.efermi, color='b')
                    label_bar(f'{i-self.nbnds_excl}', emax, rect[0])

        ax.set(ylabel='Energy / eV')
        ax.grid()

        plt.savefig(savefig, dpi=200, bbox_inches='tight', transparent=True)

    def report_eigenval(self, erange:List[float, float]=None, separate:bool=False):
        r'''
        Print the table of the eigenvalue distribution of each bands. When `separate` is `False`, we will merge the distribution of bands when there is no global gap.

        :param erange: The energy interval you concerned about.
        :param separate: Control whether to merge the bands without the global gap or not. Default is `False`.
        '''
        print(f'EFERMI: {self.efermi: 2.6f}')
        print('--------------------------------')
        print('Band No.     EMIN        EMAX')
        print('--------------------------------')
        # merge bands without global gap
        if not separate:
            idx = self.eband_min[1:] > self.eband_max[:-1]
            idx = np.nonzero(idx)[0]

            if self.nbnds_excl:
                idx = idx[idx > (self.nbnds_excl - 1)]
                idx_min, idx_max = [self.nbnds_excl]+list(idx+1), list(idx)+[self.nbnds-1]
            else:
                idx_min, idx_max = [0]+list(idx+1), list(idx)+[self.nbnds-1]
            
            eplot_min, eplot_max = self.eband_min[idx_min], self.eband_max[idx_max]

            ymin, ymax = erange if erange else (-1e4, 1e4)

            for i, (emin, emax) in enumerate(zip(eplot_min, eplot_max)):
                if emax > ymin and emin < ymax:
                    print(f'{idx_min[i]:3d}~{idx_max[i]:3d}  {emin:+10.5f}  {emax:+10.5f}')

        # dont merge bands without global gap
        else:
            idx = np.arange(self.nbnds_excl, self.nbnds-1, self.ndeg) if self.nbnds_excl else np.arange(0, self.nbnds-1, self.ndeg)
            ymin, ymax = erange if erange else (-1e4, 1e4)
            for i in idx:
                emin, emax = self.eband_min[i], self.eband_max[i]
                if emax > ymin and emin < ymax:
                    print(f'  {i:3d}    {emin:+10.5f}  {emax:+10.5f}')
        print('--------------------------------')
    
    def count_states_most(self, erange:List[float, float]) -> int:
        r'''
        Return maximium number of states inside the energy interval defined in `erange` on all calculated k-points. Used in `dis_froz_max` and `dis_froz_min` check.

        :param erange: The energy interval you concerned about.
        '''
        emin, emax = erange
        mask = np.logical_and(self.eband_min <= emax, self.eband_max >= emin)
        return sum(mask)

    def count_states_least(self, erange:List[float, float]) -> int:
        r'''
        Return the minimum number of states inside the energy interval defined in `erange` on all calculated k-points. Used in `dis_win_max` check.

        :param erange: The energy interval you concerned about.
        '''
        emin, emax = erange
        mask = np.logical_and(emin <= self.ir_ebands[0], self.ir_ebands[0] <= emax)
        return np.sum(mask, axis=1).min()

    def suggest_win_max(self, emin:float, nwann:int=None) -> float:
        r'''
        Return lower limit of `dis_win_max` for given `dis_win_min` and number of WFs. If there is no `nwann` input, we will use `self.nwann` from `.yaml` config file instead.  The results are obtained from largest eigenvalues with its band index equal to `#WFs + i + 1`. `i` is the band index of first band with all energy larger than given `emin`.
        
        In order to avoid the insufficient sampling of k-points in Brilliouin zone, we manually add `self.eps` 

        :param emin: `dis_win_min`.
        :param nwann: Number of Wannier functions. If no input from arguments, we will use `self.nwann` from `.yaml` config file. Default is `None`.
        '''
        # TODO handle error case for no enough states when `emin` is too large
        nwann = nwann if nwann else self.nwann
        mask = self.eband_min >= emin
        idx = np.argmax(mask) + nwann - 1
        if idx >= self.nbnds:
            print(f'There is no enough states for {emin} with {nwann} WFs!')
            res = int(self.emax) + 1.
        else:
            res = self.eband_max[idx] + self.eps
        return res

    def suggest_froz_min(self, emax:float, nwann:int=None) -> float:
        r'''
        Return lower limit of `dis_froz_min` for given `dis_froz_max` and `nwann`.

        :param emax: `dis_froz_max`.
        :param nwann: Number of Wannier functions. If no input from arguments, we will use `self.nwann` from `.yaml` config file. Default is `None`.
        '''
        nwann = nwann if nwann else self.nwann
        mask_emax = self.eband_min <= emax
        idx = np.argmin(mask_emax) - nwann - 1
        res = int(self.emin) - 1. if idx < 0 else self.eband_max[idx] + self.eps
        return res

    def suggest_froz_max(self, emin:float, nwann:int=None) -> float:
        r'''
        Return upper limit of `dis_froz_max` for given `dis_froz_min` and `nwann`.

        :param emax: `dis_froz_min`.
        :param nwann: Number of Wannier functions. If no input from arguments, we will use `self.nwann` from `.yaml` config file. Default is `None`.
        '''
        nwann = nwann if nwann else self.nwann
        mask_emin = self.eband_max >= emin
        idx = np.argmax(mask_emin) + nwann
        res = int(self.emax) + 1. if idx >= self.nbnds else self.eband_min[idx] - self.eps
        return res

    def get_dis_froz_df(self, erange:List[float, float]) -> pd.DataFrame:
        r'''
        Return `pd.DataFrame` with suggested frozen window inside the given energy interval.

        :param erange: The energy interval you concerned about.
        '''
        N = self.count_states(erange)
        print(f'There are {N} states in {erange} with Fermi level at {self.efermi}.')
        emin, emax = erange
        dN = N - self.nwann
        froz_min_list, froz_max_list = [], []

        if self.nwann <= 0:
            print(f'Please input vaild number of WF, now is {self.nwann}.')
            return pd.DataFrame(columns=['dis_froz_min', 'dis_froz_max'])

        elif dN > 0:
            print('Suggest dis_froz_min & dis_froz_max as following:')
            print(f'nwann: {self.nwann}    degenercy: {self.ndeg}    Fermi: {self.efermi:12.6f}')

            for i in range(1, dN + 1, self.ndeg):
                # First get froz_max for nwann = nwann_input + i, then get froz_min for nwann = nwann_input and froz_max
                froz_max = self.suggest_froz_max(emin, nwann=self.nwann+i)
                froz_min = self.suggest_froz_min(froz_max)
                froz_max_list.append(froz_max)
                froz_min_list.append(froz_min)

            # number of missing states between `emin` and lowest `dis_froz_min`
            num_missing = self.count_states((emin, min(froz_min_list)))
            df = pd.DataFrame(zip(froz_min_list, froz_max_list), columns=['dis_froz_min', 'dis_froz_max'])

            if num_missing > 0:
                print(f'\nWANRING: There are states between given `emin`: {emin} and lowest `dis_froz_min`: {min(froz_min_list)}. Please carefully treat the suggestion of dis_froz_min / dis_froz_max and check energy range of each bands again. This situation usually happens in no-SOC system with many denegeracy point. But we still want to give you some useful energy window information with states less than number of WFs.\n')
                for i in range(self.nwann - num_missing + 1, self.nwann, 1):
                    new = pd.DataFrame({"dis_froz_min" : emin,
                                        "dis_froz_max" : self.suggest_froz_max(emin, nwann=i)}, index=[1])
                    df = df.append(new, ignore_index=True)

            return df
        else:
            return pd.DataFrame(columns=['dis_froz_min', 'dis_froz_max'])

    def edit_win(self, dis_dict:Dict[str: float]):
        r'''
        Edit `self._sys.win` file with input dis windows from `dis_dict`.
        '''
        with open(self._win, 'r+') as file:
            lines = file.readlines()

        for idx, line in enumerate(lines):
            for key in dis_dict.keys():
                m = re.match(r'^( |\t)*VARIABLE'.replace('VARIABLE', key), line)
                # match and replace
                if m:
                    lines[idx] = '{0:<17s} = {1}\n'.format(m.group(0), dis_dict[key])

        with open(self._win, 'w+') as file:
            file.writelines(lines)

    def evaluate(self, mode:str='AbAk') -> float:
        r"""
        Calculate the difference between VASP band and Wannier90 band.

        self.dEs     : summation of band difference including weight factor from kernel function over all k-points.
        self.dEs_max : maximium of band difference includeing weight factor from kernel function over all k-points.

        :param mode: 'AbAk', 'AbMk', 'MbAk', 'MbMk'. The abberation of `A` and `M` means `Average` and `Maximum`. And the abberation of `b` and `k` means bands and k-points. During the test, `AbAk` gives the best result and it becomes the default mode. It might need to switch to other evaluation in some special cases.
        """
        vkk, vee = self._parse_dat(f'{self.config.vasp_bnd}')
        wkk, wee = self._parse_dat(f'{self.config.w90_bnd}')
        kernel = self.config.kernel
        nbnds, _ = wee.shape                    # num of bands in wannier90
        w2v_ratio = np.max(vkk) / np.max(wkk)   # Theoretically, it should be 2 * pi or 1

        # Use eigenvalues in first k-point to align W90 data and VASP data
        ve, we = vee[:, 0], wee[:, 0]
        Nv, Nw = len(ve), len(we)
        diff = np.array([np.sum(np.abs(ve[i:i+Nw] - we)) for i in range(Nv-Nw+1)])
        nbnds_excl = np.argmin(diff)

        dEs, wgts, dEs_max, wgts_max = [], [], [], []

        # mask of VASP data
        diff_vkk = vkk[1:] - vkk[:-1]
        vmask = [True] + list(np.logical_not(diff_vkk < 1e-7))
        # mask of W90 data
        diff_wkk = wkk[1:] - wkk[:-1]
        wmask = list(np.logical_not(diff_wkk < 1e-7)) + [True]

        for i in range(nbnds):
            # Interpolation need to remove the duplicates
            # REF: https://stackoverflow.com/questions/12054060/scipys-splrep-splev-for-python-interpolation-returns-nan

            tck = interpolate.splrep(wkk[wmask], wee[i][wmask])
            fit_wee_i = interpolate.splev(vkk[vmask] / w2v_ratio, tck)
            vee_i = vee[i + nbnds_excl][vmask]

            dEi = fit_wee_i - vee_i
            # Using filter to smooth the data and remove
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
            dEi =  np.abs(savgol_filter(dEi, 9, 2))

            # ! AVERAGE DISTANCE
            dEs.append(np.sum(kernel(vee[i + nbnds_excl][vmask]) * dEi) / len(dEi) * 1000)
            wgts.append(np.sum(kernel(vee[i + nbnds_excl][vmask])) / len(dEi))
            # ! MAX DISTANCE
            dEs_max.append(np.max(kernel(vee_i) * dEi) * 1000)
            wgts_max.append(np.max(kernel(vee_i)))

        # dEs     store sum over k with weighted dEi
        # dEs_max store maximum over k with weighted dEi
        self.dEs    , self.wgts     = np.array(dEs)    , np.array(wgts)
        self.dEs_max, self.wgts_max = np.array(dEs_max), np.array(wgts_max)

        if mode == 'AbAk':      # Average over bands. Average over kpoints
            return np.sum(self.dEs) / np.sum(self.wgts)
        elif mode == 'AbMk':    # Average over bands. Choose maximum over kpoints
            return np.sum(self.dEs_max) / np.sum(self.wgts_max)
        elif mode == 'MbAk':    # Choose maximum over bands. Average over kpoints
            return max(self.dEs)
        elif mode == 'MbMk':    # Choose maximum over bands and kpoints
            return max(self.dEs_max)
        else:
            raise ValueError(f'Unrecognized mode: {mode}! Please check your input.')

    def total_spread(self) -> ArrayLike:
        r'''
        Return spread message from `.wout` file.
        '''
        conv_str = os.popen(f'grep CONV {self._sys}.wout').read().split('\n')
        if len(conv_str) > 4:
            conv = conv_str[3:-1]
            spread = np.array([eval(c.split()[3]) for c in conv])
            return spread
        else:
            return None

    def opt_dis_dict(self, opt_input:ArrayLike) -> Dict[str: float]:
        r"""
        generate optimized `dis windows` dict from input array

        :param opt_input: List of all energy values to be optimized with default order in `self.opt_ini_dis`.
        """
        if len(opt_input) != len(self.config.opt_ini_dis):
            raise ValueError(f'There are only {len(self.config.opt_ini_dis)} parameters but {len(opt_input)} input!')
        return OrderedDict([(k, v) for v, k in zip(opt_input, self.config.opt_ini_dis.keys())])

    def total_dis_dict(self, opt_input:ArrayLike) -> Dict[str: float]:
        r'''
        Get full `dis_windows` with value with replacing `None` to `\pm self.inf`

        :param opt_input: List of all energy values to be optimized with default order in `self.opt_ini_dis`.
        '''
        if len(opt_input) != len(self.config.opt_ini_dis):
            raise ValueError(f'There are only {len(self.config.opt_ini_dis)} parameters but {len(opt_input)} input!')
        d = copy.copy(self.config.ini_dis)
        for k in self.config.ini_dis.keys():
            if d[k] == None and k[-3:] == 'max': d[k] =  self.inf
            if d[k] == None and k[-3:] == 'min': d[k] = -self.inf
        for i, k in enumerate(self.config.opt_ini_dis.keys()):
            d[k] = opt_input[i]
        return d

    def check_total_dis_dict(self, opt_input:ArrayLike) -> bool:
        r"""
        Check input value. When number of states inside frozen windows is greater than number of WFs or number of states inside dis windows is less than number of WFs, the function will return `False` with WARNING.

        :param opt_input: List of all energy values to be optimized with default order in `self.opt_ini_dis`.
        """
        dis = self.total_dis_dict(opt_input)
        n_froz = self.count_states_most((dis['dis_froz_min'], dis['dis_froz_max']))
        n_win  = self.count_states_least((dis['dis_win_min'], dis['dis_win_max']))
        if n_froz > self.nwann:
            logger.warning(f"There are {n_froz:3d} states inside frozen window but {self.nwann:3d} WFs given! \nFix `dis windows` as following".ljust(80, ' '))
            return False
        if n_win < self.nwann:
            logger.warning(f"There are {n_win:3d} states inside dis window but {self.nwann:3d} WFs given!\nFix `dis windows` as following".ljust(80, ' '))
            return False
        return True

    def show_total_dis_dict(self, opt_input:ArrayLike):
        r'''
        Display full dis energy windows with value and range via `logging.Logger`.

        :param opt_input: List of all energy values to be optimized with default order in `self.opt_ini_dis`.
        '''
        opt_input_dict = self.opt_dis_dict(opt_input)
        logger.info(f"+{'ENERGY-WINDOWS'.center(78, '-')}+")
        logger.info("|                Energy Windows     Value      Range                           |")
        logger.info("|                --------------     -----      -----                           |")
        for idx, key in enumerate(self.config.ini_dis.keys()):
            v = self.config.ini_dis[key]
            if v != None:
                fix = ' ' if self.config.opt_dis[key] else 'F'
                if key in opt_input_dict.keys():
                    v = opt_input_dict[key]
                logger.info("|              {0} {1:<12s}      {2: <+9.4f}   {3:<8} to {4:<18}  |".format(fix, key, v, self.config.lbs[idx], self.config.ubs[idx]))
        logger.info(f"+{''.center(78, '-')}+")

    def rational_opt_input(self, opt_input:ArrayLike) -> Dict[str: float]:
        r'''
        Return rationalized full dis windows dict from input. `Wannier90` requires at most #WFs of states inside the frozen energy window and at least #WFs of states inside the dis energy window. We treat the input `dis_win_min` as beginning to generate `dis_froz_max`, `dis_froz_min` and `dis_win_max` which satisifying the requirement of `Wannier90`.
        
        :param opt_input: List of all energy values to be optimized with default order in `self.opt_ini_dis`.
        '''
        dis = self.total_dis_dict(opt_input)
        emin = dis['dis_win_min']
        emax = self.suggest_froz_max(emin)
        dis['dis_froz_max'] = emax
        dis['dis_froz_min'] = self.suggest_froz_min(emax)
        dis['dis_win_max']  = self.suggest_win_max(emin)
        return [dis[k] for k in self.config.opt_ini_dis.keys()]

    def show_spread(self):
        r'''
        Display spreading via `loggind.Logger`.
        '''
        spread = self.total_spread()
        len_unit_spread = max(spread / 20)
        len_spread = (spread / len_unit_spread).astype(int)
        N2 = (len(spread) + 1)// 2
        logger.info(f"Spread (Ang^2) in `{self._sys}.wout`:".ljust(80, ' '))
        logger.info(f"+-----+{'-'*32}++-----+{'-'*32}+")
        logger.info(f"|   i |{' Spread':32s}||   i |{' Spread':32s}|")
        logger.info(f"+-----+{'-'*32}++-----+{'-'*32}+")
        for i in range(N2):
            il, ir = i, i + N2
            sl = f'{self._block * len_spread[il]} {spread[il]:.1f}'
            sr = f'{self._block * len_spread[ir]} {spread[ir]:.1f}' if ir < len(spread) else ''
            logger.info(f'|{il+1:4d} |{sl:<32s}||{ir+1:4d} |{sr:<32s}|')
        logger.info(f"+-----+{'-'*32}++-----+{'-'*32}+")

    def show_dEs(self):
        r'''
        Display the maximium and average band difference via `loggind.Logger`.

        Notice: The maximuim value of difference has multiplied with kernel factor. It will have some influence when using Gaussian function as kernel function. The average value of difference has removed the kernel factor. And we only counts the value over non-zero k-points with window function.
        '''
        dEs     = self.dEs /(self.wgts+1e-6)  # Correct dEs with dividing weight
        dEs_max = self.dEs_max
        len_unit_max = max(dEs_max / 22)
        len_unit_ave = max(dEs     / 22)
        len_dEs_max = (dEs_max / len_unit_max).astype(int)
        len_dEs_ave = (dEs     / len_unit_ave).astype(int)
        logger.info(f"+-----+{'-'*32}++-----+{'-'*32}+")
        logger.info(f"|   i |{' MAX DIFF (meV)':32s}||   i |{' AVERAGE DIFF (meV)':32s}|")
        logger.info(f"+-----+{'-'*32}++-----+{'-'*32}+")
        for i in range(self.nwann):
            max_str_i = f'{self._block * len_dEs_max[i]} {dEs_max[i]:.2f}'
            ave_str_i = f'{self._block * len_dEs_ave[i]} {dEs[i]:.2f}'
            logger.info(f'|{i+1:4d} |{max_str_i:<32s}||{i+1:4d} |{ave_str_i:<32s}|')
        logger.info(f"+-----+{'-'*32}++-----+{'-'*32}+")

    def _parse_dat(self, datfile:str) -> List[ArrayLike, ArrayLike]:
        r'''
        return the k-distance array and eigenvalue array from `.dat` file.
        '''
        data = pd.read_csv(datfile, header=None, sep=r'\s+', comment='#')
        k_idxs = np.where(data[0]==np.min(data[0]))[0]
        nk = k_idxs[1] - k_idxs[0]
        kk = np.array(data[0])[:nk]
        EE = np.array(data[1]).reshape(-1, nk)
        return kk, EE
