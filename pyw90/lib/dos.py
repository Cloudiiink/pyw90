from math import sqrt
import pandas as pd
import numpy as np
from numpy.typing import ArrayLike
from typing import Tuple
import warnings
from collections import Counter
from functools import reduce
import os

# Calculate the percentage of given sites and given orbitals inside the energy interval
from scipy.interpolate import splrep, splev
from scipy.integrate import fixed_quad

from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.core import structure

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# pyw90
from pyw90.utility.utility import num2str, str2num, get_id_from_specie
from pyw90.utility.utility import bc

def check_orb_name(func):
    def __wrapper__(*args, **kwargs):
        if args[1].strip() not in args[0].orbital_names:
            raise ValueError(f'Unknown orbital name: {args[1]}')
        return func(*args, **kwargs)
    return __wrapper__

def check_orb_type_name(func):
    def __wrapper__(*args, **kwargs):
        if args[1].strip() not in args[0].orbital_type_names:
            raise ValueError(f'Unknown orbital type name: {args[1]}')
        return func(*args, **kwargs)
    return __wrapper__

class DOS():
    r"""
    Class offer orbitals and DOS (density of states) manipulation.
    """
    # Orbital Notation for VASP since 5.4.4 also have some difference compare to
    # list `orbital_names` below. Here we changed some of the notations 
    # (x2-y2 -> dx2) for comparision between `VASP` and `Wannier90`.

    # See https://www.vasp.at/forum/viewtopic.php?t=18044
    # https://www.vasp.at/forum/viewtopic.php?f=4&t=6924

    # TODO Need to check for f-electron
    # v5.4.4        v5.4.1      pymatgen
    # x2-y2    <-->   dx2   <-->   dx2
    # fy3x2    <-->   f1    <-->   f_3
    # fxyz     <-->   f2    <-->   f_2
    # fyz2     <-->   f3    <-->   f_1
    # fz3      <-->   f4    <-->   f0
    # fxz2     <-->   f5    <-->   f1
    # fzx2     <-->   f6    <-->   f2
    # fx3      <-->   f7    <-->   f3

    # pymatgen
    pymat_orbital_names = ['s', 
                          'py', 'pz', 'px', 
                          'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                          'f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']
    vasp_orbital_names  = ['s', 
                          'py', 'pz', 'px', 
                          'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2',
                          'fy3x2', 'fxyz', 'fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3']
    orbital_type_names  = ['s', 'p', 'd', 'f']
    # ordered according to `mr` value
    w90_orbital_names   = ['s',
                          'pz', 'px', 'py',
                          'dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy',
                          'fz3', 'fxz2', 'fyz2', 'fz(x2-y2)', 'fxyz', 'fx(x2-3y2)', 'fy(3x2-y2)']
    orbital_names       = pymat_orbital_names
    vasp2w90            = {'s'    : 's',
                          'py'    : 'py',
                          'pz'    : 'pz',
                          'px'    : 'px', 
                          'dxy'   : 'dxy',
                          'dyz'   : 'dyz',
                          'dz2'   : 'dz2',
                          'dxz'   : 'dxz',
                          'dx2'   : 'dx2-y2',
                          'fy3x2' : 'fy(3x2-y2)', 
                          'fxyz'  : 'fxyz',
                          'fyz2'  : 'fyz2', 
                          'fz3'   : 'fz3',
                          'fxz2'  : 'fxz2',
                          'fzx2'  : 'fz(x2-y2)' , 
                          'fx3'   : 'fx(x2-3y2)'}
    pymat2w90           = dict([(k, v) 
                              for k, v in zip(pymat_orbital_names, vasp2w90.values())])

    @classmethod
    @check_orb_name
    def get_mr(cls, name):
        r"""
        Return mr quantum number for given Orbital name (e.g. px)
        """
        name = cls.vasp2w90[name]
        value = cls.w90_orbital_names.index(name) + 1
        return value - int(sqrt(value))**2
    
    @classmethod
    @check_orb_name
    def get_l(cls, name):
        r"""
        Return l quantum number for given Orbital name (e.g. px)
        """
        return cls.orbital_type_names.index(name[0])

    @classmethod
    def is_orbital_type(cls, name):
        r"""
        Return whether `name` is Orbital_Type (s, p , d, f) or not.
        """
        return name in cls.orbital_type_names

    @classmethod
    def is_orbital(cls, name):
        r"""
        Return whether `name` is Orbital (e.g. px, dz2 )or not.

        **Note:** Using `pymatgen` notation
        """
        return name in cls.orbital_names

    @classmethod
    @check_orb_name
    def get_orbital(cls, name):
        r"""
        Return `Orbital` from given name. The indices are basically the order in
        which the orbitals are reported in VASP and has no special meaning.

        - `s   = 0`
        - `py  = 1`  `pz  =  2`  `px  =  3`
        - `dxy = 4`  `dyz =  5`  `dz2 =  6`  `dxz =  7`  `dx2 = 8`
        - `f_3 = 9`  `f_2 = 10`  `f_1 = 11`  `f0  = 12`  `f1 = 13`  `f2 = 14`  `f3 = 15`
        """
        return Orbital(cls.orbital_names.index(name))

    @classmethod
    @check_orb_type_name
    def get_orbital_type(cls, name):
        r"""
        Return `Orbital` from given name. 

        s = 0, p = 1, d = 2, f = 3
        """
        return OrbitalType(cls.orbital_type_names.index(name))

    @classmethod
    @check_orb_type_name
    def get_orbitals_from_type(cls, name):
        indices = [i for i, s in enumerate(cls.orbital_names) if s[0]==name]
        return [Orbital(i) for i in indices]

    @classmethod
    def get_w90_orbs_str(cls, orb_list: list[str]) -> str:
        r"""
        Convert `orb_list` (e.g., ['s', 'p', 'dx2'] in `pymatgen` notation) to format used in `Wannier90`. (e.g. l=0;l=1;dx2-y2)
        """
        l_orbs_str = []
        for orb in orb_list:
            if cls.is_orbital_type(orb):
                l_orbs_str.append(f'l={cls.get_orbital_type(orb).value}')
            else:
                l_orbs_str.append(cls.pymat2w90[orb])
        return ';'.join(l_orbs_str)

    @staticmethod
    def get_projections_selected_str(res_df: pd.DataFrame, struc: structure,
                                  connect: str='-', sep :str='|') -> str:
        r"""
        Convert `res_df` with simplified projection information to formatted string for `--extra`` argument to input for pyw90.

        e.g. `'Ga,0,1-3;As,1,1-3|5'`
        """
        res = []
        for spec, site, orbs in zip(res_df.species, res_df.site, res_df.orb):
        # we only merge all sites having the same orbital selection with -1
        # so the `site` column should be `integer`
            orb_res = []
            for orb in orbs:
                orb = orb.strip()
                if DOS.is_orbital_type(orb):
                    orb_res += DOS.get_orbitals_from_type(orb)
                elif DOS.is_orbital(orb):
                    orb_res.append(DOS.get_orbital(orb))
                else:
                    raise ValueError(f'Unknown Orbital: {orb}')
            orb_str = num2str([orb.value for orb in orb_res],
                            connect=connect, sep=sep)

            if site == -1:
                site = get_id_from_specie(struc, spec)
                site = num2str(site, connect=connect, sep=sep)

            res.append(f'{spec},{site},{orb_str}')
        return ';'.join(res)

    @classmethod
    def get_projections_w90_str(cls, res_df: pd.DataFrame, struc: structure) -> str:
        r"""
        Return the projection card for `Wannier90`

        ### Parameters
        - `res_df`: Obtained from `dos_analysis_df` with **simplified** projection given by: `species`, `site` and `orb`
        - `struc`: Crystal structure
        
        ### Note
        We use `Wannier90`==1.2 to interface with `VASP`. So the format for projection card is very limited. 
        You can use site (with fractional coordinates or Cartesian coordinates) to define the location of single projection.
        You can also use species to define. (For more information, please check the user guide of `Wannier90`)

        `Wannier90` also implements the selected column of the density matrix (SCDM) method to automaticaly generate
        the initial projection (See documentation of `auto_projections` block). Here we do not use this feature.
        """
        l_w90_str = []
        for spec, site, orbs in zip(res_df.species, res_df.site, res_df.orb):
            if site == -1:
                l_w90_str.append(f"{spec}:{cls.get_w90_orbs_str(orbs)}")
            else:
                site_string = ','.join([f"{i:.6f}" for i in struc[site].frac_coords])
                # ? `f=` for fraction coordinates
                l_w90_str.append(f"f={site_string}:{cls.get_w90_orbs_str(orbs)}")
        return '\n'.join(l_w90_str)

    @classmethod
    def parse_projections_from_selected(cls, s: str) -> set[str]:
        r"""
        Convert input string to `set` with `key_string` format (`spec`_`struc_id`_`orb_name`).
        We  use `1-4|6` to represent `1,2,3,4,6` both in site_id and orbital_id.

        **Note**: `pymatgen` notation used in returned `key_string`.
        """

        if len(s) == 0:
            return []

        res = []
        for l in s.split(';'):
            spec, pos, orbs = l.split(',')
            pos, orb = str2num(pos), str2num(orbs) 
            for p in pos:
                for o in orb:
                    res.append(f'{spec}_{p}_{cls.orbital_names[o]}')
        return set(res)

    @staticmethod
    def get_dos_integral(e: ArrayLike, dos: ArrayLike,
                         window: Tuple[float, float]) -> float:
        r"""
        Return the integral of `dos` value inside the `window`.

        The left value of `window` is required to be greater than the minimium of `e` 
        and the right value of `window` is required to be less than the maximium of `e`.
        
        ### Note

        Here we use interpolation to calculate the integral. 
        If the input integral variable is outside the whole energy range,
        We will calculate the energy having definition and treat the DOS outside as 0.
        """
        # ! Throw the `warning` at `get_dos_df` method
        # if min(e) > window[0] or max(e) < window[-1]:
        #     warnings.warn(f"Since some region of the {window} is outside the input variable `e` ranging from {min(e)} to {max(e)}, the integral from interpolation might not be very convincing.")
        left  = max(min(e), window[0])
        right = min(max(e), window[-1])
        dos_cubic_fitting = splrep(e, dos, s=0)
        dos_func = lambda x: splev(x, dos_cubic_fitting, der=0)
        res = fixed_quad(dos_func, left, right)
        return res[0]

    @classmethod
    def get_dos_df(cls, dos_data_total: CompleteDos, left: float, right: float, 
                   spin: Spin=Spin.up) -> pd.DataFrame:
        r"""
        Return `pandas.DataFrame` with columns: `species`, `structure_id`, `orb_id`, `orb_name`, `key_string` and `dos`.

        The value of `dos` has been normalized with the maximum is adjusted to 1. 
        """
        structure = dos_data_total.structure
        num_sites = structure.num_sites
        ee = dos_data_total.energies

        # TODO: use `PROCAR` to generate a result, although it might not be accurate enough.
        if left < min(ee) or right > max(ee):
            warnings.warn(f'CHECK YOUR INPUT! The energies in `EIGENVAL` is ranged from {ee.min()} to {ee.max()}, ' \
                          f'which does not include all the energy ranged from {left} to {right} you input.')

        lspec, lstruc, lorbid, lorbname, lstr, ldos = [list() for _ in range(6)]
        for i in range(num_sites):
            for j in range(len(cls.orbital_names)):
                site, orb = structure[i], Orbital(j)
                dos = dos_data_total.get_site_orbital_dos(site, orb).densities
                if len(dos) > 0:
                    dis = cls.get_dos_integral(ee, dos[spin], [left, right])
                    # how to formulate the key string for every entries
                    key = site.species_string + '_' + str(i) + '_' + cls.orbital_names[j]
                    lspec.append(site.species_string)
                    lstruc.append(i)
                    lorbid.append(j)
                    lorbname.append(cls.orbital_names[j])
                    lstr.append(key)
                    ldos.append(dis)
        dos_df = pd.DataFrame(zip(lspec, lstruc, lorbid, lorbname, lstr, ldos),
                            columns=["species", "structure_id", "orb_id", "orb_name", "key_string", "dos"])
        dos_df['dos'] = dos_df['dos'] / max(dos_df['dos'])      # Renormalized
        return dos_df

    @classmethod
    def get_dos_analysis_df(cls, dos_df: pd.DataFrame, lb: float=0.1) -> Tuple[int, pd.DataFrame]:
        r"""
        Analyze the DOS distribution and recommend the projection according to `lb` (threshold for selecting)

        ### Parameters
        - `dos_df` : `pd.DataFrame` which containg the DOS contribution from each sites and each orbitals. Obtained from `get_dos_df`.
        - `lb` : The threshold (or lower bound) for selecting the projection.

        ### Return
        - number of orbitals selected.
        - `pd.DataFrame` with simplified projection information. Columns: `species`, `site`, `orb`. Usually the value of `site` is the `structure_id`.
        But we will use -1 when all site with same species have the same projections. See example as following:

        ```
        species site  orb
        0      Ga   -1  [p]
        1      As   -1  [p]
        ```
        """
        threshold = lb * dos_df['dos'].max()
        df = dos_df[dos_df['dos'] > threshold]
        n_orb = len(df)
        res_df = pd.DataFrame(columns=["species", "site", "orb"])       # Store the recommendation results

        # Simplify the results using spdf electron configuration
        # Iterate to check all the `species` in the material
        for species in df.species.unique():
            species_sites = df[df["species"]==species].structure_id.unique()

            # Simplify the orbitals using `spd` electron configuration
            # e.g. convert [s, px, py, pz, dx2] to [s, p, dx2]
            for idx in species_sites:
                site_df = df[df["structure_id"]==idx]
                site_orbs = site_df.orb_name.unique()
                for orb in cls.orbital_type_names:
                    orb_list = [i.name for i in cls.get_orbitals_from_type(orb)]
                    if set(orb_list).issubset(set(site_orbs)):
                        site_orbs = list(set(site_orbs) - set(orb_list)) + [orb]
                new = pd.DataFrame({"species": species,
                                    "site": idx,
                                    "orb": [site_orbs]}, index=[1])
                res_df = pd.concat([res_df, new], ignore_index=True)

        # display(res_df)

        # merging the sites that have the same projections
        for species in res_df.species.unique():
            # get structure id
            species_site_id = dos_df[dos_df["species"]==species].structure_id.unique()
            species_site = res_df[res_df["species"]==species].site.unique()
            species_site_orbs = res_df[res_df["species"]==species].orb
            if len(species_site_id) == len(species_site):
                ref_orb = species_site_orbs.to_list()[0]
                is_orb_united = reduce(lambda x, y: x and y, [Counter(ref_orb)==Counter(orb) for orb in species_site_orbs])
                if is_orb_united:
                    # delete all the entries of the `speices` and Add data after merging, use `-1` to represent the site
                    res_df = res_df[res_df['species']!=species]
                    # Sorted the `orb` column
                    sorted_orb = sorted(ref_orb, key=lambda s: cls.orbital_type_names.index(s[0]))
                    new = pd.DataFrame({"species": species,
                                        "site": -1,
                                        "orb": [list(sorted_orb)]})
                    res_df = pd.concat([res_df, new], ignore_index=True)

        

        # display(res_df)
        return n_orb, res_df

    @classmethod
    def get_dos_given_site_orb(cls, dos_data_total: CompleteDos, erange: Tuple[float, float], 
                               key_string: str, spin: Spin=Spin.up) -> float:
        r"""
        Return DOS contribution of `erange` for site and orbital from `key_string` (e.g. `C_1_pz`) via interpolation.
        """
        structure = dos_data_total.structure
        _, site_id, orb = key_string.split('_')
        site = structure[int(site_id)]
        orb = DOS.get_orbital(orb)
        dos = dos_data_total.get_site_orbital_dos(site, orb)
        dis = cls.get_dos_integral(dos.energies, dos.densities[spin], erange)
        return dis

    @classmethod
    def get_dos_from_selected(cls, dos_data_total: CompleteDos, erange: Tuple[float, float], 
                              selected: set[str]) -> float:
        r"""
        Return DOS contribution from selected sites and orbitals inside energy interval interpolation.
        """
        dis = [cls.get_dos_given_site_orb(dos_data_total, erange, key) for key in selected]
        return np.sum(dis)

    @staticmethod
    def plot_dos_df(dos_df: pd.DataFrame, lb: float=1, selected: set=set(),
                    colors: Tuple[str, str]=['brown', 'orange'],
                    savefig: str='dos_analysis.pdf'):
        r"""
        Bar plot of DOS distribution.

        ### Parameters: 
        - `dos_df`: `pandas.DataFrame` from `gen_dos_df` which containing the DOS contribution of each sites and orbitals.
        - `lb`: Percentage of lower bound. Since there are too many data to show, here we only present entries that the contribution is greater than `lb`%.
        - `selected`: Selected entries to show (in key_string format: {species}_{structure id}_{orb_name} )
        - `colors`: Plot with selected orbitals in the first color and not selected orbitals in the last color. The default value is ['brown', 'orange'].
        - `savefig`: Saved figure.
        """
        threshold = lb / 100 * dos_df['dos'].max()

        df = dos_df[dos_df['dos'] > threshold]
        df = df.sort_values(by='dos', ascending=False)
        # color = [colors[0] if row['key_string'] in selected else colors[-1]
        #              for _, row in df.iterrows()]
        color = [colors[0] if row in selected else colors[-1] for row in df['key_string']]
        mask  = (np.array(color) == colors[0])
        bc.cprint(bc.BLUE, f'Plotted with selected orbitals in `{colors[0]}` and non-selected orbitals in `{colors[-1]}`.')
        bc.cprint(bc.BLUE, f'Plotted with a total of {len(df)} orbitals, {sum(mask)} of which are selected.')
        bc.cprint(bc.BLUE, f"Figure should be stored at {os.path.abspath(savefig)}")
        
        import matplotlib.patheffects as pe

        # dynamic figsize and font support
        plt.rcParams['font.family'] = "Open Sans"
        ax = df.plot.barh(x="key_string", y="dos", color=color, figsize=(8, len(df)/32*8))
        ax.set_axisbelow(True)
        label = ax.bar_label(ax.containers[0], fmt='%.3f', padding=5)
        for l in label:
            l.set_path_effects([pe.withStroke(linewidth=8, foreground="w")])
        plt.grid(axis='x')

        # save figure in local folder not in source folder
        # plt.savefig(savefig, dpi=200, bbox_inches='tight')
        # plt.savefig(savefig, dpi=200, bbox_inches='tight', transparent=True)
        plt.savefig(savefig, bbox_inches='tight')

    # TODO: DOS plot
    @staticmethod
    def plot_dos_selected():
        r"""
        Line plot for given projections from `--extra` argument.
        """
        pass

if __name__ == '__main__':
    print(DOS.get_mr('px'), DOS.get_l('px'))
    print(DOS.get_orbitals_from_type('f'))
    print(DOS.get_w90_orbs_str(['s', 'p', 'dx2']))