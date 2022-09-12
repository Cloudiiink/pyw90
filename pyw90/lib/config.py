import yaml, os, re
import numpy as np
from collections import OrderedDict
from scipy.optimize import Bounds, LinearConstraint

from utility.utility import gaussian, unit, parse_kernel
from utility.utility import _replace_str_boolean, _replace_str_none

import logging
logger = logging.getLogger(__name__)

class Config():
    r'''
    config
    '''
    def __init__(self, yaml_file:str=None):
        r'''
        Parse input parameters from `.yaml` file. If there is no specific `.yaml` file given, class will treat the first `.yaml` file as config file.
        '''
        if not yaml_file:
            yaml_file_list = [f for f in os.listdir() if f.endswith('.yaml')]
            if len(yaml_file_list) != 0:
                yaml_file = yaml_file_list[0]
            else:
                raise ValueError('There is no `.yaml` file in this folder!')
        with open(yaml_file,'r', encoding='utf-8') as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        
        self.win      = data['winfile']
        self.run      = data['runfile']
        self.vasp_bnd = data['vasp_band_file']
        self.w90_bnd  = data['wann_band_file']
        self.job_name = data['jobname']
        self.usr_name = data['username']
        self.efermi   = data['efermi']
        self.nwann    = data['nwann']
        # Minimization Related parameters
        self.method   = data['method']
        self.tol      = float(data['tol'])
        self.maxiter  = int(data['maxiter'])
        self.niter    = 0

        kernel_str, mid, width = data['kernel']
        self.kernel_str = f'{kernel_str}, mid: {mid}, width: {width}'
        self.kernel = parse_kernel(kernel_str, mid, width)

        self.dis_keys = ['dis_win_max', 'dis_froz_max', 'dis_froz_min', 'dis_win_min']
        self.ini_dis = OrderedDict([(key, _replace_str_none(data['ini_dis'])[key]) 
                        for key in self.dis_keys])
        self.opt_dis = OrderedDict([(key, data['opt_dis'][key]) 
                        for key in self.dis_keys])
        self.opt_dis = _replace_str_boolean(self.opt_dis)
        
        bounds = _replace_str_none(data['bounds'])
        
        # replace the boundary using +/- np.inf
        self.lbs = np.array([b[0] if b[0] else -np.inf for b in bounds])
        self.ubs = np.array([b[1] if b[1] else +np.inf for b in bounds])
        self.mask = list(self.opt_dis.values())
        self.bounds = Bounds(lb=self.lbs, ub=self.ubs, keep_feasible=True)
        self.opt_bounds = Bounds(lb=self.lbs[self.mask], ub=self.ubs[self.mask],
                                 keep_feasible=True)
        # dis windows dict with keys and values to be optimized
        self.opt_ini_dis = OrderedDict([(k, self.ini_dis[k]) 
                            for i, k in enumerate(self.dis_keys) if self.mask[i] ])
                            
        # TODO `LinearConstraint` is used in `COBYLA` method, but Constraint option `keep_feasible` is ignored by this method now (scipy==1.9.1).
        self.opt_constraint = LinearConstraint(np.diag([1] * len(self.opt_ini_dis)),
                                               self.lbs[self.mask],
                                               self.ubs[self.mask])
        # self.opt_constraint = LinearConstraint(np.diag([1] * len(self.opt_ini_dis)),
        #                                        self.lbs[self.mask],
        #                                        self.ubs[self.mask],
        #                                        keep_feasible=True)

        self.num_print_check = int(data['num_print_check'])
        self.check_time = int(data['check_time'])
        self.display = _replace_str_boolean(data['display'])

        self.check_input()

    def check_input(self):
        # check files exist or not
        for file in [self.win, self.run, self.vasp_bnd]:
            if not os.path.isfile(file):
                raise FileNotFoundError('File {0} Not Found!'.format(file))

        # TODO : add PBS job queue system support
        # check job name
        pattern = r'^\s*#\s*SBATCH -(\-job\-name|J)\s*=?(.+)\n'
        with open(self.run) as f:
            lines = f.readlines()
        for line in lines:
            m = re.match(pattern, line)
            if m:
                run_file_job_name = m.group(2).strip()
                break
        if self.job_name != run_file_job_name:
            raise ValueError('The jobname is inconsistent! In batch file the jobname is {0}. But it\'s {1} in settings.yaml file instead.'.format(run_file_job_name, self.job_name))

        # check user name
        from getpass import getuser
        local_usr_name = getuser()
        if self.usr_name != local_usr_name:
            raise ValueError('The username is inconsistent! Your username is {0}. But it\'s {1} in settings.yaml file instead.'.format(local_usr_name, self.usr_name))

        # check Fermi level
        if 'vasprun.xml' in os.listdir():
            efermi_str = os.popen(f'grep fermi vasprun.xml').read().strip()
            m = re.match('.+ ((\-|\+)?\d+(\.\d+)?) .+', efermi_str)
            efermi = float(m.groups()[0])
            if np.abs(efermi - self.efermi) < 1e-6:
                raise ValueError('The Fermi level {0} searched from `vasprun.xml` is inconsistant with your input {1} in `.yaml`.'.format(efermi, self.efermi))

    def info(self):
        r'''
        Visualize the input through `logging.Logger`.
        '''
        logger.info(f"+{'FILE-LOCATION'.center(78, '-')}+")
        logger.info(f"|{' Wannier90 Input File':>35s}   :   {self.win:<35s} |")
        logger.info(f"|{'   Slurm Batch Script':>35s}   :   {self.run:<35s} |")
        logger.info(f"|{'   DFT Band Structure':>35s}   :   {self.vasp_bnd:<35s} |")
        logger.info(f"|{'   W90 Band Structure':>35s}   :   {self.w90_bnd:<35s} |")
        logger.info(f"+{''.center(78, '-')}+")

        logger.info(f"+{'SYSTEM-INFORMATION'.center(78, '-')}+")
        logger.info(f"|{'              Jobname':>35s}   :   {self.job_name:<35s} |")
        logger.info(f"|{'             Username':>35s}   :   {self.usr_name:<35s} |")
        logger.info(f"|{'         Fermi Energy':>35s}   :   {self.efermi:<35f} |")
        logger.info(f"|{'        number of WFs':>35s}   :   {self.nwann:<35d} |")
        logger.info(f"|{'               Method':>35s}   :   {self.method:<35s} |")
        logger.info(f"|{'Tolerance of minimize':>35s}   :   {self.tol:<35f} |")
        logger.info(f"|{' Max Iter of minimize':>35s}   :   {self.maxiter:<35d} |")
        logger.info(f"|{'               Kernel':>35s}   :   {self.kernel_str:<35s} |")
        logger.info(f"+{''.center(78, '-')}+")

        logger.info(f"+{'INITIAL-ENERGY-WINDOWS'.center(78, '-')}+")
        logger.info("|                Energy Windows     Value      Range                           |")
        logger.info("|                --------------     -----      -----                           |")
        for idx, key in enumerate(self.ini_dis.keys()):
            v = self.ini_dis[key]
            if v != None:
                fix = ' ' if self.opt_dis[key] else 'F'
                logger.info("|              {0} {1:<12s}      {2: <+9.4f}   {3:<8} to {4:<18}  |".format(fix, key, v, self.lbs[idx], self.ubs[idx]))
        logger.info(f"+{''.center(78, '-')}+")

        logger.info(f"+{'OUTPUT'.center(78, '-')}+")
        logger.info(f"|{'Num Check Print Once':>35s}   :   {self.num_print_check:<35d} |")
        logger.info(f"|{'  Check Period (sec)':>35s}   :   {self.check_time:<35d} |")
        logger.info(f"|{'   Display Spreading':>35s}   :   {str(self.display['spread']):<35s} |")
        logger.info(f"|{'   Display band diff':>35s}   :   {str(self.display['diff']):<35s} |")
        logger.info(f"+{''.center(78, '-')}+")