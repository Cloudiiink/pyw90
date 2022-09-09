import time, os
from scipy.optimize import minimize
import numpy as np
from numpy.typing import ArrayLike

from lib.config import Config
from lib.job import Job
from lib.w90 import W90

import logging
logger = logging.getLogger(__name__)

strtime = time.strftime(r'%Y%m%d%H%M%S', time.localtime())
logfilename = 'log_autow90_{0}.log'.format(strtime)
logging.basicConfig(level=logging.INFO,     # loggingDEBUG
                    filename=logfilename,
                    datefmt='%H:%M:%S',
                    format='%(message)s ! %(asctime)s/%(levelname)s/%(lineno)d')
logger = logging.getLogger(__name__)

def main_task(opt_input:ArrayLike, w90:W90) -> float:
    w90.config.niter += 1
    logger.info(f"\n{''.center(80, '-')}")
    logger.info(f"Round {w90.config.niter} with input `dis_windows`:".ljust(80, ' '))
    w90.show_total_dis_dict(opt_input, logger)

    if not w90.check_total_dis_dict(opt_input, logger):
        opt_input = w90.rational_opt_input(opt_input)
        w90.show_total_dis_dict(opt_input, logger)

    opt_input_dict = w90.opt_dis_dict(opt_input)
    w90.edit_win(opt_input_dict)

    # TODO: run job on local machine (or without queueing system) through `wannier90.x`. Check outputs in `wannier90.wout` file or `wannier90_band.dat` file?
    job = Job(w90.config)
    job.submit()
    logger.info("Successfully submit the task!".ljust(80, ' '))
    time.sleep(5)
    job.check_run_until_stop(logger)
    logger.info('Job <{0}> is done.'.format(w90.config.job_name).ljust(80, ' '))

    mtime = time.ctime(os.path.getmtime(w90.config.w90_bnd))
    logger.info(f'`{w90.config.w90_bnd}` was last modified at {mtime}.'.ljust(80, ' '))
    logger.info('Evaluating the quality ...'.ljust(80, ' '))

    res = w90.evaluate()
    if w90.config.display['spread']: w90.show_spread(logger)
    if w90.config.display['diff']:   w90.show_dEs(logger)
    logger.info('Final Value: {0} meV (with Average of all bands all kpoints)'.format(res).ljust(80, ' '))

    return res

if __name__ == '__main__':
    logger.info("""
                    ──█▀▀█─█  █ ▀▀█▀▀ █▀▀█ ░█  ░█ ▄▀▀▄ █▀▀█
                     ░█▄▄█ █  █   █ ──█──█─░█░█░█ ▀▄▄█ █▄▀█
                     ░█ ░█ ─▀▀▀───▀─  ▀▀▀▀ ░█▄▀▄█──▄▄▀ █▄▄█
                                                                                """)

    config = Config(yaml_file="auto_w90_input.yaml")
    config.info(logger)
    w90 = W90(config)

    logger.info(' Calculation Start! '.center(80, '='))

    opt_ini_input = np.array([v for v in config.opt_ini_dis.values()])

    if not w90.check_total_dis_dict(opt_ini_input, logger):
        opt_ini_input = w90.rational_opt_input(opt_ini_input)
        w90.show_total_dis_dict(opt_ini_input, logger)

    if config.method[0].lower() == 'n':
        res = minimize(main_task, opt_ini_input, args=w90, method='Nelder-Mead', tol=config.tol, bounds=config.opt_bounds, options={'maxiter': config.maxiter})
    elif config.method[0].lower() == 'p':
        res = minimize(main_task, opt_ini_input, args=w90, method='Powell', tol=config.tol, bounds=config.opt_bounds, options={'maxiter': config.maxiter})
    elif config.method[0].lower() == 'c':
        res = minimize(main_task, opt_ini_input, args=w90, method='COBYLA', tol=config.tol, constraints=config.opt_constraint, options={'maxiter': config.maxiter})
    else:
        res = minimize(main_task, opt_ini_input, args=w90, method=config.method, tol=config.tol, bounds=config.opt_bounds, constraints=config.opt_constraint, options={'maxiter': config.maxiter})
        # raise ValueError(f'Unsupprot method input: {config.method}')

    logger.info('Final Quality: {0}'.format(res.fun).ljust(80, ' '))
    logger.info("""
                    ──░█▀▀▀─░█▄ ░█ ░█▀▀▄  
                      ░█▀▀▀ ░█░█░█─░█─░█──
                     ─░█▄▄▄─░█─ ▀█ ░█▄▄▀─ 
                                                                                """)