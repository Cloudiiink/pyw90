#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) En Wang (Cloudiiink) <wangenzj@outlook.com>, SF10, IOP/CAS.
# Distributed under the terms of GPLv3.Please see the LICENSE file that should have been included as part of this package.
# For more information, please refer to https://github.com/Cloudiiink/pyw90.
# All rights reserved.

from __future__ import annotations
import argparse
import time, os
from scipy.optimize import minimize
import numpy as np
from numpy.typing import ArrayLike
import argparse

# pyw90
from pyw90.lib.config import Config
from pyw90.lib.job import Job
from pyw90.lib.w90 import W90

import logging

def get_args():
    r'''
    CML Parser
    '''
    parser = argparse.ArgumentParser(description='Auto-Wannier90-Fit. Using minimize method to choose the most suitable energy windows.')
    parser.add_argument('--path', action='store', type=str, default='.',
                        help='Default: .')
    parser.add_argument('--environ', action='store', type=str, default=None,
                        help='Full Environment Variables from `os.environ`. Default: None')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Change config level for logging to `DEBUG`. Defulat is `INFO`.')
    args = parser.parse_args()
    return args

def main_task(opt_input:ArrayLike, w90:W90) -> float:
    w90.config.niter += 1
    logger.info(f"\n{''.center(80, '-')}")
    logger.info(f"Round {w90.config.niter} with input `dis_windows`:".ljust(80, ' '))
    w90.show_total_dis_dict(opt_input)

    if not w90.check_total_dis_dict(opt_input):
        opt_input = w90.rational_opt_input(opt_input)
        w90.show_total_dis_dict(opt_input)

    input_dict = w90.total_dis_dict(opt_input)
    w90.edit_win(input_dict)

    # TODO: run job on local machine (or without queueing system) through `wannier90.x`. Check outputs in `wannier90.wout` file or `wannier90_band.dat` file?
    job = Job(w90.config)
    job.submit()
    logger.info("Successfully submit the task!".ljust(80, ' '))
    time.sleep(5)
    job.check_run_until_stop()
    logger.info('Job <{0}> is done.'.format(w90.config.job_name).ljust(80, ' '))

    w90_bnd = os.path.join(w90.path, w90.config.w90_bnd)
    mtime = time.ctime(os.path.getmtime(w90_bnd))
    logger.info(f'{os.path.abspath(w90_bnd)}'.ljust(80, ' '))
    logger.info(f'    was last modified at {mtime}.'.ljust(80, ' '))
    logger.info('Evaluating the quality ...'.ljust(80, ' '))

    res = w90.evaluate()
    if w90.config.display['spread']: w90.show_spread()
    if w90.config.display['diff']:   w90.show_dEs()
    logger.info('Final Value: {0} meV (with Average of all bands all kpoints)'.format(res).ljust(80, ' '))

    return res

if __name__ == '__main__':
    args = get_args()

    strtime = time.strftime(r'%Y%m%d%H%M%S', time.localtime())
    logfilename = 'log_autow90_{0}.log'.format(strtime)
    level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=level,
                        filename=os.path.join(args.path, logfilename),
                        datefmt='%H:%M:%S',
                        format='%(message)s ! %(asctime)s/%(levelname)s/%(module)s/%(lineno)d')
    logger = logging.getLogger(__name__)

    logger.info("""
                    ──█▀▀█─█  █ ▀▀█▀▀ █▀▀█ ░█  ░█ ▄▀▀▄ █▀▀█
                     ░█▄▄█ █  █   █ ──█──█─░█░█░█ ▀▄▄█ █▄▀█
                     ░█ ░█ ─▀▀▀───▀─  ▀▀▀▀ ░█▄▀▄█──▄▄▀ █▄▄█
                     
    A module inside `pyw90` for automatically optimizing dis energy windows.
    For more information, please refer to https://github.com/Cloudiiink/pyw90
                                                                                """)

    config = Config(yaml_file= os.path.join(args.path, "auto_w90_input.yaml"),
                    environ=eval(args.environ))
    config.info()
    w90 = W90(config, path=args.path)

    logger.info(' Calculation Start! '.center(80, '='))

    opt_ini_input = np.array([v for v in config.opt_ini_dis.values()])

    if not w90.check_total_dis_dict(opt_ini_input):
        opt_ini_input = w90.rational_opt_input(opt_ini_input)
        w90.show_total_dis_dict(opt_ini_input)

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