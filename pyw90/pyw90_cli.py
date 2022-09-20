#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) En Wang (Cloudiiink) <wangenzj@outlook.com>, SF10, IOP/CAS.
# Distributed under the terms of GPLv3.Please see the LICENSE file that should have been included as part of this package.
# For more information, please refer to https://github.com/Cloudiiink/pyw90.
# All rights reserved.

import subprocess
import signal
import os, sys
import argparse
import shutil

# pyw90
from pyw90.lib.w90 import W90
from pyw90.lib.job import Job
from pyw90.lib.config import Config
from pyw90.utility.utility import get_efermi, show_all_fonts, parse_kernel
from pyw90.utility.utility import bc
from pyw90.pre_w90_tool import main_features

def get_args():
    r'''
    CML Parser
    '''
    parser = argparse.ArgumentParser(description='python Command-line toolbox for VASP and Wannier90 interface with utility. You can also use -h to display the submenu help message. e.g. pyw90 pre -h')
    subparsers  = parser.add_subparsers(help='Main features')
    parser_eig  = subparsers.add_parser('eig' , help='Show distribution of eigenvalues.')
    parser_pre  = subparsers.add_parser('pre' , help='(Pre)-analysis before `Wannier90` Interpolation.')
    parser_auto = subparsers.add_parser('auto', help='(Auto Wannier90 Fit) Using minimize method to choose the most suitable dis energy windows.')
    parser_cmp  = subparsers.add_parser('cmp' , help='(Comparison) Show difference between VASP bands and Wannier90 bands via plotting and report. `bnd.dat` for VASP band data in `p4vasp` format and `wannier90_band.dat`, `wannier90_band.labelinfo.dat`, and `wannier90.wout` are required for plotting and analysis.')
    
    # show distribution of eigenvalues
    parser_eig.add_argument('mode', help='Mode: report, plot, count, suggest')
    parser_eig.add_argument('-e', dest='erange', action='store', type=float,
                            default=None, nargs=2,
                            help='Energy range.')
    parser_eig.add_argument('--config', action='store_true', default=False,
                            help='Read input from config file `auto_w90_input.yaml` directly. Default: False')
    parser_eig.add_argument('--path', default='.',
                            help='The path of working dir. Default: .')
    parser_eig.add_argument('-i', dest='eig', action='store', type=str,
                            default='EIGENVAL',
                            help='Select wannier90.eig file or EIGENVAL file. Default: EIGENVAL')
    parser_eig.add_argument('--rm-fermi', action='store_true', default=False,
                            help="Whether or not the input `erange` has removed the Fermi energy is indicated by this flag. Default: False")
    parser_eig.add_argument('--efermi', dest='efermi', action='store',
                            default=None,
                            help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser_eig.add_argument('-w', dest='nwann', action='store', type=int,
                            default=0,
                            help='Number of Wannier Functions. Default: 0')
    parser_eig.add_argument('-n', dest='nbnds_excl', action='store', type=int,
                            default=0,
                            help='Number of bands excluded below the bands from `Wannier90`.')
    parser_eig.add_argument('--deg', dest='ndeg', action='store', type=int,
                            default=1,
                            help='Number of degeneracy. Default: 1')
    parser_eig.add_argument('--separate', default=False, action="store_true",
                            help='Calculate bands not separately.')
    parser_eig.add_argument('--eps', action='store', type=float, default=4e-3,
                            help="Tolerance for dis energy window suggestion. Default: 0.004")
    parser_eig.set_defaults(func=eig)

    # pre-process of VASP data
    parser_pre.add_argument('mode', help='Mode: kpath, band, template, dos')
    parser_pre.add_argument('--path', default='.',
                            help='The path of working dir. Default: .')
    parser_pre.add_argument('-e', dest='erange', action='store', type=float,
                            default=[-1e3, 1e3], nargs=2,
                            help='Energy range.')
    parser_pre.add_argument('--lb', action='store', type=float, default=0.1,
                            help='Lower bound for selected orbital / max single orbital. default: 0.1')
    parser_pre.add_argument('--rm-fermi', action='store_true', default=False,
                            help="Whether or not the input `erange` has removed the Fermi energy is indicated by this flag. Default: False")
    parser_pre.add_argument('--extra', action='store', type=str, default='',
                            help='Extra input. In `template` mode and within extra input (basic, wann, band), we can choose one of the detailed parts to print.' \
                                 'In `dos` mode and within extra input (`species`, `structure_id`, `orbital_id` list separated by ;), ' \
                                 'we can treat input as projections for `Wannier90` input to suggest dis frozen energy. Details can be found in the document.')
    parser_pre.add_argument('--spin-down', action='store_true', default=False,
                            help="Specify the spin channel to `Spin.down`. Without this argument, the default one is `Spin.up`.")
    parser_pre.add_argument('--plot', default=False, action="store_true",
                            help='plot the dos distribution')
    parser_pre.add_argument('--eps', action='store', type=float, default=4e-3,
                            help="Tolerance for dis energy window suggestion. Default: 0.004")
    parser_pre.add_argument('--deg', action='store', type=int, default=1,
                            help='Degeneracy of bands. Default: 1')
    parser_pre.set_defaults(func=pre)

    # auto-wannier90-fit
    parser_auto.add_argument('mode', 
                             help='Mode: run, term(inate), input. Only first character is recognized.')
    parser_auto.add_argument('--path', default='.',
                             help='The path of working dir. Default: .')
    parser_auto.add_argument('--pid', default=None, type=int,
                             help='PID to terminate.')
    parser_auto.set_defaults(func=auto)

    # compare VASP & Wannier90 result
    parser_cmp.add_argument('name',
                            help='name of system')
    parser_cmp.add_argument('--config', action='store_true', default=False,
                            help='Read input from config file `auto_w90_input.yaml` directly. Default: False')
    parser_cmp.add_argument('--path', default='.',
                            help='The path of working dir. Default: .')
    parser_cmp.add_argument('--efermi', default=None, 
                            help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser_cmp.add_argument('--vasp', dest='vasp', default='bnd.dat',
                            help="Path of VASP band file in `p4vasp` format. Default: bnd.dat")
    parser_cmp.add_argument('--seedname', dest='seedname', default='wannier90',
                            help="Seedname of Wannier90 input. Default: wannier90")
    parser_cmp.add_argument('--ylim', default=None, nargs=2, type=float,
                            help="Energy bound for plot. Since the Fermi level has been shift to 0 during the plotting, please mind your input. Default: [E_w90.min - 1, E_w90.max + 1]")
    parser_cmp.add_argument('--kernel', default='unit,0,1',
                            help="Formatted input: type, middle, width (Defalut: unit,0,1). Kernel functions are classified into two types: `unit` and `gaussian`. **The Fermi level is not subtracted from eigenvalues when using the kernel function to evaluate difference**.")
    parser_cmp.add_argument('--show-fonts', default=False, action="store_true",
                            help="Show all availabel font families can be used in `rcParams`")
    parser_cmp.add_argument('--fontfamily', default='Open Sans',
                            help="Set font family manually. Default: Open Sans")
    parser_cmp.add_argument('--fontsize', default=18, type=int,
                            help="Set font size manually. Default: 18")
    parser_cmp.add_argument("--no-spread", default=False, action="store_true",
                            help="Don't plot spreading")
    parser_cmp.add_argument("--no-quality", default=False, action="store_true",
                            help="Don't show quality of fitting")
    parser_cmp.add_argument("--quiet", default=False, action="store_true",
                            help="Equal to --no-spreading --no-quality")
    parser_cmp.set_defaults(func=cmp)

    args = parser.parse_args()
    return args

def auto(args):
    r'''
    (Auto Wannier90 Fit) Using minimize method to choose the most suitable energy windows.
    '''
    path    = os.path.dirname(os.path.realpath(__file__))
    environ = dict(os.environ)

    if args.mode.lower()[0] == 'r':  # run
        log  = os.path.join(args.path, 'auto_w90_output.txt')
        p = subprocess.Popen([sys.executable, os.path.join(path, 'auto_w90_fit.py'), '--path', args.path, "--environ", str(environ)],
                               env    = environ,
                               stdin  = subprocess.DEVNULL,
                               stdout = open(log, 'w'),
                               stderr = open(log, 'w'),
                               start_new_session=True)
        bc.cprint(bc.BLUE, f'Auto-W90-Fit run with PID: {p.pid}')

    elif args.mode.lower()[0] == 'i':   # input
        input_path  = os.path.join(path, 'auto_w90_input.yaml')
        target_path = os.path.abspath(args.path)
        target_file = os.path.join(target_path, 'auto_w90_input.yaml')
        bc.cprint(bc.BLUE, f'Create example input file for `auto` menu at input folder')
        bc.cprint(bc.BLUE, f'    {target_path}\n')
        if os.path.exists(target_file):
            bc.cprint(bc.RED, 'There are already input file for pyw90. Do you want to overwrite it?')
            bc.cprint(bc.RED, 'Input Y(es) / N(o)')
            usr_input = input()
            if usr_input.lower()[0] == 'y':
                shutil.copy(input_path, args.path)
                return
            else:
                return

        shutil.copy(input_path, args.path)

    elif args.mode.lower()[0] == 't':   # terminate
        bc.cprint(bc.BLUE, f'Kill the job with PID {args.pid} as listed')

        from getpass import getuser
        local_usr_name = getuser()

        ps = os.popen(f'ps aux | grep {local_usr_name}').read()
        print(ps)
        all_pid = [int(s.split()[1]) for s in ps.strip().split('\n')]
        bc.cprint(bc.RED, 'Input Y(es) / N(o) (Or input the pid listed above you want to terminate): ')
        bc.cprint(bc.RED, 'IMPORTANT: If you run task local, input the `auto_w90_fit.py` and `wannier90.x` separate with comma (e.g. 2101, 2121)')
        bc.cprint(bc.RED, '           Otherwise `auto_w90_fit.py` will repeatedly submit new tasks.')
        usr_input = input()
        if usr_input.lower()[0] == 'y':
            os.killpg(os.getpgid(args.pid), signal.SIGTERM)  # Send the signal to all the process groups
            config = Config(yaml_file='auto_w90_input.yaml')
            jobs = Job(config)
            jobs.cancel()
            return
        elif usr_input[0].isdigit():
            usr_input_pid = [int(i.strip()) for i in usr_input.split(',')]
            for upid in usr_input_pid:
                if upid in all_pid:
                    os.killpg(os.getpgid(upid), signal.SIGTERM)  # Send the signal to all the process groups
            config = Config(yaml_file=os.path.join(args.path, 'auto_w90_input.yaml'))
            jobs = Job(config)
            jobs.cancel()
        else:
            return

def cmp(args):
    r'''
    (Comparison) Show difference between VASP bands and Wannier90 bands via plotting and report.
    '''
    if args.show_fonts:
        show_all_fonts()
        return

    if args.config:
        bc.cprint(bc.BLUE, f'Reading Data and config file from {os.path.relpath(args.path)}')
        config = Config(yaml_file=os.path.join(args.path, 'auto_w90_input.yaml'))
        w90 = W90(config=config, path=args.path)
        kernel = config.kernel
    else:
        bc.cprint(bc.BLUE, f'Reading Data from {os.path.relpath(args.path)}')
        efermi = get_efermi(args)
        w90 = W90(path=args.path, efermi=efermi, vasp_bnd=args.vasp, seedname=args.seedname)
        l = args.kernel.split(',')
        kernel_str, mid, width = l[0], float(l[1]), float(l[2])
        kernel = parse_kernel(kernel_str, mid, width)
    
    bc.cprint(bc.BLUE, f"The output figure should be stored at")
    bc.cprint(bc.BLUE, f"    {os.path.join(args.path, args.name+'_VASP_W90_cmp.png')}")
    w90.plot_cmp_vasp_w90(args.name, ylim=args.ylim,
                          font=args.fontfamily, size=args.fontsize)

    if not args.no_quality and not args.quiet:
        res = w90.evaluate(kernel=kernel)
        bc.cprint(bc.BLUE, f'Final Quality {res:.6f} meV with dEs for each bands')
        w90.show_dEs(update=True, terminal=True)

    if not args.no_spread and not args.quiet:
        print()
        bc.cprint(bc.BLUE, f'Spreading for each iteration is displayed below')
        w90.show_spread(update=True, terminal=True)

def pre(args):
    r'''
    (Pre)-analysis before `Wannier90` Interpolation.
    '''
    main_features(args)

def eig(args):
    r'''
    Show distribution of eigenvalues.
    '''
    efermi = get_efermi(args)
    if args.config:
        config = Config(yaml_file='auto_w90_input.yaml')
        w90 = W90(config=config, 
                  eig=args.eig, 
                  path=args.path, 
                  nbnds_excl=args.nbnds_excl,
                  eps=args.eps)
    else:
        w90 = W90(eig=args.eig,
                  path=args.path,
                  efermi=efermi, 
                  nbnds_excl=args.nbnds_excl, 
                  nwann=args.nwann, 
                  ndeg=args.ndeg,
                  eps=args.eps)
    if args.rm_fermi:
        erange = min(args.erange) + efermi, max(args.erange) + efermi
    else:
        erange = min(args.erange), max(args.erange)
    print(f"Calculated Energy Range: {erange[0]}, {erange[1]} with Fermi level {efermi:.6f}")

    if args.mode[0].lower() == 'p': # plot
        w90.plot_eigenval(erange=erange, separate=args.separate,
                          savefig=os.path.join(os.getcwd(), 'eigenval_dis.png'))
    elif args.mode[0].lower() == 'r': # report
        w90.report_eigenval(erange=erange, separate=args.separate)
    elif args.mode[0].lower() == 'c': # count
        # Count how many states inside the energy interval
        bc.cprint(bc.BLUE, 'For `dis_win_min` and `dis_win_max` settings:')
        print(f'    {w90.count_states_least(erange)} states at least in {erange}.')
        bc.cprint(bc.BLUE, 'For `dis_froz_min` and `dis_froz_max` settings:')
        print(f'    {w90.count_states_most(erange)} states at most in {erange} (At some kpoints).')
    elif args.mode[0].lower() == 's': # suggest
        # suggest dis_win_min and dis_win_max
        Nmost  = w90.count_states_most(erange)
        Nleast = w90.count_states_least(erange)
        print(f'There are at most {Nmost} states and at least {Nleast} states in {args.erange}.\n')

        bc.cprint(bc.RED, '`dis_win_min` and `dis_win_max` Table:')
        bc.cprint(bc.BLUE, '    Column `dis_win_max` shows the **lowest** dis_win_max for `dis_win_min`\n' \
                           '    Column `i+1_min` / `i_max` shows the band minimum / maximum near the gap\n' \
                           '    Column `Nleast`  / `Nmost` shows the least / most number of states inside `dis_win_min` and Fermi level.\n')
        df = w90.suggest_win_table()
        print(df)

        # suggest frozen window with given energy interval
        bc.cprint(bc.RED, '\n`dis_froz_min` and `dis_froz_max` Table:')
        bc.cprint(bc.BLUE, f'    Suggest `dis_froz_min` & `dis_froz_max` as following with #WFs = {w90.nwann} from input')
        bc.cprint(bc.BLUE,  '    Use `pre dos` menu to get suggestions based on projected density of states\n')
        df = w90.get_dis_froz_df(args.erange)
        print(df)

    else:
        bc.cprint(bc.BLUE, f'Unsupported mode: {args.mode}')

# MAIN
def main_cli():
    bc.cprint(bc.RED, '                                   \n' \
                      '         █▀▀█─█  █─░█──░█ ▄▀▀▄ █▀▀█\n' \
                      '        ─█──█ █▄▄█ ░█░█░█ ▀▄▄█ █▄▀█\n' \
                      '         █▀▀▀ ▄▄▄█ ░█▄▀▄█ ─▄▄▀ █▄▄█\n' \
                      '                                   ')
    bc.cprint(bc.RED, 'A VASP and Wannier90 interfaced tool for projection analysis\n' \
                      'and fully automated dis energy window optimization. \n' \
                      '                                   \n' \
                      'For more information, please refer to https://github.com/Cloudiiink/pyw90 \n')
    args = get_args()
    args.path = os.path.join(os.getcwd(), args.path)
    args.func(args)

if __name__ == '__main__':
    main_cli()


