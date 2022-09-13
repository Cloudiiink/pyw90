import subprocess
import os
import argparse

from lib.w90 import W90
from lib.config import Config
from utility.utility import get_efermi, show_all_fonts, parse_kernel
from pre_w90_tool import main_features

def get_args():
    r'''
    CML Parser
    '''
    parser = argparse.ArgumentParser(description='python Command-line toolbox for VASP and Wannier90 interface with utility')
    subparsers = parser.add_subparsers(help='Main features')

    # auto-wannier90-fit
    parser_auto = subparsers.add_parser('auto', help='(Auto Wannier90 Fit) Using minimize method to choose the most suitable energy windows.')
    parser_auto.add_argument('--path', default='.',
                             help='The path of working dir. Please use relative path. Default: .')
    parser_auto.set_defaults(func=auto)

    # compare VASP & Wannier90 result
    parser_cmp = subparsers.add_parser('cmp', help='(Comparison) Show difference between VASP bands and Wannier90 bands via plotting and report. `bnd.dat` for VASP band data in `p4vasp` format and `wannier90_band.dat`, `wannier90_band.labelinfo.dat`, and `wannier90.wout` are required for plotting and analysis.')
    parser_cmp.add_argument('name',
                            help='name of system')
    parser_cmp.add_argument('--config', action='store_true', default=False,
                            help='Read input from config file `auto_w90_input.yaml` directly. Default: False')
    parser_cmp.add_argument('--path', default='.',
                            help='The path of working dir. Please use relative path. Default: .')
    parser_cmp.add_argument('--efermi', default=None, 
                            help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser_cmp.add_argument('--vasp', default='bnd.dat',
                            help="location of VASP band file in `p4vasp` format. Default: bnd.dat")
    parser_cmp.add_argument('--ylim', default=None, nargs=2, type=float,
                            help="Energy bound for plot. Since the Fermi level has been shift to 0 during the plotting, please mind your input. Default: [E_w90.min - 1, E_w90.max + 1]")
    parser_cmp.add_argument('--kernel', default='unit,0,1',
                            help="kernel function for evaluating diff with formatted input: type, middle, width and **the Fermi level is not subtracted from eigenvalues**. There are two type of kernel function: `unit` and `gaussian`. Defalut: unit,0,1")
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

    # pre-process of VASP data
    parser_pre = subparsers.add_parser('pre', help='(Pre)-analysis before `Wannier90` Interpolation.')
    parser_pre.add_argument('mode', help='Mode: kpath, band, template, dos')
    parser_pre.add_argument('--path', default='.',
                            help='The path of working dir. Please use relative path. Default: .')
    parser_pre.add_argument('--deg', action='store', type=int, default=1,
                            help='Degeneracy of bands. Default: 1')
    parser_pre.add_argument('--sub-fermi', action='store_true', default=False,
                            help="Flag for whether the input `erange` has subtract the Fermi energy or not. Default: False")
    parser_pre.add_argument('--lb', action='store', type=float, default=0.1,
                            help='Lower bound for selected orbital / max single orbital. default: 0.1')
    parser_pre.add_argument('-e', dest='erange', action='store', type=float,
                            default=None, nargs=2,
                            help='Energy range.')
    parser_pre.add_argument('--plot', default=False, action="store_true",
                            help='plot the dos distribution')
    parser_pre.add_argument('--extra', action='store', type=str, default='',
                            help='Extra input. In `template` mode and within extra input (basic, wann, band), we can choose one of the detailed parts to print. In `dos` mode and within extra input (`species`, `structure_id`, `orbital_id` list separated by ;), we can decide the projections for `Wannier90` input to analyze. See example in document for details.')
    parser_pre.set_defaults(func=pre)

    # show distribution of eigenvalues
    parser_eig = subparsers.add_parser('eig', help='Show distribution of eigenvalues.')
    parser_eig.add_argument('mode', help='Mode: report, plot, count, suggest')
    parser_eig.add_argument('--config', action='store_true', default=False,
                            help='Read input from config file `auto_w90_input.yaml` directly. Default: False')
    parser_eig.add_argument('--path', default='.',
                            help='The path of working dir. Please use relative path. Default: .')
    parser_eig.add_argument('-i', dest='eig', action='store', type=str,
                            default='EIGENVAL',
                            help='Select wannier90.eig file or EIGENVAL file. Default: EIGENVAL')
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
    parser_eig.add_argument('-e', dest='erange', action='store', type=float,
                            default=None, nargs=2,
                            help='Energy range.')
    parser_eig.add_argument('--separate', default=False, action="store_true",
                            help='Calculate bands not separately.')
    parser_eig.set_defaults(func=eig)

    args = parser.parse_args()
    return args

def auto(args):
    r'''
    (Auto Wannier90 Fit) Using minimize method to choose the most suitable energy windows.
    '''
    log  = os.path.join(args.path, 'auto_w90_output.txt')
    subprocess.Popen(['python', 'auto_w90_fit.py'],
                     stdin  = subprocess.DEVNULL,
                     stdout = open(log, 'w'),
                     stderr = open(log, 'w'),
                     start_new_session=True)

def cmp(args):
    r'''
    (Comparison) Show difference between VASP bands and Wannier90 bands via plotting and report.
    '''
    if args.show_fonts:
        show_all_fonts()
        return

    if args.config:
        config = Config(yaml_file='auto_w90_input.yaml')
        w90 = W90(config=config, path=args.path)
        kernel = config.kernel
    else:
        efermi = get_efermi(args)
        w90 = W90(path=args.path, efermi=efermi)
        l = args.kernel.split(',')
        kernel_str, mid, width = l[0], float(l[1]), float(l[2])
        kernel = parse_kernel(kernel_str, mid, width)
    
    print(f'Reading Data from {os.path.relpath(args.path)}')

    w90.plot_cmp_vasp_w90(args.name, ylim=args.ylim,
                          font=args.fontfamily, size=args.fontsize)

    if not args.no_quality and not args.quiet:
        res = w90.evaluate(kernel=kernel)
        print(f'Final Quality: {res}')
        w90.show_dEs(update=True, terminal=True)

    if not args.no_spread and not args.quiet:
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
    if args.config:
        config = Config(yaml_file='auto_w90_input.yaml')
        w90 = W90(config=config, eig=args.eig, path=args.path, nbnds_excl=args.nbnds_excl)
    else:
        w90 = W90(eig=args.eig,
                  path=args.path,
                  efermi=get_efermi(args), 
                  nbnds_excl=args.nbnds_excl, 
                  nwann=args.nwann, 
                  ndeg=args.ndeg)

    if args.mode[0].lower() == 'p': # plot
        w90.plot_eigenval(erange=args.erange, separate=args.separate,
                          savefig=os.path.join(os.getcwd(), 'eigenval_dis.png'))
    elif args.mode[0].lower() == 'r': # report
        w90.report_eigenval(erange=args.erange, separate=args.separate)
    elif args.mode[0].lower() == 'c': # count
        # Count how many states inside the energy interval
        print('For `dis_win_min` and `dis_win_max` settings:')
        print(f'    {w90.count_states_least(args.erange)} states at least in {args.erange}.')
        print('For `dis_froz_min` and `dis_froz_max` settings:')
        print(f'    {w90.count_states_most(args.erange)} states at most in {args.erange} (At some kpoints).')
    elif args.mode[0].lower() == 's': # suggest
        # suggest frozen window with given energy interval
        print(f'`dis_froz_min` and `dis_froz_max` Table:')
        df = w90.get_dis_froz_df(args.erange)
        if len(df) > 0:
            print(df)
        # dis_windows require energy window containing states larger than number of target WFs. This will also generate some constraints for dis_windows
        dis_win_max = w90.suggest_win_max(args.erange[0])
        print(f'\nLowest `dis_win_max` for {args.erange[0]}: {dis_win_max}')

    else:
        print(f'Unsupported mode: {args.mode}')

# MAIN

if __name__ == '__main__':
    args = get_args()
    args.path = os.path.join(os.getcwd(), args.path)
    args.func(args)

