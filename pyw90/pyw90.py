import subprocess
import os
import argparse

from lib.w90 import W90

def get_args():
    r'''
    CML Parser
    '''
    parser = argparse.ArgumentParser(description='python Command-line toolbox for VASP and Wannier90 interface with utility: Using minimize method to choose the most suitable energy windows.')
    subparsers = parser.add_subparsers(help='Main features')

    # auto-wannier90-fit
    parser_auto = subparsers.add_parser('auto', help='Auto Wannier90 Fit')
    parser_auto.add_argument('--path', default='.',
                             help='Default: .')
    parser_auto.set_defaults(func=auto)

    # compare VASP & Wannier90 result
    parser_cmp = subparsers.add_parser('cmp', help='Comparison between VASP band and Wannier90 band. `bnd.dat` for VASP band data in p4vasp format and `wannier90_band.dat`, `wannier90_band.labelinfo.dat`, and `wannier90.wout` are required for plotting and analysis.')
    parser_cmp.add_argument('name',
                            help='name of system')
    parser_cmp.add_argument('--efermi', default=None, 
                            help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser_cmp.add_argument('--path', default='.',
                            help='Default: .')
    parser_cmp.add_argument('--vasp', default='bnd.dat',
                            help="location of VASP band file in p4vasp format. Default: bnd.dat")
    parser_cmp.add_argument('--ylim', default=None, nargs=2, type=float,
                            help="Energy bound for plot. Default: [E_w90.min - 1, E_w90.max + 1]")
    parser_cmp.add_argument('--kernel', default='unit,2,5',
                            help="kernel function for evaluating diff: type, middle, width. There are two type of kernel function: `unit` and `gaussian`. Defalut: unit,2,5")
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
    parser_pre = subparsers.add_parser('pre', help='Pre analysis before Wannier90 Interpolation.')
    parser_pre.add_argument('mode', help='Mode: kpath, band, template, dos')
    parser_pre.add_argument('--path', default='.',
                            help='Default: .')
    parser_pre.add_argument('--no-soc', action='store_true', default=False,
                            help='with SOC or not. Default: False')
    parser_pre.add_argument('--pick', action='store', type=float, default=0.1,
                            help='Minimum of selected orbital / max single orbital. default: 0.1')
    parser_pre.add_argument('-e', dest='erange', action='store', type=float,
                            default=None, nargs=2,
                            help='Energy range.')
    parser_pre.add_argument('--plot', default=False, action="store_true",
                            help='plot the dos distribution')
    parser_pre.add_argument('--extra', action='store', type=str, default='',
                            help='Extra input.')
    parser_pre.set_defaults(func=pre)

    # show distribution of eigenvalues
    parser_dis = subparsers.add_parser('dis', help='Command-line Tool for W90 energy windows.')
    parser_dis.add_argument('mode', help='Mode: report, plot, count, suggest')
    parser_dis.add_argument('-i', dest='eig', action='store', type=str,
                            default='EIGENVAL',
                            help='Select wannier90.eig file or EIGENVAL file. Default: EIGENVAL')
    parser_dis.add_argument('--path', default='.',
                            help='Default: .')
    parser_dis.add_argument('--efermi', dest='efermi', action='store',
                            default=None,
                            help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser_dis.add_argument('-w', dest='nwann', action='store', type=int,
                            default=0,
                            help='Number of Wannier Functions. Default: 0')
    parser_dis.add_argument('-n', dest='nbnds_excl', action='store', type=int,
                            default=0,
                            help='Number of bands excluded')
    parser_dis.add_argument('-d', dest='ndeg', action='store', type=int,
                            default=2,
                            help='Number of degeneracy')
    parser_dis.add_argument('-e', dest='erange', action='store', type=float,
                            default=None, nargs=2,
                            help='Energy range.')
    parser_dis.add_argument('--separate', default=False, action="store_true",
                            help='Calculate bands not separately.')
    parser_dis.set_defaults(func=dis)

    args = parser.parse_args()
    return args

def auto(args):
    path = os.getcwd() + '/' + args.path
    log  = f'{path}/auto_w90_output.txt'
    print(log)
    # comment for test
    # subprocess.Popen(['python', 'auto_w90_fit.py'],
    #                 stdin  = subprocess.DEVNULL,
    #                 stdout = open(log, 'w'),
    #                 stderr = open(log, 'w'),
    #                 start_new_session=True)

def cmp(args):
    name = args.name
    print(name)
    print(args.path)
    # output_figure = f'{args.path}/{name}_VASP_W90_cmp.png'
    # efermi = get_efermi(args)
    
    # print(f'Reading Data from {args.path}/{args.vasp}')
    # vkk, vee = parse_dat(f'{args.path}/{args.vasp}')
    # wkk, wee = parse_dat(f'{args.path}/wannier90_band.dat')

    # plot_cmp_vasp_w90(vkk, vee, wkk, wee, 
    #                   ylim=args.ylim,
    #                   efermi=efermi,
    #                   font=args.fontfamily, size=args.fontsize)

    # if not args.no_quality and not args.quiet:
    #     print('Evaluating Band Quality:')
    #     l = args.kernel.split(',')
    #     kernel, mid, width = l[0], float(l[1]), float(l[2])
    #     dEs, wgts = evaluate_cmp_vasp_w90(vkk, vee, wkk, wee,
    #                                       kernel=kernel, mid=mid, width=width)
    #     show_vasp_w90_diff(dEs, wgts)

    # if not args.no_spread and not args.quiet:
    #     logger.info('Show spreading convergence:')
    #     show_spreading(args.path)

    # if args.show_fonts:
    #     show_all_fonts()

def pre(args):
    pass

def dis(args):
    w90 = W90(eig=args.eig,
              path=args.path,
              efermi=0, #get_efermi(args), 
              nbnds_excl=args.nbnds_excl, 
              nwann=args.nwann, 
              ndeg=args.ndeg)

    if args.mode[0].lower() == 'p': # plot
        w90.plot_eigenval(erange=args.erange, separate=args.separate)
    elif args.mode[0].lower() == 'r': # report
        w90.report_eigenval(erange=args.erange, separate=args.separate)
    elif args.mode[0].lower() == 'c': # count
        # Count how many states inside the energy interval
        print(f'There are {w90.count_states(args.erange)} states in {args.erange}.')
    elif args.mode[0].lower() == 's': # suggest
        # suggest frozen window with given energy interval
        print(f'`dis_froz_min` and `dis_froz_max` Table:')
        df = w90.get_dis_froz_df(args.erange, eps=4e-3)
        if len(df) > 0:
            print(df)
        # dis_windows require energy window containing states larger than number of target WFs. This will also generate some constraint for dis_windows
        dis_win_max = w90.suggest_win_max(args.erange[0])
        print(f'\nLowest `dis_win_max` for {args.erange[0]}: {dis_win_max}')

    else:
        print(f'Unsupported mode: {args.mode}')

if __name__ == '__main__':
    args = get_args()
    args.func(args)

