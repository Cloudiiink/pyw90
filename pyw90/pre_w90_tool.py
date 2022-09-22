#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) En Wang (Cloudiiink) <wangenzj@outlook.com>, SF10, IOP/CAS.
# Distributed under the terms of GPLv3.Please see the LICENSE file that should have been included as part of this package.
# For more information, please refer to https://github.com/Cloudiiink/pyw90.
# All rights reserved.

import argparse
import os
from os.path import abspath, relpath

# pymatgen
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin

# pyw90
from pyw90.utility.utility import get_efermi
from pyw90.utility.utility import bc
from pyw90.lib.w90 import W90
from pyw90.lib.dos import DOS

def kpath(kpoints: str, delimeter: str=None):
    r"""
    Generate `Wannier90` `Kpoint_Path` bolck from `KPOINTS`

    **Note**: `delimeter` is used to identify kpoint, such as `0.0  0.0  0.0  ! Gamma`.
    """
    delimeter = delimeter if delimeter else '!'
    with open(kpoints, 'r') as f:
        lines = f.readlines()
    lines    = [l.strip() for l in lines if len(l.strip())!=0 ][4:]
    labels   = [s.split(delimeter)[-1].strip() for s in lines]
    nums_str = [s.split(delimeter)[0].strip() for s in lines]
    nums     = [[eval(n) for n in s.split()] for s in nums_str]
    nice_string = lambda p, label: f"{label:5s} {p[0]: 02.5f}    {p[1]: 02.5f}    {p[2]: 02.5f}"
    for i in range(0, len(labels), 2):
        s1 = nice_string(nums[i], labels[i])
        s2 = nice_string(nums[i+1], labels[i+1])
        print(f'{s1}  {s2}')

def export_vasp_band(path: str):
    r"""
    Export `VASP` band with `p4vasp` format in `path` folder.
    """
    def export2dat(kk, EE, filename):
        lines = []
        nbnds, nk = EE.shape
        for i in range(nbnds):
            for j in range(nk):
                lines.append(f'{kk[j]: 12.6f} {EE[i, j]: 12.6f}\n')
            lines.append('\n')
        with open(filename, 'w') as f:
            f.writelines(lines)

    run = Vasprun(os.path.join(path, 'vasprun.xml'))
    bands = run.get_band_structure(os.path.join(path, 'KPOINTS'), line_mode=True)
    nspin = len(bands.bands.keys())
    kk = bands.distance
    if nspin == 1:
        export2dat(kk, bands.bands[Spin.up], os.path.join(path, 'bnd.dat'))
        os.popen(f"ln -s {relpath(os.path.join(path, 'bnd.dat'))} .")
    elif nspin == 2:
        export2dat(kk, bands.bands[Spin.up], os.path.join(path, 'bnd_up.dat'))
        export2dat(kk, bands.bands[Spin.down], os.path.join(path, 'bnd_down.dat'))
        os.popen(f"ln -s {relpath(os.path.join(path, 'bnd_up.dat'))} .")
        os.popen(f"ln -s {relpath(os.path.join(path, 'bnd_down.dat'))} .")

def get_args():
    '''
    CML parser.
    '''
    parser = argparse.ArgumentParser(description="Pre analysis before Wannier90 Interpolation. dos: pyw90 pre dos --plot -e -4 7 --extra 'Bi,4-7,0-3;F,8-23,1-3'", add_help=True)

    parser.add_argument('mode', help='Mode: kpath, band, template, dos. Only the first character is recognized.')
    parser.add_argument('--path', default='.',
                        help='Default: .')
    parser.add_argument('-e', dest='erange', action='store', type=float,
                        default=[-1e3, 1e3], nargs=2,
                        help='Energy range.')
    parser.add_argument('--rm-fermi', action='store_true', default=False,
                        help="Whether or not the input `erange` has removed the Fermi energy is indicated by this flag. Default: False")
    parser.add_argument('--deg', action='store', type=int, default=1,
                        help='Degeneracy of bands. Default: 1')
    parser.add_argument('--lb', action='store', type=float, default=0.1,
                        help='Lower bound for selected orbital / max single orbital. default: 0.1')
    parser.add_argument('--spin-down', action='store_true', default=False,
                        help="Specify the spin channel to `Spin.down`. Without this argument, the default one is `Spin.up`.")
    parser.add_argument('--plot', default=False, action="store_true",
                        help='plot the dos distribution')
    parser.add_argument('--extra', action='store', type=str, default='',
                        help='Extra input.')
    parser.add_argument('--eps', action='store', type=float, default=4e-3,
                        help="Tolerance for dis energy window suggestion. Default: 0.004")
    return parser.parse_args()

def main_features(args):
    r'''
    Main features integrated
    '''
    path = os.path.join(os.getcwd(), args.path)
    if args.mode[0].lower() == 'k': # generate kpath
        if 'KPOINTS' in os.listdir(path):
            bc.cprint(bc.BLUE, f'Generate Kpoint_Path from\n    {abspath(os.path.join(path, "KPOINTS"))}\n')
            kpath(os.path.join(path, "KPOINTS"), delimeter=args.extra)
        else:
            bc.cprint(bc.BLUE, f'There is no KPOINTS file in\n   {path}\nPlease check your input!')
    elif args.mode[0].lower() == 't':   # print W90 parameter template
        if len(args.extra) == 0:
            W90.win_template('basic')
            W90.win_template('wann')
            W90.win_template('band')
        else:
            W90.win_template(args.extra)
    elif args.mode[0].lower() == 'b':   # generate p4vasp-like band data file
        bc.cprint(bc.BLUE, f"Generate `p4vasp` format `bnd.dat` file from {abspath(os.path.join(path, 'vasprun.xml'))}")
        export_vasp_band(path)
        print()
        bc.cprint(bc.BLUE, f'Export band data to `bnd.dat` (or `bnd_up(down).dat` for spinful system).')
        bc.cprint(bc.BLUE, f'Files should be stored at {abspath(path)}')
        bc.cprint(bc.BLUE, f"We also create symbolic link to current folder for `bnd.dat`.")
    elif args.mode[0].lower() == 'd':   # DOS Analysis
        efermi = get_efermi(args, from_file=True)
        spin   = Spin.up if not args.spin_down else Spin.down
        lb     = args.lb
        left, right = min(args.erange), max(args.erange)
    
        vasprun = Vasprun(os.path.join(path, "vasprun.xml"))
        dos_data_total = vasprun.complete_dos       # get dos data

        bc.cprint(bc.BLUE, f"\nReading vasprun.xml file from\n    `{abspath(os.path.join(path, 'vasprun.xml'))}` \n" \
                            "for DOS analysis...")
        print(f"    Fermi level : {efermi:.5f} eV")
        print(f"    DOS Gap     : {dos_data_total.get_gap(spin=spin):.5f} eV\n")

        if args.rm_fermi:
            left, right = left + efermi, right + efermi
            
        erange = left, right
        bc.cprint(bc.BLUE, f"Calculated DOS Energy Range: {left}, {right}\n")
        
        structure = dos_data_total.structure
        dos_df = DOS.get_dos_df(dos_data_total, left, right, spin=spin)
        
        if len(args.extra) == 0:
            print()
            print(dos_df.sort_values('dos', ascending=False))
            print()
            norb, simple_res_df = DOS.get_dos_analysis_df(dos_df, lb=lb)
            nwann = int(norb * args.deg)

            bc.cprint(bc.RED, f"Based on your input, set the lower selection bound to {lb} and {norb} orbitals are selected.")
            bc.cprint(bc.RED, f"Number of WFs selected: {nwann} (with degeneracy {args.deg})\n")
            bc.cprint(bc.RED,  "Orbitals Selected: ")
            print(simple_res_df)
            print()
            bc.cprint(bc.RED,  "Wannier90 Projection:")
            print(DOS.get_projections_w90_str(simple_res_df, structure))
            print()
            bc.cprint(bc.RED, "pyw90 --extra input:")
            selected_str = DOS.get_projections_selected_str(simple_res_df, structure)
            print(selected_str)
            print()

            if args.plot:
                selected = DOS.parse_projections_from_selected(selected_str)
                savefig  = os.path.join(os.getcwd(), 'dos_analysis.pdf')
                DOS.plot_dos_df(dos_df, selected=selected, savefig=savefig)

        else:
            selected = DOS.parse_projections_from_selected(args.extra)
            nwann = int(args.deg * len(selected))

            w90 = W90(eig='EIGENVAL',
                      path=path,
                      efermi=efermi, # get_efermi(args, from_file=True),
                      nbnds_excl=0,
                      nwann=nwann,
                      ndeg=1,
                      eps=args.eps)

            dis_froz_df = w90.get_dis_froz_df(erange)
            dis_tdos_l, dis_pdos_l, percent_l = [], [], []
            for fmin, fmax in zip(dis_froz_df['dis_froz_min'], dis_froz_df['dis_froz_max']):
                dis_pdos = DOS.get_dos_from_selected(dos_data_total, (fmin, fmax), selected)
                dis_tdos = DOS.get_dos_integral(vasprun.tdos.energies, vasprun.tdos.densities[spin], (fmin, fmax))
                percent  = dis_pdos / dis_tdos
                dis_tdos_l.append(dis_tdos)
                dis_pdos_l.append(dis_pdos)
                percent_l.append(percent)

            dis_froz_dos_df = dis_froz_df.assign(pdos=dis_pdos_l, tdos=dis_tdos_l, percent=percent_l)
            dis_froz_dos_df = dis_froz_dos_df.sort_values('percent', ascending=False)
            # N = len(dis_froz_dos_df)
            print()
            bc.cprint(bc.RED,  "Dis frozen window table")
            bc.cprint(bc.BLUE, "The table is sorted according to the percentage of pdos/tdos.")
            bc.cprint(bc.BLUE, "If you want to see `dis_win_min(max)` suggestion, please use `pyw90 eig suggest` menu.")
            print(dis_froz_dos_df)

            dis_win_max = w90.suggest_win_max(erange[0])
            bc.cprint(bc.RED, f'\nLowest `dis_win_max` for {erange[0]}: {dis_win_max}')
            print()

            if args.plot:
                # print(f'Use best dis frozen window in above table to regenerate dos_df to plot ...')
                # left  = dis_froz_dos_df['dis_froz_min'][0]
                # right = dis_froz_dos_df['dis_froz_max'][0]
                dos_df  = DOS.get_dos_df(dos_data_total, left, right)
                savefig = os.path.join(os.getcwd(), 'dos_analysis_selected.pdf')
                DOS.plot_dos_df(dos_df, selected=selected,
                                savefig=savefig)

if __name__ == "__main__":
    args = get_args()
    main_features(args)