#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) En Wang (Cloudiiink) <wangenzj@outlook.com>, SF10, IOP/CAS.
# Distributed under the terms of GPLv3.Please see the LICENSE file that should have been included as part of this package.
# For more information, please refer to https://github.com/Cloudiiink/pyw90.
# All rights reserved.

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
import argparse
import os
from os.path import abspath, relpath
import warnings
from typing import Tuple

# pymatgen
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.core import structure

# Calculate the percentage of given sites and given orbitals inside the energy interval
# from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev
from scipy.integrate import fixed_quad

from collections import Counter
from functools import reduce

import matplotlib
matplotlib.use('Agg')

# pyw90
from pyw90.utility.utility import get_efermi
from pyw90.utility.utility import num2str, str2num, get_id_from_specie
from pyw90.lib.w90 import W90
from pyw90.lib.dos import DOS

def dos_distribute(e: ArrayLike, dos: ArrayLike,
                   window: Tuple[float, float]) -> float:
    r"""
    Return the integral of `dos` value inside the `window`.

    The left value of `window` is required to be greater than the minimium of `e` and the right value of `window` is required to be less than the maximium of `e`.
    """
    if min(e) > window[0] or max(e) < window[-1]:
        warnings.warn(f"Since some region of the {window} is outside the input variable `e` ranging from {min(e)} to {max(e)}, the intergal from interpolation might not be very convincing.")
    left  = max(min(e), window[0])
    right = min(max(e), window[-1])
    dos_cubic_fitting = splrep(e, dos, s=0)
    dos_func = lambda x: splev(x, dos_cubic_fitting, der=0)
    res = fixed_quad(dos_func, left, right)
    return res[0]

def gen_dos_df(dos_data_total: CompleteDos, left: float, right: float, spin: Spin=Spin.up,
               orb_names: list[str]=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']) -> pd.DataFrame:
    r"""
    Return `pandas.DataFrame` with columns: `species`, `structure_id`, `orb_id`, `orb_name`, `key_string` and `dos`.

    The value of `dos` has been normalized with the maximum is adjusted to 1. 
    """
    structure = dos_data_total.structure
    num_sites = structure.num_sites
    ee = dos_data_total.energies

    # TODO: use `PROCAR` to generate a result, although it might not be accurate enough.
    if left < min(ee) or right > max(ee):
        warnings.warn(f'CHECK YOUR INPUT! The energies in `EIGENVAL` is ranged from {ee.min()} to {ee.max()}' \
                       ' which not include all the energy ranged from {left} to {right} you input.')

    lspec, lstruc, lorbid, lorbname, lstr, ldos = [list() for _ in range(6)]
    for i in range(num_sites):
        for j in range(len(orb_names)):
            site, orb = structure[i], Orbital(j)
            dos = dos_data_total.get_site_orbital_dos(site, orb)
            dis = dos_distribute(ee, dos.densities[spin], [left, right])
            key = site.species_string + '_' + str(i) + '_' + orb_names[j]       # how to formulate the key string for every entries
            lspec.append(site.species_string)
            lstruc.append(i)
            lorbid.append(j)
            lorbname.append(orb_names[j])
            lstr.append(key)
            ldos.append(dis)
    dos_df = pd.DataFrame(zip(lspec, lstruc, lorbid, lorbname, lstr, ldos),
                          columns=["species", "structure_id", "orb_id", "orb_name", "key_string", "dos"])
    dos_df['dos'] = dos_df['dos'] / max(dos_df['dos'])      # Renormalized
    return dos_df

def plot_dos_dis(dos_df: pd.DataFrame, lb: float=1, selected: set=set(),
                 colors: Tuple[str, str]=['brown', 'orange'],
                 savefig: str='dos_analysis.png'):
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
    print(f'Plot with selected orbitals in `{colors[0]}` and not selected orbitals in `{colors[-1]}`.')
    print(f'{len(selected)} orbitals selected')
    mask = (np.array(color) == colors[0])
    print(f'Plot {len(df)} orbitals with {sum(mask)} selected.')
    
    ax = df.plot.barh(x="key_string", y="dos", color=color) #, figsize=(8, 20))
    ax.set_axisbelow(True)
    plt.grid(axis='x')

    # save figure in local folder not in source folder
    plt.savefig(savefig, dpi=200, bbox_inches='tight')
    # plt.savefig(savefig, dpi=200, bbox_inches='tight', transparent=True)

def dos_analysis_df(dos_df: pd.DataFrame, lb: float=0.1) -> Tuple[int, pd.DataFrame]:
    r"""
    Analyze the DOS distribution and recommend the projection according to `lb` (threshold for selecting)

    ### Parameters
    - `dos_df` : `pd.DataFrame` which containg the DOS contribution from each sites and each orbitals. Obtained from `gen_dos_df`.
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

    # Simplify the results using spd electron configuration
    d_orb = ['dxy', 'dyz', 'dxz', 'dx2', 'dz2']
    p_orb = ['px', 'py', 'pz']
    s_orb = ['s']

    # Iterate to check all the `species` in the material
    for species in df.species.unique():
        species_sites = df[df["species"]==species].structure_id.unique()

        # Simplify the orbitals using `spd` electron configuration
        # e.g. convert [s, px, py, pz, dx2] to [s, p, dx2]
        for idx in species_sites:
            site_df = df[df["structure_id"]==idx]
            site_orbs = site_df.orb_name.unique()
            for orb in ["s", "p", "d"]:
                if orb == "s":
                    orb_list = s_orb
                elif orb == "p":
                    orb_list = p_orb
                elif orb == "d":
                    orb_list = d_orb
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
                new = pd.DataFrame({"species": species,
                                    "site": -1,
                                    "orb": [list(ref_orb)]})
                res_df = pd.concat([res_df, new], ignore_index=True)

    # display(res_df)
    return n_orb, res_df

def orb_string(orb_list: list[str]) -> str:
    r"""
    Convert `orb_list` (e.g., ['s', 'p', 'dz2']) to format used in `Wannier90`. (e.g. l=0;l=1)
    """
    orb_string_list = []
    for o in orb_list:
        if o == "s":
            orb_string_list.append('l=0')
        elif o == "p":
            orb_string_list.append('l=1')
        elif o == "d":
            orb_string_list.append('l=2')
        else:
            orb_string_list.append(o)
    return ';'.join([s for s in orb_string_list])

def w90_string(res_df: pd.DataFrame, struc: structure) -> str:
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
    w90_string_list = []
    for _, row in res_df.iterrows():
        if row['site'] == -1:
            w90_string_list.append(f"{row['species']}:{orb_string(row['orb'])}")
        else:
            site_string = ','.join([f"{i:.6f}" for i in struc[row['site']].frac_coords])
            w90_string_list.append(f"f={site_string}:{orb_string(row['orb'])}") # !!! f=... ? or c=...?
    return '\n'.join([s for s in w90_string_list])

def extra_string(res_df: pd.DataFrame, struc: structure,
                 connect: str='-', sep :str='|') -> str:
    r"""
    Convert `res_df` with simplified projection information to formatted string for `--extra`` argument to input for pyw90.

    e.g. 'Ga,0,1-3;As,1,1-3|5'
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

def dos_given_site_orb(dos_data_total: pd.DataFrame, erange: Tuple[float, float], 
                       key_string: str, spin: Spin=Spin.up,
                       orb_names: list[str]=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']) -> float:
    r"""
    Return DOS contribution of `erange` for site and orbital from `key_string` (e.g. `C_1_pz`) via interpolation.
    """
    structure = dos_data_total.structure
    _, site_id, orb_id = key_string.split('_')
    site = structure[int(site_id)]
    orb = Orbital(orb_names.index(orb_id))
    dos = dos_data_total.get_site_orbital_dos(site, orb)
    dis = dos_distribute(dos.energies, dos.densities[spin], erange)
    return dis

def dos_given_selected(dos_data_total: pd.DataFrame, erange: Tuple[float, float], selected: set[str]) -> float:
    r"""
    Return DOS contribution from selected sites and orbitals inside energy interval interpolation.
    """
    dis = [dos_given_site_orb(dos_data_total, erange, key) for key in selected]
    return np.sum(dis)

def select_str2list(s: str, orb_names: list[str]=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']) -> set[str]:
    r"""
    Convert input string to `set` with `key_string` format.
    
    The default deliminater is `;` and we can use `1-4|6` to represent `1,2,3,4,6` both in site_id and orbital_id.
    """

    if len(s) == 0:
        return []

    res = []
    for l in s.split(';'):
        spec, pos, orb = l.split(',')
        pos, orb = str2num(pos), str2num(orb) 
        for p in pos:
            for o in orb:
                res.append(f'{spec}_{p}_{orb_names[o]}')
    return set(res)

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

    print(f"Generate `p4vasp` format bnd.dat file from {abspath(os.path.join(path, 'vasprun.xml'))}")
    run = Vasprun(os.path.join(path, 'vasprun.xml'))
    bands = run.get_band_structure(os.path.join(path, 'KPOINTS'), line_mode=True)
    nspin = len(bands.bands.keys())
    kk = bands.distance
    if nspin == 1:
        print(f'export band data to `bnd.dat` in folder\n    {abspath(path)}')
        export2dat(kk, bands.bands[Spin.up], os.path.join(path, 'bnd.dat'))
        os.popen(f"ln -s {relpath(os.path.join(path, 'bnd.dat'))} .")
    elif nspin == 2:
        print('NSPIN = 2')
        print(f'export band data to `bnd_up.dat` and `bnd_down.dat` separatly in folder\n    {abspath(path)}')
        export2dat(kk, bands.bands[Spin.up], os.path.join(path, 'bnd_up.dat'))
        export2dat(kk, bands.bands[Spin.down], os.path.join(path, 'bnd_down.dat'))
        os.popen(f"ln -s {relpath(os.path.join(path, 'bnd_up.dat'))} .")
        os.popen(f"ln -s {relpath(os.path.join(path, 'bnd_down.dat'))} .")
 
def template(flag:str):
    r"""
    Print template for `Wannier90`.

    `flag`: `wann`, `band` and `basic`
    """
    if flag == 'wann':
        print('###############\n' \
              '#     W90     #\n' \
              '###############\n' \
              'dis_win_max       =  0.0\n' \
              'dis_froz_max      =  0.0\n' \
              'dis_froz_min      = -5.0\n' \
              'dis_win_min       = -5.0\n' \
              '\n'                         \
              'num_iter          = 1000\n' \
              'num_print_cycles  =   40\n' \
              'dis_num_iter      = 5000\n' \
              'dis_mix_ratio     =  1.0\n')
    elif flag == 'band':
        print('###############\n' \
              '#  Band Plot  #\n' \
              '###############   # restart = plot\n' \
              'write_hr          = true\n' \
              'bands_plot        = true\n' \
              'bands_num_points  = 151 \n' \
              'bands_plot_format = gnuplot\n' \
              '\n' \
              'Begin Kpoint_Path\n' \
              'End Kpoint_Path\n')
    elif flag == 'basic':
        print('num_wann  = 8\n' \
              '# num_bands = 15\n' \
              '\n' \
              'begin projections\n' \
              'Ga:l=0;l=1\n' \
              'As:l=0;l=1\n' \
              'end projections\n' \
              '\n' \
              '# spinors = .true.\n')
    else:
        print('Unknown flag: ', flag)

def get_args():
    '''
    CML parser.
    '''
    parser = argparse.ArgumentParser(description="Pre analysis before Wannier90 Interpolation. dos: pyw90 pre dos --plot -e -4 7 --extra 'Bi,4-7,0-3;F,8-23,1-3'", add_help=True)

    parser.add_argument('mode', help='Mode: kpath, band, template, dos')
    parser.add_argument('--path', default='.',
                        help='Default: .')
    parser.add_argument('--sub-fermi', action='store_true', default=False,
                        help="Flag for whether the input `erange` has subtract the Fermi energy or not. Default: False")
    parser.add_argument('--deg', action='store', type=int, default=1,
                        help='Degeneracy of bands. Default: 1')
    parser.add_argument('--lb', action='store', type=float, default=0.1,
                        help='Lower bound for selected orbital / max single orbital. default: 0.1')
    parser.add_argument('--spin-down', action='store_true', default=False,
                        help="Specify the spin channel to `Spin.down`. Without this argument, the default one is `Spin.up`.")
    parser.add_argument('-e', dest='erange', action='store', type=float,
                        default=[-1e3, 1e3], nargs=2,
                        help='Energy range.')
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
            print(f'Generate Kpoint_Path from\n    {abspath(os.path.join(path, "KPOINTS"))}\n')
            kpath(os.path.join(path, "KPOINTS"), delimeter=args.extra)
        else:
            print(f'There is no KPOINTS file in\n   {path}\nPlease check your input!')
    elif args.mode[0].lower() == 't':   # print W90 parameter template
        if len(args.extra) == 0:
            template('basic')
            template('wann')
            template('band')
        else:
            template(args.extra)
    elif args.mode[0].lower() == 'b':   # generate p4vasp-like band data file
        export_vasp_band(path)
    elif args.mode[0].lower() == 'd':   # DOS Analysis
        efermi = get_efermi(args, from_file=True)
        spin   = Spin.down if args.spin_down else Spin.up
        lb     = args.lb
        left, right = min(args.erange), max(args.erange)
    
        vasprun = Vasprun(os.path.join(path, "vasprun.xml"))
        dos_data_total = vasprun.complete_dos       # get dos data

        print(f"\nReading vasprun.xml file from\n    `{abspath(os.path.join(path, 'vasprun.xml'))}` \n" \
               "for DOS analysis...")
        print(f"    Fermi level : {efermi:.5f} eV")
        print(f"    DOS Gap     : {dos_data_total.get_gap():.5f} eV")

        if args.sub_fermi:
            left, right = left + efermi, right + efermi
            print(f"\nCalculated DOS Energy Range: {left}, {right}")
        else:
            print(f"\nCalculated DOS Energy Range: {left}, {right}")
        erange = left, right
        
        structure = dos_data_total.structure
        dos_df = gen_dos_df(dos_data_total, left, right, spin=spin)
        
        if len(args.extra) == 0:
            print()
            print(dos_df.sort_values('dos', ascending=False))
            norb, simple_res_df = dos_analysis_df(dos_df, lb=lb)
            nwann = int(norb * args.deg)
            print(f"\nNumber of Selected Orbitals: {norb}")
            print(f"\nNumber of Selected WFs: {nwann}")
            print("\nSelected Orbitals: ")
            print(simple_res_df)
            print("\nWannier90 Projection:")
            print(w90_string(simple_res_df, structure))
            print("\npyw90 --extra input:")
            print(extra_string(simple_res_df, structure))

            if args.plot:
                selected = select_str2list(args.extra)
                plot_dos_dis(dos_df, selected=selected, savefig=os.path.join(os.getcwd(), 'dos_analysis.png'))

        else:
            selected = select_str2list(args.extra)
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
                dis_pdos = dos_given_selected(dos_data_total, (fmin, fmax), selected)
                dis_tdos = dos_distribute(vasprun.tdos.energies, vasprun.tdos.densities[spin], (fmin, fmax))
                percent  = dis_pdos / dis_tdos
                dis_tdos_l.append(dis_tdos)
                dis_pdos_l.append(dis_pdos)
                percent_l.append(percent)

            dis_froz_dos_df = dis_froz_df.assign(pdos=dis_pdos_l, tdos=dis_tdos_l, percent=percent_l)
            dis_froz_dos_df = dis_froz_dos_df.sort_values('percent', ascending=False)
            N = len(dis_froz_dos_df)
            print()
            print(dis_froz_dos_df)

            dis_win_max = w90.suggest_win_max(erange[0])
            print(f'\nLowest `dis_win_max` for {erange[0]}: {dis_win_max}')

            if args.plot:
                # print(f'Use best dis frozen window in above table to regenerate dos_df to plot ...')
                # left  = dis_froz_dos_df['dis_froz_min'][0]
                # right = dis_froz_dos_df['dis_froz_max'][0]
                dos_df = gen_dos_df(dos_data_total, left, right)
                plot_dos_dis(dos_df, selected=selected,
                             savefig=os.path.join(os.getcwd(), 'dos_analysis_selected.png'))

if __name__ == "__main__":
    args = get_args()
    main_features(args)