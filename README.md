# pyw90

A `VASP` and `Wannier90` interfaced tool for projection analysis and fully automated dis energy window optimization.

## Key features

1. Show distribution of eigenvalues.
2. Pre-analysis before `Wannier90` interpolation with projection and dis energy window recommendation.
3. Auto Wannier90 Fit. Using minimize method to choose the most suitable dis energy windows. 
4. Comparison. Show difference between `VASP` bands and `Wannier90` bands via plotting and report.

## Installation via `pip`

`pyw90` is installed using `pip` by

```bash
pip install pyw90
# update the package
pip install pyw90 --upgrade
```

The `pyw90` dependencies are listed as following:

- python (>= 3.8, < 3.12)
- pymatgen
- scipy (>= 1.8)
- numpy (>= 1.20.1)
- python-yaml (PyYAML)
- pandas
- matplotlib (>=3.4)

## Usage

You can use `pyw90 -h` to see help message. You can also use `-h` in submenu (e.g. `pyw90 auto -h`).

```
usage: pyw90 [-h] {eig,pre,auto,cmp} ...

python Command-line toolbox for VASP and Wannier90 interface with utility. You
can also use -h to display the submenu help message. e.g. pyw90 pre -h

positional arguments:
  {eig,pre,auto,cmp}  Main features
    eig               Show distribution of eigenvalues.
    pre               (Pre)-analysis before `Wannier90` Interpolation.
    auto              (Auto Wannier90 Fit) Using minimize method to choose the
                      most suitable dis energy windows.
    cmp               (Comparison) Show difference between VASP bands and
                      Wannier90 bands via plotting and report. `bnd.dat` for
                      VASP band data in `p4vasp` format and
                      `wannier90_band.dat`, `wannier90_band.labelinfo.dat`, and
                      `wannier90.wout` are required for plotting and analysis.

optional arguments:
  -h, --help          show this help message and exit
```

### 1. `eig` Menu

This menu allows you have knowledge to the full distribution of eigenvalues. This may help you to decide the `exclude_bands` block in `Wannier90` and the `dis_win_min` and `dis_win_max`.

The help message of `eig` menu is listed as below.

```
usage: pyw90 eig [-h] [-e ERANGE ERANGE] [--config] [--path PATH] [-i EIG]
                 [--rm-fermi] [--efermi EFERMI] [-w NWANN] [-n NBNDS_EXCL]
                 [--deg NDEG] [--separate] [--eps EPS]
                 mode

positional arguments:
  mode              Mode: report, plot, count, suggest

optional arguments:
  -h, --help        show this help message and exit
  -e ERANGE ERANGE  Energy range.
  --config          Read input from config file `auto_w90_input.yaml` directly.
                    Default: False
  --path PATH       The path of working dir. Default: .
  -i EIG            Select wannier90.eig file or EIGENVAL file. Default:
                    EIGENVAL
  --rm-fermi        Whether or not the input `erange` has removed the Fermi
                    energy is indicated by this flag. Default: False
  --efermi EFERMI   Fermi level. Default value is generated from `vasprun.xml`.
  -w NWANN          Number of Wannier Functions. Default: 0
  -n NBNDS_EXCL     Number of bands excluded below the bands from `Wannier90`.
  --deg NDEG        Number of degeneracy. Default: 1
  --separate        Calculate bands not separately.
  --eps EPS         Tolerance for dis energy window suggestion. Default: 0.004
```

The `report` mode allow to directly print the results to terminal. But in `plot` mode, program will generate `dos_analysis.pdf` of bands.

The `count` can calculate how many states at least and at most for each k-point inside the giving `erange`. And the `suggest` mode are used to for dis energy window suggestion. 

Program reads the data from `EIGENVAL` file. You can specify the `EIGENVAL` location by `--path` argument or direct input the absolute or relative path to `EIGENVAL` by `-i` argument. `-e` arguments can be used to control the output range in `report` mode and `plot` mode. It can also be used to calculate the states numbrt inside the energe interval and the suggestion of `dis_win_min` and `dis_win_max` values.  `-w` argument reads the input as the number of Wannier functions (WFs). 

We directly reads in the data from `EIGENVAL` without modification. So the Fermi level for materail usually is non-zero value. If you input the `--rm-fermi` argument, the program will be shift and treat fermi level as 0. This will affects the `-e` argument. The default Fermi level is read from the `vasprun.xml` file in `--path` folder. You can also specify other number you want become Fermi level via `--efermi` argument.

Other arguments:

- `-n NBNDS_EXCL`. Number of bands exclude from the bottom bands. These bands will not displayed in `report` mode and `plot` mode.
- `--deg NDEG`. Number of degenercy of each bands. If you want to consider the extra spin freedom degree for number of Wannier functions, please input `--deg 2`. Default value is 1.
- `--eps EPS`. Since we only do the finite samplign in Brilliouin zone and the numerical error, there might some missing states when we count the number of states inside the energy interval. This might cause error in `Wannier90`. So we substract or add `eps` to the minimum and maximum of each bands to make sure the input is much reliable.

**One of the key feature of `eig` menu is the dis energy suggestion.** The outer window suggestion requires `-e`  input. The `dis_win_min` values generates from the midpoint of global gap below the Fermi level. And based on the `dis_win_min`, `dis_win_max` is also generated which is the lower limit with given number of WFs. `Wannier90` restricts that the number of states should be at least number of WFs (#WFs for simplicity). 

The inner window suggestion can be generated with `-w` input (the default value of #WFs is 0). Based on all recommended `dis_win_min`, we also treat them as `dis_froz_min` to generate possible `dis_froz_max`. 

The basic procedure of generate `dis_froz_max` is count how many states at least over all k-points in Brilliouin zone. We assume there are $N>N_{\text WF}$ states at least. Then program will generate upper bound $E_{\text{fmax}}$ for at least $N_i = N_{\text{WF}} + i$ states inside the given `dis_win_min`. And then generate the lower bound of `dis_froz_min` by counting down. This procedure ensures the generated dis frozen energy window satisfy the requirement of `Wannier90`. But this simple solution might ignore some valuable frozen energy window due to the degenerate points. Considering the mis-sampling in BZ and the numerical errpr, once we use procedure above to generate a degenerate `dis_froz_max`, add `eps` to output will add an additional states inside the frozen window. So the `dis_froz_min` will skip one state and become discontinous.

To avoiding situation above, we add additional check whether there are skipped states between the input $E_{\text{fmin}}$ and the `dis_froz_min`. If there exist $\Delta N$ skipped states, we will go through to check all situation of with $N'_{\text{WF}}=N_{\text{WF}}-1, \cdots, N_{\text{WF}}-\Delta N$ to generate possible dis frozen energy windows.

### 2. `pre` Menu

This menu offers some basic input of `Wannier90`. The help message of `pre` menu is listed as below.

```
positional arguments:
  mode              Mode: kpath, band, template, dos

optional arguments:
  -h, --help        show this help message and exit
  --path PATH       The path of working dir. Default: .
  -e ERANGE ERANGE  Energy range.
  --lb LB           Lower bound for selected orbital / max single orbital.
                    default: 0.1
  --rm-fermi        Whether or not the input `erange` has removed the Fermi
                    energy is indicated by this flag. Default: False
  --extra EXTRA     Extra input. In `template` mode and within extra input
                    (basic, wann, band), we can choose one of the detailed parts
                    to print.In `dos` mode and within extra input (`species`,
                    `structure_id`, `orbital_id` list separated by ;), we can
                    treat input as projections for `Wannier90` input to suggest
                    dis frozen energy. Details can be found in the document.
  --spin-down       Specify the spin channel to `Spin.down`. Without this
                    argument, the default one is `Spin.up`.
  --plot            plot the dos distribution
  --eps EPS         Tolerance for dis energy window suggestion. Default: 0.004
  --deg DEG         Degeneracy of bands. Default: 1
```

The `kpath`, `band` and `template` modes each provide different necessary input sections for `Wannier90` respectively.

- `kpath` : Converting `VASP` input file `KPOINTS` with line mode to `Wannier90` input for `Kpoint_Path` block. 
- `band` : Converting `VASP` band structure to a simpler k-E `.dat` file in `p4vasp` format. Default file name is `bnd.dat`. If the system is spinful, it will exported to `bnd_up.dat` and `bnd_down.dat` for different spin channel separatly. 
- `template` : Print basic input for `Wannier90` including basic file structure of `seedname.win`, dis energy window block and band structure calculating block.

**The `dos` mode is the most important feature of `pre` menu.** From the first-principle calculation in VASP with `LORBIT=11`, `vasprun.xml` are written with decomposed projection information.
> For further details, see [LORBIT - Vaspwiki](https://www.vasp.at/wiki/index.php/LORBIT).

For accurate and complete present the results from `VASP`, `Wannier90` needs to select the most representive projection as the initial guess for `projection` block. Here we use projected density of states (pdos) and total density of states (tdos) to analyze and give the suggestion of `dis_froz_min(max)` based on the ratio of pdos / tdos. 

**In `kpath` and `band` mode**, you need to input the path to the folder where `VASP` calculation results located via `--path`. The default path is the current folder. **In `template` mode**, default settings will print all stored input template for `Wannier90`. You can also input key words `basic`, `wann` and `band` to print the basic structure of `seedname.win`, dis energy window block and band structure calculating block separately.

In `dos` mode, you also need to input location of `vasprun.xml` via `--path`. We can calculate the pdos distribution of each orbital at each atom inside desired energy interval from `-e` argument. The `--rm-fermi` argument is used if you want to assuming the input energy interval has already remove the Fermi level.

$$
\text{pDOS} = \int_{E_1}^{E_2} \text{dos}_{\text{atom, orb}}(E)\,\mathrm{dE}
$$

The result of pdos is presented as table with columns including: `species`, `structure_id`, `orb_id`, `orb_name`, `key_string` and `dos`. The `structure_id` is obtained from `vasprun.xml` (also same as listed in `POSCAR`) with index started from $0$. The `orb_id` and `orb_name` use the notation in pymatgen (same as `VASP` in old version). 

> For further details, see note at the end of the section.

The `key_string` is formatted with `species`, `structure_id` and `orb_name` (e.g. `Ga_0_px`). The dos column is calculated from above formula. Since we only need to compare the pdos relatively, the maximum of this column is normalized to $1$.

The `lb` argument is used to as lower bound to select the most reprentitive projection. Only when the results in `dos` column are larger than `lb`, the orbitals and sites will be selected as the `projection` block input for `Wannier90`. We will print the total number of selected orbitals. And the selection information will also be simplified with merging the situation in which orbitals locate in one site and all sites of same species have the same orbitals. The results display as a table with `species`, `site` and `orb` columns. The default number of `site` column is `structure_id`. When all sites of same species have been selected, the value of `site` column will be $-1$.

The selected information will be converted to the format for `projection` block in `Wannier90`. These projection information will also be printed as `key_string` format for `pyw90` input. (The represention of orbital is `orb_id`. Default delimeter for keys inside one species is comma and the delimeter between species is ;. )

With `--extra` input in `dos` mode and `-e` given energy interval, we will also present the dis frozen window suggestion based on pdos. The result will be presented with columns including `dis_froz_min`, `dis_froz_max`, `N`, `pdos`, `tdos` and `percent`. The final column of table is the percentage of pdos/tdos. The full table is also sorted with `percent` in descending order. The lower bound of `dis_win_max` is also generated.

**Note:**

The notation of orbitals has changed a lot among different software and different version. `pyw90` use the same notation as `pymatgen` (2022.8.23). All orbitals are listed as below:

```bash
# pymatgen 2022.8.23
s 
py   pz   px
dxy  dyz  dz2  dxz  dx2
f_3  f_2  f_1  f0   f1   f2  f3

# VASP 5.4.4
s
py     pz    px 
dxy    dyz   dz2   dxz  x2-y2
fy3x2  fxyz  fyz2  fz3  fxz2   fzx2  fx3

# Wannier90
s
pz   px    py
dz2  dxz   dyz   dx2-y2     dxy
fz3  fxz2  fyz2  fz(x2-y2)  fxyz  fx(x2-3y2)  fy(3x2-y2)
```

### 3. `auto` Menu

This menu allows you do fully automated dis energy optimization and help message of `auto` menu is listed as below.

```
usage: pyw90 auto [-h] [--path PATH] [--pid PID] mode

positional arguments:
  mode         Mode: run, term(inate), input. Only first character is
               recognized.

optional arguments:
  -h, --help   show this help message and exit
  --path PATH  The path of working dir. Default: .
  --pid PID    PID to terminate.
```

In `input` mode, program will generate the configuration file `auto_w90_input.yaml` in foloar `--path` which contains the all necessary message for running `Wannier90`. Some the parameters is deplicated with other menu, we also add `--config` argument in other menus to directly get information from configuration.

`run` mode offers two method to automated run `Wannier90` job: via SLURM job system or direct run. The quality of `Wannier90` result is evaulated from

$$
\Delta=\frac{1}{C} \frac{1}{N_{\mathbf{k}}} \sum_{i=1}^{M} \sum_{\mathbf{k}}\left|\varepsilon_{i, \mathbf{k}}^{\mathrm{DFT}}-\varepsilon_{i, \mathbf{k}}^{\mathrm{TB}}\right|{\color{red} f\left(\varepsilon_{i, \mathbf{k}}^{\mathrm{DFT}}\right)}.
$$

$f(\cdot)$ is the kernel function (unit function or gaussian function). $i$ and $\bf k$ represent band index and kpoint separately. $C$ is the normalization constant obtained from

$$
C= \frac{1}{N_{\mathbf{k}}} \sum_{i=1}^{M} \sum_{\mathbf{k}}{\color{red} f\left(\varepsilon_{i, \mathbf{k}}^{\mathrm{DFT}}\right)}
$$

With initial input dis energy window and suitable minimization method, we will get the final result after some iterations. You can terminate your task in `term` mode. Please follow the instrcution message to kill all process. 

⚠ Notice. The server prohibits heavy computing task on the login node generally. So we recommmend you to run the job via slurm system. You can also login to the compute node to run the task locally. Since the `pyw90` will create a daemon process to monitoring the status of your `Wannier90` task, be careful to terminate this daemon process. Actually this process is really easy to find. Because of the necessity to pass the full environment, the process for `pyw90 auto run` usually become extremely long. If you have any question, please ask the administrator of your server.

**Configuration file**

The configuration file is written in YAML format allowing `#` to comment. The contents are case-sensitive and they use indentation to indicate hierarchical relationships. Tabs are not allowed for indentation with only spaces allowed.

| key              | description  |
|------------------|---|
|`winfile`         | File for `Wannier90` input. Default: `wannier90.win` | 
|`runfile`         | Job file for SLURM job system to submit the `Wannier90` task  | 
|`vasp_band_file`  | Band structure file of `VASP` in `p4vasp` format. This file can be generate via `pyw90 pre band`. Default: `bnd.dat`  | 
|`wann_band_file`  | Band structure results of `Wannier90`. Default: `wannier90_band.dat`  | 
|`jobname`         | Name of your task. | 
|`username`        | Your username. `jobname` and `username` is used to identify which is the automated `Wannier90` job from `squeue` result. Default `username` is obtained from `getpass.getuser`.| 
|`local`           | Boolean. If it's `True`, the `Wannier90` task will run via the input from `localrun` | 
|`localrun`        | Command to run `Wannier90`. Only is used when `local: True`  | 
|`efermi`          | Fermi level of material. | 
|`kernel`          | Kernel funtion for evaluating the quality of `Wannier90` band structure result. Formatted input: type, middle, width (e.g. unit,0,1). Kernel functions are classified into two types: `unit` and `gaussian`. **The Fermi level is not subtracted from eigenvalues when using the kernel function to evaluate difference.**  | 
|`method`          | The solving method of minimization. This parameter is treated as the input for `scipy.optimize.minimize`. **We recommend you to use solver allow no gradient information, such as Nelder-Mead, Powell and COBYLA.** (See [scipy.optimize.minimize — SciPy v1.9.1 Manual](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)). Default: COBYLA  | 
|`tol`             | Tolerance for termination. Default: 0.02. | 
|`maxiter`         | Maximum number of iterations to perform. Both `tol` and `maxiter` are also passed to `scipy.optimize.minimize`. Default: 100.| 
|`ini_dis`         | Initial dis energy windows dictionary. Input each value with indent. If some parameters is not needed, please set to `None` or `~`. | 
|`opt_dis`         | Boolean dictionary control which parameters to optimize. `pyw90` will only optimized the key is `True`.  | 
|`bounds`          | Boundary of energy window parameters. Please using `None` or `+/- inf` to set the boundary which is unbounded. Since there is also parameters correction during each iteration of `Wannier90`, you can set a really rough bound.| 
|`num_print_check` | How many iteration to print check message once. Default: 1  |
|`check_time`      | How many seconds to check the job status once. The job status is examined via the result from `squeue` in SLURM system or `pid` when `local` is `True`. Default: 30 |
|`display`         | Boolean dictionary control display the band difference (`diff`) and total spreading of Wannier functions (`spread`). Default: both `True`.   |


**Before running your task, remenber to modify the value of each key to your own.**

### 4. `cmp` Menu

```
usage: pyw90 cmp [-h] [--config] [--path PATH] [--efermi EFERMI] [--vasp VASP]
                 [--seedname SEEDNAME] [--ylim YLIM YLIM] [--kernel KERNEL]
                 [--show-fonts] [--fontfamily FONTFAMILY] [--fontsize FONTSIZE]
                 [--no-spread] [--no-quality] [--quiet]
                 name

positional arguments:
  name                  name of system

optional arguments:
  -h, --help            show this help message and exit
  --config              Read input from config file `auto_w90_input.yaml`
                        directly. Default: False
  --path PATH           The path of working dir. Default: .
  --efermi EFERMI       Fermi level. Default value is generated from
                        `vasprun.xml`.
  --vasp VASP           Path of VASP band file in `p4vasp` format. Default:
                        bnd.dat
  --seedname SEEDNAME   Seedname of Wannier90 input. Default: wannier90
  --ylim YLIM YLIM      Energy bound for plot. Since the Fermi level has been
                        shift to 0 during the plotting, please mind your input.
                        Default: [E_w90.min - 1, E_w90.max + 1]
  --kernel KERNEL       Formatted input: type, middle, width (Defalut:
                        unit,0,1). Kernel functions are classified into two
                        types: `unit` and `gaussian`. **The Fermi level is not
                        subtracted from eigenvalues when using the kernel
                        function to evaluate difference**.
  --show-fonts          Show all availabel font families can be used in
                        `rcParams`
  --fontfamily FONTFAMILY
                        Set font family manually. Default: Open Sans
  --fontsize FONTSIZE   Set font size manually. Default: 18
  --no-spread           Don't plot spreading
  --no-quality          Don't show quality of fitting
  --quiet               Equal to --no-spreading --no-quality
```


## Example

