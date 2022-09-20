# pyw90

A `VASP` and `Wannier90` interfaced tool for projection analysis and fully automated dis energy window optimization.

## Key features

1. Show distribution of eigenvalues.
2. Pre-analysis before `Wannier90` interpolation with projection and dis energy window recommendation
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

## Example

