import numpy as np
from numpy.typing import ArrayLike
from typing import Callable
import os, re, copy

def gaussian(x:ArrayLike, mid:float= 0, width:float=3, normalized:bool=False) -> ArrayLike:
    r'''
    Return values of Gaussian function with given input. Used as kernel function for band fitting quality evaluation. 

    :param mid: $\mu$. center of distribution
    :param width: $\sigma$. Standard deviation of Gaussian function.
    :param normalized: Default: False
    '''
    mu, sig = mid, width
    if normalized:
        return 1/(np.sqrt(2*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2)/2)
    else:
        return np.exp(-np.power((x - mu)/sig, 2)/2)

def unit(x:ArrayLike, mid:float=0, width:float=3) -> ArrayLike:
    r'''
    Return values of window function with given input. Used as kernel function for band fitting quality evaluation.

    :param mid: Center of the window.
    :param width: Half length of the interval with non-zero value.
    '''
    res = np.abs(x-mid) < width
    return res.astype(float)

def parse_kernel(kernel_str:str, mid:float, width:float) -> Callable[[ArrayLike],ArrayLike]:
    r'''
    Return the lambda function with given parameters.
    '''
    if kernel_str[0].lower() == 'u':
        return lambda x: unit(x, mid=mid, width=width)
    elif kernel_str[0].lower() == 'g':
        return lambda x: gaussian(x, mid=mid, width=width)
    else:
        raise ValueError(f"Unrecognized parameters: {kernel_str}")

def _replace_str_none(l):
    r'''
    Replace unused values in List or dict to `None`
    '''
    res = copy.copy(l)
    none_type = ['none', 'None', 'n', 'N', 'inf', '-inf', '+inf']
    is_none_check = lambda x: True if x in none_type else False
    if isinstance(l, str) and is_none_check(l):
        res = None
    elif isinstance(l, list):
        for idx, ele in enumerate(res):
            if isinstance(ele, list):
                res[idx] = _replace_str_none(ele)
            else:
                if is_none_check(ele):
                    res[idx] = None
    elif isinstance(l, dict):
        for key in res.keys():
            res[key] = _replace_str_none(res[key])
    return res

def _replace_str_boolean(d):
    r'''
    Replace string values to boolean
    '''
    for key in d.keys():
        if d[key][0].lower() == 't':
            d[key] = True
        elif d[key][0].lower() == 'f':
            d[key] = False
        else:
            raise ValueError(f'Wrong Value in {d}')
    return d

def show_all_fonts():
    r'''
    Print the list of all the fonts currently available for Matplotlib
    '''
    # Ref: https://stackoverflow.com/questions/8753835/how-to-get-a-list-of-all-the-fonts-currently-available-for-matplotlib
    from matplotlib import font_manager as fm
    fpaths = fm.findSystemFonts()
    family_name = set([fm.get_font(i).family_name for i in fpaths])
    print(sorted(family_name))

def get_efermi(args, from_file=False) -> float:
    r'''
    Return the Fermi level from `vasprun.xml` in `args.path` folder.
    
    If `from_file` is False, function will return `args.efermi` (Fermi level from the argument and not `None`) directly.
    '''
    if not from_file and args.efermi:
        efermi = float(args.efermi)
    else:
        # The input args.path is the relative path, converting to absolute path of the directory.
        vasprun = os.path.join(os.getcwd(), args.path, 'vasprun.xml')
        efermi_str = os.popen(f'grep fermi {vasprun}').read().strip()
        m = re.match('.+ ((\-|\+)?\d+(\.\d+)?) .+', efermi_str)
        efermi = float(m.groups()[0])
    return efermi