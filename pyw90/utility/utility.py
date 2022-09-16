import numpy as np
from numpy.typing import ArrayLike
from typing import Callable
import os, re, copy

from pymatgen.core import structure

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

def num2str(l: list[int], connect: str='-', sep :str='|') -> str:
    r"""
    Replace non negative continous number in list to formatted string.
    
    e.g. [1, 2, 3, 4, 5, 7, 8, 9, 11] -> '1-5|7-9|11'
    """
    l = sorted(set(l))
    res = []
    left = right = None
    for v in l + [None]:
        if right != None and v == right + 1:
            right += 1
        else:
            if left != None:
                s = str(left) if left == right else f'{left}{connect}{right}'
                res.append(s)
            left = right = v

    return sep.join(res)

def str2num(s: str, connect: str='-', sep: str='|') -> list[int]:
    r"""
    Replace formatted string to list
    
    e.g.'1-5|7-9|11' -> [1, 2, 3, 4, 5, 7, 8, 9, 11]
    """
    segs = s.strip().split(sep)
    res = []
    for seg in segs:
        if connect in seg:
            lseg = seg.strip().split(connect)
            left, right = int(lseg[0]), int(lseg[-1])
            res += list(range(left, right+1))
        else:
            res.append(int(seg))
    return res

def get_id_from_specie(struc: structure, specie: str) -> list[int]:
    r"""
    Return `struc` structure indices with `spcies`
    """
    return [i for i in range(len(struc)) if struc[i].specie.name == specie]

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

    for i, name in enumerate(sorted(family_name)):
        print('{:25s}'.format(name), end='')
        if i % 4 == 3: print()

def get_efermi(args, from_file=False) -> float:
    r'''
    Return the Fermi level from `vasprun.xml` in `args.path` folder.
    
    If `from_file` is False, function will return `args.efermi` (Fermi level from the argument and not `None`) directly.
    '''
    if not from_file and args.efermi != None:
        efermi = float(args.efermi)
    else:
        # The input args.path is the relative path, converting to absolute path of the directory.
        vasprun = os.path.join(os.getcwd(), args.path, 'vasprun.xml')
        efermi_str = os.popen(f'grep fermi {vasprun}').read().strip()
        m = re.match('.+ ((\-|\+)?\d+(\.\d+)?) .+', efermi_str)
        efermi = float(m.groups()[0])
    return efermi

if __name__ == '__main__':
    print(num2str([1, 2, 3, 4, 5, 7, 8, 9, 11]))
    print(str2num('1-5|7-9|11'))