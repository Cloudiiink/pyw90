from math import sqrt
from pymatgen.electronic_structure.core import Orbital, OrbitalType

def check_orb_name(func):
    def __wrapper__(*args, **kwargs):
        if args[1].strip() not in args[0].orbital_names:
            raise ValueError(f'Unknown orbital name: {args[1]}')
        return func(*args, **kwargs)
    return __wrapper__

def check_orb_type_name(func):
    def __wrapper__(*args, **kwargs):
        if args[1].strip() not in args[0].orbital_type_names:
            raise ValueError(f'Unknown orbital type name: {args[1]}')
        return func(*args, **kwargs)
    return __wrapper__

class DOS():
    r"""
    Class offer orbitals and DOS (density of states) manipulation.
    """
    
    # Orbital Notation for VASP since 5.4.4 also have some difference compare to
    # list `orbital_names` below. Here we changed some of the notations 
    # (x2-y2 -> dx2) for comparision between `VASP` and `Wannier90`.

    # See https://www.vasp.at/forum/viewtopic.php?t=18044
    # https://www.vasp.at/forum/viewtopic.php?f=4&t=6924

    # TODO Need to check for f-electron
    # v5.4.4        v5.4.1      pymatgen
    # x2-y2    <-->   dx2   <-->   dx2
    # fy3x2    <-->   f1    <-->   f_3
    # fxyz     <-->   f2    <-->   f_2
    # fyz2     <-->   f3    <-->   f_1
    # fz3      <-->   f4    <-->   f0
    # fxz2     <-->   f5    <-->   f1
    # fzx2     <-->   f6    <-->   f2
    # fx3      <-->   f7    <-->   f3

    # pymatgen
    pymat_orbital_names = ['s', 
                          'py', 'pz', 'px', 
                          'dxy', 'dyz', 'dz2', 'dxz', 'dx2',
                          'f_3', 'f_2', 'f_1', 'f0', 'f1', 'f2', 'f3']
    vasp_orbital_names  = ['s', 
                          'py', 'pz', 'px', 
                          'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2',
                          'fy3x2', 'fxyz', 'fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3']
    orbital_type_names  = ['s', 'p', 'd', 'f']
    # ordered according to `mr` value
    w90_orbital_names   = ['s',
                          'pz', 'px', 'py',
                          'dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy',
                          'fz3', 'fxz2', 'fyz2', 'fz(x2-y2)', 'fxyz', 'fx(x2-3y2)', 'fy(3x2-y2)']
    orbital_names       = pymat_orbital_names
    vasp2w90            = {'s'    : 's',
                          'py'    : 'py',
                          'pz'    : 'pz',
                          'px'    : 'px', 
                          'dxy'   : 'dxy',
                          'dyz'   : 'dyz',
                          'dz2'   : 'dz2',
                          'dxz'   : 'dxz',
                          'dx2'   : 'dx2-y2',
                          'fy3x2' : 'fy(3x2-y2)', 
                          'fxyz'  : 'fxyz',
                          'fyz2'  : 'fyz2', 
                          'fz3'   : 'fz3',
                          'fxz2'  : 'fxz2',
                          'fzx2'  : 'fz(x2-y2)' , 
                          'fx3'   : 'fx(x2-3y2)'}
    pymat2w90           = dict([(k, v) 
                              for k, v in zip(pymat_orbital_names, vasp2w90.values())])

    @classmethod
    @check_orb_name
    def get_mr(cls, name):
        name = cls.vasp2w90[name]
        value = cls.w90_orbital_names.index(name) + 1
        return value - int(sqrt(value))**2
    
    @classmethod
    @check_orb_name
    def get_l(cls, name):
        return cls.orbital_type_names.index(name[0])

    @classmethod
    def is_orbital_type(cls, name):
        return name in cls.orbital_type_names

    @classmethod
    def is_orbital(cls, name):
        return name in cls.orbital_names

    @classmethod
    @check_orb_name
    def get_orbital(cls, name):
        r"""
        Return `Orbital` with given name. The indices are basically the order in
        which the orbitals are reported in VASP and has no special meaning.

        - `s   = 0`
        - `py  = 1`  `pz  =  2`  `px  =  3`
        - `dxy = 4`  `dyz =  5`  `dz2 =  6`  `dxz =  7`  `dx2 = 8`
        - `f_3 = 9`  `f_2 = 10`  `f_1 = 11`  `f0  = 12`  `f1 = 13`  `f2 = 14`  `f3 = 15`
        """
        return Orbital(cls.orbital_names.index(name))

    @check_orb_type_name
    def get_orbital_type(cls, name):
        r"""
        Return `Orbital` with given name. 

        s = 0, p = 1, d = 2, f = 3
        """
        return OrbitalType(cls.orbital_type_names.index(name))

    @classmethod
    @check_orb_type_name
    def get_orbitals_from_type(cls, name):
        indices = [i for i, s in enumerate(cls.orbital_names) if s[0]==name]
        return [Orbital(i) for i in indices]

if __name__ == '__main__':
    print(DOS.get_mr('px'), DOS.get_l('px'))
    # DOS.get_orbitals_from_type('f')