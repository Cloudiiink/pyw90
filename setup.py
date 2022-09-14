#-*- encoding: UTF-8 -*-
from setuptools import setup, find_packages

# REF https://www.jianshu.com/p/eb27d5cb5e1d

VERSION = '0.1.0'

setup(name='pyW90',
      version=VERSION,
      python_requires='>=3.8',
      description="A tool interfaced to VASP and Wannier90 with projection analysis and automatically dis energy window optimization.",
      long_description='Command-line tool for constructing tight-binding model interfacing with VASP and Wannier90. Key features includes: 1. Show distribution of eigenvalues. 2. Pre-analysis before `Wannier90` interpolation with projection and dis energy window recommendation 3. Auto Wannier90 Fit. Using minimize method to choose the most suitable dis energy windows. 4. (Comparison) Show difference between VASP bands and Wannier90 bands via plotting and report. ',
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='python wannier90 VASP tight-binding solid-state-physics',
      author='Cloudiiink',
      author_email='wangenzj@outlook.com',
      url='https://github.com/cloudiiink/pyw90',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=True,
      install_requires=[
        'pymatgen',
        'scipy',
        'numpy',
        'PyYAML',
        'pandas'
      ],
      entry_points={
        'console_scripts':[
            'pyw90 = pyw90.pyw90:main'
        ]
      },
)