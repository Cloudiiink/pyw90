#-*- encoding: UTF-8 -*-
from setuptools import setup, find_packages

# REF https://www.jianshu.com/p/eb27d5cb5e1d

VERSION = '0.1.0'

setup(name='pyW90',
      version=VERSION,
      description="python Command-line toolbox for VASP and Wannier90 interface with utility",
      long_description='To be filled',
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='python wannier90 VASP DFT terminal',
      author='Cloudiiink',
      author_email='wangenzj@outlook.com',
      url='https://github.com/cloudiiink/pyw90',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=True,
      install_requires=[
        'requests',
      ],
      entry_points={
        'console_scripts':[
            'pyw90 = pyw90.pyw90:main'
        ]
      },
)