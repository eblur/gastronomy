#! /usr/bin/env python

from setuptools import setup

import glob

##-------------------------------------------------
## Package setup

setup(name='gastronomy',
      version='0.0',
      description='Library for handling molecules and minerals optical properties',
      author='Lia Corrales',
      author_email='liac@umich.edu',
      url='https://github.com/eblur/gastronomy',
      packages=['gastronomy'],
      package_data={'gastronomy': ['tables/*']}
)
