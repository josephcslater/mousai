#!/usr/bin/env python

# from distutils.core import setup
import os
import sys
from setuptools import setup

if sys.version_info < (3, 5):
    sys.exit('Sorry, Python < 3.5 is not supported.')
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    """You want a docstring? Tell me why."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


with open('harmonic_balance/__init__.py', 'rb') as fid:
    for line in fid:
        line = line.decode('utf-8')
        if line.startswith('__version__'):
            version = line.strip().split()[-1][1:-1]
            break


setup(name='harmonic_balance',
      version=version,
      description='Harmonic Balance solvers.',
      author='Joseph C. Slater',
      author_email='joseph.c.slater@gmail.com',
      url='https://github.com/josephcslater/harmonic_balance',
      download_url='https://github.com/josephcslater/harmonic_balance',
      packages=['harmonic_balance'],
      package_data={'harmonic_balance': ['../README.md'], '': ['README.md']},
      long_description=read('README.md'),
      keywords=['harmonic balance', 'numerical solvers',
                'vibration', 'oscillation'],
      install_requires=['numpy', 'scipy', 'matplotlib'])

# https://docs.python.org/3/distutils/setupscript.html#additional-meta-data
