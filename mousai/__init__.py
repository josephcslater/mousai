# -*- coding: utf-8 -*-
# Copyright (c) Joseph C. Slater developers.
# Distributed under the terms of the BSD-3-Clause License.

"""
Mousai is a set of implementations of the Harmonic Balance Method
that can be coupled with external codes through connectivity modules.

Documentation is available in the docstrings and online at
https://github.com/josephcslater/mousai.
Joseph C. Slater

.. code-block:: python

    >>> import mousai as ms
"""
# import sys

__title__ = 'mousai'
__version__ = '0.2'
__author__ = 'Joseph C. Slater'
__license__ = 'BSD-3-Clause'
__copyright__ = 'Copyright 2017 Joseph C. Slater'
__name__ = 'mousai'
__package__ = 'mousai'

import sys
import matplotlib as mpl

if 'pytest' in sys.argv[0]:
    print('Setting backend to agg to run tests')
    mpl.use('agg')

from .har_bal import *

# __all__ = ['har_bal', '__version__']

# When everything is imported from mousai, define "everything below"
'''__all__ += ["hb_so",
            "harmonic_deriv",
            "solmf",
            "duff_osc",
            "hb_so_err"]
'''
