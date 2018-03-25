"""
Implementations of the Harmonic Balance Method.

Mousai provides tools for harmonic balance solutions.

Solvers can be coupled with external codes through connectivity modules.

Documentation is available in the docstrings and online at
https://github.com/josephcslater/mousai.
Joseph C. Slater

.. code-block:: python

    >>> import mousai as ms

"""

import sys
import matplotlib as mpl

from .har_bal import *

# Copyright (c) Joseph C. Slater developers.

# Distributed under the terms of the BSD-3-Clause License.


__title__ = 'mousai'
# version may have no more than numerical digits after decimal point.
# 1.11 is actually a higher release than 1.2 (confusing)
__version__ = '0.2.5'
__author__ = 'Joseph C. Slater'
__license__ = 'BSD-3-Clause'
__copyright__ = 'Copyright 2017 Joseph C. Slater'
# __name__ = 'mousai'
# __package__ = 'mousai'


if 'pytest' in sys.argv[0]:
    print('Setting backend to agg to run tests')
    mpl.use('agg')



# __all__ = ['har_bal', '__version__']

# When everything is imported from mousai, define "everything below"
'''__all__ += ["hb_time",
            "harmonic_deriv",
            "solmf",
            "duff_osc",
            "hb_time_err"]
'''
