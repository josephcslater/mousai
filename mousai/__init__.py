"""
Implementations of the Harmonic Balance Method.

Mousai provides tools for harmonic balance solutions.

Solvers can be coupled with external codes through connectivity modules.

Documentation is available in the docstrings and online at
https://github.com/josephcslater/mousai.
Joseph C. Slater

.. code-block:: python

    >>> import mousai as ms

Copyright (c) Joseph C. Slater

Distributed under the terms of the BSD-3-Clause License.
"""

__title__ = 'mousai'
# 1) Version may have other than numerical digits after decimal point.
# 2) 1.11 is actually a higher release than 1.2 (confusing).
# 3) Let's just increment with single digits.

# Moved to Semantic Versioning. 0.3 was last non-semantic release.
__version__ = '0.3.1'
__author__ = 'Joseph C. Slater'
__license__ = 'BSD-3-Clause'
__copyright__ = 'Copyright 2017 Joseph C. Slater'
# __name__ = 'mousai'
# __package__ = 'mousai'

import sys
import matplotlib as mpl

from .har_bal import *

# Accomodate using Travis-ci and potential matlotlib results.
if 'pytest' in sys.argv[0]:
    print('Setting backend to agg to run tests')
    mpl.use('agg')

# This would have value for using
# from mousai import *
# That's bad practice and I'd rather not support it
# __all__ = ['har_bal', '__version__']
