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

from .har_bal import *

__title__ = 'mousai'
__version__ = '0.1a'
__author__ = 'Joseph C. Slater'
__license__ = 'BSD-3-Clause'
__copyright__ = 'Copyright 2017 Joseph C. Slater'
__all__ = ['har_bal', '__version__']
