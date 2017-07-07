Mousai
======

.. image:: https://badge.fury.io/py/mousai.png/
    :target: http://badge.fury.io/py/mousai

.. image:: https://travis-ci.org/josephcslater/mousai.svg?branch=master
    :target: https://travis-ci.org/josephcslater/mousai

Python-based Harmonic Balance solvers and relevant tools.

Can solve first-order and second-order ordinary differential equations written in state-space form (solved for acceleration for second-order). All you need to provide is the name of a python function. The function accepts:
1. The states (displacements and derivatives for second order forms)
2. Time (interpreted as time within the period)
3. Fundamental frequency of the harmonic representation

The function must return the state derivatives (acceleration for second order form). It is expected that in most cases this will simply be a python wrapper function to call an external finite element, electro-magnetic, molecular dynamic, computational fluid dynamic, or other simulation code.

No traditional numerical integration is performed in harmonic balance, so issues of stability, reliability, etc. are separate concerns. An inability to determine a solution is reported if failure occurs A good initial guess is always helpful, sometimes essential, and may simply have the form of a nearby solution. Sweeping through frequencies or increasing amplitudes of excitation is the best way to ensure swift convergence.

Quasi-linear models can (not yet implemented) obtain good low amplitude solutions that can be used as starting guesses for solutions (appropriate scaling may also help).

The manual is still under development, but the Tutorial provides two examples of solutions to the Duffing Oscillator, one in first order form, the other in second order form.

The Reference section of the manual illustrates supporting code for building out more time-refined solutions, finding velocities from displacements, and specifics on function calls. 

Please see the `manual <https://josephcslater.github.io/mousai/>`__ for usage and installation.
