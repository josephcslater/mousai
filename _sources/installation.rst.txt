Installation
------------

Easy Installation
_________________

If you aren't familiar at all with Python, please see  `Installing Python <https://github.com/josephcslater/mousai/blob/master/docs/Installing_Python.rst>`_.

Installation is made easy with ``pip`` (or ``pip3``), with releases as we have time while we try
to create a full first release. Much of it works already, but we certainly need
issue reports (on `github <http://github.com/josephcslater/mousai>`_).

To install::

  pip install --user mousai

where ``--user`` isn't necessary if you are using a locally installed version of Python such as `Anaconda <https://www.continuum.io/downloads>`_.

To run, I recommend you open a `Jupyter`_ notebook by using ``jupyter notebook`` and then type::

  import mousai as ms

For examples, see the `example ipynb notebooks <https://github.com/josephcslater/mousai/tree/master/docs/tutorial>`_. Some of these have interactive capabilities that are only apparent when you load them with `Jupyter`_ instead of just looking at them on github.

Installation of current development version
___________________________________________

The usage documentation is far behind the current code, while the reference is way ahead of the released code due to the `autodoc <http://www.sphinx-doc.org/en/stable/ext/autodoc.html>`_ capability of `Sphinx <http://www.sphinx-doc.org/en/stable/>`_. Especially as of early 2017, the code is in rapid development. So is the documentation. Releases to `pypi <https://pypi.python.org/pypi>`_ are far behind current status as stopping to deploy would cost more time that it is worth. We have the objective of releasing a first non-beta version at the end of May, but even this cannot be promised.

If you wish to stay current with the software, fork the repository, then clone it to your desktop.

Then, in the root directory of ``mousai``, type::

  pip install -e .

Now if sync/update your clone, you will have the current version each time you do so.

The added benefit is that you can also develop and contribute by submitting pull requests. Please don't edit inside the ``mousai`` directories unless you intend to submit a pull request because your edits may be overwritten or lost.

That should be it. Please note issues on the `issues tab <https://github.com/josephcslater/mousai/issues>`_ on github.

.. _Jupyter: jupyter.org
