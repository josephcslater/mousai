Notes for Developers
--------------------

To develop, it's best to NOT use ``pip install mousai`` as you will not know which mousai you are running versus testing.

The easiest way is to:

1. Fork Mousai
2. Clone to your hard drive/computer
3. To work in `developer mode <https://packaging.python.org/distributing/#working-in-development-mode>`_, at the top level directory inside ``mousai`` at a terminal (or anaconda) prompt, type::

    $ pip install -e .

Submitting a pull request
-------------------------

If you make minor changes, before submitting a pull request (most likely via github), please::

    $ git pull origin
    
Then run tests (``pytest`` above), then push to your github repository. 

After that, please submit your pull request. 

Testing Docstrings
------------------

Typing ``pytest`` at the top of the directory structure will run all unit tests in docstrings.

Making Distributions
--------------------

To make a distribution. (not completely valid- pulled from Vibration Toolbox)

1) Edit the version number in ``mousai/__init__.py``
2) Use the Makefile, ``make release``

The ``conf.py`` file for the documentation pulls the version from ``__init__.py``

To make a wheel file to test before deployment::

  >>> make wheel

Testing a release
-----------------

To test before release::

  >>> pip install --force-reinstall --upgrade --no-deps dist/mousai-0.5b9-py3-none-any.whl

See ``create_distro.rst`` for explicit ``pypi`` commands that may not be necessary.

See `twine notes <https://packaging.python.org/distributing/#working-in-development-mode>`_ on modern pypi connectivity.

Once travis-ci checking is done, release checking will be unecessary. 
Still needs work on travis-ci.

Debugging/writing
-----------------

Please avoid using `print` and instead use `logging <https://inventwithpython.com/blog/2012/04/06/stop-using-print-for-debugging-a-5-minute-quickstart-guide-to-pythons-logging-module/>`_. 

See the top of `har_bal.py` for the `logging` commands. The `logging.debug` should be used in place of a `print` statement.

Uncommenting the appropriately commentes line will turn off all debuging output. Please do this before commits to the main branch. 
