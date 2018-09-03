"""Harmonic balance solvers and other related tools."""
import warnings
# import logging
import numpy as np
import scipy as sp
import scipy.fftpack as fftp
import scipy.linalg as la
from scipy.optimize import newton_krylov, anderson, broyden1, broyden2, \
    excitingmixing, linearmixing, diagbroyden


# logging.basicConfig(level=print)
# Use `logging.debug` in place of print.
# for instance logging.debug(pformat('e {} X {}'.format(e,X)))

# This will output only info and warnings
# logging.basicConfig(level=logging.INFO)

# This will output only warnings
# logging.basicConfig(level=logging.WARNING)


def hb_time(sdfunc, x0=None, omega=1, method='newton_krylov', num_harmonics=1,
            num_variables=None, eqform='second_order', params={}, realify=True,
            **kwargs):
    r"""Harmonic balance solver for first and second order ODEs.

    Obtains the solution of a first-order and second-order differential
    equation under the presumption that the solution is harmonic using an
    algebraic time method.

    Returns `t` (time), `x` (displacement), `v` (velocity), and `a`
    (acceleration) response of a first- or second- order linear ordinary
    differential equation defined by
    :math:`\ddot{\mathbf{x}}=f(\mathbf{x},\mathbf{v},\omega)` or
    :math:`\dot{\mathbf{x}}=f(\mathbf{x},\omega)`.

    For the state space form, the function `sdfunc` should have the form::

        def duff_osc_ss(x, params):  # params is a dictionary of parameters
            omega = params['omega']  # `omega` will be put into the dictionary
                                     # for you
            t = params['cur_time']   # The time value is available as
                                     # `cur_time` in the dictionary
            xdot = np.array([[x[1]],[-x[0]-.1*x[0]**3-.1*x[1]+1*sin(omega*t)]])
            return xdot

    In a state space form solution, the function must take the states and the
    `params` dictionary. This dictionary should be used to obtain the
    prescribed response frequency and the current time. These plus any other
    parameters are used to calculate the state derivatives which are returned
    by the function.

    For the second order form the function `sdfunc` should have the form::

        def duff_osc(x, v, params):  # params is a dictionary of parameters
            omega = params['omega']  # `omega` will be put into the dictionary
                                     # for you
            t = params['cur_time']   # The time value is available as
                                     # `cur_time` in the dictionary
            return np.array([[-x-.1*x**3-.2*v+sin(omega*t)]])

    In a second-order form solution the function must take the states and the
    `params` dictionary. This dictionary should be used to obtain the
    prescribed response frequency and the current time. These plus any other
    parameters are used to calculate the state derivatives which are returned
    by the function.

    Parameters
    ----------
    sdfunc : function
        For `eqform='first_order'`, name of function that returns **column
        vector** first derivative given `x`, `omega` and \*\*kwargs. This is
        *NOT* a string.

        :math:`\dot{\mathbf{x}}=f(\mathbf{x},\omega)`

        For `eqform='second_order'`, name of function that returns **column
        vector** second derivative given `x`, `v`, `omega` and \*\*kwargs. This
        is *NOT* a string.

        :math:`\ddot{\mathbf{x}}=f(\mathbf{x},\mathbf{v},\omega)`
    x0 : array_like, somewhat optional
        n x m array where n is the number of equations and m is the number of
        values representing the repeating solution.
        It is required that :math:`m = 1 + 2 num_{harmonics}`. (we will
        generalize allowable default values later.)
    omega : float
        assumed fundamental response frequency in radians per second.
    method : str, optional
        Name of optimization method to be used.
    num_harmonics : int, optional
        Number of harmonics to presume. The omega = 0 constant term is always
        presumed to exist. Minimum (and default) is 1. If num_harmonics*2+1
        exceeds the number of columns of `x0` then `x0` will be expanded, using
        Fourier analaysis, to include additional harmonics with the starting
        presumption of zero values.
    num_variables : int, somewhat optional
        Number of states for a state space model, or number of generalized
        dispacements for a second order form.
        If `x0` is defined, num_variables is inferred. An error will result if
        both `x0` and num_variables are left out of the function call.
        `num_variables` must be defined if `x0` is not.
    eqform : str, optional
        `second_order` or `first_order`. (second order is default)
    params : dict, optional
        Dictionary of parameters needed by sdfunc.
    realify : boolean, optional
        Force the returned results to be real.
    other : any
        Other keyword arguments available to nonlinear solvers in
        `scipy.optimize.nonlin
        <https://docs.scipy.org/doc/scipy/reference/optimize.nonlin.html>`_.
        See `Notes`.

    Returns
    -------
    t, x, e, amps, phases : array_like
        time, displacement history (time steps along columns), errors,
    amps : float array
        amplitudes of displacement (primary harmonic) in column vector format.
    phases : float array
        amplitudes of displacement (primary harmonic) in column vector format.

    Examples
    --------
    >>> import mousai as ms
    >>> t, x, e, amps, phases = ms.hb_time(ms.duff_osc,
    ...                                    np.array([[0,1,-1]]),
    ...                                    omega = 0.7)

    Notes
    -----
    .. seealso::

       ``hb_freq``

    This method is not reliable for a low number of harmonics.

    Calls a linear algebra function from
    `scipy.optimize.nonlin
    <https://docs.scipy.org/doc/scipy/reference/optimize.nonlin.html>`_ with
    `newton_krylov` as the default.

    Evaluates the differential equation/s at evenly spaced points in time. Each
    point in time yields a single equation. One harmonic plus the constant term
    results in 3 points in time over the cycle.

    Solver should gently "walk" solution up to get to nonlinearities for hard
    nonlinearities.

    Algorithm:
        1. calls `hb_err` with `x` as the variable to solve for.
        2. `hb_err` uses a Fourier representation of `x` to obtain
           velocities (after an inverse FFT) then calls `sdfunc` to determine
           accelerations.
        3. Accelerations are also obtained using a Fourier representation of x
        4. Error in the accelerations (or state derivatives) are the functional
           error used by the nonlinear algebraic solver
           (default `newton_krylov`) to be minimized by the solver.

    Options to the nonlinear solvers can be passed in by \*\*kwargs (keyward
    arguments) identical to those available to the nonlinear solver.

    """
    # Initial conditions exist?
    if x0 is None:
        if num_variables is not None:
            x0 = np.zeros((num_variables, 1 + num_harmonics * 2))
        else:
            print('Error: Must either define number of variables or initial\
                  guess for x.')
            return
    elif num_harmonics is None:
        num_harmonics = int((x0.shape[1] - 1) / 2)
    elif 1 + 2 * num_harmonics > x0.shape[1]:
        x_freq = fftp.fft(x0)
        x_zeros = np.zeros((x0.shape[0], 1 + num_harmonics * 2 - x0.shape[1]))
        x_freq = np.insert(x_freq, [x0.shape[1] - x0.shape[1] // 2], x_zeros,
                           axis=1)

        x0 = fftp.ifft(x_freq) * (1 + num_harmonics * 2) / x0.shape[1]
        x0 = np.real(x0)
    if isinstance(sdfunc, str):
        sdfunc = globals()[sdfunc]
        print("sdfunc is expected to be a function name, not a string")
    params['function'] = sdfunc  # function that returns SO derivative
    time = np.linspace(0, 2 * np.pi / omega, num=x0.shape[1], endpoint=False)
    params['time'] = time
    params['omega'] = omega
    params['n_har'] = num_harmonics

    def hb_err(x):
        r"""Array (vector) of hamonic balance second order algebraic errors.

        Given a set of second order equations
        :math:`\ddot{x} = f(x, \dot{x}, \omega, t)`
        calculate the error :math:`E = \ddot{x} - f(x, \dot{x}, \omega, t)`
        presuming that :math:`x` can be represented as a Fourier series, and
        thus :math:`\dot{x}` and :math:`\ddot{x}` can be obtained from the
        Fourier series representation of :math:`x`.

        Parameters
        ----------
        x : array_like
            x is an :math:`n \\times m` by 1 array of presumed displacements.
            It must be a "list" array (not a linear algebra vector). Here
            :math:`n` is the number of displacements and :math:`m` is the
            number of times per cycle at which the displacement is guessed
            (minimum of 3)

        **kwargs : string, float, variable
            **kwargs is a packed set of keyword arguments with 3 required
            arguments.
                1. `function`: a string name of the function which returned
                the numerically calculated acceleration.

                2. `omega`: which is the defined fundamental harmonic
                at which the is desired.

                3. `n_har`: an integer representing the number of harmonics.
                Note that `m` above is equal to 1 + 2 * `n_har`.

        Returns
        -------
        e : array_like
            2d array of numerical error of presumed solution(s) `x`.

        Notes
        -----
        `function` and `omega` are not separately defined arguments so as to
        enable algebraic solver functions to call `hb_time_err` cleanly.

        The algorithm is as follows:
            1. The velocity and accelerations are calculated in the same shape
               as `x` as `vel` and `accel`.
            3. Each column of `x` and `v` are sent with `t`, `omega`, and other
               `**kwargs** to `function` one at a time with the results
               agregated into the columns of `accel_num`.
            4. The difference between `accel_num` and `accel` is reshaped to be
               :math:`n \\times m` by 1 and returned as the vector error used
               by the numerical algebraic equation solver.

        """
        nonlocal params  # Will stay out of global/conflicts
        n_har = params['n_har']
        omega = params['omega']
        time = params['time']
        m = 1 + 2 * n_har
        vel = harmonic_deriv(omega, x)
        if eqform == 'second_order':
            accel = harmonic_deriv(omega, vel)
            accel_from_deriv = np.zeros_like(accel)

            # Should subtract in place below to save memory for large problems
            for i in np.arange(m):
                # This should enable t to be used for current time in loops
                # might be able to be commented out, left as example
                t = time[i]
                params['cur_time'] = time[i]  # loops
                # Note that everything in params can be accessed within
                # `function`.
                accel_from_deriv[:, i] = params['function'](x[:, i], vel[:, i],
                                                            params)[:, 0]
            e = accel_from_deriv - accel
        elif eqform == 'first_order':

            vel_from_deriv = np.zeros_like(vel)
            # Should subtract in place below to save memory for large problems
            for i in np.arange(m):
                # This should enable t to be used for current time in loops
                t = time[i]
                params['cur_time'] = time[i]
                # Note that everything in params can be accessed within
                # `function`.
                vel_from_deriv[:, i] =\
                    params['function'](x[:, i], params)[:, 0]

            e = vel_from_deriv - vel
        else:
            print('eqform cannot have a value of {}', eqform)
            return 0, 0, 0, 0, 0
        return e

    try:
        x = globals()[method](hb_err, x0, **kwargs)
    except:
        x = x0  # np.full([x0.shape[0],x0.shape[1]],np.nan)
        amps = np.full([x0.shape[0], ], np.nan)
        phases = np.full([x0.shape[0], ], np.nan)
        e = hb_err(x)  # np.full([x0.shape[0],x0.shape[1]],np.nan)
        raise
    else:
        xhar = fftp.fft(x) * 2 / len(time)
        amps = np.absolute(xhar[:, 1])
        phases = np.angle(xhar[:, 1])
        e = hb_err(x)

    if realify is True:
        x = np.real(x)
    else:
        print('x was real')
    return time, x, e, amps, phases


def hb_freq(sdfunc, x0=None, omega=1, method='newton_krylov', num_harmonics=1,
            num_variables=None, mask_constant=True, eqform='second_order',
            params={}, realify=True, num_time_steps=51, **kwargs):
    r"""Harmonic balance solver for first and second order ODEs.

    Obtains the solution of a first-order and second-order differential
    equation under the presumption that the solution is harmonic using an
    algebraic time method.

    Returns `t` (time), `x` (displacement), `v` (velocity), and `a`
    (acceleration) response of a first or second order linear ordinary
    differential equation defined by
    :math:`\ddot{\mathbf{x}}=f(\mathbf{x},\mathbf{v},\omega)` or
    :math:`\dot{\mathbf{x}}=f(\mathbf{x},\omega)`.

    For the state space form, the function `sdfunc` should have the form::

        def duff_osc_ss(x, params):  # params is a dictionary of parameters
            omega = params['omega']  # `omega` will be put into the dictionary
                                     # for you
            t = params['cur_time']   # The time value is available as
                                     # `cur_time` in the dictionary
            x_dot = np.array([[x[1]],
                              [-x[0]-.1*x[0]**3-.1*x[1]+1*sin(omega*t)]])
            return xdot

    In a state space form solution, the function must take the states and the
    `params` dictionary. This dictionary should be used to obtain the
    prescribed response frequency and the current time. These plus any other
    parameters are used to calculate the state derivatives which are returned
    by the function.

    For the second order form the function `sdfunc` should have the form::

        def duff_osc(x, v, params):  # params is a dictionary of parameters
            omega = params['omega']  # `omega` will be put into the dictionary
                                     # for you
            t = params['cur_time']   # The time value is available as
                                     # `cur_time` in the dictionary
            return np.array([[-x-.1*x**3-.2*v+sin(omega*t)]])

    In a second-order form solution the function must take the states and the
    `params` dictionary. This dictionary should be used to obtain the
    prescribed response frequency and the current time. These plus any other
    parameters are used to calculate the state derivatives which are returned
    by the function.

    Parameters
    ----------
    sdfunc : function
        For `eqform='first_order'`, name of function that returns **column
        vector** first derivative given `x`, `omega` and \*\*kwargs. This is
        *NOT* a string.

        :math:`\dot{\mathbf{x}}=f(\mathbf{x},\omega)`

        For `eqform='second_order'`, name of function that returns **column
        vector** second derivative given `x`, `v`, `omega` and \*\*kwargs. This
        is *NOT* a string.

        :math:`\ddot{\mathbf{x}}=f(\mathbf{x},\mathbf{v},\omega)`
    x0 : array_like, somewhat optional
        n x m array where n is the number of equations and m is the number of
        values representing the repeating solution.
        It is required that :math:`m = 1 + 2 num_{harmonics}`. (we will
        generalize allowable default values later.)
    omega : float
        assumed fundamental response frequency in radians per second.
    method : str, optional
        Name of optimization method to be used.
    num_harmonics : int, optional
        Number of harmonics to presume. The `omega` = 0 constant term is always
        presumed to exist. Minimum (and default) is 1. If num_harmonics*2+1
        exceeds the number of columns of `x0` then `x0` will be expanded, using
        Fourier analaysis, to include additional harmonics with the starting
        presumption of zero values.
    num_variables : int, somewhat optional
        Number of states for a state space model, or number of generalized
        dispacements for a second order form.
        If `x0` is defined, num_variables is inferred. An error will result if
        both `x0` and num_variables are left out of the function call.
        `num_variables` must be defined if `x0` is not.
    eqform : str, optional
        `second_order` or `first_order`. (`second order` is default)
    params : dict, optional
        Dictionary of parameters needed by sdfunc.
    realify : boolean, optional
        Force the returned results to be real.
    mask_constant : boolean, optional
        Force the constant term of the series representation to be zero.
    num_time_steps : int, default = 51
        number of time steps to use in time histories for derivative
        calculations.
    other : any
        Other keyword arguments available to nonlinear solvers in
        `scipy.optimize.nonlin
        <https://docs.scipy.org/doc/scipy/reference/optimize.nonlin.html>`_.
        See Notes.

    Returns
    -------
    t, x, e, amps, phases : array_like
        time, displacement history (time steps along columns), errors,
    amps : float array
        amplitudes of displacement (primary harmonic) in column vector format.
    phases : float array
        amplitudes of displacement (primary harmonic) in column vector format.

    Examples
    --------
    >>> import mousai as ms
    >>> t, x, e, amps, phases = ms.hb_freq(ms.duff_osc,
    ...                                    np.array([[0,1,-1]]),
    ...                                    omega = 0.7)

    Notes
    -----
    .. seealso::

       `hb_time`

    Calls a linear algebra function from
    `scipy.optimize.nonlin
    <https://docs.scipy.org/doc/scipy/reference/optimize.nonlin.html>`_ with
    `newton_krylov` as the default.

    Evaluates the differential equation/s at evenly spaced points in time
    defined by the user (default 51). Uses error in FFT of derivative
    (acceeration or state equations) calculated based on:

    1. governing equations
    2. derivative of `x` (second derivative for state method)

    Solver should gently "walk" solution up to get to nonlinearities for hard
    nonlinearities.

    Algorithm:
        1. calls `hb_time_err` with x as the variable to solve for.
        2. `hb_time_err` uses a Fourier representation of x to obtain
           velocities (after an inverse FFT) then calls `sdfunc` to determine
           accelerations.
        3. Accelerations are also obtained using a Fourier representation of x
        4. Error in the accelerations (or state derivatives) are the functional
           error used by the nonlinear algebraic solver
           (default `newton_krylov`) to be minimized by the solver.

    Options to the nonlinear solvers can be passed in by \*\*kwargs.

    """
    # Initial conditions exist?
    if x0 is None:
        if num_variables is not None:
            x0 = np.zeros((num_variables, 1 + num_harmonics * 2))
            x0 = x0 + np.random.randn(*x0.shape)
        else:
            print('Error: Must either define number of variables or initial\
                  guess for x.')
            return
    elif num_harmonics is None:
        num_harmonics = int((x0.shape[1] - 1) / 2)
    elif 1 + 2 * num_harmonics > x0.shape[1]:
        x_freq = fftp.fft(x0)
        x_zeros = np.zeros((x0.shape[0], 1 + num_harmonics * 2 - x0.shape[1]))
        x_freq = np.insert(x_freq, [x0.shape[1] - x0.shape[1] // 2], x_zeros,
                           axis=1)

        x0 = fftp.ifft(x_freq) * (1 + num_harmonics * 2) / x0.shape[1]
        x0 = np.real(x0)
    if isinstance(sdfunc, str):
        sdfunc = globals()[sdfunc]
        print("sdfunc is expected to be a function name, not a string")
    params['function'] = sdfunc  # function that returns SO derivative
    time = np.linspace(0, 2 * np.pi / omega, num=x0.shape[1], endpoint=False)
    params['time'] = time
    params['omega'] = omega
    params['n_har'] = num_harmonics

    X0 = fftp.rfft(x0)
    if mask_constant is True:
        X0 = X0[:, 1:]

    params['mask_constant'] = mask_constant

    def hb_err(X):
        """Return errors in equation eval versus derivative calculation."""
        # r"""Array (vector) of hamonic balance second order algebraic errors.
        #
        # Given a set of second order equations
        # :math:`\ddot{x} = f(x, \dot{x}, \omega, t)`
        # calculate the error :math:`E = \mathcal{F}(\ddot{x}
        # - \mathcal{F}\left(f(x, \dot{x}, \omega, t)\right)`
        # presuming that :math:`x` can be represented as a Fourier series, and
        # thus :math:`\dot{x}` and :math:`\ddot{x}` can be obtained from the
        # Fourier series representation of :math:`x` and :math:`\mathcal{F}(x)`
        # represents the Fourier series of :math:`x(t)`
        #
        # Parameters
        # ----------
        # X : float array
        #     X is an :math:`n \\times m` by 1 array of sp.fft.rfft
        #     fft coefficients lacking the constant (first) element.
        #     Here :math:`n` is the number of displacements and :math:`m` 2
        #     times the number of harmonics to be solved for.
        #
        # **kwargs : string, float, variable
        #     **kwargs is a packed set of keyword arguments with 3 required
        #     arguments.
        #         1. `function`: a string name of the function which returned
        #         the numerically calculated acceleration.
        #
        #         2. `omega`: which is the defined fundamental harmonic
        #         at which the is desired.
        #
        #         3. `n_har`: an integer representing the number of harmonics.
        #         Note that `m` above is equal to 2 * `n_har`.
        #
        # Returns
        # -------
        # e : float array
        #     2d array of numerical errors of presumed solution(s) `X`. Error
        #     between first (or second) derivative via Fourier analysis and via
        #     solution of the governing equation.
        #
        # Notes
        # -----
        # `function` and `omega` are not separately defined arguments so as to
        # enable algebraic solver functions to call `hb_err` cleanly.
        #
        # The algorithm is as follows:
        #     1. X is prepended with a zero vector (to represent the constant
        #        value)
        #     2. `x` is calculated via an inverse `numpy.fft.rfft`
        #     1. The velocity and accelerations are calculated in the same
        #        shape as `x` as `vel` and `accel`.
        #     3. Each column of `x` and `v` are sent with `t`, `omega`, and
        #        other `**kwargs** to `function` one at a time with the results
        #        agregated into the columns of `accel_num`.
        #     4. The rfft is taken of `accel_num` and `accel`.
        #     5. The first column is stripped out of both `accel_num_freq and
        #        `accel_freq`.

        # """
        nonlocal params  # Will stay out of global/conflicts
        omega = params['omega']
        time = params['time']
        mask_constant = params['mask_constant']
        if mask_constant is True:
            X = np.hstack((np.zeros_like(X[:, 0]).reshape(-1, 1), X))

        x = fftp.irfft(X)
        time_e, x = time_history(time, x, num_time_points=num_time_steps)

        vel = harmonic_deriv(omega, x)

        m = num_time_steps

        if eqform == 'second_order':
            accel = harmonic_deriv(omega, vel)
            accel_from_deriv = np.zeros_like(accel)

            # Should subtract in place below to save memory for large problems
            for i in np.arange(m):
                # This should enable t to be used for current time in loops
                # might be able to be commented out, left as example
                # t = time_e[i]
                params['cur_time'] = time_e[i]  # loops
                # Note that everything in params can be accessed within
                # `function`.
                accel_from_deriv[:, i] = params['function'](x[:, i], vel[:, i],
                                                            params)[:, 0]
            e = (accel_from_deriv - accel)  # /np.max(np.abs(accel))

            states = accel

        elif eqform == 'first_order':

            vel_from_deriv = np.zeros_like(vel)
            # Should subtract in place below to save memory for large problems
            for i in np.arange(m):
                # This should enable t to be used for current time in loops
                # t = time_e[i]
                params['cur_time'] = time_e[i]
                # Note that everything in params can be accessed within
                # `function`.
                vel_from_deriv[:, i] =\
                    params['function'](x[:, i], params)[:, 0]

            e = (vel_from_deriv - vel)  # /np.max(np.abs(vel))

            states = vel
        else:
            print('eqform cannot have a value of {}', eqform)
            return 0, 0, 0, 0, 0

        states_fft = fftp.rfft(states)

        e_fft = fftp.rfft(e)

        states_fft_condensed = condense_rfft(states_fft, num_harmonics)

        e = condense_rfft(e_fft, num_harmonics)

        if mask_constant is True:
            e = e[:, 1:]

        e = e / np.max(np.abs(states_fft_condensed))
        return e

    try:
        X = globals()[method](hb_err, X0, **kwargs)
        e = hb_err(X)
        if mask_constant is True:
            X = np.hstack((np.zeros_like(X[:, 0]).reshape(-1, 1), X))
        xhar = rfft_to_fft(X) * 2 / len(time)
        amps = np.sqrt(X[:, 1]**2+X[:, 2]**2)*2/X.shape[1]
        phases = np.arctan2(X[:, 1], -X[:, 2])
    except:  # Catches and raises errors- needs actual error listed.
        print(
            'Excepted- search failed for omega = {:6.4f} rad/s.'.format(omega))
        print("""What ever error this is, please put into har_bal
               after the excepts (2 of them)""")
        X = X0
        print(mask_constant)
        e = hb_err(X)
        if mask_constant is True:
            X = np.hstack((np.zeros_like(X[:, 0]).reshape(-1, 1), X))
        amps = np.sqrt(X[:, 1]**2+X[:, 2]**2)*2/X.shape[1]
        phases = np.arctan2(X[:, 1], -X[:, 2])

        raise

    x = fftp.irfft(X)

    if realify is True:
        x = np.real(x)
    else:
        print('x was real')
    return time, x, e, amps, phases


def hb_so(sdfunc, **kwargs):
    """Deprecated function name. Use hb_time."""
    message = 'hb_so is deprecated. Please use hb_time or an alternative.'
    warnings.warn(message, DeprecationWarning)
    return hb_time(sdfunc, kwargs)


def harmonic_deriv(omega, r):
    r"""Return derivative of a harmonic function using frequency methods.

    Parameters
    ----------
    omega: float
        Fundamendal frequency, in rad/sec, of repeating signal
    r: array_like
        | Array of rows of time histories to take the derivative of.
        | The 1 axis (each row) corresponds to a time history.
        | The length of the time histories *must be an odd integer*.

    Returns
    -------
    s: array_like
        Function derivatives.
        The 1 axis (each row) corresponds to a time history.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from mousai import *
    >>> import scipy as sp
    >>> from scipy import pi, sin, cos
    >>> f = 2
    >>> omega = 2.*pi * f
    >>> numsteps = 11
    >>> t = sp.arange(0,1/omega*2*pi,1/omega*2*pi/numsteps)
    >>> x = sp.array([sin(omega*t)])
    >>> v = sp.array([omega*cos(omega*t)])
    >>> states = sp.append(x,v,axis = 0)
    >>> state_derives = harmonic_deriv(omega,states)
    >>> plt.plot(t,states.T,t,state_derives.T,'x')
    [<matplotlib.line...]

    """
    s = np.zeros_like(r)
    for i in np.arange(r.shape[0]):
        s[i, :] = fftp.diff(r[i, :]) * omega
    return np.real(s)


def solmf(x, v, M, C, K, F):
    r"""Return acceleration of second order linear matrix system.

    Parameters
    ----------
    x, v, F : array_like
        :math:`n\times 1` arrays of current displacement, velocity, and Force.
    M, C, K : array_like
        Mass, damping, and stiffness matrices.

    Returns
    -------
    a : array_like
        :math:`n\\times 1` acceleration vector

    Examples
    --------
    >>> import numpy as np
    >>> M = np.array([[2,0],[0,1]])
    >>> K = np.array([[2,-1],[-1,3]])
    >>> C = 0.01 * M + 0.01 * K
    >>> x = np.array([[1],[0]])
    >>> v = np.array([[0],[10]])
    >>> F = v * 0.1
    >>> a = solmf(x, v, M, C, K, F)
    >>> print(a)
        [[-0.95]
         [ 1.6 ]]

    """
    return -la.solve(M, C @ v + K @ x - F)


def duff_osc(x, v, params):
    """Duffing oscillator acceleration."""
    omega = params['omega']
    t = params['cur_time']
    acceleration = np.array([[-x - .1 * x**3. - 0.2 * v + np.sin(omega * t)]])
    return acceleration


def time_history(t, x, num_time_points=200, realify=True):
    r"""Generate refined time history from harmonic balance solution.

    Harmonic balance solutions presume a limited number of harmonics in the
    solution. The result is that the time history is usually a very limited
    number of values. Plotting these results implies that the solution isn't
    actually a continuous one. This function fills in the gaps using the
    harmonics obtained in the solution.

    Parameters
    ----------
    t: array_like
        1 x m array where m is the number of
        values representing the repeating solution.
    x: array_like
        n x m array where m is the number of equations and m is the number of
        values representing the repeating solution.
    realify: boolean
        Force the returned results to be real.
    num_time_points: int
        number of points desired in the "smooth" time history.

    Returns
    -------
    t: array_like
        1 x num_time_points array.
    x: array_like
        n x num_time_points array.

    Examples
    --------
    >>> import numpy as np
    >>> import mousai as ms
    >>> x = np.array([[-0.34996499,  1.36053998, -1.11828552]])
    >>> t = np.array([0.        , 2.991993  , 5.98398601])
    >>> t_full, x_full = ms.time_history(t, x, num_time_points=300)

    Notes
    -----
    The implication of this function is that the higher harmonics that
    were not determined in the solution are zero. This is indeed the assumption
    made when setting up the harmonic balance solution. Whether this is a valid
    assumption is something that the user must judge when obtaining the
    solution.

    """
    dt = t[1]
    t_length = t.size
    t = np.linspace(0, t_length * dt, num_time_points, endpoint=False)
    x_freq = fftp.fft(x)
    x_zeros = np.zeros((x.shape[0], t.size - x.shape[1]))
    x_freq = np.insert(x_freq, [t_length - t_length // 2], x_zeros, axis=1)

    x = fftp.ifft(x_freq) * num_time_points / t_length
    if realify is True:
        x = np.real(x)
    else:
        print('x was real')
    return t, x


def condense_fft(X_full, num_harmonics):
    """Create equivalent amplitude reduced-size FFT from longer FFT."""
    X_red = (np.hstack((X_full[:, 0:(num_harmonics + 1)],
                        X_full[:, -1:-(num_harmonics + 1):-1]))
             * (2 * num_harmonics + 1) / X_full[0, :].size)
    return X_red


def condense_rfft(X_full, num_harmonics):
    """Return real fft with fewer harmonics."""
    X_len = X_full.shape[1]
    X_red = X_full[:, :(num_harmonics) * 2 + 1] / \
        X_len * (1 + 2 * num_harmonics)
    return X_red


def expand_rfft(X, num_harmonics):
    """Return real fft with mor harmonics."""
    X_len = X.shape[1]
    cur_num_harmonics = (X_len - 1) / 2
    X_expanded = np.hstack((X / X_len * (1 + 2 * num_harmonics),
                            np.zeros((X.shape[0],
                                      int(2 * (num_harmonics
                                               - cur_num_harmonics))))
                            ))
    return X_expanded


def rfft_to_fft(X_real):
    """Switch from SciPy real fft form to complex fft form."""
    X = fftp.fft(fftp.irfft(X_real))
    return X


def fft_to_rfft(X):
    """Switch from complex form fft form to SciPy rfft form."""
    X_real = fftp.rfft(np.real(fftp.ifft(X)))
    return X_real


def time_history_r(t, x, num_time_points=200, realify=True):
    r"""Generate refined time history from harmonic balance solution.

    Harmonic balance solutions presume a limited number of harmonics in the
    solution. The result is that the time history is usually a very limited
    number of values. Plotting these results implies that the solution isn't
    actually a continuous one. This function fills in the gaps using the
    harmonics obtained in the solution.

    Parameters
    ----------
    t: array_like
        1 x m array where m is the number of
        values representing the repeating solution.
    x: array_like
        n x m array where m is the number of equations and m is the number of
        values representing the repeating solution.
    realify: boolean
        Force the returned results to be real.
    num_time_points: int
        number of points desired in the "smooth" time history.

    Returns
    -------
    t: array_like
        1 x num_time_points array.
    x: array_like
        n x num_time_points array.

    Examples
    --------
    >>> import numpy as np
    >>> import mousai as ms
    >>> x = np.array([[-0.34996499,  1.36053998, -1.11828552]])
    >>> t = np.array([0.        , 2.991993  , 5.98398601])
    >>> t_full, x_full = ms.time_history(t, x, num_time_points=300)

    Notes
    -----
    The implication of this function is that the higher harmonics that
    were not determined in the solution are zero. This is indeed the assumption
    made when setting up the harmonic balance solution. Whether this is a valid
    assumption is something that the user must judge when obtaining the
    solution.

    """
    dt = t[1]
    t_length = t.size
    t = sp.linspace(0, t_length * dt, num_time_points, endpoint=False)
    x_freq = fftp.fft(x)
    x_zeros = sp.zeros((x.shape[0], t.size - x.shape[1]))
    x_freq = sp.insert(x_freq, [t_length - t_length // 2], x_zeros, axis=1)
    # print(x_freq)
    # x_freq = np.hstack((x_freq, x_zeros))
    # print(x_freq)
    x = fftp.ifft(x_freq) * num_time_points / t_length
    if realify is True:
        x = np.real(x)
    else:
        print('x was real')
    return t, x
