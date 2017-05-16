import scipy as sp
import scipy.fftpack as fftp
import scipy.linalg as la
# from scipy import pi, sin,cos
# import matplotlib.pyplot as plt

__all__ = ["mousai_so",
           "harmonic_deriv",
           "somf"]


def mousai_so(sdfunc, x0, omega, method, *kwargs, num_harmonics=1):
    r"""Harmonic balance solver for second order ODEs.

    Obtains the solution of a second order differential equation under the
    presumption that the solution is harmonic.

    Returns t (time), x (displacement), v (velocity), and a (acceletation)
    response of a second order linear ordinary differential
    equation defined by
    :math:`\ddot{\mathbf{x}}=f(\mathbf{x},\mathbf{v},\omega)`.

    Parameters
    ----------
    sdfunc: str
        name of function that returns second derivative given omega and
        \*kwargs
        :math:`\ddot{\mathbf{x}}=f(\mathbf{x},\mathbf{v},\omega)`
    omega:  float
        assumed fundamental response frequency.
    num_harmonics: int
        number of harmonics to presume. Constant term is always presumed.
    x0: ndarray
        n x m array where n is the number of equations and m is the number of
        values representing the repeating solution.
        It is required that :math:`m = 1 + 2  num_{harmonics}`.
    method: str
        Name of optimization method to be used.

    Returns
    -------
    t, x, v, a : ndarrays
    """

    '''
    a) define a function that returns errors in time domain as vector
    b) define function to obtain velocity and accelerations from displacement
    and frequencies.
    c)
    '''

    return


def harmonic_deriv(omega, r):
    """Derivative of a harmonic function using frequency methods.

    Returns the derivatives of a harmonic function

    Parameters
    ----------
    omega: float
        Fundamendal frequency, in rad/sec, of repeating signal
    r: array
        | Array of rows of time histories to take the derivative of.
        | The 1 axis (each row) corresponds to a time history.
        | The length of the time histories *must be an odd integer*.

    Returns
    -------
    s: array
        array of function derivatives.
        The 1 axis (each row) corresponds to a time history.

    Notes
    -----
    At this time, the length of the time histories *must be an odd integer*.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from har_bal import *
    >>> import scipy as sp
    >>> from scipy import pi,sin,cos
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

    n = r.shape[1]
    omega_half = -sp.arange((n-1)/2+1) * omega * 2j/(n-2)
    omega_whole = sp.append(sp.conj(omega_half[-1:0:-1]), omega_half)
    r_freq = fftp.fft(r)
    s_freq = r_freq * omega_whole
    s = fftp.ifft(s_freq)
    return sp.real(s)


def somf(x, v, M, C, K, F):
    """Acceleration of second order linear matrix system.

    Parameters
    ----------
    x, v, F : arrays
        :math:`n\\times 1` arrays of current displacement, velocity, and Force.
    M, C, K : arrays
        Mass, damping, and stiffness matrices.

    Returns
    -------
    a : array
        :math:`n\\times 1` acceleration vector

    Examples
    --------
    >>> import scipy as sp
    >>> M = sp.array([[2,0],[0,1]])
    >>> K = sp.array([[2,-1],[-1,3]])
    >>> C = 0.01 * M + 0.01 * K
    >>> x = sp.array([[1],[0]])
    >>> v = sp.array([[0],[10]])
    >>> F = v * 0.1
    >>> a = somf(x, v, M, C, K, F)
    >>> print(a)
        [[-0.95]
         [ 1.6 ]]
    """

    return -la.solve(M, C @ v + K @ x - F)


def hb_so_err(x, **kwargs):
    """Array (vector) of hamonic balance second order algebraic errors.

    Given a set of second order equations
    :math:`\ddot{x} = f(x, \dot{x}, \omega, t)`
    calculate the error :math:`E = \ddot{x} - f(x, \dot{x}, \omega, t)`
    presuming that :math:`x` can be represented as a Fourier series, and thus
    :math:`\dot{x}` and :math:`\ddot{x}` can be obtained from the Fourier
    series representation of :math:`x`.

    Parameters
    ----------
    x : array
        x is an :math:`n \\times m` by 1 array of presumed displacements. It
        must be a "list" array (not a linear algebra vector). Here
        :math:`n` is the number of displacements and :math:`m` is the number of
        times per cycle at which the displacement is guessed (minimum of 3)

    **kwargs : string, float, variable
        **kwargs is a packed set of keyword arguments with 3 required
        arguments.
            1. `function`: a string name of the
            function which returned the numerically calculated acceleration.

            2. `omega`: which is the defined fundamental harmonic
            at which the is desired.

            3. `n_har`: an integer representing the number of harmonics. Note
            that `m` above is equal to 1 + 2 * `n_har`.

    Returns
    -------
    e : array
        2d numpy vector array of numerical error of presumed solution(s) `x`

    Notes
    -----
    `function` and `omega` are not separately defined arguments so as to enable
    algebraic solver functions to call `hb_so_err` cleanly.

    The algorithm is as follows:
        1. The vector `x` is reshaped into an :math:`n` by :math:`m` array
        2. The velocity and accelerations are calculated in the same shape as
           `x` as `vel` and `accel`.
        3. Each column of `x` and `v` are sent with `t`, `omega`, and other
           `**kwargs** are sent to `function` one at a time with the results
           agregated into the columns of `accel_num`.
        4. The difference between `accel_num` and `accel` is reshaped to be
           :math:`n \\times m` by 1 and returned as the vector error used by
           the numerical algebraic equation solver.
    """

    n_har = kwargs['n_har']
    omega = kwargs['omega']
    function = kwargs['function']
    m = 1 + 2 * n_har
    vel = harmonic_deriv(omega, x)
    accel = harmonic_deriv(omega, vel)
    accel_num = sp.zeros_like(accel)

    for i in sp.arange(m):
        accel_num[:, i] = globals()[function](x[:, i], vel[:, i], kwargs)

    e = accel_num - accel

    return e

"""
optimizer will actually be solver from scipy.optimize
from scipy.optimize
Say:
newton_krylov
broyden1
https://docs.scipy.org/doc/scipy/reference/optimize.nonlin.html

"""


if __name__ == "__main__":
    """Run doctests.

    python (name of this file)  -v
    will test all of the examples in the help.

    Leaving off -v will run the tests without any output. Success will return
    nothing.

    See the doctest section of the python manual.
    https://docs.python.org/3.5/library/doctest.html
    """

    # import mousai as ms

    # doctest.run_docstring_examples(frfest,globals(),optionflags=doctest.ELLIPSIS)
    # doctest.run_docstring_examples(asd,globals(),optionflags=doctest.ELLIPSIS)
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS |
                    doctest.NORMALIZE_WHITESPACE)
