import numpy as np
from scipy.optimize import minimize
import pylupnt as pnt


def qllo_opt(a0, e0):
    """
    Optimize the initial conditions for a quasi-LLO.

    Parameters
    ----------
    a0 : float
        Initial semi-major axis [km].
    e0 : float
        Initial eccentricity [-].

    Returns
    -------
    oe : numpy.ndarray
        Optimized orbital elements [a, e, i, Omega, w, M] in op frame.
    """

    # Initial guess and bounds
    x0 = np.array([a0, e0, np.pi / 2, 0, np.pi / 2])
    bounds = [
        (0, None),
        (0, 0.6),
        (85 * np.pi / 180, 95 * np.pi / 180),
        (0, 2 * np.pi),
        (80 * np.pi / 180, 100 * np.pi / 180),
    ]

    # Optimization
    result = minimize(
        objfunc,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints={"type": "ineq", "fun": nonlcon},
    )

    if not result.success:
        print("Warning: Optimization did not converge!")

    x = result.x

    # Output dictionary
    a = x[0]
    e = x[1]
    i = x[2]
    Omega = x[3]
    w = x[4]
    M = 0.0

    oe = np.array([a, e, i, Omega, w, M])

    return oe


def objfunc(x):
    # Objective function to minimize
    dot = EMOs(x)
    return dot[0] ** 2 + dot[1] ** 2 + dot[2] ** 2


def nonlcon(x):
    # Nonlinear equality constraint
    a = x[0]
    e = x[1]
    return [a * (1 - e) - 1738]  # Constrains the periapsis distance


def EMOs(x):
    # Parameters from Singh et al. (2020)
    k = 0.98784941573006
    J2 = 2.0323e-4
    n = 2.64907232701554e-06
    mu = 4902.8000661638
    RM = 1738

    # Orbital elements
    a = x[0]
    e = x[1]
    i = x[2]
    Omega = x[3]
    omega = x[4]

    # Equations for e_dot, i_dot, and omega_dot
    e_dot = (
        (15 * k * n**2 * a ** (3 / 2) * e * np.sqrt(1 - e**2) / (8 * np.sqrt(mu)))
        * np.sin(i) ** 2
        * np.sin(2 * omega)
    )

    i_dot = (
        -(15 * k * n**2 * a ** (3 / 2) * e**2 / (16 * np.sqrt(mu) * np.sqrt(1 - e**2)))
        * np.sin(2 * i)
        * np.sin(2 * omega)
    )

    omega_dot = (3 * J2 * np.sqrt(mu) * RM**2) / (
        4 * a ** (7 / 2) * (1 - e**2) ** 2
    ) * (5 * np.cos(i) ** 2 - 1) + (3 * k * n**2 * a ** (3 / 2)) / (
        8 * np.sqrt(mu) * np.sqrt(1 - e**2)
    ) * (
        (5 * np.cos(i) ** 2 - 1 + e**2)
        + 5 * (1 - e**2 - np.cos(i) ** 2) * np.cos(2 * omega)
    )

    return np.array([e_dot, i_dot, omega_dot])
