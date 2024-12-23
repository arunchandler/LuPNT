import numpy as np
import cvxpy as cp
from scipy.interpolate import CubicSpline


# linear fit
def linear_angle_fit(t_array, df_true_val, n=50, solver="ECOS"):
    # get the three angles
    psi_ang = CubicSpline(t_array, df_true_val["psi"])
    theta_ang = CubicSpline(t_array, df_true_val["theta"])
    phi_ang = CubicSpline(t_array, df_true_val["phi"])
    # interpolation points
    t_interpolation = np.linspace(t_array[0], t_array[-1], n + 1)
    # true values interpolated
    psi_interp = psi_ang(t_interpolation)
    theta_interp = theta_ang(t_interpolation)
    phi_interp = phi_ang(t_interpolation)
    # System of equations Ta = x
    B_matrix = np.zeros((3 * n, 6))
    y_matrix = np.zeros((3 * n, 1))
    # build design matrix
    for j in range(n):  # row
        B_matrix[j, 0] = 1
        B_matrix[j, 1] = t_interpolation[j]
        B_matrix[j + n, 2] = 1
        B_matrix[j + n, 3] = t_interpolation[j]
        B_matrix[j + 2 * n, 4] = 1
        B_matrix[j + 2 * n, 5] = t_interpolation[j]

        y_matrix[j] = psi_interp[j]
        y_matrix[j + n] = theta_interp[j]
        y_matrix[j + 2 * n] = phi_interp[j]

    # solve optimization problem
    a_vec = cp.Variable([6, 1])

    cost = cp.sum_squares(B_matrix @ a_vec - y_matrix)
    # set up constraints
    const = []
    const.append(B_matrix[0, :] @ a_vec == y_matrix[0])  # initial point
    const.append(B_matrix[n, :] @ a_vec == y_matrix[n])  # initial point
    const.append(B_matrix[2 * n, :] @ a_vec == y_matrix[2 * n])  # initial point

    prob = cp.Problem(cp.Minimize(cost), const)
    prob.solve()
    a_soln = a_vec.value

    if prob.status == "optimal":
        return a_soln.flatten()
    else:
        return a_soln.flatten()


# quadratic fit
def quad_angle_fit(t_array, df_true_val, n=50, solver="ECOS"):
    # get the three angles
    psi_ang = CubicSpline(t_array, df_true_val["psi"])
    theta_ang = CubicSpline(t_array, df_true_val["theta"])
    phi_ang = CubicSpline(t_array, df_true_val["phi"])
    # interpolation points
    t_interpolation = np.linspace(t_array[0], t_array[-1], n + 1)
    # true values interpolated
    psi_interp = psi_ang(t_interpolation)
    theta_interp = theta_ang(t_interpolation)
    phi_interp = phi_ang(t_interpolation)
    # System of equations Ta = x
    B_matrix = np.zeros((3 * n, 9))
    y_matrix = np.zeros((3 * n, 1))
    # build design matrix
    for j in range(n):  # row
        B_matrix[j, 0] = 1
        B_matrix[j, 1] = t_interpolation[j]
        B_matrix[j, 2] = t_interpolation[j] ** 2
        B_matrix[j + n, 3] = 1
        B_matrix[j + n, 4] = t_interpolation[j]
        B_matrix[j, 5] = t_interpolation[j] ** 2
        B_matrix[j + 2 * n, 6] = 1
        B_matrix[j + 2 * n, 7] = t_interpolation[j]
        B_matrix[j, 8] = t_interpolation[j] ** 2

        y_matrix[j] = psi_interp[j]
        y_matrix[j + n] = theta_interp[j]
        y_matrix[j + 2 * n] = phi_interp[j]

    # solve optimization problem
    a_vec = cp.Variable([9, 1])

    cost = cp.sum_squares(B_matrix @ a_vec - y_matrix)
    # set up constraints
    const = []
    const.append(B_matrix[0, :] @ a_vec == y_matrix[0])  # initial point
    const.append(B_matrix[n, :] @ a_vec == y_matrix[n])  # initial point
    const.append(B_matrix[2 * n, :] @ a_vec == y_matrix[2 * n])  # initial point

    prob = cp.Problem(cp.Minimize(cost), const)
    prob.solve()
    a_soln = a_vec.value

    if prob.status == "optimal":
        return a_soln.flatten()
    else:
        return a_soln.flatten()


def angles_from_coeffs(time, angles_coeffs):
    if len(angles_coeffs) == 6:
        # retrieve angles
        psi_dot = angles_coeffs[1]
        psi = angles_coeffs[0] + psi_dot * time

        theta_dot = angles_coeffs[3]
        theta = angles_coeffs[2] + theta_dot * time

        phi_dot = angles_coeffs[5]
        phi = angles_coeffs[4] + phi_dot * time

    elif len(angles_coeffs) == 9:
        # retrieve angles
        psi_ddot = angles_coeffs[2]
        psi_dot = angles_coeffs[1] + psi_ddot * time
        psi = angles_coeffs[0] + psi_dot * time + psi_ddot * (time**2)

        theta_ddot = angles_coeffs[5]
        theta_dot = angles_coeffs[4] + theta_ddot * time
        theta = angles_coeffs[3] + theta_dot * time + theta_ddot * (time**2)

        phi_ddot = angles_coeffs[8]
        phi_dot = angles_coeffs[7] + phi_ddot * time
        phi = angles_coeffs[6] + phi_dot * time + phi_ddot * (time**2)

    angles = [phi, theta, psi, phi_dot, theta_dot, psi_dot]

    return angles


def rotX(theta):
    return np.array(
        [
            [1, 0, 0],
            [0, np.cos(theta), np.sin(theta)],
            [0, -np.sin(theta), np.cos(theta)],
        ]
    )


def rotY(theta):
    return np.array(
        [
            [np.cos(theta), 0, -np.sin(theta)],
            [0, 1, 0],
            [np.sin(theta), 0, np.cos(theta)],
        ]
    )


def rotZ(theta):
    return np.array(
        [
            [np.cos(theta), np.sin(theta), 0],
            [-np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
    )


def rot_mat_mci2pa(angles):
    """
    Returns the rotation matrix from MCI to PA frame
    """
    phi = angles[0]
    theta = angles[1]
    psi = angles[2]
    phi_dot = angles[3]
    theta_dot = angles[4]
    psi_dot = angles[5]

    cpsi = np.cos(psi)
    spsi = np.sin(psi)

    mat = np.array(
        [
            [-psi_dot * spsi, psi_dot * cpsi, 0],
            [-psi_dot * cpsi, -psi_dot * spsi, 0],
            [0, 0, 0],
        ]
    )
    R_mi2pa = rotZ(psi) @ rotX(theta) @ rotZ(phi)
    R_mi2pa_dot = mat @ rotX(theta) @ rotZ(phi)

    return R_mi2pa, R_mi2pa_dot
