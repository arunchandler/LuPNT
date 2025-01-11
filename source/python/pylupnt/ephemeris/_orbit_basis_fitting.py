try:
    from ._basis import *
except ImportError:
    from _basis import *
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def cheby_interp_points(n):
    """
    Generate Chebyshev interpolation points
    """
    # x_j = [np.cos( np.pi * (j+0.5)/(n+1) ) for j in range(n+1)]
    x_j = [np.cos(np.pi * j / (n)) for j in range(n + 1)]
    # make sure the output is -1 to 1
    return x_j[::-1]


def fit_basis_function(t_array,n_array,f_true_val,fdot_true_val,t_interval, basis, # Chebyshev, Polynomial
        solver="ECOS", add_sise_constraint=True, res_pos=13.43 / 4, res_vel=1.2 / 4, alpha=None, 
        fit_velocity=True, cheby_interp=True, end_points=True, tol=1e-7):
    """
    Fit the basis functions to the given data with the sise constraint

    Args:
        t_array (np.array): time vector
        n_array (np.array): array of polynomial orders
        f_true_val (np.array): true function values
        fdot_true_val (np.array): true function derivative values
        t_interval (np.array): time interval
        basis (Basis): basis functions
        solver (str): solver for the convex optimization problem
        add_sise_constraint (bool): add the sise constraint
        res_pos (float): position error requirement
        res_vel (float): velocity error requirement
        alpha (np.array): array of L-norm fitting coefficients
        cheby_interp (bool): use Chebyshev interpolation points
        end_points (bool): add end points constraints
        tol (float): tolerance for the solver
    """
    lent = t_array.size
    truncation_err = np.zeros((2, len(n_array)))
    lenn = len(n_array)

    status_list = []

    if alpha is None:
        alpha = np.zeros(lenn)

    # estimated values
    f_est_val = np.zeros((lenn, lent))
    fdot_est_val = np.zeros((lenn, lent))

    # coefficients
    coeff_list = np.zeros((lenn, np.max(n_array) + 1))

    # construct Cubic Spline Interporation (used to compute interporaltion point values)
    cs_f = CubicSpline(t_array, f_true_val)
    cs_fdot = CubicSpline(t_array, fdot_true_val)

    for k in range(len(n_array)):
        n = n_array[k]
        # print("n = {}".format(n))
        # interpolation points
        if cheby_interp:
            t_interpolation = cheby_interp_points(n)
        else:
            t_interpolation = np.linspace(-1, 1, n + 1)

        # true values interpolated
        f_true_interp = cs_f(t_interpolation)
        f_dot_true_interp = cs_fdot(t_interpolation)

        # System of equations Ta = x
        t_matrix = np.ones((2 * (n + 1), n + 1))
        x_matrix = np.ones((2 * (n + 1)))

        # build chebyshev polynomial matrix
        for j in range(n + 1):  # row
            for i in range(n + 1):  # column
                t_matrix[j, i] = basis.generate_pol(i, t_interpolation[j])
                t_matrix[j + n + 1, i] = basis.generate_der_pol(
                    i, t_interpolation[j], t_interval
                )
            x_matrix[j] = f_true_interp[j]
            x_matrix[j + n + 1] = f_dot_true_interp[j]

        # solve convex optimization problem ------------------------------------------
        ax = cp.Variable(n + 1)
        if fit_velocity:    
            cost = cp.sum_squares(t_matrix @ ax - x_matrix) + alpha[k] * cp.norm1(ax)
        else:
            cost = cp.sum_squares(t_matrix[0 : n + 1, :] @ ax - x_matrix[0 : n + 1]) + alpha[k] * cp.norm1(ax)

        # set up constraints
        const = []
        if add_sise_constraint:
            const.append(
                cp.norm_inf(t_matrix[0 : n + 1, :] @ ax - x_matrix[0 : n + 1]) <= res_pos
            )  # position error
            if fit_velocity:
                const.append(
                    cp.norm_inf(t_matrix[n + 1 : 2 * (n + 1), :] @ ax - x_matrix[n + 1 :])
                    <= res_vel
                )  # velocity error

        if end_points:
            const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
            const.append(t_matrix[n, :] @ ax == x_matrix[n])  # final point
            if fit_velocity:
                const.append(t_matrix[n + 1, :] @ ax == x_matrix[n + 1])  # initial velocity
                const.append(t_matrix[-1, :] @ ax == x_matrix[-1])  # final velocity

        prob = cp.Problem(cp.Minimize(cost), const)

        try:
            prob.solve(solver=solver, feastol=tol, max_iters=100000)
            # print("status:", prob.status)
            # print("optimal value", prob.value)
            a_x = ax.value

            status_list.append(prob.status)

            if prob.status == "optimal":
                # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                # calculate function approximations, maintain the same coefficients
                f_est_val[k, :] = basis.get_f_approx(a_x, n, t_array)
                fdot_est_val[k, :] = basis.get_f_approx_der(a_x, n, t_array, t_interval)

                # error as the L2 norm
                err = f_true_val - f_est_val[k, :]
                err_dot = fdot_true_val - fdot_est_val[k, :]

                truncation_err[0, k] = np.linalg.norm(err)
                truncation_err[1, k] = np.linalg.norm(err_dot)

                coeff_list[k, 0 : n + 1] = a_x

            if prob.status == "optimal_inaccurate":
                # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                # calculate function approximations, maintain the same coefficients
                f_est_val[k, :] = basis.get_f_approx(a_x, n, t_array)
                fdot_est_val[k, :] = basis.get_f_approx_der(a_x, n, t_array, t_interval)

                # error as the L2 norm
                err = f_true_val - f_est_val[k, :]
                err_dot = fdot_true_val - fdot_est_val[k, :]

                truncation_err[0, k] = np.linalg.norm(err)
                truncation_err[1, k] = np.linalg.norm(err_dot)

                coeff_list[k, 0 : n + 1] = a_x

        except cp.error.SolverError:
            f_est_val[k, :] = np.zeros(len(t_array))
            fdot_est_val[k, :] = np.zeros(len(t_array))

            # error as the L2 norm
            err = f_true_val - f_est_val[k, :]
            err_dot = fdot_true_val - fdot_est_val[k, :]

            truncation_err[0, k] = np.linalg.norm(err)
            truncation_err[1, k] = np.linalg.norm(err_dot)

            coeff_list[k, 0 : n + 1] = np.zeros(n + 1)

            status_list.append("infeasible")

    return coeff_list, truncation_err, f_est_val, fdot_est_val, status_list


def fit_basis_orbit(df_interp, t_interval, n_array, basis, add_sise_constraint=True, 
                    fit_velocity=True, cheby_interp=True, plot=False):
    """
    Fit the basis functions to the given data with the sise constraint

    """
      # evaluation points
    t_array = np.array(df_interp["time"])
    # specify alpha array if L-norm fitting
    alpha = None
    solver = "ECOS"

    # requirements (scaled since this is a single variable fitting)
    pos_req = 0.01343 / 6
    vel_req = (1.2e-6) / 6

    # fitting
    pos_labels = ["x", "y", "z"]
    vel_labels = ["v_x", "v_y", "v_z"]
    for i in range(3):  # x, y, z
        f_true_val = np.array(df_interp[pos_labels[i]])
        fdot_true_val = np.array(df_interp[vel_labels[i]])

        a_list, truncation_err, f_est_val, fdot_est_val, status = (
            fit_basis_function(t_array, n_array, f_true_val, fdot_true_val,
                t_interval, basis, solver="ECOS", add_sise_constraint=add_sise_constraint,
                res_pos=pos_req, res_vel=vel_req, alpha=None, fit_velocity=fit_velocity, 
                cheby_interp=cheby_interp, end_points=False, tol=1e-8,
            )
        )

        if i == 0:
            ax_list = a_list
            fx_est_val = f_est_val
            fdotx_est_val = fdot_est_val
            x_status = status
        elif i == 1:
            ay_list = a_list
            fy_est_val = f_est_val
            fdoty_est_val = fdot_est_val
            y_status = status
        else:
            az_list = a_list
            fz_est_val = f_est_val
            fdotz_est_val = fdot_est_val
            z_status = status

    infeas_status = []
    for s in range(len(x_status)):
        if (
            x_status[s] == "infeasible"
            or y_status[s] == "infeasible"
            or z_status[s] == "infeasible"
        ):
            infeas_status.append(False)
        else:
            infeas_status.append(True)

    return (
        ax_list,
        ay_list,
        az_list,
        fx_est_val,
        fy_est_val,
        fz_est_val,
        fdotx_est_val,
        fdoty_est_val,
        fdotz_est_val,
        infeas_status,
    )
