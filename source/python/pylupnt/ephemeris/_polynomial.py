import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd


def generate_polynomial_basis(k, t):
    """Generates k-th order polynomial evaluation at time t

    Args:
        k (int): polynomial order
        t (float): time value

    Returns:
        float: component of order k at time t
    """
    return t ** (k)


def generate_der_polynomial_basis(k, t, t_interval):
    """Generates k-th order polynomial derivative evaluation at time t

    Args:
        k (int): polynomial order
        t (float): time value
        t_interval (float): total time interval

    Returns:
        float: derivative component of order k at time t
    """
    if k == 0:
        return 0
    else:
        return k * (t ** (k - 1)) * (2 / t_interval)


def get_f_approx_poly(a, k, t_array):
    """Builds approximation

    Args:
        a (float): coefficients obtained from fitting
        k (int): polynomial order
        t_array (numpy.array): time span

    Returns:
        list: approximation at each time step
    """
    f_approx = np.zeros(len(t_array))

    for t in range(len(t_array)):
        for k_val in range(k + 1):
            f_approx[t] += a[k_val] * generate_polynomial_basis(k_val, t_array[t])
    return f_approx


def get_fdot_approx_poly(a, k, t_array, t_interval):
    """Builds approximation of derivative

    Args:
        a (float): coefficients obtained from fitting
        k (int): polynomial order
        t_array (numpy.array): time span
        t_interval (float): total time interval

    Returns:
        list: approximation at each time step
    """
    f_approx = np.zeros(len(t_array))
    for t in range(len(t_array)):
        for k_val in range(k + 1):
            f_approx[t] += a[k_val] * generate_der_polynomial_basis(
                k_val, t_array[t], t_interval
            )
    return f_approx


def cheby_interp_points(n):
    # x_j = [np.cos( np.pi * (j+0.5)/(n+1) ) for j in range(n+1)]
    x_j = [np.cos(np.pi * j / (n)) for j in range(n + 1)]
    # make sure the output is -1 to 1
    return x_j[::-1]


##FITTING FUNCTIONS
# Original fitting function
def fit_polynomial(
    t_array,
    k_array,
    f_true_val,
    fdot_true_val,
    t_interval,
    solver="ECOS",
    alpha=None,
    weight=1,
    f_normalize=1,
    fdot_normalize=1,
    cheby_interp=True,
    more_points=True,
):

    lent = t_array.size
    truncation_err = np.zeros((2, len(k_array)))
    lenn = len(k_array)

    if alpha is None:
        alpha = np.zeros(lenn)

    # estimated values
    f_est_val = np.zeros((lenn, lent))
    fdot_est_val = np.zeros((lenn, lent))

    # construct Cubic Spline Interporation (used to compute interporaltion point values)
    cs_f = CubicSpline(t_array, f_true_val)
    cs_fdot = CubicSpline(t_array, fdot_true_val)

    max_err = np.zeros(len(alpha))
    max_err_dot = np.zeros(len(alpha))

    for p in range(len(k_array)):
        k = k_array[p]
        print("Polynomial order (k) = {}".format(k))

        if more_points:
            num_points = 4 * k
        else:
            num_points = k

        # interpolation points (there are k+1 points just like the number of coefficient to solve for)
        if cheby_interp:
            t_interpolation = cheby_interp_points(num_points)
            # t_interpolation = [10*(1 + t)/2 for t in t_interpolation]
        else:
            t_interpolation = np.linspace(t_array[0], t_array[-1], num_points + 1)

        # print(t_interpolation)

        # true values interpolated
        f_true_interp = cs_f(t_interpolation)
        f_dot_true_interp = cs_fdot(t_interpolation)

        # System of equations Xa = y
        # X is a (2n x (k+1)) matrix, where n is the number of interpolation points and k is the polynomial order
        # the current system makes it so that there are as many interpolation points as the number of coefficients to solve for (k+1)

        X_matrix = np.ones((2 * (num_points + 1), k + 1))
        y_matrix = np.ones((2 * (num_points + 1)))

        # build matrices
        for j in range(num_points + 1):  # row
            for i in range(k + 1):  # column
                X_matrix[j, i] = generate_polynomial_basis(i, t_interpolation[j])
                X_matrix[j + num_points + 1, i] = generate_der_polynomial_basis(
                    i, t_interpolation[j], t_interval
                )

            y_matrix[j] = f_true_interp[j]
            y_matrix[j + num_points + 1] = f_dot_true_interp[j]

        # print(X_matrix)
        # print(y_matrix)

        # solve convex optimization problem
        ax = cp.Variable(k + 1)
        cost = cp.sum_squares(X_matrix @ ax - y_matrix) + alpha[p] * cp.norm1(ax)
        const = []
        const.append(X_matrix[0, :] @ ax == y_matrix[0])  # initial point
        const.append(
            X_matrix[num_points, :] @ ax == y_matrix[num_points]
        )  # final point
        const.append(
            X_matrix[num_points + 1, :] @ ax == y_matrix[num_points + 1]
        )  # initial velocity
        const.append(X_matrix[-1, :] @ ax == y_matrix[-1])  # final velocity
        prob = cp.Problem(cp.Minimize(cost), const)
        # prob = cp.Problem(cp.Minimize(cost))
        prob.solve(verbose=False, solver=cp.ECOS)
        print("status:", prob.status)
        print("optimal value", prob.value)
        a_x = ax.value

        if prob.status == "optimal":
            # calculate function approximations, maintain the same coefficients
            f_est_val[p, :] = get_f_approx_poly(a_x, k, t_array)
            fdot_est_val[p, :] = get_fdot_approx_poly(a_x, k, t_array, t_interval)

            # error as the L2 norm
            err = f_true_val - f_est_val[p, :]
            err_dot = fdot_true_val - fdot_est_val[p, :]

            truncation_err[0, p] = np.linalg.norm(err)
            truncation_err[1, p] = np.linalg.norm(err_dot)

            max_err[p] = np.max(np.absolute(err))
            max_err_dot[p] = np.max(np.absolute(err_dot))

    return a_x, truncation_err, f_est_val, fdot_est_val


# position or velocity only
def singular_fit_polynomial(
    t_array,
    n_array,
    f_true_val,
    fdot_true_val,
    solver="ECOS",
    pos=True,
    cheby_interp=True,
    alpha=None,
):
    # coefficients
    coeff_list = np.zeros((lenn, np.max(n_array) + 1))

    if pos:
        # fit to position only
        lent = t_array.size
        truncation_err = np.zeros((1, len(n_array)))
        lenn = len(n_array)

        sise_tot = np.zeros((lenn, lent))

        if alpha is None:
            alpha = np.zeros(lenn)

        # estimated values
        f_est_val = np.zeros((lenn, lent))

        # construct Cubic Spline Interporation (used to compute interporaltion point values)
        cs_f = CubicSpline(t_array, f_true_val)

        for k in range(len(n_array)):
            n = n_array[k]
            # print("n = {}".format(n))
            # interpolation points

            num_points = 4 * n

            if cheby_interp:
                t_interpolation = cheby_interp_points(num_points)
            else:
                t_interpolation = np.linspace(-1, 1, num_points + 1)

            # true values interpolated
            f_true_interp = cs_f(t_interpolation)

            # System of equations Ta = x
            t_matrix = np.ones(((num_points + 1), n + 1))
            x_matrix = np.ones(((num_points + 1)))

            # build chebyshev polynomial matrix
            for j in range(num_points + 1):  # row
                for i in range(n + 1):  # column
                    t_matrix[j, i] = generate_polynomial_basis(i, t_interpolation[j])
                x_matrix[j] = f_true_interp[j]

            # solve convex optimization problem
            ax = cp.Variable(n + 1)
            cost = cp.sum_squares(t_matrix @ ax - x_matrix) + alpha[k] * cp.norm1(ax)
            const = []
            const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
            const.append(
                t_matrix[num_points, :] @ ax == x_matrix[num_points]
            )  # final point
            prob = cp.Problem(cp.Minimize(cost), const)
            prob.solve()
            # print("status:", prob.status)
            # print("optimal value", prob.value)
            a_x = ax.value

            if prob.status == "optimal":
                # calculate function approximations, maintain the same coefficients
                f_est_val[k, :] = get_f_approx_poly(a_x, n, t_array)

                # error as the L2 norm
                err = f_true_val - f_est_val[k, :]
                truncation_err[0, k] = np.linalg.norm(err) * 1000

                sise_tot[k, :] = [np.linalg.norm(error) * 1000 for error in err]
                coeff_list[k, 0 : n + 1] = a_x

        return coeff_list, truncation_err, f_est_val

    else:
        # fit to position only
        lent = t_array.size
        truncation_err = np.zeros((1, len(n_array)))
        lenn = len(n_array)

        sise_tot = np.zeros((lenn, lent))

        if alpha is None:
            alpha = np.zeros(lenn)

        # estimated values
        fdot_est_val = np.zeros((lenn, lent))

        # construct Cubic Spline Interporation (used to compute interporaltion point values)
        cs_f = CubicSpline(t_array, fdot_true_val)

        for k in range(len(n_array)):
            n = n_array[k]
            # print("n = {}".format(n))
            # interpolation points
            num_points = 4 * n

            if cheby_interp:
                t_interpolation = cheby_interp_points(num_points)
            else:
                t_interpolation = np.linspace(-1, 1, num_points + 1)

            # true values interpolated
            fdot_true_interp = cs_f(t_interpolation)

            # System of equations Ta = x
            t_matrix = np.ones(((num_points + 1), n + 1))
            x_matrix = np.ones(((num_points + 1)))

            # build chebyshev polynomial matrix
            for j in range(num_points + 1):  # row
                for i in range(n + 1):  # column
                    t_matrix[j, i] = generate_polynomial_basis(i, t_interpolation[j])
                x_matrix[j] = fdot_true_interp[j]

            # solve convex optimization problem
            ax = cp.Variable(n + 1)
            cost = cp.sum_squares(t_matrix @ ax - x_matrix) + alpha[k] * cp.norm1(ax)
            const = []
            const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
            const.append(
                t_matrix[num_points, :] @ ax == x_matrix[num_points]
            )  # final point
            prob = cp.Problem(cp.Minimize(cost), const)
            prob.solve()
            # print("status:", prob.status)
            # print("optimal value", prob.value)
            a_x = ax.value

            if prob.status == "optimal":
                # calculate function approximations, maintain the same coefficients
                fdot_est_val[k, :] = get_f_approx_poly(a_x, n, t_array)

                # error as the L2 norm
                err = fdot_true_val - fdot_est_val[k, :]
                truncation_err[0, k] = np.linalg.norm(err) * 1000000

                sise_tot[k, :] = [np.linalg.norm(error) * 1000000 for error in err]
                coeff_list[k, 0 : n + 1] = a_x

        return coeff_list, truncation_err, fdot_est_val


# SISE constrained polynomial fitting
# def sise_fit_polynomial(t_array, n_array, f_true_val, fdot_true_val, t_interval, interp_points, solver="ECOS",
#                         res_pos = 13.43/4, res_vel = 1.2/4, alpha = None, cheby_interp = True, end_points = True, tol = 1e-6, sample_test = False):
#     lent = t_array.size
#     truncation_err = np.zeros((2,len(n_array)))
#     lenn = len(n_array)

#     status_list = []

#     if alpha is None:
#         alpha = np.zeros(lenn)

#     # estimated values
#     f_est_val = np.zeros((lenn, lent))
#     fdot_est_val = np.zeros((lenn, lent))

#     #coefficients (polynomial order +1 )
#     coeff_list = np.zeros((lenn, np.max(n_array)+1))

#     # construct Cubic Spline Interporation (used to compute interporaltion point values)
#     cs_f = CubicSpline(t_array, f_true_val)
#     cs_fdot = CubicSpline(t_array, fdot_true_val)

#     error_tot = {key:{} for key in n_array}

#     for k in range(len(n_array)):
#         # print(f'n = {n_array[k]}')
#         n = n_array[k]
#         #interpolation points

#         #let's try something out
#         if sample_test:
#             n_samples_list = [n, 5*n, 10*n, 25*n, 50*n, 100*n]
#             mse_dict = {key : { k:[] for k in n_samples_list} for key in ['coeffs', 'mse_pos', 'mse_vel', 'status']}
#         else:
#             n_samples_list = [n]
#             mse_dict = {key : { k:[] for k in n_samples_list} for key in ['coeffs', 'mse_pos', 'mse_vel', 'status']}
#         # n_samples = [n]

#         for n_sample in n_samples_list:
#             if interp_points:
#                 t_interpolation = interp_points[k]
#             else:
#                 if cheby_interp:
#                     t_interpolation = cheby_interp_points(n_sample)
#                 else:
#                     t_interpolation = np.linspace(-1,1,n_sample+1)

#             #true values interpolated
#             f_true_interp = cs_f(t_interpolation)
#             f_dot_true_interp = cs_fdot(t_interpolation)

#             # System of equations Ta = x
#             t_matrix = np.ones((2*(n_sample+1), n+1))
#             x_matrix = np.ones((2*(n_sample+1)))

#             #build matrices
#             for j in range(n_sample+1):               # row
#                 for i in range(n+1):           # column
#                     t_matrix[j,i] = generate_polynomial_basis(i,t_interpolation[j])
#                     t_matrix[j+n_sample+1,i] = generate_der_polynomial_basis(i,t_interpolation[j],t_interval)
#                 x_matrix[j] = f_true_interp[j]
#                 x_matrix[j+n_sample+1] = f_dot_true_interp[j]

#             #solve convex optimization problem
#             ax = cp.Variable(n+1)
#             cost = cp.sum_squares(t_matrix @ ax - x_matrix)
#             #set up constraints
#             const = []
#             const.append(cp.norm_inf(t_matrix[0 : n_sample+1, :] @ ax - x_matrix[0:n_sample+1]) <= res_pos)              # position error
#             const.append(cp.norm_inf(t_matrix[n_sample+1 :, :] @ ax - x_matrix[n_sample+1:]) <= res_vel)        # velocity error

#             if end_points:
#                 const.append(t_matrix[0, :] @ ax == x_matrix[0])      # initial point
#                 const.append(t_matrix[n, :] @ ax == x_matrix[n])  # final point
#                 const.append(t_matrix[n+1, :] @ ax == x_matrix[n+1])  # initial velocity
#                 const.append(t_matrix[-1, :] @ ax == x_matrix[-1])    # final velocity

#             prob = cp.Problem(cp.Minimize(cost), const)

#             try:
#                 prob.solve(solver = solver, feastol=tol, max_iters=100000)
#                 # print("status:", prob.status)
#                 # print("optimal value", prob.value)

#                 a_x = ax.value

#                 if prob.status == "optimal" or prob.status == "optimal_inaccurate":
#                     f_approx = get_f_approx_poly(a_x,n,t_array)
#                     f_dot_approx = get_fdot_approx_poly(a_x,n,t_array,t_interval)

#                     mse = [(f_approx[i] - f_true_val[i])**2 for i in range(len(t_array))]
#                     # print(f'The positional MSE is {mse} m')
#                     mse_dict['mse_pos'][n_sample] = mse
#                     mse = [(f_dot_approx[i] - fdot_true_val[i])**2 for i in range(len(t_array))]
#                     # print(f'The velocity MSE is {mse} mm/s')
#                     mse_dict['mse_vel'][n_sample] = mse
#                     mse_dict['coeffs'][n_sample] = a_x.flatten()
#                     mse_dict['status'][n_sample] = prob.status

#                 else :
#                     mse_dict['mse_pos'][n_sample] = len(t_array)*[np.nan]
#                     mse_dict['mse_vel'][n_sample] = len(t_array)*[np.nan]
#                     mse_dict['coeffs'][n_sample] = (n+1)*[0]
#                     mse_dict['status'][n_sample] = prob.status

#                 if n_sample == n:
#                     status_list.append(prob.status)
#                     if prob.status == "optimal":
#                         # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
#                         # calculate function approximations, maintain the same coefficients
#                         f_est_val[k, :] = get_f_approx_poly(a_x,n,t_array)
#                         fdot_est_val[k, :]  = get_fdot_approx_poly(a_x,n,t_array,t_interval)

#                         #error as the L2 norm
#                         err = f_true_val - f_est_val[k, :]
#                         err_dot = fdot_true_val - fdot_est_val[k, :]

#                         truncation_err[0,k] = np.linalg.norm(err)
#                         truncation_err[1,k] = np.linalg.norm(err_dot)

#                         coeff_list[k,0:n+1] = a_x

#                     if prob.status == 'optimal_inaccurate':
#                         # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
#                         #calculate function approximations, maintain the same coefficients
#                         f_est_val[k, :] = get_f_approx_poly(a_x,n,t_array)
#                         fdot_est_val[k, :]  = get_fdot_approx_poly(a_x,n,t_array,t_interval)

#                         #error as the L2 norm
#                         err = f_true_val - f_est_val[k, :]
#                         err_dot = fdot_true_val - fdot_est_val[k, :]

#                         truncation_err[0,k] = np.linalg.norm(err)
#                         truncation_err[1,k] = np.linalg.norm(err_dot)

#                         coeff_list[k,0:n+1] = a_x

#                     # fig = plt.figure()
#                     # plt.plot(t_array, f_est_val[k, :])
#                     # plt.scatter(t_interpolation, f_true_interp)
#                     # plt.plot(t_array, f_true_val)
#                     # plt.scatter(cheby_interp_points(n), cs_f(cheby_interp_points(n)) )

#             except cp.error.SolverError:
#                 f_est_val[k, :] = np.zeros(len(t_array))
#                 fdot_est_val[k, :]  =  np.zeros(len(t_array))

#                 #error as the L2 norm
#                 err = f_true_val - f_est_val[k, :]
#                 err_dot = fdot_true_val - fdot_est_val[k, :]

#                 truncation_err[0,k] = np.linalg.norm(err)
#                 truncation_err[1,k] = np.linalg.norm(err_dot)

#                 coeff_list[k,0:n+1] = np.zeros(n+1)
#                 mse_dict['mse_pos'][n_sample] = len(t_array)*[np.nan]
#                 mse_dict['mse_vel'][n_sample] = len(t_array)*[np.nan]
#                 mse_dict['coeffs'][n_sample] = (n+1)*[0]
#                 mse_dict['status'][n_sample] =' infeasible'


#                 status_list.append('infeasible')

#         error_tot[n] = mse_dict

#     return coeff_list, truncation_err, f_est_val, fdot_est_val, status_list, error_tot

# NEW THINGS UPCOMING --> Chance constraint as a MC approximation? but it's not convex, right?...


def sise_fit_polynomial(
    t_array,
    n_array,
    f_true_val,
    fdot_true_val,
    t_interval,
    interp_points,
    solver="ECOS",
    res_pos=13.43 / 4,
    res_vel=1.2 / 4,
    alpha=None,
    cheby_interp=True,
    end_points=True,
    tol=1e-6,
    sample_test=False,
    plot_figure=False,
):
    lent = t_array.size
    truncation_err = np.zeros((2, len(n_array)))
    lenn = len(n_array)

    status_list = []

    if alpha is None:
        alpha = np.zeros(lenn)

    # estimated values
    f_est_val = np.zeros((lenn, lent))
    fdot_est_val = np.zeros((lenn, lent))

    # coefficients (polynomial order +1 )
    coeff_list = np.zeros((lenn, np.max(n_array) + 1))

    # construct Cubic Spline Interporation (used to compute interporaltion point values)
    cs_f = CubicSpline(t_array, f_true_val)
    cs_fdot = CubicSpline(t_array, fdot_true_val)

    error_tot = {key: {} for key in n_array}

    for k in range(len(n_array)):
        # print(f'n = {n_array[k]}')
        n = n_array[k]
        # interpolation points

        # let's try something out
        if sample_test:
            n_samples_list = [n, 5 * n, 10 * n, 25 * n, 50 * n, 100 * n]
            mse_dict = {
                key: {k: [] for k in n_samples_list}
                for key in ["coeffs", "mse_pos", "mse_vel", "status"]
            }
        else:
            n_samples_list = [n]
            mse_dict = {
                key: {k: [] for k in n_samples_list}
                for key in ["coeffs", "mse_pos", "mse_vel", "status"]
            }
        # n_samples = [n]

        for n_sample in n_samples_list:
            if interp_points:
                t_interpolation = interp_points[k]
            else:
                if cheby_interp:
                    t_interpolation = cheby_interp_points(n_sample)
                else:
                    t_interpolation = np.linspace(-1, 1, n_sample + 1)

            # true values interpolated
            f_true_interp = cs_f(t_interpolation)
            f_dot_true_interp = cs_fdot(t_interpolation)

            # System of equations Ta = x
            t_matrix = np.ones((2 * (n_sample + 1), n + 1))
            x_matrix = np.ones((2 * (n_sample + 1)))

            # build matrices
            for j in range(n_sample + 1):  # row
                for i in range(n + 1):  # column
                    t_matrix[j, i] = generate_polynomial_basis(i, t_interpolation[j])
                    t_matrix[j + n_sample + 1, i] = generate_der_polynomial_basis(
                        i, t_interpolation[j], t_interval
                    )
                x_matrix[j] = f_true_interp[j]
                x_matrix[j + n_sample + 1] = f_dot_true_interp[j]

            # solve convex optimization problem
            ax = cp.Variable(n + 1)
            cost = cp.sum_squares(t_matrix @ ax - x_matrix)
            # set up constraints
            const = []
            # const.append(cp.norm_inf(t_matrix[0 : n_sample+1, :] @ ax - x_matrix[0:n_sample+1]) <= res_pos)              # position error
            # const.append(cp.norm_inf(t_matrix[n_sample+1 :, :] @ ax - x_matrix[n_sample+1:]) <= res_vel)        # velocity error

            if end_points:
                const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
                const.append(t_matrix[n, :] @ ax == x_matrix[n])  # final point
                const.append(
                    t_matrix[n + 1, :] @ ax == x_matrix[n + 1]
                )  # initial velocity
                const.append(t_matrix[-1, :] @ ax == x_matrix[-1])  # final velocity

            prob = cp.Problem(cp.Minimize(cost), const)

            try:
                prob.solve(solver=solver, feastol=tol, max_iters=100000)
                # print("status:", prob.status)
                # print("optimal value", prob.value)

                a_x = ax.value

                if prob.status == "optimal" or prob.status == "optimal_inaccurate":
                    f_approx = get_f_approx_poly(a_x, n, t_array)
                    f_dot_approx = get_fdot_approx_poly(a_x, n, t_array, t_interval)

                    mse = [
                        (f_approx[i] - f_true_val[i]) ** 2 for i in range(len(t_array))
                    ]
                    # print(f'The positional MSE is {mse} m')
                    mse_dict["mse_pos"][n_sample] = mse
                    mse = [
                        (f_dot_approx[i] - fdot_true_val[i]) ** 2
                        for i in range(len(t_array))
                    ]
                    # print(f'The velocity MSE is {mse} mm/s')
                    mse_dict["mse_vel"][n_sample] = mse
                    mse_dict["coeffs"][n_sample] = a_x.flatten()
                    mse_dict["status"][n_sample] = prob.status

                else:
                    mse_dict["mse_pos"][n_sample] = len(t_array) * [np.nan]
                    mse_dict["mse_vel"][n_sample] = len(t_array) * [np.nan]
                    mse_dict["coeffs"][n_sample] = (n + 1) * [0]
                    mse_dict["status"][n_sample] = prob.status

                if n_sample == n:
                    status_list.append(prob.status)
                    if prob.status == "optimal":
                        # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                        # calculate function approximations, maintain the same coefficients
                        f_est_val[k, :] = get_f_approx_poly(a_x, n, t_array)
                        fdot_est_val[k, :] = get_fdot_approx_poly(
                            a_x, n, t_array, t_interval
                        )

                        # error as the L2 norm
                        err = f_true_val - f_est_val[k, :]
                        err_dot = fdot_true_val - fdot_est_val[k, :]

                        truncation_err[0, k] = np.linalg.norm(err)
                        truncation_err[1, k] = np.linalg.norm(err_dot)

                        coeff_list[k, 0 : n + 1] = a_x

                    if prob.status == "optimal_inaccurate":
                        # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                        # calculate function approximations, maintain the same coefficients
                        f_est_val[k, :] = get_f_approx_poly(a_x, n, t_array)
                        fdot_est_val[k, :] = get_fdot_approx_poly(
                            a_x, n, t_array, t_interval
                        )

                        # error as the L2 norm
                        err = f_true_val - f_est_val[k, :]
                        err_dot = fdot_true_val - fdot_est_val[k, :]

                        truncation_err[0, k] = np.linalg.norm(err)
                        truncation_err[1, k] = np.linalg.norm(err_dot)

                        coeff_list[k, 0 : n + 1] = a_x

                if plot_figure:
                    fig, axs = plt.subplots(2, 1, figsize=(8, 5), dpi=500)
                    # fig = plt.figure(figsize=(10,5))
                    axs[0].scatter(t_interpolation, f_true_interp)
                    axs[0].plot(t_array, f_true_val)
                    axs[0].plot(t_array, f_est_val[k, :], "--")
                    axs[0].set_ylabel("z-position [km]")

                    axs[1].scatter(t_interpolation, f_dot_true_interp)
                    axs[1].plot(t_array, fdot_true_val)
                    axs[1].plot(t_array, fdot_est_val[k, :], "--")
                    axs[1].set_xlabel("Normalized Time [-]")
                    axs[1].set_ylabel("z-velocity [km/s]")
                    # plt.scatter(cheby_interp_points(n), cs_f(cheby_interp_points(n)) )

            except cp.error.SolverError:
                f_est_val[k, :] = np.zeros(len(t_array))
                fdot_est_val[k, :] = np.zeros(len(t_array))

                # error as the L2 norm
                err = f_true_val - f_est_val[k, :]
                err_dot = fdot_true_val - fdot_est_val[k, :]

                truncation_err[0, k] = np.linalg.norm(err)
                truncation_err[1, k] = np.linalg.norm(err_dot)

                coeff_list[k, 0 : n + 1] = np.zeros(n + 1)
                mse_dict["mse_pos"][n_sample] = len(t_array) * [np.nan]
                mse_dict["mse_vel"][n_sample] = len(t_array) * [np.nan]
                mse_dict["coeffs"][n_sample] = (n + 1) * [0]
                mse_dict["status"][n_sample] = " infeasible"

                status_list.append("infeasible")

        error_tot[n] = mse_dict

    return coeff_list, truncation_err, f_est_val, fdot_est_val, status_list, error_tot


def solve_polynomial(
    df_interp, t_interval, n_array, solve_type, df_OE_interp=None, sample_test=False
):
    # evaluation points
    t_array = np.array(df_interp["time"])
    # specify alpha array if L-norm fitting
    alpha = None
    solver = "ECOS"

    if solve_type == "original":
        print("not working")
        solve_type = "sise_constrained"
        # f_true_val = np.array(df_interp['x'])
        # fdot_true_val = np.array(df_interp['v_x'])
        # # normalization parameters
        # weight_pos = 1
        # weight_vel = 1
        # f_c = 1
        # fdot_c = 1
        # # f_c = np.max(abs(f_true_val))
        # # fdot_c = np.max(abs(fdot_true_val))
        # #apply to data
        # f_true_val = f_true_val/f_c * weight_pos
        # fdot_true_val = fdot_true_val/fdot_c *weight_vel

        # #fitting
        # a_x, truncation_err, f_est_val, fdot_est_val \
        #     = fit_polynomial(t_array, n_array, f_true_val, fdot_true_val, t_interval, solver=solver, alpha = alpha, weight = 1,
        #             f_normalize=f_c/weight_pos, fdot_normalize=fdot_c/weight_vel, cheby_interp = True, more_points = False)

    elif solve_type == "MLLS":
        print("not working")
        solve_type = "sise_constrained"
        # n = 20
        # pos_req = (0.01343/5)
        # vel_req =  ((1.2e-6)/5)

        # # the true value
        # f_true_val = np.array(df_interp['x'])
        # fdot_true_val = np.array(df_interp['v_x'])

        # a_x_MLLS, f_est_val, fdot_est_val, coeff_num = MLLS_chebyshev(t_array, n, f_true_val, fdot_true_val, t_interval, solver="ECOS", maximum_iter=1000000000, set_tol = 1e-25,
        #             res_pos = pos_req, res_vel = vel_req,alpha = None, weight = 1, f_normalize=1, fdot_normalize=1, cheby_interp = True, end_points = True, tol = 1e-50,
        #             plot = False, more_points = False, low = 50, high= 52)

        # f_true_val = np.array(df_interp['y'])
        # fdot_true_val = np.array(df_interp['v_y'])

        # a_y_MLLS, f_est_val, fdot_est_val, coeff_num = MLLS_chebyshev(t_array, n, f_true_val, fdot_true_val, t_interval, solver="ECOS", maximum_iter=1000000000, set_tol = 1e-25,
        #     res_pos = pos_req, res_vel = vel_req,alpha = None, weight = 1, f_normalize=1, fdot_normalize=1, cheby_interp = True, end_points = True, tol = 1e-50,
        #     plot = False, more_points = False, low = 50, high= 52)

        # # the true value
        # f_true_val = np.array(df_interp['z'])
        # fdot_true_val = np.array(df_interp['v_z'])

        # a_z_MLLS, f_est_val, fdot_est_val, coeff_num = MLLS_chebyshev(t_array, n, f_true_val, fdot_true_val, t_interval, solver="ECOS", maximum_iter=1000000000, set_tol = 1e-25,
        #     res_pos = pos_req, res_vel = vel_req,alpha = None, weight = 1, f_normalize=1, fdot_normalize=1, cheby_interp = True, end_points = True, tol = 1e-50,
        #     plot = False, more_points = False, low = 50, high= 52)

        # return a_x_MLLS, a_y_MLLS, a_z_MLLS

    elif solve_type == "individual":

        # the true value
        f_true_val = np.array(df_interp["x"])
        fdot_true_val = np.array(df_interp["v_x"])

        # change pos argument to fit either position or velocity
        ax_list, truncation_err, f_est_val = singular_fit_polynomial(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=True
        )
        axdot_list, truncation_err_dot, fdot_est_val = singular_fit_polynomial(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=False
        )

        # the true value
        f_true_val = np.array(df_interp["y"])
        fdot_true_val = np.array(df_interp["v_y"])

        # change pos argument to fit either position or velocity
        ay_list, truncation_err, f_est_val = singular_fit_polynomial(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=True
        )
        aydot_list, truncation_err_dot, fdot_est_val = singular_fit_polynomial(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=False
        )

        # the true value
        f_true_val = np.array(df_interp["z"])
        fdot_true_val = np.array(df_interp["v_z"])

        # change pos argument to fit either position or velocity
        az_list, truncation_err, f_est_val = singular_fit_polynomial(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=True
        )
        azdot_list, truncation_err_dot, fdot_est_val = singular_fit_polynomial(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=False
        )

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
        )

    elif solve_type == "sise_constrained":
        # requirements (scaled since this is a single variable fitting)
        pos_req = 0.01343 / 6
        vel_req = (1.2e-6) / 6
        # fitting

        # interp_points = []
        # #true anomaly based interpolation points
        # if df_OE_interp['f'].iloc[-1] < df_OE_interp['f'].iloc[0]:
        #     anom_f = df_OE_interp['f'].iloc[-1] + 2*np.pi
        # else:
        #     anom_f = df_OE_interp['f'].iloc[-1]

        # for n in n_array:
        #     #just have to sample a number of evenly spaced points
        #     true_anom_points = np.linspace(df_OE_interp['f'].iloc[0], anom_f ,n+1)%(2*np.pi)
        #     # construct Cubic Spline Interporation (used to compute interporaltion point values)

        #     f_true = df_OE_interp['f']
        #     t_list = []
        #     for f in true_anom_points:
        #         t = df_OE_interp.iloc[(f_true-f).abs().argsort()[:1]].iloc[0]['time']
        #         t_list.append(2*t/(t_interval) - 1)

        #     interp_points.append(t_list)

        interp_points = None
        # the true value
        f_true_val = np.array(df_interp["x"])
        fdot_true_val = np.array(df_interp["v_x"])

        ax_list, truncation_err, fx_est_val, fdotx_est_val, x_status, x_mse = (
            sise_fit_polynomial(
                t_array,
                n_array,
                f_true_val,
                fdot_true_val,
                t_interval,
                interp_points,
                solver="ECOS",
                res_pos=pos_req,
                res_vel=vel_req,
                alpha=None,
                cheby_interp=False,
                end_points=False,
                tol=1e-8,
                sample_test=sample_test,
            )
        )

        # the true value
        f_true_val = np.array(df_interp["y"])
        fdot_true_val = np.array(df_interp["v_y"])

        ay_list, truncation_err, fy_est_val, fdoty_est_val, y_status, y_mse = (
            sise_fit_polynomial(
                t_array,
                n_array,
                f_true_val,
                fdot_true_val,
                t_interval,
                interp_points,
                solver="ECOS",
                res_pos=pos_req,
                res_vel=vel_req,
                alpha=None,
                cheby_interp=False,
                end_points=False,
                tol=1e-8,
                sample_test=sample_test,
            )
        )

        # the true value
        f_true_val = np.array(df_interp["z"])
        fdot_true_val = np.array(df_interp["v_z"])

        az_list, truncation_err, fz_est_val, fdotz_est_val, z_status, z_mse = (
            sise_fit_polynomial(
                t_array,
                n_array,
                f_true_val,
                fdot_true_val,
                t_interval,
                interp_points,
                solver="ECOS",
                res_pos=pos_req,
                res_vel=vel_req,
                alpha=None,
                cheby_interp=False,
                end_points=False,
                tol=1e-8,
                sample_test=sample_test,
            )
        )

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

        if sample_test:

            fig, axs = plt.subplots(1, 2)
            fig.set_size_inches(13, 5)

            for n_idx, n in enumerate(n_array):
                x_mse_n = x_mse[n]["mse_pos"]
                y_mse_n = y_mse[n]["mse_pos"]
                z_mse_n = z_mse[n]["mse_pos"]

                n_samples_list = [n, 5 * n, 10 * n, 25 * n, 50 * n, 100 * n]

                n_sise_array = np.ones((len(n_samples_list), 6)) * np.nan

                for k_idx, k in enumerate(n_samples_list):
                    mci_sise = np.array(
                        [
                            np.sqrt(x_mse_n[k][i] + y_mse_n[k][i] + z_mse_n[k][i])
                            * 1000
                            for i in range(len(t_array))
                        ]
                    )
                    pos_err = [
                        mci_sise.mean() - 3 * mci_sise.std(),
                        mci_sise.mean(),
                        mci_sise.mean() + 3 * mci_sise.std(),
                    ]
                    n_sise_array[k_idx, 0:3] = pos_err

                x_mse_n = x_mse[n]["mse_vel"]
                y_mse_n = y_mse[n]["mse_vel"]
                z_mse_n = z_mse[n]["mse_vel"]

                for k_idx, k in enumerate(n_samples_list):
                    mci_sise = np.array(
                        [
                            np.sqrt(x_mse_n[k][i] + y_mse_n[k][i] + z_mse_n[k][i])
                            * 1000000
                            for i in range(len(t_array))
                        ]
                    )
                    vel_err = [
                        mci_sise.mean() - 3 * mci_sise.std(),
                        mci_sise.mean(),
                        mci_sise.mean() + 3 * mci_sise.std(),
                    ]
                    n_sise_array[k_idx, 3:6] = vel_err

                n_factors = [1, 5, 10, 25, 50, 100]

                # axs[0].plot(n_factors, np.log10(n_sise_array[:,1]), '.-', label = f'{n}')
                axs[0].errorbar(
                    n_factors,
                    np.log10(n_sise_array[:, 1]),
                    linestyle="-",
                    yerr=np.abs(np.log10(n_sise_array[:, 2] - n_sise_array[:, 0])),
                    elinewidth=0.5,
                    label=f"n = {n}",
                    fmt=".",
                    capsize=4,
                )
                # axs[1].plot(n_factors, np.log10(n_sise_array[:,4]), '.-', label = f'{n}')
                axs[1].errorbar(
                    n_factors,
                    np.log10(n_sise_array[:, 4]),
                    linestyle="-",
                    yerr=np.abs(np.log10(n_sise_array[:, 5] - n_sise_array[:, 3])),
                    elinewidth=0.5,
                    label=f"n = {n}",
                    fmt=".",
                    capsize=4,
                )

            axs[0].plot(
                [1, 150],
                [np.log10(13.34)] * 2,
                "--",
                color="black",
                label="Requirement",
            )
            axs[1].plot(
                [1, 150], [np.log10(1.2)] * 2, "--", color="black", label="Requirement"
            )
            # print(f'Order is {n}')
            # print(n_sise_array)

            axs[0].set_xlabel(r"$m$ factor | $m \times n$ samples")
            axs[1].set_xlabel(r"$m$ factor | $m \times n$ samples")
            axs[0].set_ylabel("log($SISE_{pos}$) in MCI")
            axs[1].set_ylabel("log($SISE_{vel}$) in MCI")
            plt.legend()

            axs[0].set_xlim(0, 110)
            axs[1].set_xlim(0, 110)

        # TODO: make it so that we output MCI SISE for different sampled points and can calculate mean and 3sigma dev across true anomalies?
        # I think it's good enough to remain at one true anomaly

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

    else:
        print("No valid fitting algorithm specified")
