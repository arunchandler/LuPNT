import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


def generate_chebyshev_pol(n, t):
    """Generates Chebyshev coefficients recursively

    Args:
        n (int): polynomial order
        t (float): time value

    Returns:
        float: coefficient of order n at time t
    """
    if n == 0:
        return 1
    if n == 1:
        return t
    if n >= 2:
        return 2 * t * (generate_chebyshev_pol(n - 1, t)) - (
            generate_chebyshev_pol(n - 2, t)
        )


def generate_der_chebyshev_pol(n, t, t_interval):
    """Generates derivative of Chebyshev coefficients recursively

    Args:
        n (int): polynomial order
        t (float): time value

    Returns:
        float: derivative of coefficient of order n at time t
    """
    if n == 0:
        return 0
    if n == 1:
        return 1 * (2 / t_interval)
    if n >= 2:
        return (
            2 * t * (generate_der_chebyshev_pol(n - 1, t, t_interval))
            + 2 * (generate_chebyshev_pol(n - 1, t)) * (2 / t_interval)
            - (generate_der_chebyshev_pol(n - 2, t, t_interval))
        )


def get_f_approx(a, n, t_array):
    """Builds approximation

    Args:
        a (float): coefficients obtained
        n (int): chebyshev coefficients highest order
        t_array (numpy.array): time span

    Returns:
        list: approximation at each time step
    """
    f_approx = np.zeros(len(t_array))

    # TODO: let's make this a lot faster?

    for t in range(len(t_array)):
        for n_val in range(n + 1):
            f_approx[t] += a[n_val] * generate_chebyshev_pol(n_val, t_array[t])

    return f_approx


def get_fdot_approx(a, n, t_array, t_interval):
    """Builds approximation of derivative

    Args:
        a (float): coefficients obtained (from approximatio the original function)
        n (int): chebyshev coefficients highest order
        t_array (numpy.array): time span

    Returns:
        list: approximation of the derivative at each time step
    """
    f_approx = np.zeros(len(t_array))
    for t in range(len(t_array)):
        for n_val in range(n + 1):
            f_approx[t] += a[n_val] * generate_der_chebyshev_pol(
                n_val, t_array[t], t_interval
            )

    return f_approx


def get_vel_coeffs(n, a):
    """Calculate derivative approximation coefficients based on function approximation coefficients

    Args:
        n (int): polynomial order
        a (list): list of approximation coefficients

    Returns:
        float: derivative coefficients not using derivative of Cheybeshev polynomials
    """
    if n == 0:
        return a[1] + get_vel_coeffs(2, a) / 2

    if n >= 1 and n <= len(a) - 3:
        return 2 * (n + 1) * (a[n + 1]) - get_vel_coeffs(n + 2, a)

    if n == len(a) - 2:
        return 2 * (len(a) - 1) * a[-1]

    if n == len(a) - 1:  # zero indexed
        return 0


def cheby_interp_points(n):
    # x_j = [np.cos( np.pi * (j+0.5)/(n+1) ) for j in range(n+1)]
    x_j = [np.cos(np.pi * j / (n)) for j in range(n + 1)]
    # make sure the output is -1 to 1
    return x_j[::-1]


##FITTING FUNCTION
# Original fitting function
def fit_chebyshev(
    t_array,
    n_array,
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
    plot=False,
):
    lent = t_array.size
    truncation_err = np.zeros((2, len(n_array)))
    lenn = len(n_array)

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

    for k in range(len(n_array)):
        n = n_array[k]
        # print("n = {}".format(n))
        # interpolation points

        if more_points:
            num_points = 50
        else:
            num_points = n

        if cheby_interp:
            t_interpolation = cheby_interp_points(num_points)
        else:
            t_interpolation = np.linspace(-1, 1, num_points + 1)

        # true values interpolated
        f_true_interp = cs_f(t_interpolation)
        f_dot_true_interp = cs_fdot(t_interpolation)

        # System of equations Ta = x
        t_matrix = np.ones((2 * (num_points + 1), n + 1))
        x_matrix = np.ones((2 * (num_points + 1)))

        # build chebyshev polynomial matrix
        for j in range(num_points + 1):  # row
            for i in range(n + 1):  # column
                t_matrix[j, i] = generate_chebyshev_pol(i, t_interpolation[j])
                t_matrix[j + num_points + 1, i] = generate_der_chebyshev_pol(
                    i, t_interpolation[j], t_interval
                )
            x_matrix[j] = f_true_interp[j]
            x_matrix[j + num_points + 1] = f_dot_true_interp[j]

        # scale f and fdot
        # x_max = max(abs(x_matrix[:n+1]))
        # xdot_max = max(abs(x_matrix[n+1:]))
        # x_matrix[:n+1] = x_matrix[:n+1]/x_max
        # x_matrix[n+1:] = x_matrix[n+1:]/xdot_max * weight
        # t_matrix[:n+1, :] = t_matrix[:n+1, :]/x_max
        # t_matrix[n+1:, :] = t_matrix[n+1:, :]/xdot_max * weight

        # solve convex optimization problem
        ax = cp.Variable(n + 1)
        cost = cp.sum_squares(t_matrix @ ax - x_matrix) + alpha[k] * cp.norm1(ax)
        const = []
        const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
        const.append(
            t_matrix[num_points, :] @ ax == x_matrix[num_points]
        )  # final point
        const.append(
            t_matrix[num_points + 1, :] @ ax == x_matrix[num_points + 1]
        )  # initial velocity
        const.append(t_matrix[-1, :] @ ax == x_matrix[-1])  # final velocity
        prob = cp.Problem(cp.Minimize(cost), const)
        # prob = cp.Problem(cp.Minimize(cost))
        # print(prob)
        prob.solve()
        print("status:", prob.status)
        print("optimal value", prob.value)
        a_x = ax.value

        if prob.status == "optimal":
            # calculate function approximations, maintain the same coefficients
            f_est_val[k, :] = get_f_approx(a_x, n, t_array)
            fdot_est_val[k, :] = get_fdot_approx(a_x, n, t_array, t_interval)

            # error as the L2 norm
            err = f_true_val - f_est_val[k, :]
            err_dot = fdot_true_val - fdot_est_val[k, :]

            truncation_err[0, k] = np.linalg.norm(err)
            truncation_err[1, k] = np.linalg.norm(err_dot)

            max_err[k] = np.max(np.absolute(err))
            max_err_dot[k] = np.max(np.absolute(err_dot))

    if plot:
        # plot values
        fig, [ax1, ax2] = plt.subplots(2, 1)
        # f(x)
        ax1.grid()
        ax1.set_xlabel("t")
        ax1.set_ylabel("f(x) [km]")
        for k in range(lenn):
            ax1.plot(
                t_array,
                f_normalize * f_est_val[k, :],
                "--",
                label="n={}".format(n_array[k]),
            )
        ax1.plot(t_array, f_normalize * f_true_val, label="True")
        ax1.legend()

        # fdot(x)
        ax2.grid()
        ax2.set_xlabel("t")
        ax2.set_ylabel("$\dot{f}$(x) [km/s]")
        for k in range(lenn):
            ax2.plot(
                t_array,
                fdot_normalize * fdot_est_val[k, :],
                "--",
                label="n={}".format(n_array[k]),
            )
        ax2.plot(t_array, fdot_normalize * fdot_true_val, label="True")
        ax2.legend()

        plt.show()

        # plot errors
        plt.figure()
        plt.grid()
        plt.xlabel("Number of Interp Points")
        plt.ylabel("Error")
        plt.scatter(n_array, f_normalize * truncation_err[0, :], label="f")
        plt.scatter(n_array, fdot_normalize * truncation_err[1, :], label="fdot")
        plt.yscale("log")
        plt.legend()

        plt.show()

        # print(max_err)
        # print(max_err_dot)

    return a_x, truncation_err, f_est_val, fdot_est_val


# single position or velocity chebyshev fitting
def singular_fit_chebyshev(
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
                    t_matrix[j, i] = generate_chebyshev_pol(i, t_interpolation[j])
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
                f_est_val[k, :] = get_f_approx(a_x, n, t_array)

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
                    t_matrix[j, i] = generate_chebyshev_pol(i, t_interpolation[j])
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
                fdot_est_val[k, :] = get_f_approx(a_x, n, t_array)

                # error as the L2 norm
                err = fdot_true_val - fdot_est_val[k, :]
                truncation_err[0, k] = np.linalg.norm(err) * 1000000

                sise_tot[k, :] = [np.linalg.norm(error) * 1000000 for error in err]
                coeff_list[k, 0 : n + 1] = a_x

        return coeff_list, truncation_err, fdot_est_val


# SISE constrained chebyshev fitting
def sise_fit_chebyshev(
    t_array,
    n_array,
    f_true_val,
    fdot_true_val,
    t_interval,
    solver="ECOS",
    res_pos=13.43 / 4,
    res_vel=1.2 / 4,
    alpha=None,
    cheby_interp=True,
    end_points=True,
    tol=1e-7,
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
                t_matrix[j, i] = generate_chebyshev_pol(i, t_interpolation[j])
                t_matrix[j + n + 1, i] = generate_der_chebyshev_pol(
                    i, t_interpolation[j], t_interval
                )
            x_matrix[j] = f_true_interp[j]
            x_matrix[j + n + 1] = f_dot_true_interp[j]

        # solve convex optimization problem
        ax = cp.Variable(n + 1)
        cost = cp.sum_squares(t_matrix @ ax - x_matrix) + alpha[k] * cp.norm1(ax)
        # set up constraints
        const = []
        const.append(
            cp.norm_inf(t_matrix[0 : n + 1, :] @ ax - x_matrix[0 : n + 1]) <= res_pos
        )  # position error
        const.append(
            cp.norm_inf(t_matrix[n + 1 : 2 * (n + 1), :] @ ax - x_matrix[n + 1 :])
            <= res_vel
        )  # velocity error

        if end_points:
            const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
            const.append(t_matrix[n, :] @ ax == x_matrix[n])  # final point
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
                f_est_val[k, :] = get_f_approx(a_x, n, t_array)
                fdot_est_val[k, :] = get_fdot_approx(a_x, n, t_array, t_interval)

                # error as the L2 norm
                err = f_true_val - f_est_val[k, :]
                err_dot = fdot_true_val - fdot_est_val[k, :]

                truncation_err[0, k] = np.linalg.norm(err)
                truncation_err[1, k] = np.linalg.norm(err_dot)

                coeff_list[k, 0 : n + 1] = a_x

            if prob.status == "optimal_inaccurate":
                # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                # calculate function approximations, maintain the same coefficients
                f_est_val[k, :] = get_f_approx(a_x, n, t_array)
                fdot_est_val[k, :] = get_fdot_approx(a_x, n, t_array, t_interval)

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


def MLLS_chebyshev(
    t_array,
    n,
    f_true_val,
    fdot_true_val,
    t_interval,
    solver="ECOS",
    set_tol=1e-10,
    maximum_iter=1e6,
    res_pos=(0.01343 / 3),
    res_vel=((1.2e-6) / 3),
    alpha=None,
    weight=1,
    f_normalize=1,
    fdot_normalize=1,
    cheby_interp=True,
    end_points=True,
    tol=1e-7,
    plot=False,
    more_points=True,
    low=12,
    high=17,
):

    lent = t_array.size
    # construct Cubic Spline Interporation (used to compute interporaltion point values)
    cs_f = CubicSpline(t_array, f_true_val)
    cs_fdot = CubicSpline(t_array, fdot_true_val)

    if more_points:
        num_points = 100
    else:
        num_points = n

    print("Maximum n = {}".format(n))
    # interpolation points
    if cheby_interp:
        t_interpolation = cheby_interp_points(num_points)
    else:
        t_interpolation = np.linspace(-1, 1, num_points + 1)

    # true values interpolated
    f_true_interp = cs_f(t_interpolation)
    f_dot_true_interp = cs_fdot(t_interpolation)

    # System of equations Ta = x
    t_matrix = np.ones((2 * (num_points + 1), n + 1))
    x_matrix = np.ones((2 * (num_points + 1)))

    # build chebyshev polynomial matrix
    for j in range(num_points + 1):  # row
        for i in range(n + 1):  # column
            t_matrix[j, i] = generate_chebyshev_pol(i, t_interpolation[j])
            t_matrix[j + num_points + 1, i] = generate_der_chebyshev_pol(
                i, t_interpolation[j], t_interval
            )
        x_matrix[j] = f_true_interp[j]
        x_matrix[j + num_points + 1] = f_dot_true_interp[j]

    # solve convex optimization problem
    ax = cp.Variable(n + 1)
    # mean square error
    mse = cp.sum_squares(t_matrix @ ax - x_matrix)

    # set up constraints
    const = []
    # need a better combined tollerance
    const.append(mse <= set_tol * n)
    const.append(
        cp.norm_inf(t_matrix[0 : num_points + 1, :] @ ax - x_matrix[0 : num_points + 1])
        <= res_pos
    )  # position error
    const.append(
        cp.norm_inf(t_matrix[num_points + 1 :, :] @ ax - x_matrix[num_points + 1 :])
        <= res_vel
    )  # velocity error
    # # end point constraints
    if end_points:
        const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
        const.append(
            t_matrix[num_points, :] @ ax == x_matrix[num_points]
        )  # final point
        const.append(
            t_matrix[num_points + 1, :] @ ax == x_matrix[num_points + 1]
        )  # initial velocity
        const.append(t_matrix[-1, :] @ ax == x_matrix[-1])  # final velocity

    # print(x_matrix[num_points+1])
    # print(x_matrix[-1])
    # print(f_dot_true_interp)

    prob = cp.Problem(cp.Minimize(cp.length(ax)), const)
    prob.solve(
        qcp=True, solver=cp.ECOS, max_iters=maximum_iter, low=low, high=high, eps=tol
    )
    # prob.solve(qcp = True, solver = cp.ECOS,  feastol = 1e-8)#, solver.max_iter = 100000)

    # print(const)

    print("status:", prob.status)
    print("Found a solution, with length: ", prob.value)
    print("Mean Square Error: ", mse.value)
    print("Mean Square Error threshold: ", set_tol * n)
    print("Coefficients found: ", ax.value)
    a_x = ax.value[0 : int(prob.value)]
    lenn = len(a_x)

    # print(np.linalg.norm(t_matrix @ a_x - x_matrix)**2)

    if prob.status == "optimal":
        n = int(prob.value) - 1
        print(n)
        # estimated values
        f_est_val = get_f_approx(a_x, n, t_array)
        fdot_est_val = get_fdot_approx(a_x, n, t_array, t_interval)

        # error as the L2 norm
        err = f_true_val - f_est_val
        err_dot = fdot_true_val - fdot_est_val

        # print(f'Max error: {np.max(np.absolute(err))}')
        # print(f'Max dot error: {np.max(np.absolute(err_dot))}')

        # truncation_err[0,k] = np.linalg.norm(err)
        # truncation_err[1,k] = np.linalg.norm(err_dot)

        # print(np.absolute(t_matrix[n+1 : 2*(n+1), :] @ a_x - x_matrix[n+1:]))

    if plot:
        # plot values
        fig, [ax1, ax2] = plt.subplots(2, 1)
        # f(x)
        ax1.grid()
        ax1.set_xlabel("t")
        ax1.set_ylabel("f(x)")
        for k in range(lenn):
            ax1.plot(
                t_array, f_normalize * f_est_val[:], "--", label="n={}".format(lenn)
            )
        ax1.plot(t_array, f_normalize * f_true_val, label="True")
        ax1.legend()

        # fdot(x)
        ax2.grid()
        ax2.set_xlabel("t")
        ax2.set_ylabel("fdot(x)")
        for k in range(lenn):
            ax2.plot(
                t_array,
                fdot_normalize * fdot_est_val[:],
                "--",
                label="n={}".format(lenn),
            )
        ax2.plot(t_array, fdot_normalize * fdot_true_val, label="True")
        ax2.legend()

        plt.show()

        # # plot errors
        # plt.figure()
        # plt.grid()
        # plt.xlabel('Number of Interp Points')
        # plt.ylabel('Error')
        # plt.scatter(n_array, f_normalize*truncation_err[0,:], label='f')
        # plt.scatter(n_array, fdot_normalize*truncation_err[1,:], label='fdot')
        # plt.yscale("log")
        # plt.legend()

        # plt.show()

    return a_x, f_est_val, fdot_est_val, lenn


def solve_chebyshev(df_interp, t_interval, n_array, solve_type, plot=False):
    # evaluation points
    t_array = np.array(df_interp["time"])
    # specify alpha array if L-norm fitting
    alpha = None
    solver = "ECOS"

    if solve_type == "original":
        print("not working")
        f_true_val = np.array(df_interp["x"])
        fdot_true_val = np.array(df_interp["v_x"])
        # normalization parameters
        weight_pos = 1
        weight_vel = 1
        f_c = 1
        fdot_c = 1
        # f_c = np.max(abs(f_true_val))
        # fdot_c = np.max(abs(fdot_true_val))
        # apply to data
        f_true_val = f_true_val / f_c * weight_pos
        fdot_true_val = fdot_true_val / fdot_c * weight_vel

        # fitting
        a_x, truncation_err, f_est_val, fdot_est_val = fit_chebyshev(
            t_array,
            n_array,
            f_true_val,
            fdot_true_val,
            t_interval,
            solver=solver,
            alpha=alpha,
            weight=1,
            f_normalize=f_c / weight_pos,
            fdot_normalize=fdot_c / weight_vel,
            cheby_interp=True,
            more_points=False,
        )

    elif solve_type == "MLLS":
        n = 20
        pos_req = 0.01343 / 5
        vel_req = (1.2e-6) / 5

        # the true value
        f_true_val = np.array(df_interp["x"])
        fdot_true_val = np.array(df_interp["v_x"])

        a_x_MLLS, f_est_val, fdot_est_val, coeff_num = MLLS_chebyshev(
            t_array,
            n,
            f_true_val,
            fdot_true_val,
            t_interval,
            solver="ECOS",
            maximum_iter=1000000000,
            set_tol=1e-25,
            res_pos=pos_req,
            res_vel=vel_req,
            alpha=None,
            weight=1,
            f_normalize=1,
            fdot_normalize=1,
            cheby_interp=True,
            end_points=True,
            tol=1e-50,
            plot=False,
            more_points=False,
            low=50,
            high=52,
        )

        f_true_val = np.array(df_interp["y"])
        fdot_true_val = np.array(df_interp["v_y"])

        a_y_MLLS, f_est_val, fdot_est_val, coeff_num = MLLS_chebyshev(
            t_array,
            n,
            f_true_val,
            fdot_true_val,
            t_interval,
            solver="ECOS",
            maximum_iter=1000000000,
            set_tol=1e-25,
            res_pos=pos_req,
            res_vel=vel_req,
            alpha=None,
            weight=1,
            f_normalize=1,
            fdot_normalize=1,
            cheby_interp=True,
            end_points=True,
            tol=1e-50,
            plot=False,
            more_points=False,
            low=50,
            high=52,
        )

        # the true value
        f_true_val = np.array(df_interp["z"])
        fdot_true_val = np.array(df_interp["v_z"])

        a_z_MLLS, f_est_val, fdot_est_val, coeff_num = MLLS_chebyshev(
            t_array,
            n,
            f_true_val,
            fdot_true_val,
            t_interval,
            solver="ECOS",
            maximum_iter=1000000000,
            set_tol=1e-25,
            res_pos=pos_req,
            res_vel=vel_req,
            alpha=None,
            weight=1,
            f_normalize=1,
            fdot_normalize=1,
            cheby_interp=True,
            end_points=True,
            tol=1e-50,
            plot=False,
            more_points=False,
            low=50,
            high=52,
        )

        return a_x_MLLS, a_y_MLLS, a_z_MLLS

    elif solve_type == "individual":

        # the true value
        f_true_val = np.array(df_interp["x"])
        fdot_true_val = np.array(df_interp["v_x"])

        # change pos argument to fit either position or velocity
        ax_list, truncation_err, f_est_val = singular_fit_chebyshev(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=True
        )
        axdot_list, truncation_err_dot, fdot_est_val = singular_fit_chebyshev(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=False
        )

        # the true value
        f_true_val = np.array(df_interp["y"])
        fdot_true_val = np.array(df_interp["v_y"])

        # change pos argument to fit either position or velocity
        ay_list, truncation_err, f_est_val = singular_fit_chebyshev(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=True
        )
        aydot_list, truncation_err_dot, fdot_est_val = singular_fit_chebyshev(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=False
        )

        # the true value
        f_true_val = np.array(df_interp["z"])
        fdot_true_val = np.array(df_interp["v_z"])

        # change pos argument to fit either position or velocity
        az_list, truncation_err, f_est_val = singular_fit_chebyshev(
            t_array, n_array, f_true_val, fdot_true_val, solver="ECOS", pos=True
        )
        azdot_list, truncation_err_dot, fdot_est_val = singular_fit_chebyshev(
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
        pos_labels = ["x", "y", "z"]
        vel_labels = ["v_x", "v_y", "v_z"]
        for i in range(3):  # x, y, z
            f_true_val = np.array(df_interp[pos_labels[i]])
            fdot_true_val = np.array(df_interp[vel_labels[i]])

            a_list, truncation_err, f_est_val, fdot_est_val, status = (
                sise_fit_chebyshev(
                    t_array,
                    n_array,
                    f_true_val,
                    fdot_true_val,
                    t_interval,
                    solver="ECOS",
                    res_pos=pos_req,
                    res_vel=vel_req,
                    alpha=None,
                    cheby_interp=True,
                    end_points=False,
                    tol=1e-8,
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

    else:
        print("No valid fitting algorithm specified")
