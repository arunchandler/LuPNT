import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


def generate_sinusoidal_bases(n, t):
    """Generates sinusoidal bases
    Args:
        n (int): 2*n+1 bases to generate
        t (float): time of evaluation

    Returns:
        float: 2n+1 coefficients at time t
    """

    # note, we shall consider a normalized time interval
    bases = np.ones(2 * n + 1)
    bases[0] = 1 / 2

    for k in range(0, n):
        # bases[2*i+1] = np.sin((2*np.pi*t*(i+1))/t_interval)
        # bases[2*i+2] = np.cos((2*np.pi*t*(i+1))/t_interval)
        bases[2 * k + 1] = np.sin((k + 1) * t)
        bases[2 * k + 2] = np.cos((k + 1) * t)

    return bases


def generate_der_sinusoidal_bases(n, t, t_interval):
    """Generates derivative of sinusoidal bases
    Args:
        n (int): 2*n+1 bases to generate
        t (float): time of evaluation
        t_interval (float): total time interval

    Returns:
        float: 2n+1 coefficients at time t
    """

    bases = np.zeros(2 * n + 1)
    bases[0] = 0

    for k in range(0, n):
        # bases[2*i+1] = np.cos((2*np.pi*t*(i+1))/t_interval) * ((2*np.pi*(i+1))/t_interval)
        # bases[2*i+2] = -np.sin((2*np.pi*t*(i+1))/t_interval) * ((2*np.pi*(i+1))/t_interval)
        bases[2 * k + 1] = np.cos((k + 1) * t) * ((k + 1)) * (2 * np.pi / t_interval)
        bases[2 * k + 2] = -np.sin((k + 1) * t) * ((k + 1)) * (2 * np.pi / t_interval)

    return bases


def get_f_approx_sin(theta, n, t_array):
    """Builds approximation

    Args:
        theta (float): coefficients obtained
        n (int): fourier series order
        t_array (numpy.array): time span (it's normalized)

    Returns:
        list: approximation at each time step
    """
    f_approx = np.zeros(len(t_array))

    for jj in range(len(t_array)):
        t = t_array[jj]
        sin_base = generate_sinusoidal_bases(n, t)  # a 2n+1 list of bases
        for n_val in range(2 * n + 1):  # again, there are 2n+1 values
            f_approx[jj] += theta[n_val] * sin_base[n_val]

    return f_approx


def get_fdot_approx_sin(theta, n, t_array, t_interval):
    """Builds approximation of derivative

    Args:
        theta (float): coefficients obtained
        n (int): fourier series order
        t_array (numpy.array): time span (it's normalized)
        t_interval (float): total time interval

    Returns:
        list: approximation of the derivative at each time step
    """

    f_approx = np.zeros(len(t_array))

    for jj in range(len(t_array)):
        t = t_array[jj]
        sin_base = generate_der_sinusoidal_bases(
            n, t, t_interval
        )  # a 2n+1 list of bases
        for n_val in range(2 * n + 1):  # again, there are 2n+1 values
            f_approx[jj] += theta[n_val] * sin_base[n_val]

    return f_approx


def cheby_interp_points(n):
    # mindful here, we need 2n+1 points
    x_j = [np.cos(np.pi * j / (2 * n + 1)) for j in range(2 * n + 2)]
    # make sure the output is -1 to 1
    return x_j[::-1]


##FITTING FUNCTION
def fit_sinusoidal(
    t_array,
    n_array,
    f_true_val,
    fdot_true_val,
    t_interval,
    interp_points,
    solver="ECOS",
    res_pos=13.43 / 6,
    res_vel=1.2 / 6,
    alpha=None,
    cheby_interp=True,
    end_points=False,
    tol=1e-6,
    sample_test=True,
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

    # coefficients - remeber, n givens 2*n+1 coefficients
    coeff_list = np.zeros((lenn, 2 * np.max(n_array) + 1))

    # construct Cubic Spline Interporation (used to compute interporaltion point values)
    cs_f = CubicSpline(t_array, f_true_val)
    cs_fdot = CubicSpline(t_array, fdot_true_val)

    error_tot = {key: {} for key in n_array}

    for k in range(len(n_array)):
        # print(f'n = {n_array[k]}')
        n = n_array[k]
        n_samples = (2 * n) + 1

        # let's try something out
        if sample_test:
            n_samples_list = [n, 5 * n, 25 * n, 100 * n]
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

        for n_more in n_samples_list:
            # interpolation points
            if interp_points:
                t_interpolation = interp_points[k]
            else:
                if cheby_interp:
                    # passing n will output 2n_more+1
                    t_interpolation = cheby_interp_points(n_more)
                    t_interpolation = [np.pi * (c + 1) for c in t_interpolation]
                else:
                    t_interpolation = np.linspace(0, 2 * np.pi, int(2(n_more) + 1))

            # true values interpolated
            f_true_interp = cs_f(t_interpolation)
            f_dot_true_interp = cs_fdot(t_interpolation)

            # System of equations Ta = x
            t_matrix = np.ones((2 * (2 * n_more + 1), 2 * n + 1))
            x_matrix = np.ones((2 * (2 * n_more + 1)))

            # # System of equations Ta = x
            # t_matrix = np.ones((2*n_more+1, 2*n+1))
            # x_matrix = np.ones((2*n_more+1))

            # build matrices
            for j in range((2 * n_more + 1)):  # row
                # get the polynomial bases [2n+1]
                pos_bases = generate_sinusoidal_bases(n, t_interpolation[j])
                vel_bases = generate_der_sinusoidal_bases(
                    n, t_interpolation[j], t_interval
                )
                for i in range((2 * n + 1)):  # column
                    t_matrix[j, i] = pos_bases[i]
                    t_matrix[j + (2 * n_more + 1), i] = vel_bases[i]

                x_matrix[j] = f_true_interp[j]
                x_matrix[j + (2 * n_more + 1)] = f_dot_true_interp[j]

            # solve convex optimization problem
            ax = cp.Variable(2 * n + 1)
            cost = cp.sum_squares(t_matrix @ ax - x_matrix)
            const = []
            # const.append(cp.norm_inf(t_matrix[0 : (2*n_more+1), :] @ ax - x_matrix[0:(2*n_more+1)]) <= res_pos)              # position error
            # const.append(cp.norm_inf(t_matrix[(2*n_more+1) : , :] @ ax - x_matrix[(2*n_more+1):]) <= res_vel)        # velocity error

            if end_points:
                const.append(t_matrix[0, :] @ ax == x_matrix[0])  # initial point
                const.append(t_matrix[2 * n, :] @ ax == x_matrix[2 * n])  # final point
                const.append(
                    t_matrix[2 * n + 1, :] @ ax == x_matrix[2 * n + 1]
                )  # initial velocity
                const.append(t_matrix[-1, :] @ ax == x_matrix[-1])  # final velocity

            prob = cp.Problem(cp.Minimize(cost), const)

            try:
                prob.solve(solver=solver, feastol=tol, max_iters=100000)
                # print("status:", prob.status)
                # print("optimal value", prob.value)

                a_x = ax.value
                status_list.append(prob.status)

                if prob.status == "optimal" or prob.status == "optimal_inaccurate":
                    f_approx = get_f_approx_sin(a_x, n, t_array)
                    f_dot_approx = get_fdot_approx_sin(a_x, n, t_array, t_interval)

                    mse = [
                        (f_approx[i] - f_true_val[i]) ** 2 for i in range(len(t_array))
                    ]
                    # print(f'The positional MSE is {mse} m')
                    mse_dict["mse_pos"][n_more] = mse
                    mse = [
                        (f_dot_approx[i] - fdot_true_val[i]) ** 2
                        for i in range(len(t_array))
                    ]
                    # print(f'The velocity MSE is {mse} mm/s')
                    mse_dict["mse_vel"][n_more] = mse
                    mse_dict["coeffs"][n_more] = a_x.flatten()
                    mse_dict["status"][n_more] = prob.status

                else:
                    mse_dict["mse_pos"][n_more] = len(t_array) * [np.nan]
                    mse_dict["mse_vel"][n_more] = len(t_array) * [np.nan]
                    mse_dict["coeffs"][n_more] = (2 * n + 1) * [0]
                    mse_dict["status"][n_more] = prob.status

                if prob.status == "optimal":
                    # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                    # calculate function approximations, maintain the same coefficients
                    f_est_val[k, :] = get_f_approx_sin(a_x, n, t_array)
                    fdot_est_val[k, :] = get_fdot_approx_sin(
                        a_x, n, t_array, t_interval
                    )

                    # error as the L2 norm
                    err = f_true_val - f_est_val[k, :]
                    err_dot = fdot_true_val - fdot_est_val[k, :]

                    truncation_err[0, k] = np.linalg.norm(err)
                    truncation_err[1, k] = np.linalg.norm(err_dot)

                    coeff_list[k, 0 : (2 * n + 1)] = a_x

                if prob.status == "optimal_inaccurate":
                    # print("The norm of the residual is ", (cp.norm(t_matrix @ ax - x_matrix, p=2).value)**2)
                    # calculate function approximations, maintain the same coefficients
                    f_est_val[k, :] = get_f_approx_sin(a_x, n, t_array)
                    fdot_est_val[k, :] = get_fdot_approx_sin(
                        a_x, n, t_array, t_interval
                    )

                    # error as the L2 norm
                    err = f_true_val - f_est_val[k, :]
                    err_dot = fdot_true_val - fdot_est_val[k, :]

                    truncation_err[0, k] = np.linalg.norm(err)
                    truncation_err[1, k] = np.linalg.norm(err_dot)

                    coeff_list[k, 0 : (2 * n + 1)] = a_x

                # fig = plt.figure()
                # plt.plot(t_array, f_true_val)
                # plt.plot(t_array, f_est_val[k, :])
                # plt.scatter(t_interpolation, f_true_interp)

                # fig = plt.figure()
                # plt.plot(t_array, fdot_true_val)
                # plt.plot(t_array, fdot_est_val[k, :])
                # plt.scatter(t_interpolation, f_dot_true_interp)

            except cp.error.SolverError:
                f_est_val[k, :] = np.zeros(len(t_array))
                fdot_est_val[k, :] = np.zeros(len(t_array))

                # error as the L2 norm
                err = f_true_val - f_est_val[k, :]
                err_dot = fdot_true_val - fdot_est_val[k, :]

                truncation_err[0, k] = np.linalg.norm(err)
                truncation_err[1, k] = np.linalg.norm(err_dot)

                coeff_list[k, 0 : (2 * n + 1)] = np.zeros((2 * n + 1))
                mse_dict["mse_pos"][n_more] = len(t_array) * [np.nan]
                mse_dict["mse_vel"][n_more] = len(t_array) * [np.nan]
                mse_dict["coeffs"][n_more] = (2 * n + 1) * [0]
                mse_dict["status"][n_more] = " infeasible"

                status_list.append("infeasible")
        error_tot[n] = mse_dict
    return coeff_list, truncation_err, f_est_val, fdot_est_val, status_list, error_tot


def solve_sin(
    df_interp,
    t_interval,
    n_array,
    solve_type,
    df_OE_interp,
    sample_test=True,
    plot_figure=False,
):
    # evaluation points
    t_array = np.array(df_interp["time"])
    # specify alpha array if L-norm fitting
    alpha = None
    solver = "ECOS"

    # interp_points = []
    # #true anomaly based interpolation points
    # if df_OE_interp['f'].iloc[-1] < df_OE_interp['f'].iloc[0]:
    #     anom_f = df_OE_interp['f'].iloc[-1] + 2*np.pi
    # else:
    #     anom_f = df_OE_interp['f'].iloc[-1]

    # for n in n_array:
    #     #just have to sample a number of evenly spaced points
    #     true_anom_points = np.linspace(df_OE_interp['f'].iloc[0], anom_f ,(2*n+1))%(2*np.pi)
    #     # construct Cubic Spline Interporation (used to compute interporaltion point values)

    #     f_true = df_OE_interp['f']
    #     t_list = []
    #     for f in true_anom_points:
    #         t = df_OE_interp.iloc[(f_true-f).abs().argsort()[:1]].iloc[0]['time']
    #         t_list.append(2*t/(t_interval) - 1)

    #     interp_points.append(t_list)

    interp_points = None
    # requirements (scaled since this is a single variable fitting)
    pos_req = 0.01343 / 6
    vel_req = (1.2e-6) / 6

    # the true value

    f_true_val = np.array(df_interp["x"])
    # f_true_val = [2 * val/(np.max(f_true) - np.min(f_true)) for val in f_true]

    fdot_true_val = np.array(df_interp["v_x"])
    # fdot_true_val = [2 * val/(np.max(fdot_true) - np.min(fdot_true)) for val in fdot_true]

    ax_list, truncation_err, fx_est_val, fdotx_est_val, x_status, x_mse = (
        fit_sinusoidal(
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
            cheby_interp=True,
            end_points=False,
            tol=1e-8,
            sample_test=sample_test,
        )
    )

    # fig = plt.figure()
    # plt.scatter( n_array, np.log10(truncation_err[0, :]))
    # plt.scatter( n_array, np.log10(truncation_err[1, :]))
    # plt.legend(); plt.grid();plt.title("Trucation error");

    # the true value
    f_true_val = np.array(df_interp["y"])
    fdot_true_val = np.array(df_interp["v_y"])
    ay_list, truncation_err, fy_est_val, fdoty_est_val, y_status, y_mse = (
        fit_sinusoidal(
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
            cheby_interp=True,
            end_points=False,
            tol=1e-8,
            sample_test=sample_test,
        )
    )

    # the true value
    f_true_val = np.array(df_interp["z"])
    fdot_true_val = np.array(df_interp["v_z"])
    az_list, truncation_err, fz_est_val, fdotz_est_val, z_status, z_mse = (
        fit_sinusoidal(
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
            cheby_interp=True,
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

    if sample_test and plot_figure:

        fig, axs = plt.subplots(1, 2)
        fig.set_size_inches(13, 5)

        for n_idx, n in enumerate(n_array):
            x_mse_n = x_mse[n]["mse_pos"]
            y_mse_n = y_mse[n]["mse_pos"]
            z_mse_n = z_mse[n]["mse_pos"]

            n_samples_list = [n, 5 * n, 25 * n, 100 * n]

            n_sise_array = np.ones((len(n_samples_list), 6)) * np.nan

            for k_idx, k in enumerate(n_samples_list):
                mci_sise = np.array(
                    [
                        np.sqrt(x_mse_n[k][i] + y_mse_n[k][i] + z_mse_n[k][i]) * 1000
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
                        np.sqrt(x_mse_n[k][i] + y_mse_n[k][i] + z_mse_n[k][i]) * 1000000
                        for i in range(len(t_array))
                    ]
                )
                vel_err = [
                    mci_sise.mean() - 3 * mci_sise.std(),
                    mci_sise.mean(),
                    mci_sise.mean() + 3 * mci_sise.std(),
                ]
                n_sise_array[k_idx, 3:6] = vel_err

            n_factors = [1, 5, 25, 100]

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
            [1, 150], [np.log10(13.34)] * 2, "--", color="black", label="Requirement"
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
