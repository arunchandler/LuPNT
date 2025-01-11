import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import pylupnt as pnt
import pylupnt.ephemeris as eph
from tqdm import tqdm


def basis_fit_eval_ephemeris(
    df_dict,
    t_intervals,
    tinv_eval,
    anom_num,
    T_orbit,
    n_array,
    add_sise_constraint,
    basis,
    angle_model,
    fit_velocity=True,
    cheby_interp=True,
    plot_fig=False,
):
    """
    Fit and evaluate the ephemeris for different time intervals and starting conditions

    Parameters
    ----------
    df_dict : dict
        Dictionary with the dataframes for the different datasets

    t_intervals : list
        List with the different time intervals to evaluate

    tinv_eval : float
        Time interval for SISE constraint evaluation

    anom_num : int
        Number of true anomalies to sample at

    T_orbit : float
        Orbital period of the satellite [s]

    n_array : list
        List with the different polynomial degrees to evaluate

    add_sise_constraint : bool
        Flag to add the SISE constraint to the optimization problem

    basis : str
        Basis to use for the polynomial fitting (cheby or polynomial)

    angle_model : str
        Type of model to use for the angles fitting (linear or quadratic)

    test : bool
        Flag to plot the results of the fitting

    Returns
    -------
    sise_dict : dict
        Dictionary with the SISE values for the different time intervals and starting conditions

    """
    df_MCI = df_dict["MCI"]
    df_MCMF = df_dict["MCMF"]
    df_angles = df_dict["angles"]
    df_OE = df_dict["OE"]

    sise_dict = [
        {key: [] for key in ["f0", "t0", "n", "sise_pos", "sise_vel", "coefficients"]}
        for t_idx in range(len(t_intervals))
    ]

    # uniformly sample from mean anomaly to obtain different starting points -----------
    # figure out orbital periods

    # loacre all the true anomalies within an orbital period
    f_true = df_OE.loc[(df_OE["time"] < T_orbit)].reset_index(drop=True)["f"]
    # these are the true anomalies we hope to sample at
    f_list = np.linspace(0, 2 * np.pi, anom_num)[:-1]
    # to collect inital time and true anomaly (precisely)
    t0_list = []
    f0_list = []

    # plot f_true
    for f in f_list:
        closest_row = df_OE.iloc[(f_true - f).abs().argsort()[:1]].iloc[0]
        t0_list.append(closest_row["time"])
        f0_list.append(closest_row["f"])

    # loop over ---------------------------------------------------------------
    # 1. different time intervals
    # 2. different starting points (true anomalies)
    # 3. different polynomial degrees (n_array, inside function)
    for t_idx in range(len(t_intervals)):
        # for each time interval, need to gather data across different starting conditions
        print(
            "----------------------------------------------------------------------------------------------------------"
        )
        print(
            "{0}/{1} | Basis: {2} - Approx interval: {3:4.3f} hours".format(
                t_idx, len(t_intervals), basis, t_intervals[t_idx] / 3600
            )
        )

        for start_idx in tqdm(range(len(t0_list))):
            # print('True Anom: {0:4.2f} deg'.format(np.rad2deg(f0_list[start_idx])))
            # fix untis for starting time and true anomaly
            t_start = t0_list[start_idx]
            f_start = np.rad2deg(f0_list[start_idx])
            t_interval = t_intervals[t_idx]

            time_fit, df_interp, df_OE_interp, df_MCMF_interp, df_angles_interp = (
                tshift_interp_df(
                    t_start, t_interval, tinv_eval, df_MCI, df_MCMF, df_angles, df_OE
                )
            )
            t0 = df_interp["time"].iloc[0]
            tf = t_interval + t0

            # orbit fitting ------------------------------------------------------------
            if basis == "cheby":
                basis_inst = eph.Chebyshev()
            elif basis == "polynomial":
                basis_inst = eph.Polynomial()
            elif basis == "legendre":
                basis_inst = eph.Legendre()
            elif basis == "fourier":
                basis_inst = eph.Fourier()
            else:
                # error
                print("Basis not recognized, default to Chebyshev")
                basis_inst = eph.Chebyshev()

            # normalize time
            df_interp["time"] = df_interp["time"].apply(
                lambda x: (2 / (tf - t0)) * x - 1
            )
            t_array = np.array(df_interp["time"])

            # get coefficient lists
            (ax_list, ay_list, az_list, fx_est_val, fy_est_val, fz_est_val,
             fdotx_est_val, fdoty_est_val, fdotz_est_val, status
            ) = eph.fit_basis_orbit(df_interp, t_interval, n_array, basis_inst, add_sise_constraint,
                                    cheby_interp=cheby_interp, fit_velocity=fit_velocity, plot=True)

            # fix back time
            df_interp["time"] = df_interp["time"].apply(
                lambda x: (x + 1) / (2 / (tf - t0))
                )

            # angles fitting ------------------------------------------------------------
            # only keep points from which we will evaluate the approximation
            t_array = np.array(df_angles_interp["time"])  # evaluation points
            lent = t_array.shape[0]

            # return coefficients in this order [psi coeffs, theta coeffs, phi coeffs]
            if angle_model == "linear":
                angle_coeffs = pnt.ephemeris.linear_angle_fit(t_array, df_angles_interp)
            elif angle_model == "quadratic":
                angle_coeffs = pnt.ephemeris.quad_angle_fit(t_array, df_angles_interp)
            else:
                print(
                    "Did not specify an available angle fitting model, default to linear"
                )
                angle_coeffs = pnt.ephemeris.linear_angle_fit(t_array, df_angles_interp)

            # evaluation ---------------------------------------------------------------
            for k in range(len(n_array)):

                # store in dict ------------------------------------------------------------
                fit_dict = {
                    "time_fit": time_fit,
                    "ax_list": ax_list[k, :],  # list of coefficients for x
                    "ay_list": ay_list[k, :],
                    "az_list": az_list[k, :],
                    "fx_est_val": fx_est_val[k, :],
                    "fy_est_val": fy_est_val[k, :],
                    "fz_est_val": fz_est_val[k, :],
                    "fdotx_est_val": fdotx_est_val[k, :],
                    "fdoty_est_val": fdoty_est_val[k, :],
                    "fdotz_est_val": fdotz_est_val[k, :],
                    "angle_coeffs": angle_coeffs,
                    "status": status,
                }

                if not status[k]:
                    sise_dict[t_idx]["f0"].append(f_start)
                    sise_dict[t_idx]["t0"].append(t_start)
                    sise_dict[t_idx]["n"].append(n_array[k])
                    sise_dict[t_idx]["sise_pos"].append(np.nan)
                    sise_dict[t_idx]["sise_vel"].append(np.nan)
                    ax_n = ax_list[k].tolist()
                    ay_n = ay_list[k].tolist()
                    az_n = az_list[k].tolist()
                    sise_dict[t_idx]["coefficients"].append(
                        [ax_n, ay_n, az_n, angle_coeffs.flatten().tolist()]
                    )
                    continue

                else:
                    f_df, fdot_df, angles_store = rot_mci2pa(time_fit, fit_dict)

                    # evaluate SISE in position
                    sise_n_pos = sise_evaluation(
                        df_MCMF_interp, f_df, time_fit, typ="pos"
                    )  # in MCMF
                    sise_n_abs = np.abs(sise_n_pos)
                    sigma3_err = np.percentile(sise_n_abs, 99.7)

                    sise_dict[t_idx]["f0"].append(f_start)
                    sise_dict[t_idx]["t0"].append(t_start)
                    sise_dict[t_idx]["n"].append(n_array[k])
                    sise_dict[t_idx]["sise_pos"].append(sigma3_err)

                    # evaluate SISE in velocity
                    sise_n_vel = sise_evaluation(
                        df_MCMF_interp, fdot_df, time_fit, typ="vel"
                    )  # in MCMF

                    sise_10sec = []
                    i = 0  # must consider 10 second time intervals
                    dt = time_fit[1] - time_fit[0]
                    inv_10sec = int(np.ceil(10 / (dt * tinv_eval)))

                    for i in range(lent - 10):
                        sise_10sec.append(np.max(np.abs(sise_n_vel[i : i + inv_10sec])))
                    sise_10sec = np.array(sise_10sec)

                    sigma3_err = np.percentile(sise_10sec, 99.7)
                    sise_dict[t_idx]["sise_vel"].append(sigma3_err)

                    # store coefficients
                    ax_n = ax_list[k].tolist()
                    ay_n = ay_list[k].tolist()
                    az_n = az_list[k].tolist()
                    sise_dict[t_idx]["coefficients"].append(
                        [ax_n, ay_n, az_n, angle_coeffs.flatten().tolist()]
                    )

                    # plot fit results
                    if plot_fig:
                        plot_angle_fit(t_array, angles_store, df_angles_interp)
                        plot_orbit_fit(
                            df_MCMF_interp, f_df, time_fit, sise_n_pos, typ="pos"
                        )
                        plot_orbit_fit(
                            df_MCMF_interp, fdot_df, time_fit, sise_n_vel, typ="vel"
                        )
            # end of n_array loop
        # end of starting conditions loop

        # compute the mean and std of the SISE values for this time interval
        disp_sise_result(sise_dict, t_idx)

    # end of time intervals loop
    return sise_dict


def tshift_interp_df(t_start, t_interval, tinv_eval, df_MCI, df_MCMF, df_angles, df_OE):
    """
    Function to shift the time of the dataframes and interpolate them to the same time points

    Parameters
    ----------
    t_start : float
        The time at which the interpolation should start
    t_interval : float
        The time interval for which the interpolation should be done
    tinv_eval : float
        Time interval for SISE constraint evaluation


    """
    df_MCI_mod = df_MCI[df_MCI["time"] >= t_start].reset_index(drop=True)
    df_MCI_mod["time"] = df_MCI_mod["time"] - df_MCI_mod["time"].iloc[0]

    df_MCMF_mod = df_MCMF[df_MCMF["time"] >= t_start].reset_index(drop=True)
    df_MCMF_mod["time"] = df_MCMF_mod["time"] - df_MCMF_mod["time"].iloc[0]

    df_angles_mod = df_angles[df_angles["time"] >= t_start].reset_index(drop=True)
    df_angles_mod["time"] = df_angles_mod["time"] - df_angles_mod["time"].iloc[0]

    df_OE_mod = df_OE[df_OE["time"] >= t_start].reset_index(drop=True)
    df_OE_mod["time"] = df_OE_mod["time"] - df_OE_mod["time"].iloc[0]

    # obtain end points in time and build equally spaced array
    dt = df_MCI_mod["time"].iloc[1] - df_MCI_mod["time"].iloc[0]
    tf_idx = int(t_interval / dt)
    tfit_idx = np.arange(0, tf_idx, int(tinv_eval / dt))
    time_fit = df_MCI_mod["time"].iloc[tfit_idx].values
    lent = len(time_fit)

    df_interp = df_MCI_mod.iloc[tfit_idx].reset_index(drop=True)
    df_OE_interp = df_OE_mod.iloc[tfit_idx].reset_index(drop=True)
    df_MCMF_interp = df_MCMF_mod.iloc[tfit_idx].reset_index(drop=True)
    df_angles_interp = df_angles_mod.iloc[tfit_idx].reset_index(drop=True)

    if lent != len(df_interp):
        print("lent: ", lent)
        print("len(df_interp): ", len(df_interp))
        print("Was not able to store the desired amount of points")
        return (None,)

    return time_fit, df_interp, df_OE_interp, df_MCMF_interp, df_angles_interp


def rot_mci2pa(t_array, fit_dict):
    """
    Rotate the estimated position and velocity vectors from MCI to MCMF using te estimated angles

    Parameters
    ----------
    t_array : np.array
        Array with the time points to evaluate the fit

    fit_dict : dict
        Dictionary with the coefficients and estimated values for the fit
    """
    lent = t_array.shape[0]
    angles_store = np.zeros((6, lent))
    f_approx = np.zeros((lent, 3))
    fdot_approx = np.zeros((lent, 3))

    angle_coeffs = fit_dict["angle_coeffs"]
    fx_est_val = fit_dict["fx_est_val"]
    fy_est_val = fit_dict["fy_est_val"]
    fz_est_val = fit_dict["fz_est_val"]
    fdotx_est_val = fit_dict["fdotx_est_val"]
    fdoty_est_val = fit_dict["fdoty_est_val"]
    fdotz_est_val = fit_dict["fdotz_est_val"]

    for jj in range(lent):
        angles_t = eph.angles_from_coeffs(t_array[jj], angle_coeffs)
        angles_store[:, jj] = angles_t
        R_mat, R_dot = eph.rot_mat_mci2pa(angles_t)
        # get MCMF values at t_array[jj]
        f_MCI = np.array([[fx_est_val[jj]], [fy_est_val[jj]], [fz_est_val[jj]]])
        fdot_MCI = np.array(
            [[fdotx_est_val[jj]], [fdoty_est_val[jj]], [fdotz_est_val[jj]]]
        )
        f_MCMF = R_mat @ f_MCI
        fdot_MCMF = R_dot @ f_MCI + R_mat @ fdot_MCI
        # save values
        f_approx[jj, :] = f_MCMF.T
        fdot_approx[jj, :] = fdot_MCMF.T

    # convert to dataframe
    f_df = pd.DataFrame(
        np.array([f_approx[:, 0], f_approx[:, 1], f_approx[:, 2]]).T,
        columns=["x", "y", "z"],
    )
    fdot_df = pd.DataFrame(
        np.array([fdot_approx[:, 0], fdot_approx[:, 1], fdot_approx[:, 2]]).T,
        columns=["v_x", "v_y", "v_z"],
    )

    return f_df, fdot_df, angles_store


def plot_angle_fit(t_array, angles_store, df_angles_interp):
    """
    Plot the estimated angles and the true angles

    Parameters
    ----------
    t_array : np.array
        Array with the time points to evaluate the fit

    angles_store : np.array
        Array with the estimated angles

    df_angles_interp : pd.DataFrame
        Dataframe with the true angles
    """

    # plot angles
    fig = plt.figure(figsize=(12, 6))
    matplotlib.rcParams.update({"font.size": 12})
    x = (t_array - t_array[0]) / pnt.SECS_HOUR
    plt.suptitle("Lunar Orientation Angles")

    t0 = df_angles_interp["time"].iloc[0]

    angle_labels = ["phi", "theta", "psi", "phi_dot", "theta_dot", "psi_dot"]
    plot_labels = [
        "Phi [deg]",
        "Theta [deg]",
        "Psi [deg]",
        "Phi_dot [deg/s]",
        "Theta_dot [deg/s]",
        "Psi_dot [deg/s]",
    ]

    for i in range(6):
        plt.subplot(3, 2, i + 1)
        plt.plot(x, np.mod(angles_store[i, :] * pnt.DEG, 360))  # predicted angles
        plt.plot(
            x, np.mod(df_angles_interp[angle_labels[i]] * pnt.DEG, 360)
        )  # actual angles)
        plt.xlabel("Hours past " + pnt.time2gregorian_string(t0) + " TAI")
        plt.ylabel(plot_labels[i])
        plt.grid()
        plt.xlim(x[0], x[-1])
        plt.legend(["Predicted", "Actual"])

    plt.tight_layout()
    plt.show()


def plot_orbit_fit(df_interp, df_approx, time_fit, sise_t, typ="pos"):
    fig, ax = plt.subplots(4, 1, figsize=(6, 10))
    ax[0].plot(time_fit, sise_t, "o-")
    if typ == "pos":
        ax[0].set_title("SISE (Position)")
        ax[0].set_ylabel("SISE [m]")
    else:
        ax[0].set_title("SISE (Velocity)")
        ax[0].set_ylabel("SISE [mm/s]")
    ax[0].set_xlabel("Time [hr]")
    ax[0].grid(True)

    # plot x error
    if typ == "pos":
        ax[1].plot(
            time_fit, 1000 * (df_interp["x"] - df_approx["x"]), "o-", label="True"
        )
        ax[1].set_title("X")
        ax[1].set_ylabel("X [m]")
    else:
        ax[1].plot(
            time_fit, 1e6 * (df_interp["v_x"] - df_approx["v_x"]), "o-", label="True"
        )
        ax[1].set_title("V_X")
        ax[1].set_ylabel("V_X [mm/s]")
    ax[1].set_xlabel("Time [hr]")
    ax[1].legend()
    ax[1].grid(True)

    # plot y error
    if typ == "pos":
        ax[2].plot(
            time_fit, 1000 * (df_interp["y"] - df_approx["y"]), "o-", label="True"
        )
        ax[2].set_title("Y")
        ax[2].set_ylabel("Y [m]")
    else:
        ax[2].plot(
            time_fit, 1e6 * (df_interp["v_y"] - df_approx["v_y"]), "o-", label="True"
        )
        ax[2].set_title("V_Y")
        ax[2].set_ylabel("V_Y [mm/s]")

    ax[2].set_xlabel("Time [hr]")
    ax[2].legend()
    ax[2].grid(True)

    # plot z error
    if typ == "pos":
        ax[3].plot(
            time_fit, 1000 * (df_interp["z"] - df_approx["z"]), "o-", label="True"
        )
        ax[3].set_title("Z")
        ax[3].set_ylabel("Z [m]")
    else:
        ax[3].plot(
            time_fit, 1e6 * (df_interp["v_z"] - df_approx["v_z"]), "o-", label="True"
        )
        ax[3].set_title("V_Z")
        ax[3].set_ylabel("V_Z [mm/s]")

    ax[3].set_xlabel("Time [hr]")
    ax[3].legend()
    ax[3].grid(True)

    plt.tight_layout()
    plt.show()


def sise_evaluation(df_interp, df_approx, time_fit, typ="pos"):
    """
    Evaluate the SISE for the position or velocity vectors

    Parameters
    ----------
    df_interp : pd.DataFrame
        Dataframe with the true values

    df_approx : pd.DataFrame
        Dataframe with the approximated values

    time_fit : np.array
        Array with the time points to evaluate the fit

    typ : str
        Type of evaluation (pos or vel)

    Returns
    -------
    sise_t : np.array
        Array with the SISE values for the time points
    """
    time_fit = np.array(time_fit) / 3600
    sise_t = np.ones((len(time_fit)))
    # careful!!!
    if typ == "pos":
        for i in range(len(time_fit)):
            true_vec = np.array(
                [df_interp["x"].iloc[i], df_interp["y"].iloc[i], df_interp["z"].iloc[i]]
            )
            approx_vec = np.array(
                [df_approx["x"].iloc[i], df_approx["y"].iloc[i], df_approx["z"].iloc[i]]
            )
            sise_t[i] = np.linalg.norm(true_vec - approx_vec, ord=2) * 1000  # meters

    else:
        for i in range(len(time_fit)):
            true_vec = np.array(
                [
                    df_interp["v_x"].iloc[i],
                    df_interp["v_y"].iloc[i],
                    df_interp["v_z"].iloc[i],
                ]
            )
            approx_vec = np.array(
                [
                    df_approx["v_x"].iloc[i],
                    df_approx["v_y"].iloc[i],
                    df_approx["v_z"].iloc[i],
                ]
            )
            sise_t[i] = np.linalg.norm(true_vec - approx_vec, ord=2) * (1000000)  # mm/s

    return sise_t


def disp_sise_result(sise_dict, t_idx):
    """
    compute the mean and std of the SISE values for this time interval

    Parameters
    ----------
    sise_dict : dict
        Dictionary with the SISE values for the different time intervals and starting conditions

    t_idx : int
        Index of the time interval to evaluate

    f0_list : list
        List with the different true anomalies to evaluate

    n_array : list
        List with the different polynomial degrees to evaluate
    """
    sise_pos = np.array(sise_dict[t_idx]["sise_pos"])
    sise_vel = np.array(sise_dict[t_idx]["sise_vel"])
    f0 = np.array(sise_dict[t_idx]["f0"])
    fmin = np.min(f0)
    # f0 value closest to 180
    fmax = f0[np.abs(f0 - 180).argmin()]
    nmin = np.min(sise_dict[t_idx]["n"])
    nmax = np.max(sise_dict[t_idx]["n"])

    idx_peri_nmin = np.where(
        (np.array(sise_dict[t_idx]["f0"]) == fmin)
        & (np.array(sise_dict[t_idx]["n"] == nmin))
    )[0]
    idx_peri_nmax = np.where(
        (np.array(sise_dict[t_idx]["f0"]) == fmin)
        & (np.array(sise_dict[t_idx]["n"] == nmax))
    )[0]
    idx_apo_nmin = np.where(
        (np.array(sise_dict[t_idx]["f0"]) == fmax)
        & (np.array(sise_dict[t_idx]["n"] == nmin))
    )[0]
    idx_apo_nmax = np.where(
        (np.array(sise_dict[t_idx]["f0"]) == fmax)
        & (np.array(sise_dict[t_idx]["n"] == nmax))
    )[0]

    sise_pos_peri_nmin = sise_pos[idx_peri_nmin]
    sise_vel_peri_nmin = sise_vel[idx_peri_nmin]
    sise_pos_peri_nmax = sise_pos[idx_peri_nmax]
    sise_vel_peri_nmax = sise_vel[idx_peri_nmax]
    sise_pos_apo_nmin = sise_pos[idx_apo_nmin]
    sise_vel_apo_nmin = sise_vel[idx_apo_nmin]
    sise_pos_apo_nmax = sise_pos[idx_apo_nmax]
    sise_vel_apo_nmax = sise_vel[idx_apo_nmax]

    print("Perilune (f={0:4.1f})".format(fmin))
    print("n = {}".format(nmin))
    print(
        "  SISE (pos): {0:6.3f}   SISE (vel): {1:6.3f}".format(
            sise_pos_peri_nmin.mean(), sise_vel_peri_nmin.mean()
        )
    )
    print("n = {}".format(nmax))
    print(
        "  SISE (pos): {0:6.3f}   SISE (vel): {1:6.3f}".format(
            sise_pos_peri_nmax.mean(), sise_vel_peri_nmax.mean()
        )
    )
    print(" ")
    print("Apolune (f={0:4.1f})".format(fmax))
    print("n = {}".format(nmin))
    print(
        "  SISE (pos): {0:6.3f}   SISE (vel): {1:6.3f}".format(
            sise_pos_apo_nmin.mean(), sise_vel_apo_nmin.mean()
        )
    )
    print("n = {}".format(nmax))
    print(
        "  SISE (pos): {0:6.3f}   SISE (vel): {1:6.3f}".format(
            sise_pos_apo_nmax.mean(), sise_vel_apo_nmax.mean()
        )
    )
