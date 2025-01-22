/**
 * @file keplarian_ephemeris.cc
 * @author Staford NAV LAB
 * @brief  Keplarian based ephemeris fitting
 * @version 0.1
 * @date 2025-01-12
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "lupnt/ephemeris/keplarian_ephemeris.h"

#include <matplot/matplot.h>

#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"
#include "lupnt/physics/orbit_state/anomaly.h"
#include "lupnt/physics/orbit_state/conversions.h"

using namespace matplot;

namespace lupnt {

  EphemFitResult IEphemeris::FitEphemeris(VecX t_fit, MatX fit_arc, Frame frame_arc,
                                          EphemFitOption fit_option, bool optimize) {
    bool debug = fit_option.debug;
    EphemFitResult fit_result;

    // Solving the linear problem to fit the ephemeris
    int lent = t_fit.size();

    // Convert to body fixed frame
    MatX fit_arc_pos_bf(lent, 3);
    MatX fit_arc_rv_bf(lent, 6);
    MatX fit_arc_rv_bf_rot(lent, 6);

    bool rotate_only = true;
    Frame frame_b = GetBodyFixedFrameName(body_id_);

    for (int i = 0; i < lent; i++) {
      Vec6 rv_in = fit_arc.row(i).transpose();
      Vec6 rv_fixed;
      if (rotate_only) {
        auto [R, r] = GetFrameRotationTranslation(t_fit(i), frame_arc, frame_b);
        rv_fixed.head(3) = R * rv_in.head(3);
        rv_fixed.tail(3) = R * rv_in.tail(3);
      } else {
        rv_fixed = ConvertFrame(t_fit(i), rv_in, frame_arc, frame_b);
      }
      Vec6 rv_fixed_rot = ConvertFrame(t_fit(i), rv_in, frame_arc, frame_b);
      fit_arc_pos_bf.row(i) = rv_fixed.head(3);
      fit_arc_rv_bf.row(i) = rv_fixed;
      fit_arc_rv_bf_rot.row(i) = rv_fixed_rot;
    }

    // Initial guess
    VecX x0 = GetInitialGuess(t_fit, fit_arc_rv_bf);  // Initial guess
    if (debug) {
      std::cout << "Initial Guess: " << x0.transpose() << std::endl;
    }
    Real t_ref = x0(0);

    // Parameters
    bool converged = false;
    double cost = 0;
    VecX x = x0.segment(1, ephem_size_ - 1);  // Remove the reference time

    if (optimize) {
      /************************************************************
       * Define the objective function
       * Note that x is a ephem_size - 1 vector
       ************************************************************ */
      auto objfunc = [this, t_ref, t_fit, fit_arc_pos_bf, lent](const VecX x, MatXd* J) {
        if (J != nullptr) {
          J->resize(3 * lent, ephem_size_ - 1);
          J->setZero();
        }

        VecXd r(3 * lent);
        VecX x_in = VecX::Zero(ephem_size_);
        x_in(0) = t_ref;  // Set the reference time
        x_in.segment(1, ephem_size_ - 1) = x;

        for (int i = 0; i < lent; i++) {
          VecX pos_fit = fit_arc_pos_bf.row(i);
          Vec3 pos_est;
          MatXd Jtmp(3, this->ephem_size_);
          VecXd dyi(3);

          if (J != nullptr) {
            pos_est = EphemerisToFixedFramePos(t_fit(i), x_in, &Jtmp);
            VecXd J_tmp_vec
                = Jtmp.block(0, 1, 3, ephem_size_ - 1);  // (1x3) x (3 x ephem) = (1xephem)
            J->block(3 * i, 0, 3, ephem_size_ - 1) = J_tmp_vec;
          } else {
            pos_est = EphemerisToFixedFramePos(t_fit(i), x_in, nullptr);
          }
          dyi = (pos_fit - pos_est).cast<double>();  // 3 x 1
          r.segment(3 * i, 3) = dyi;
        }

        return r;
      };

      /*******************************************
       * Run iteration using Lebenberg-Marquardt
       * https://people.duke.edu/~hpgavin/lm.pdf
       * ******************************************/

      VecXd dy(3 * lent);
      MatXd J(3 * lent, ephem_size_ - 1);

      double lambda = fit_option.lm_lambda0;  // Initial damping factor
      std::vector<double> cost_history;       // Store the cost history

      for (int iter = 0; iter < fit_option.max_iter; iter++) {
        dy = objfunc(x, &J);
        cost = dy.squaredNorm();
        cost_history.push_back(cost);

        // Convergence criteria1: cost
        if (cost < fit_option.eps_cost) {
          converged = true;
          if (debug) {
            std::cout << "Converged at iter: " << iter << " Cost: " << cost << std::endl;
          }
          break;
        }

        // Solve the least squares problem
        MatXd Jt = J.transpose();
        MatXd H = Jt * J;
        VecXd g = Jt * dy;

        // construct damping matrix
        MatXd H_damp = lambda * H.diagonal().asDiagonal();
        MatXd H_new = H + H_damp;

        // solve Hx = g
        // VecXd dx = H.inverse() * g;
        VecXd dx = H_new.ldlt().solve(g);
        VecXd dxT = dx.transpose();

        // Compute the new cost
        VecX x_new = x + dx;
        VecXd dy_new = objfunc(x_new, nullptr);
        double cost_new = dy_new.squaredNorm();

        double rho_tmp = (dxT * (H_damp * dx + g)).sum();
        double rho = (cost - cost_new) / std::abs(rho_tmp);

        // std::cout << "cost: " << cost << " cost_new: " << cost_new << " rho: " << rho << "
        // rho_tmp: " << rho_tmp << std::endl;

        if (rho > fit_option.eps_lm) {
          x = x_new;
          cost = cost_new;
          lambda = MaxD(lambda / fit_option.lm_down, fit_option.lm_min);

          if (debug) {
            std::cout << "  (Improved)     Iter: " << iter << " Cost: " << cost
                      << " Lambda: " << lambda << "  rho:" << rho << std::endl;
          }
        } else {
          lambda = MinD(lambda * fit_option.lm_up, fit_option.lm_max);
          if (debug) {
            std::cout << "  (Not-improved) Iter: " << iter << " Cost: " << cost
                      << " Lambda: " << lambda << "  rho:" << rho << std::endl;
          }
        }
      }  // end of iteration
    }  // end of optimization

    VecX x_sol = VecX::Zero(ephem_size_);
    x_sol(0) = t_ref;
    x_sol.segment(1, ephem_size_ - 1) = x;

    // evaluate the fitted orbit ------------------------------------
    MatXd pos_est(lent, 3), vel_est(lent, 3), pos_init(lent, 3), vel_init(lent, 3);

    for (int i = 0; i < lent; i++) {
      Vec6d rv_est = EphemerisToFixedFramePosVel(t_fit(i), x_sol).cast<double>();
      pos_est.row(i) = rv_est.head(3).transpose();
      vel_est.row(i) = rv_est.tail(3).transpose();

      Vec6d rv_init = EphemerisToFixedFramePosVel(t_fit(i), x0).cast<double>();
      pos_init.row(i) = rv_init.head(3).transpose();
      vel_init.row(i) = rv_init.tail(3).transpose();
    }

    fit_result.t = t_fit.cast<double>();
    fit_result.x = x_sol.cast<double>();
    fit_result.pos_true = fit_arc_pos_bf.cast<double>();
    fit_result.vel_true = fit_arc_rv_bf.block(0, 3, lent, 3).cast<double>();
    fit_result.vel_true_rot = fit_arc_rv_bf_rot.block(0, 3, lent, 3).cast<double>();
    fit_result.pos_est = pos_est;
    fit_result.vel_est = vel_est;
    fit_result.pos_init = pos_init;
    fit_result.vel_init = vel_init;
    fit_result.x0 = x0.cast<double>();
    fit_result.cost = cost;
    fit_result.success = converged;

    return fit_result;
  };

  Vec3 IEphemeris::EphemerisToFixedFramePos(Real t, VecX ephemeris, MatXd* J) {
    auto func = [=, this](const VecX ephem) { return EphemerisToFixedFramePos(t, ephem); };

    VecX xf;
    if (J == nullptr) {
      xf = func(ephemeris);
    } else {
      VecX ephemeris_tmp = ephemeris.cast<double>();
      *J = jacobian(func, wrt(ephemeris_tmp), at(ephemeris_tmp), xf);
    }

    return xf;
  };

  Vec6d IEphemeris::EvalFitError(EphemFitResult fit_result, bool is_plot, std::string figname,
                                 bool plot_init) {
    // Construct matrices
    VecXd fit_t = fit_result.t;
    int lent = fit_t.size();
    VecXd plot_t = (fit_t - fit_t[0] * VecXd::Ones(lent)) / 60.0;  // Convert to minutes

    MatXd rv_true(lent, 6), rv_est(lent, 6), rv_init(lent, 6), rtn_est(lent, 6), rtn_init(lent, 6),
        v_rot_true(lent, 3);
    rv_true.block(0, 0, lent, 3) = fit_result.pos_true;
    rv_true.block(0, 3, lent, 3) = fit_result.vel_true;
    v_rot_true = fit_result.vel_true_rot;
    rv_est.block(0, 0, lent, 3) = fit_result.pos_est;
    rv_est.block(0, 3, lent, 3) = fit_result.vel_est;
    rv_init.block(0, 0, lent, 3) = fit_result.pos_init;
    rv_init.block(0, 3, lent, 3) = fit_result.vel_init;

    // Convert to RTN frame
    for (int ti = 0; ti < lent; ti++) {
      Vec6d rv_true_tmp = rv_true.row(ti).transpose();
      Vec6d rv_est_tmp = rv_est.row(ti).transpose();
      Vec6d rv_init_tmp = rv_init.row(ti).transpose();

      // Get rotation matrix
      Vec3d r_c = rv_true_tmp.head(3);
      Vec3d v_c = rv_true_tmp.tail(3);
      Vec3d x = r_c.normalized();
      Vec3d y = (r_c.cross(v_c)).normalized();
      Vec3d z = y.cross(x);
      Mat3d R3;
      R3 << x.transpose(), z.transpose(), y.transpose();
      Mat6d R;
      R.block(0, 0, 3, 3) = R3;
      R.block(3, 3, 3, 3) = R3;

      // Convert to RTN
      rv_true_tmp.tail(3) = v_rot_true.row(ti).transpose();
      R = MatXd::Identity(6, 6);
      Vec6d rv_est_rtn = R * (rv_est_tmp - rv_true_tmp);
      Vec6d rv_init_rtn = R * (rv_init_tmp - rv_true_tmp);
      rtn_est.row(ti) = rv_est_rtn.transpose().cast<double>();
      rtn_init.row(ti) = rv_init_rtn.transpose().cast<double>();
    }

    std::vector<std::string> ylabels
        = {"Radial Position [m]",    "Tangential Position [m]",    "Normal Position [m]",
           "Radial Velocity [mm/s]", "Tangential Velocity [mm/s]", "Normal Velocity [mm/s]"};

    Vec6d rms_error;
    for (int i = 0; i < 6; i++) {
      VecXd x = rtn_est.col(i);
      rms_error(i) = sqrt(x.squaredNorm() / double(lent));
    }

    if (is_plot) {
      auto fig = figure(true);
      // set size
      fig->size(1200, 800);

      for (int i = 0; i < 6; i++) {
        fig->add_subplot(2, 3, i);
        hold(on);

        double scale = 1;
        if (i < 3)
          scale = 1000;
        else
          scale = 1e6;

        VecXd x = scale * rtn_est.col(i);
        VecXd x_init = scale * rtn_init.col(i);

        plot(plot_t, x, "r");
        if (plot_init) plot(plot_t, x_init, "b--");

        std::string rms_err_str;
        int prec = 3;
        rms_err_str = std::to_string(scale * rms_error(i));
        rms_err_str = rms_err_str.substr(0, rms_err_str.find(".") + (prec + 1));

        if (i < 3) {
          // rms_error in 2 decimal places
          title("Position Error (m): " + rms_err_str);
        } else {
          title("Velocity Error (mm/s): " + rms_err_str);
        }
        xlabel("Time [min]");
        xlim({0, plot_t[lent - 1]});
        ylabel(ylabels[i]);
        grid(true);
      }

      // legend({"True", "Estimated"});
      if (save_path_set_) {
        fig->save(save_path_ / figname);

        // save rtn_est into csv
        std::string csv_name = "rtn_est.csv";
        std::ofstream file(save_path_ / csv_name);
        file << "Time, r_R, r_T, r_N, v_R, v_T, v_N\n";
        for (int i = 0; i < lent; i++) {
          file << plot_t[i] * 60 << ", ";
          for (int j = 0; j < 6; j++) {
            file << rtn_est(i, j) << ", ";
          }
          file << "\n";
        }
        file.close();
      }
    }

    return rms_error;
  }

  /************************************************************************
   *  Keplarian Ephemeris
   *   - This set of ephemeris fits the orbit using Keplarian elements
   **************************************************************************/

  KeplarianEphemeris::KeplarianEphemeris(NaifId body_id, KepEphemType type) {
    name_ = "KeplarianEphemeris";
    body_id_ = body_id;

    switch (type) {
      case KepEphemType::GPS16:
        ephem_size_ = 16;
        use_dotr_ = false;
        use_dotu_ = false;
        use_cl_ = false;
        use_cr2_ = false;
        ephemeris_names_ = {
            "tref",                                                   // Reference Epoch
            "a",       "e",     "M0",        "w",   "i0",  "Omega0",  // Keplerian elements
            "delta_n", "dot_i", "dot_Omega",                          // Secular rates
            "Cuc",     "Cus",   "Crc",       "Crs", "Cic", "Cis",     // Harmonic coefficients
        };
        break;
      case KepEphemType::EPH18:
        ephem_size_ = 18;
        use_dotr_ = true;
        use_dotu_ = true;
        use_cl_ = false;
        use_cr2_ = false;
        ephemeris_names_ = {
            "tref",                                                   // Reference Epoch
            "a",       "e",     "M0",        "w",   "i0",  "Omega0",  // Keplerian elements
            "delta_n", "dot_i", "dot_Omega",                          // Secular rates
            "Cuc",     "Cus",   "Crc",       "Crs", "Cic", "Cis",     // Harmonic coefficients
            "dot_u",   "dot_r"  // additional Harmonic coefficients
        };
        break;
      case KepEphemType::EPH20:
        ephem_size_ = 20;
        use_dotr_ = true;
        use_dotu_ = true;
        use_cl_ = true;
        use_cr2_ = false;
        ephemeris_names_ = {
            "tref",                                                   // Reference Epoch
            "a",       "e",     "M0",        "w",   "i0",  "Omega0",  // Keplerian elements
            "delta_n", "dot_i", "dot_Omega",                          // Secular rates
            "Cuc",     "Cus",   "Crc",       "Crs", "Cic", "Cis",     // Harmonic coefficients
            "dot_u",   "dot_r", "Clc",       "Cls"  // additional Harmonic coefficients
        };
        break;
      case KepEphemType::EPH22:
        ephem_size_ = 22;
        use_dotr_ = true;
        use_dotu_ = true;
        use_cl_ = true;
        use_cr2_ = true;
        ephemeris_names_ = {
            "tref",  // Reference Epoch
            "a",       "e",     "M0",
            "w",       "i0",    "Omega0",     // Keplerian elements
            "delta_n", "dot_i", "dot_Omega",  // Secular rates
            "Cuc",     "Cus",   "Crc",
            "Crs",     "Cic",   "Cis",  // Harmonic coefficients
            "dot_u",   "dot_r", "Clc",
            "Cls",     "Crc2",  "Crs2"  // additional Harmonic coefficients
        };
        break;
      default: break;
    }  // end cases
  }

  Vec3 KeplarianEphemeris::EphemerisToFixedFramePos(Real t, VecX ephemeris) {
    auto internal_params = ComputeInternalParams(t, ephemeris, false);
    Real r_k = internal_params["r_k"];
    Real u_k = internal_params["u_k"];
    Real i_k = internal_params["i_k"];
    Real lambda_k = internal_params["lambda_k"];

    Real x_k = r_k * cos(u_k);
    Real y_k = r_k * sin(u_k);

    Real x_pos = x_k * cos(lambda_k) - y_k * cos(i_k) * sin(lambda_k);
    Real y_pos = x_k * sin(lambda_k) + y_k * cos(i_k) * cos(lambda_k);
    Real z_pos = y_k * sin(i_k);

    Vec3 X_pos;
    X_pos << x_pos, y_pos, z_pos;

    return X_pos;
  };

  Vec6 KeplarianEphemeris::EphemerisToFixedFramePosVel(Real t, VecX ephemeris) {
    auto internal_params = ComputeInternalParams(t, ephemeris, true);
    Real r_k = internal_params["r_k"];
    Real u_k = internal_params["u_k"];
    Real i_k = internal_params["i_k"];
    Real lambda_k = internal_params["lambda_k"];
    Real rdot = internal_params["rdot"];
    Real udot = internal_params["udot"];
    Real didt = internal_params["didt"];
    Real lambdadot = internal_params["lambdadot"];

    Real x_k = r_k * cos(u_k);
    Real y_k = r_k * sin(u_k);
    Real xdot_k = rdot * cos(u_k) - r_k * udot * sin(u_k);
    Real ydot_k = rdot * sin(u_k) + r_k * udot * cos(u_k);

    Real x_pos = x_k * cos(lambda_k) - y_k * cos(i_k) * sin(lambda_k);
    Real y_pos = x_k * sin(lambda_k) + y_k * cos(i_k) * cos(lambda_k);
    Real z_pos = y_k * sin(i_k);

    // velocity terms
    // Reference:
    // https://www.gps.gov/technical/icwg/meetings/2019/09/GPS-SV-velocity-and-acceleration.pdf
    Real x_vel = -x_k * lambdadot * sin(lambda_k)
                 - y_k * (lambdadot * cos(i_k) * cos(lambda_k) - didt * sin(i_k) * sin(lambda_k))
                 + xdot_k * cos(lambda_k) - ydot_k * cos(i_k) * sin(lambda_k);
    Real y_vel = x_k * lambdadot * cos(lambda_k)
                 - y_k * (lambdadot * cos(i_k) * sin(lambda_k) + didt * sin(i_k) * cos(lambda_k))
                 + xdot_k * sin(lambda_k) + ydot_k * cos(i_k) * cos(lambda_k);
    Real z_vel = y_k * didt * cos(i_k) + ydot_k * sin(i_k);

    Vec6 X_posvel;
    X_posvel << x_pos, y_pos, z_pos, x_vel, y_vel, z_vel;

    return X_posvel;
  };

  std::map<std::string, Real> KeplarianEphemeris::ComputeInternalParams(Real t, VecX ephemeris,
                                                                        bool compute_vel) {
    // Convert the keplarian state to the fixed frame state
    VecX state = VecX::Zero(3);

    // Get Constants
    double GM = GetBodyGM(body_id_);
    double omega_b = GetBodyOmega(body_id_);

    // Get the parameters from the ephemeris
    std::map<std::string, Real> ephemeris_map = GetEphemerisMap(ephemeris);

    // ephemeris_names_ = {"tref",                 // Reference Epoch
    //         "a", "e", "M0", "w", "i0", "Omega0",        // Keplerian elements
    //         "delta_n", "dot_i", "dot_Omega",            // Secular rates
    //         "Cuc",  "Cus", "Crc", "Crs", "Cic", "Cis"   // Harmonic coefficients
    //         };
    Real t_ref = ephemeris_map["tref"];
    Real a = ephemeris_map["a"];
    Real e = ephemeris_map["e"];
    Real M0 = ephemeris_map["M0"];
    Real w = ephemeris_map["w"];
    Real i0 = ephemeris_map["i0"];
    Real Omega0 = ephemeris_map["Omega0"];
    Real delta_n = ephemeris_map["delta_n"];
    Real dot_i = ephemeris_map["dot_i"];
    Real dot_Omega = ephemeris_map["dot_Omega"];
    Real Cuc = ephemeris_map["Cuc"];
    Real Cus = ephemeris_map["Cus"];
    Real Crc = ephemeris_map["Crc"];
    Real Crs = ephemeris_map["Crs"];
    Real Cic = ephemeris_map["Cic"];
    Real Cis = ephemeris_map["Cis"];

    a = a * a;  // Convert to semi-major axis

    Real t_k = t - t_ref;
    Real n = sqrt(GM / (a * a * a));
    Real M = M0 + (n + delta_n) * t_k;
    Real E_k = Mean2EccAnomaly(M, e);
    Real nu_k = Ecc2TrueAnomaly(E_k, e);
    Real Phi = Wrap2Pi(2 * (w + nu_k));
    Real u_k = w + nu_k + Cuc * cos(Phi) + Cus * sin(Phi);
    Real r_k = a * (1 - e * cos(E_k)) + Crc * cos(Phi) + Crs * sin(Phi);
    Real i_k = i0 + dot_i * t_k + Cic * cos(Phi) + Cis * sin(Phi);
    Real lambda_k = Wrap2Pi(Omega0 + (dot_Omega - omega_b) * t_k);

    // Optional parameters to add
    Real dot_u = 0, dot_r = 0, Clc = 0, Cls = 0, Crc2 = 0, Crs2 = 0;

    if (use_dotu_) {
      dot_u = ephemeris_map["dot_u"];
      u_k = u_k + dot_u * t_k;
    }
    if (use_dotr_) {
      dot_r = ephemeris_map["dot_r"];
      r_k = r_k + dot_r * t_k;
    }
    if (use_cl_) {
      Clc = ephemeris_map["Clc"];
      Cls = ephemeris_map["Cls"];
      lambda_k = Wrap2Pi(lambda_k + Clc * cos(Phi) + Cls * sin(Phi));
    }
    if (use_cr2_) {
      Crc2 = ephemeris_map["Crc2"];
      Crs2 = ephemeris_map["Crs2"];
      r_k = r_k + Crc2 * cos(2 * Phi) + Crs2 * sin(2 * Phi);
    }

    std::map<std::string, Real> internal_params;
    internal_params["r_k"] = r_k;
    internal_params["u_k"] = u_k;
    internal_params["i_k"] = i_k;
    internal_params["lambda_k"] = lambda_k;

    if (compute_vel) {
      Real Edot = (n + delta_n) / (1 - e * cos(E_k));
      Real nudot = Edot * sqrt(1 - e * e) / (1 - e * cos(E_k));
      Real didt = dot_i + 2 * nudot * (Cis * cos(Phi) - Cic * sin(Phi));
      Real udot = nudot + dot_u + 2 * nudot * (Cus * cos(Phi) - Cuc * sin(Phi));
      Real rdot = e * a * Edot * sin(E_k) + dot_r
                  + 2 * nudot
                        * (Crs * cos(Phi) - Crc * sin(Phi) + 2 * Crc2 * cos(2 * Phi)
                           - 2 * Crs2 * sin(2 * Phi));
      Real lambdadot = (dot_Omega - omega_b) + 2 * nudot * (Cls * cos(Phi) - Clc * sin(Phi));
      internal_params["rdot"] = rdot;
      internal_params["udot"] = udot;
      internal_params["didt"] = didt;
      internal_params["lambdadot"] = lambdadot;
    }

    return internal_params;
  }

  VecX KeplarianEphemeris::GetInitialGuess(VecX fit_t, MatX fit_arc) {
    // Initilaize ephemeris with zeros
    VecX ephemeris = VecX::Zero(ephem_size_);

    // Get Constants
    double GM = GetBodyGM(body_id_);
    Real GM_r = GM;
    double omega_b = GetBodyOmega(body_id_);

    // Get the reference time and state
    int ref_idx = 0;  // fit_t.size() / 2;
    Real t_ref = fit_t(ref_idx);
    Vec6 rv_ref = fit_arc.row(ref_idx).transpose();

    // Convert to classical orbital elements
    Vec6 coe = Cart2Classical(rv_ref, GM_r);

    // Use this for intialization
    ephemeris[GetEphemIndex("tref")] = t_ref;
    ephemeris[GetEphemIndex("a")] = sqrt(coe[0]);
    ephemeris[GetEphemIndex("e")] = coe[1];
    ephemeris[GetEphemIndex("i0")] = coe[2];
    ephemeris[GetEphemIndex("Omega0")] = Wrap2Pi(coe[3]);
    ephemeris[GetEphemIndex("w")] = coe[4];
    ephemeris[GetEphemIndex("M0")] = coe[5];

    // Compute Initial Guess for the perturbation parameters
    bool initialize_coeffs = true;

    if (initialize_coeffs) {
      std::map<std::string, Real> map = ComputeInitHarmCoeff(fit_t, fit_arc, t_ref, coe, false);

      // Vec5 crcs_dn;
      for (const auto& [key, value] : map) {
        ephemeris[GetEphemIndex(key)] = value;
      }
    }

    return ephemeris;
  };

  std::map<std::string, Real> KeplarianEphemeris::ComputeInitHarmCoeff(VecX fit_t, MatX fit_arc_bf,
                                                                       Real t_ref, Vec6 coe_ref,
                                                                       bool plot_result) {
    // Constants
    double GM = GetBodyGM(body_id_);
    Real GM_r = GM;
    double omega_b = GetBodyOmega(body_id_);
    Frame frame_b = GetBodyFixedFrameName(body_id_);
    bool rotate_only = true;
    int lent = fit_t.size();

    // Get the orbit elements at the two endpoints
    VecX delta_r(lent), delta_u(lent), delta_i(lent), delta_ni(lent), delta_omegai(lent);
    MatXd Phi_mat(lent, 2), Phidot_mat(lent, 3), Phi2dot_mat(lent, 5), dot_mat(lent, 1);

    // First compute delta_n
    Real a_ref = coe_ref[0];
    Real n_ref = sqrt(GM / (a_ref * a_ref * a_ref));

    for (int i = 0; i < lent; i++) {
      Vec6 rv = fit_arc_bf.row(i).transpose();
      Vec6 coe = Cart2Classical(rv, GM_r);
      Real Mi = coe[5];

      // delta_n
      Real ni = (Mi - coe_ref[5]);
      if (ni > (PI)) {
        ni = ni - 2 * PI;
      } else if (ni < -(PI)) {
        ni = ni + 2 * PI;
      }
      ni = ni / (fit_t(i) - t_ref);

      if (fit_t(i) == t_ref) {
        delta_ni(i) = 0.0;
      } else {
        delta_ni(i) = (ni - n_ref) / (fit_t(i) - t_ref);
      }
    }

    Real delta_n = delta_ni.sum() / (lent - 1);

    for (int i = 0; i < lent; i++) {
      Vec6 rv = fit_arc_bf.row(i).transpose();
      Vec6 coe = Cart2Classical(rv, GM_r);
      Real Mi = coe[5];

      Real a_ref = coe_ref[0];
      Real n_ref = sqrt(GM / (a_ref * a_ref * a_ref));

      // Copute Angles
      Real e_ref = coe_ref[1];
      Real w_ref = coe_ref[4];
      Real M_ref = coe_ref[5] + (n_ref + delta_n) * (fit_t(i) - t_ref);
      Real Ei_ref = Mean2EccAnomaly(M_ref, e_ref);
      Real nu_i_ref = Ecc2TrueAnomaly(Ei_ref, e_ref);
      Real Phi_i = Wrap2Pi(2 * (w_ref + nu_i_ref));

      Phi_mat.row(i) << cos(Phi_i).val(), sin(Phi_i).val();
      Phidot_mat.row(i) << cos(Phi_i).val(), sin(Phi_i).val(), (fit_t(i) - t_ref).val();
      Phi2dot_mat.row(i) << cos(Phi_i).val(), sin(Phi_i).val(), cos(2 * Phi_i).val(),
          sin(2 * Phi_i).val(), (fit_t(i) - t_ref).val(),
          dot_mat.row(i) << (fit_t(i) - t_ref).val();

      // delta_r
      delta_r(i) = rv.head(3).norm() - a_ref * (1 - e_ref * cos(Ei_ref));

      // delta_u
      Real ui = coe[4] + Mean2TrueAnomaly(Mi, coe[1]);  // The real u
      Real ui_ref = w_ref + nu_i_ref;
      if (ui - ui_ref > (PI)) {
        ui = ui - 2 * PI;
      } else if (ui - ui_ref < -(PI)) {
        ui = ui + 2 * PI;
      }
      delta_u(i) = Wrap2Pi(ui - ui_ref);

      // delta_i
      delta_i(i) = coe[2] - coe_ref[2];

      // dot_omega
      Real lambda = coe[3];
      Real lambda_ref = Wrap2Pi(coe_ref[3] - omega_b * (fit_t(i) - t_ref));
      delta_omegai(i) = lambda - lambda_ref;
      if (delta_omegai(i) > (PI)) {
        delta_omegai(i) = delta_omegai(i) - 2 * PI;
      } else if (delta_omegai(i) < -(PI)) {
        delta_omegai(i) = delta_omegai(i) + 2 * PI;
      }
    }

    std::map<std::string, Real> coeffs;

    // Solve equation for r
    VecXd x_r, y_r;
    std::string rstring;
    if ((use_dotr_) && (use_cr2_)) {
      VecXd x_r = Phi2dot_mat.colPivHouseholderQr().solve(delta_r.cast<double>());
      coeffs["Crc"] = x_r[0];
      coeffs["Crs"] = x_r[1];
      coeffs["Crc2"] = x_r[3];
      coeffs["Crs2"] = x_r[4];
      coeffs["dot_r"] = x_r[5];
      y_r = Phi2dot_mat * x_r;
      rstring = "Crc: " + std::to_string(coeffs["Crc"].val())
                + " Crs: " + std::to_string(coeffs["Crs"].val())
                + " Crc2: " + std::to_string(coeffs["Crc2"].val())
                + " Crs2: " + std::to_string(coeffs["Crs2"].val())
                + " dot_r: " + std::to_string(coeffs["dot_r"].val());
    } else if (use_dotr_) {
      x_r = Phidot_mat.colPivHouseholderQr().solve(delta_r.cast<double>());
      coeffs["Crc"] = x_r[0];
      coeffs["Crs"] = x_r[1];
      coeffs["dot_r"] = x_r[2];
      y_r = Phidot_mat * x_r;
      rstring = "Crc: " + std::to_string(coeffs["Crc"].val())
                + " Crs: " + std::to_string(coeffs["Crs"].val())
                + " dot_r: " + std::to_string(coeffs["dot_r"].val());
    } else {
      x_r = Phi_mat.colPivHouseholderQr().solve(delta_r.cast<double>());
      coeffs["Crc"] = x_r[0];
      coeffs["Crs"] = x_r[1];
      y_r = Phi_mat * x_r;
      rstring = "Crc: " + std::to_string(coeffs["Crc"].val())
                + " Crs: " + std::to_string(coeffs["Crs"].val());
    }

    // Solve equation for u
    VecXd x_u, y_u;
    std::string ustring;
    if (use_dotu_) {
      x_u = Phidot_mat.colPivHouseholderQr().solve(delta_u.cast<double>());
      coeffs["Cuc"] = x_u[0];
      coeffs["Cus"] = x_u[1];
      coeffs["dot_u"] = x_u[2];
      y_u = Phidot_mat * x_u;
      ustring = "1e8 x Cuc: " + std::to_string(1e8 * coeffs["Cuc"].val())
                + " Cus: " + std::to_string(1e8 * coeffs["Cus"].val())
                + " dot_u: " + std::to_string(1e8 * coeffs["dot_u"].val());
    } else {
      x_u = Phi_mat.colPivHouseholderQr().solve(delta_u.cast<double>());
      coeffs["Cuc"] = x_u[0];
      coeffs["Cus"] = x_u[1];
      y_u = Phi_mat * x_u;
      ustring = "1e8 x Cuc: " + std::to_string(1e8 * coeffs["Cuc"].val())
                + " Cus: " + std::to_string(1e8 * coeffs["Cus"].val());
    }

    // Solve equation for i
    std::string istring;
    Vec3d x_i = Phidot_mat.colPivHouseholderQr().solve(delta_i.cast<double>());
    coeffs["Cic"] = x_i[0];
    coeffs["Cis"] = x_i[1];
    coeffs["dot_i"] = x_i[2];
    VecXd y_i = Phidot_mat * x_i;
    istring = "1e9 * Cic: " + std::to_string(1e9 * coeffs["Cic"].val())
              + " Cis: " + std::to_string(1e9 * coeffs["Cis"].val())
              + " dot_i: " + std::to_string(1e9 * coeffs["dot_i"].val());

    // Solve equation for Omega
    VecXd x_omega, y_omega;
    std::string omegastring;
    if (use_cl_) {
      x_omega = Phidot_mat.colPivHouseholderQr().solve(delta_omegai.cast<double>());
      coeffs["Clc"] = x_omega[0];
      coeffs["Cls"] = x_omega[1];
      coeffs["dot_Omega"] = x_omega[2];
      y_omega = Phidot_mat * x_omega;
      omegastring = "1e8 x Clc: " + std::to_string(1e8 * coeffs["Clc"].val())
                    + " Cls: " + std::to_string(1e8 * coeffs["Cls"].val())
                    + " dot_Omega: " + std::to_string(1e8 * coeffs["dot_Omega"].val());
    } else {
      x_omega = dot_mat.colPivHouseholderQr().solve(delta_omegai.cast<double>());
      coeffs["dot_omega"] = x_omega[0];
      y_omega = dot_mat * x_omega;
      omegastring = "dot_Omega: " + std::to_string(coeffs["dot_Omega"].val());
    }

    coeffs["delta_n"] = delta_n;

    // Plot delta_u, delta_r, and delta_i
    if (plot_result) {
      auto fig = figure(true);
      fig->position({0, 0, 800, 1500});
      VecXd plot_t = (fit_t - fit_t(0) * VecX::Ones(lent)).cast<double>();

      // Convert to double
      std::vector<double> delta_rp, delta_up, delta_ip, delta_omegap, delta_rpest, delta_upest,
          delta_ipest, delta_omegaest, cos_Phi, sin_Phi, ddelta_omega;
      delta_rp = ToDoubleVec(delta_r);
      delta_up = ToDoubleVec(delta_u);
      delta_ip = ToDoubleVec(delta_i);
      delta_omegap = ToDoubleVec(delta_omegai);
      delta_rpest = ToDoubleVec(y_r);
      delta_upest = ToDoubleVec(y_u);
      delta_ipest = ToDoubleVec(y_i);
      delta_omegaest = ToDoubleVec(y_omega);
      for (int ti = 0; ti < lent; ti++) {
        cos_Phi.push_back(Phi_mat(ti, 0));
        sin_Phi.push_back(Phi_mat(ti, 1));
        ddelta_omega.push_back(delta_omegap[ti] - delta_omegaest[ti]);
      }

      // Generate Plot
      for (int i = 0; i < 5; i++) {
        fig->add_subplot(5, 1, i);
        hold(true);
        if (i == 0) {
          plot(plot_t, delta_rp, "r");
          plot(plot_t, delta_rpest, "b--");
          ylabel("dr [km]");
          title(rstring);
        } else if (i == 1) {
          plot(plot_t, delta_up, "r");
          plot(plot_t, delta_upest, "b--");
          ylabel("du [rad]");
          title(ustring);
        } else if (i == 2) {
          plot(plot_t, delta_ip, "r");
          plot(plot_t, delta_ipest, "b--");
          ylabel("di [rad]");
          title(istring);
        } else if (i == 3) {
          // plot(plot_t, delta_omegap, "r");
          // plot(plot_t, delta_omegaest, "b--");
          plot(plot_t, ddelta_omega, "r");
          ylabel("dlambda [rad]");
          title(omegastring);
        } else if (i == 4) {
          plot(plot_t, cos_Phi, "r-");
          plot(plot_t, sin_Phi, "b-");
          ylabel("sin(Phi), cos(Phi) [rad]");
        }
        xlabel("t [s]");
        grid(true);
      }

      if (save_path_set_) {
        std::string filename = "coeffs_size_" + std::to_string(ephem_size_) + ".png";
        fig->save(save_path_ / filename);

        // save t and deltas into csv
        std::string csv_name = "coeffs_delta.csv";
        std::ofstream file(save_path_ / csv_name);
        file << "Time,delta_r,delta_u,delta_i,delta_omega,delta_r_est,delta_u_est,delta_i_est,"
                "delta_omega_est\n";
        for (int i = 0; i < lent; i++) {
          file << plot_t[i] << ", " << delta_rp[i] << ", " << delta_up[i] << ", " << delta_ip[i]
               << ", " << delta_omegap[i] << ", " << delta_rpest[i] << ", " << delta_upest[i]
               << ", " << delta_ipest[i] << ", " << delta_omegaest[i] << "\n";
        }
        file.close();
      }
    }

    for (const auto& [key, value] : coeffs) {
      std::cout << key << " : " << value << std::endl;
    }

    return coeffs;
  };

}  // namespace lupnt
