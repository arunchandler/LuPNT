/**
 * @file ex_fit_ephemeris.cc
 * @author Stanford NAV LAB
 * @brief  Example of fitting ephemeris using keplarian state
 * @version 0.1
 * @date 2025-01-12
 *
 * @copyright Copyright (c) 2025
 *
 */

// lupnt includes
#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace matplot;

int main() {
  auto dyn = MakePtr<NBodyDynamics<Real>>(IntegratorType::RKF45);
  dyn->SetIntegratorParams(IntegratorParams(20, 1e-12, 1e-12));

  std::string body = "MOON";
  Frame inertial_frame;
  NaifId id;
  Real a, e, inc, Omega, w, M_init;
  double GM;
  double sim_min = 0;
  double dt = 1.0;

  KepEphemType eph_type;

  if (body == "MOON") {
    inertial_frame = Frame::MOON_CI;
    id = NaifId::MOON;
    GM = GM_MOON;

    a = 6541.4;
    e = 0.6;
    inc = 65.5 * RAD;
    Omega = 0.0 * RAD;
    w = 90.0 * RAD;
    M_init = 0.0 * RAD;

    sim_min = 30.0;
    dt = 10.0;

    eph_type = KepEphemType::EPH20;

    dyn->AddBody(BodyT<Real>::Moon(20, 20));
    dyn->AddBody(BodyT<Real>::Earth());
    dyn->AddBody(BodyT<Real>::Sun());
  } else if (body == "EARTH") {
    inertial_frame = Frame::GCRF;
    id = NaifId::EARTH;
    GM = GM_EARTH;

    dyn->AddBody(BodyT<Real>::Earth(20, 20));
    dyn->AddBody(BodyT<Real>::Moon());
    dyn->AddBody(BodyT<Real>::Sun());

    a = 26560;
    e = 0.01;
    inc = 55.0 * RAD;
    Omega = 0.0 * RAD;
    w = 0.0 * RAD;
    M_init = 0.0 * RAD;

    eph_type = KepEphemType::EPH20;

    sim_min = 30;
    dt = 10.0;
  }

  // Initial state;
  dyn->SetFrame(inertial_frame);
  dyn->SetTimeStep(1.0);

  // Time
  Real et0_utc = Gregorian2Time(2025, 1, 1, 12, 0, 0).val();  // in UTC
  double et0 = UTC2TAI(et0_utc).val();                        // in TAI
  double etf = et0 + sim_min * 60;                            // Final time [10 min]
  int lent = (etf - et0) / dt + 1;

  // Initial state
  ClassicalOE coe({a, e, inc, Omega, w, M_init}, inertial_frame);
  CartesianOrbitState cart0 = Classical2Cart(coe, GM);
  VecX x0 = cart0.GetVec();
  dyn->SetTimeStep(1.0);
  MatXd Phi = MatXd::Zero(x0.size(), x0.size());

  MatX fit_x(lent, 6);
  VecX fit_t(lent);

  fit_t(0) = et0;
  fit_x.row(0) = x0.transpose();

  for (int i = 1; i < lent; i++) {
    Real t_prev = et0 + (i - 1) * dt;
    Real t_curr = et0 + i * dt;
    fit_t(i) = t_curr;
    VecX x_next = dyn->Propagate(x0, t_prev, t_curr, nullptr);
    fit_x.row(i) = x_next.transpose();
    x0 = x_next;
  }

  //
  std::cout << "Propagated Fit Orbit" << std::endl;

  // Create the ephemeris object
  auto ephemeris = MakePtr<KeplarianEphemeris>(id, eph_type);

  // Fit the ephemeris
  EphemFitOption fit_option;
  fit_option.debug = true;
  fit_option.max_iter = 20;
  fit_option.lm_down = 9.0;
  fit_option.lm_up = 11.0;
  fit_option.lm_min = 1e-12;
  fit_option.lm_max = 1e20;
  fit_option.eps_cost = 1e-11;
  fit_option.eps_lm = 1e-1;
  fit_option.lm_lambda0 = 1e10;
  EphemFitResult fit_result
      = ephemeris->FitEphemeris(fit_t, fit_x, inertial_frame, fit_option, true);

  // Print the result
  std::cout << "Fit Result: " << fit_result.x.transpose() << std::endl;
  std::cout << "Cost: " << fit_result.cost << std::endl;

  // Saving Result =========================================================================
  std::string filename = "ex_ephemeris_fit/" + body;
  filename
      += "/sim_" + std::to_string(int(sim_min)) + "_M0_" + std::to_string(int(M_init.val() * DEG));
  filename += "/EPH" + std::to_string(ephemeris->GetEphemerisSize());
  auto output_path = GetOutputPath(filename);
  std::filesystem::create_directories(output_path);
  ephemeris->SetSavePath(output_path);

  // Plot the harmonics terms -----------------------------------------------------------------
  bool is_plot = true;

  int ref_idx = 0;
  MatX rv_b(lent, 6);
  rv_b.block(0, 0, lent, 3) = fit_result.pos_true;
  rv_b.block(0, 3, lent, 3) = fit_result.vel_true;
  Vec6 rv_ref = rv_b.row(ref_idx).transpose();
  Vec6 coe_ref = Cart2Classical(rv_ref, GM);

  auto params = ephemeris->ComputeInitHarmCoeff(fit_t, rv_b, fit_t[ref_idx], coe_ref, is_plot);

  // Plot the true and estimated position -------------------------------------
  bool plot_init = true;
  Vec6d rms_error;
  rms_error = ephemeris->EvalFitError(fit_result, is_plot, "true_estimated_position.png", false);
  rms_error
      = ephemeris->EvalFitError(fit_result, is_plot, "true_estimated_position_init.png", true);

  std::vector<std::string> rms_labels
      = {"RMS Radial: ", "RMS Tangential: ", "RMS Normal: ", "RMS V_R: ", "RMS V_T: ", "RMS V_N: "};
  std::cout << " " << std::endl;
  Vec6d scale = {1000, 1000, 1000, 1e6, 1e6, 1e6};

  for (int i = 0; i < 6; i++) {
    std::cout << rms_labels[i] << scale[i] * rms_error(i) << std::endl;
  }
  std::cout << "RMS Total Position [m]: " << 1000 * rms_error.head(3).norm() << std::endl;
  std::cout << "RMS Total Velocity [mm/s]: " << 1e6 * rms_error.tail(3).norm() << std::endl;

  return 0;
}