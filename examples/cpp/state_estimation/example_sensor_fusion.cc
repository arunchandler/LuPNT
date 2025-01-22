#include <lupnt/lupnt.h>

using namespace lupnt;
using std::cout;
using std::endl;
using std::format;
using std::left;
using std::right;
using std::setw;

void PrintProgressHeader() {
  std::string hline = std::string(80, '-');
  std::string sep = " | ";

  cout << endl << hline << endl;
  cout << left << setw(12) << "Time [min]";
  auto labels = {"Pos Err (3σ) [m]", "Vel Err (3σ) [mm/s]", "Clk Bias Err (3σ) [m]", "Satellites"};
  for (auto label : labels) cout << left << setw(24) << label;
  cout << endl << hline << endl;
}

VecXd ComputeEstimationErrors(const Ptr<Spacecraft> sat, EKF* ekf) {
  auto x_est = ekf->GetState();
  auto x_true = sat->GetStateVec();

  double r_err = M_KM * (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();   // [m]
  double v_err = MM_KM * (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();  // [mm/s]
  double b_err = M_KM * C * abs((x_true(6) - x_est(6)).val());                       // [m]
  double d_err = MM_KM * C * abs((x_true(7) - x_est(7)).val());                      // [mm/s]

  VecXd x_err(4);
  x_err << r_err, v_err, b_err, d_err;
  return x_err;
}

void PrintProgress(Real t, Real tf, Vec4d x_err, MatXd P, int num_sat, int num_used_meas) {
  cout.precision(5);

  int N_sigma = 3;
  double sigma_r = M_KM * sqrt(P(0, 0) + P(1, 1) + P(2, 2));   // [m]
  double sigma_v = MM_KM * sqrt(P(3, 3) + P(4, 4) + P(5, 5));  // [mm/s]
  double sigma_b = M_KM * C * sqrt(P(6, 6));                   // [m]
  std::vector<double> sigmas = {sigma_r, sigma_v, sigma_b};

  double t_min = t.val() / SECS_MINUTE, tf_min = tf.val() / SECS_MINUTE;
  cout << left << setw(12) << std::format("{:.0f}/{:.0f}", t_min, tf_min);
  for (int i = 0; i < 3; i++) {
    cout << right << setw(12) << format("{:.2f} ", x_err(i));
    cout << right << setw(12) << format("({:.2f})", N_sigma * sigmas[i]);
  }
  cout << right << setw(20) << format("{:2d} ({:2d} meas)", num_sat, num_used_meas);
  cout << endl;
};

void PrintEstimationStatistics(VecXd N_meas, MatXd error_mat, double data_ratio = 1.0) {
  int n_time = N_meas.size();
  Vec4d rms, means, stds, p68, p95, p99;

  int i_start = (int)((1.0 - data_ratio) * n_time);
  int i_end = n_time - 1;
  int n_range = i_end - i_start;

  VecXd N_meas_range(n_range);
  MatXd error_mat_range(4, n_range);

  N_meas_range = N_meas.segment(i_start, n_range);
  error_mat_range = error_mat.block(0, i_start, 4, n_range);

  // compute statistics ----------------------------------------------------
  // rms
  for (int i = 0; i < 4; i++) {
    rms(i) = RootMeanSquare(error_mat_range.row(i));
    means(i) = error_mat_range.row(i).mean();
    stds(i) = Std(error_mat_range.row(i));
    p68(i) = Percentile(error_mat_range.row(i), 0.68);
    p95(i) = Percentile(error_mat_range.row(i), 0.95);
    p99(i) = Percentile(error_mat_range.row(i), 0.99);
  }

  cout << std::format("\n\nSimulation Statistics (Last {:.2f}%)\n", data_ratio * 100);
  std::string hline = std::string(100, '-');
  std::string sep = " | ";
  int w1 = 12, w2 = 16;
  cout << hline << endl;
  cout << left << setw(w1) << "Statistics" << sep;
  auto labels = {"Position [m]", "Velocity [mm/s]", "Clk Bias [m]", "Clk Drift [mm/s]"};
  for (auto label : labels) cout << left << setw(w2) << label << sep;
  cout << endl << hline << endl << left << setw(w1) << "RMS" << sep;
  for (int i = 0; i < 4; i++) cout << left << setw(w2) << format("{:.2f}", rms(i)) << sep;
  cout << endl << left << setw(w1) << "Mean (Std)" << sep;
  for (int i = 0; i < 4; i++)
    cout << left << setw(w2) << format("{:2.2f} ({:2.2f})", means(i), stds(i)) << sep;
  cout << endl << left << setw(w1) << "68%" << sep;
  for (int i = 0; i < 4; i++) cout << left << setw(w2) << format("{:2.2f}", p68(i)) << sep;
  cout << endl << left << setw(w1) << "95%" << sep;
  for (int i = 0; i < 4; i++) cout << left << setw(w2) << format("{:2.2f}", p95(i)) << sep;
  cout << endl << left << setw(w1) << "99%" << sep;
  for (int i = 0; i < 4; i++) cout << left << setw(w2) << format("{:2.2f}", p99(i)) << sep;
  cout << endl << hline << endl;
}

class MyApp : public Application {
protected:
  double epoch0_;  // start epoch in TAI
  double epoch_;   // curent epoch
  double t_;       // Current time [s]

  Ptr<IFilter> filter_;
  FilterDynamicsFunction dynamics_func_;
  FilterMeasurementFunction meas_func_;
  JointState state_vec_;

public:
  double GetInitialEpoch() { return epoch0_; };
  double GetCurrentEpoch() { return epoch_; };
  double GetCurrrentTime() { return t_; };

  void Step(Real t_end) override;
};

int main() {
  // Time
  double t0_utc = Gregorian2Time(2025, 1, 1, 12, 0, 0).val();
  double t0_tai = ConvertTime(t0_utc, Time::UTC, Time::TAI).val();
  double dt_prop = 1.0;
  double dt_sim = 10 * SECS_MINUTE;
  double dt_print = dt_sim;
  double dt_save = dt_sim;

  // GNSS constellation
  bool use_galileo = true;
  bool use_qzss = true;
  std::string gps_tle = "gps_2025_01_01";
  std::string galileo_tle = "galileo_2025_01_01";
  std::string qzss_tle = "qzss_2025_01_01";

  // Simulation seed
  int seed = 10;
  std::srand(seed);

  // Initial State -------------------------------------------------------------
  const int Nsc = 2;      // Number of spacecraft
  Real a = 5740;          // [km] Semi-Mjor axis
  Real e = 0.58;          // [-] Eccentricity
  Real i = 54.856 * RAD;  // [rad] Inclination
  Real Omega = 0 * RAD;   // [rad] Right ascension of the ascending node
  Real w = 86.322 * RAD;  // [rad] Argument of periapsis
  Real M = 0 * RAD;       // [rad] Mean anomaly
  Real clk_bias = 0.0;    // [km]
  Real clk_drift = 0.0;   // [km/s]

  int N_orbit = 3;
  double period = GetOrbitalPeriod(a, GM_MOON).val();
  double tf = N_orbit * period;
  double tf_tai = t0_tai + tf;
  int Nt = int((tf_tai - t0_tai) / dt_sim) + 1;

  auto clk_model = ClockModel::kMiniRafs;

  // Measurements
  bool use_range = true;
  bool use_range_rate = false;

  double sis_ure_std = 5.0e-3;       // [km]
  double sis_ure_rate_std = 5.0e-6;  // [km/s]

  // Estimation Parameters ----------------------------------------------------
  int state_size = 8;                   // [km, km/s, s, s/s] r, v, bias, drift
  double sigma_r = 1.0 / sqrt(3);       // [km] Position
  double sigma_v = sigma_r * 1e-2;      // [km/s] Velocity
  double sigma_b = 1.0 / C;             // [s] Clock bias
  double sigma_d = sigma_b * 1e-3;      // [s/s] Clock drift
  double sigma_a = std::pow(10, -7.7);  // [km/s^2] Acceleration

  bool use_adaptive_proc = false;
  double alpha_Q = 0.9;  // Q_k = alpha * Q_k-1 + (1 - alpha) * (K dy dy^T K^T)

  double outlier_thresh = 3.0;

  std::vector<GnssMeasurementType> meas_types;
  if (use_range) meas_types.push_back(GnssMeasurementType::PR);
  if (use_range_rate) meas_types.push_back(GnssMeasurementType::PRR);
  int meas_type_num = meas_types.size();

  // Orbit Dynamics --------------------------------------------------------
  IntegratorParams iparams;
  iparams.abstol = 1e-12;
  iparams.reltol = 1e-12;

  auto dyn_true = MakePtr<NBodyDynamics<double>>(IntegratorType::RKF45);
  dyn_true->SetIntegratorParams(iparams);
  dyn_true->SetFrame(Frame::MOON_CI);
  dyn_true->AddBody(Body::Moon(100, 100));
  dyn_true->AddBody(Body::Earth());
  dyn_true->AddBody(Body::Sun());
  dyn_true->AddBody(Body::Mars());
  dyn_true->AddBody(Body::Venus());
  dyn_true->AddBody(Body::Jupiter());
  dyn_true->AddBody(Body::Saturn());
  dyn_true->AddBody(Body::Uranus());
  dyn_true->SetTimeStep(dt_prop);

  auto dyn_est = MakePtr<NBodyDynamics<Real>>(IntegratorType::RKF45);
  dyn_est->SetIntegratorParams(iparams);
  dyn_est->SetFrame(Frame::MOON_CI);
  dyn_est->AddBody(BodyT<Real>::Moon(18, 18));
  dyn_est->AddBody(BodyT<Real>::Earth());
  dyn_est->AddBody(BodyT<Real>::Sun());
  dyn_est->SetTimeStep(dt_prop);

  auto dyn_gnss = MakePtr<CartesianTwoBodyDynamics>(GM_EARTH);
  dyn_gnss->SetTimeStep(dt_prop);

  // Clock Dynamics --------------------------------------------------------
  auto dyn_clk_true = MakePtr<ClockDynamics>(clk_model);
  auto dyn_clk_est = MakePtr<ClockDynamics>(clk_model);
  dyn_clk_true->SetNoise(true);
  dyn_clk_est->SetNoise(false);

  // GNSS Constellation ----------------------------------------------------
  auto gnss_channel = MakePtr<GnssChannel>();
  GnssConstellation gnss_const;
  gnss_const.SetChannel(gnss_channel);
  gnss_const.SetDynamics(dyn_gnss);
  gnss_const.SetEpoch(t0_tai);
  gnss_const.LoadTleFile(GnssType::GPS, gps_tle);

  std::vector<std::string> signals = {"L1"};
  std::vector<std::string> gnss_types = {"GPS"};

  if (use_galileo) {
    gnss_const.LoadTleFile(GnssType::GALILEO, galileo_tle);
    signals.push_back("E1");
    gnss_types.push_back("GALILEO");
  }
  if (use_qzss) {
    gnss_const.LoadTleFile(GnssType::QZSS, qzss_tle);
    gnss_types.push_back("QZSS");
  }

  // Moon spacecraft --------------------------------------------------------
  ClassicalOE coe_mop({a, e, i, Omega, w, M}, Frame::MOON_OP);
  CartesianOrbitState cart_mop = Classical2Cart(coe_mop, GM_MOON);
  CartesianOrbitState cart_mci = ConvertOrbitStateFrame(cart_mop, t0_tai, Frame::MOON_CI);
  auto rv_state_mci = MakePtr<CartesianOrbitState>(cart_mci);

  Vec2 clock_vec{clk_bias, clk_drift};  // [s, s/s]
  auto clk_state = MakePtr<ClockState>(clock_vec);

  auto moon_sat = MakePtr<Spacecraft>();
  auto receiver = MakePtr<GnssReceiver>("moongpsr");

  moon_sat->AddDevice(receiver);
  moon_sat->SetDynamics(dyn_true);
  moon_sat->SetClock(*clk_state.get());
  moon_sat->SetOrbitState(rv_state_mci);
  moon_sat->SetEpoch(t0_tai);
  moon_sat->SetBodyId(NaifId::MOON);
  moon_sat->SetClockDynamics(dyn_clk_true);

  receiver->SetAgent(moon_sat);
  receiver->SetReceiverAttitudeMode("PZ_EarthPoint");
  receiver->SetChannel(gnss_channel);
  gnss_channel->AddReceiver(receiver);

  // Initial covariance -----------------------------------------------------
  MatXd P0 = InitialCovariancePosVelClock(sigma_r, sigma_v, sigma_b, sigma_d);

  // Joint state and dynamics ------------------------------------------------
  JointState joint_state;
  auto proc_noise_rv = MakePtr<FilterProcessNoiseFunction>(ProcessNoiseLinearPosVel(sigma_a));
  auto proc_noise_clk = MakePtr<FilterProcessNoiseFunction>(ProcessNoiseClock(clk_model, 2));
  joint_state.PushBackStateAndDynamics(rv_state_mci, dyn_est, proc_noise_rv);
  joint_state.PushBackStateAndDynamics(clk_state, dyn_clk_est, proc_noise_clk);
  FilterDynamicsFunction joint_dynamics = joint_state.GetFilterDynamicsFunction();

  /*********************************************
   * Define Measurement function
   * *******************************************/
  FilterMeasurementFunction meas_func_pos_clk
      = [moon_sat, receiver, state_size, meas_types, signals, sis_ure_std, sis_ure_rate_std,
         use_range, use_range_rate](const VecX x, MatXd* H, MatXd* R) -> VecX {
    // Measurements
    std::string freq = "L1";
    double epoch = moon_sat->GetEpoch().val();
    auto measall = receiver->GetMeasurement(epoch);  // measurements of all frequencies
    auto meas_L1 = measall.ExtractSignal(signals);   // measurements of L1
    auto meas = meas_L1.ApplyIonoMask();             // apply ionosphere mask
    int sat_num = meas.GetTrackedSignalNum();        // number of tracked GPS satellites
    int mtot = sat_num * meas_types.size();          // total number of measurements

    // Predict measurements
    *R = MatXd::Zero(mtot, mtot);
    VecX x_N = VecX::Zero(mtot);  // a dummy variable for carrier phase
    Frame frame_in = Frame::MOON_CI;
    VecX z = meas.GetPredictedGnssMeasurement(epoch, x.head(6), x.tail(2), x_N, meas_types,
                                              frame_in, H);  // Jacobian with autodiff

    int n_meas_sat = int(z.size() / meas_types.size());
    VecXd noise_std_vec = meas.GetGnssNoiseStdVec(meas_types).cast<double>();

    int z_idx = 0;
    if (use_range) {
      for (int idx = 0; idx < n_meas_sat; idx++) {
        (*R)(z_idx, z_idx) = std::pow(noise_std_vec(z_idx), 2) + std::pow(sis_ure_std, 2);
        z_idx++;
      }
    }
    if (use_range_rate) {
      for (int idx = 0; idx < n_meas_sat; idx++) {
        (*R)(z_idx, z_idx) = std::pow(noise_std_vec(z_idx), 2) + std::pow(sis_ure_rate_std, 2);
        z_idx++;
      }
    }

    // scaling the measurement noise
    double scale = 1000;
    // R->diagonal().array() *= scale;

    return z;
  };

  /*********************************************
   * Define Process Noise function
   * *******************************************/
  FilterProcessNoiseFunction proc_noise_func = joint_state.GetFilterProcessNoiseFunction();

  /*************************************
   * EKF Setup
   * ***********************************/
  EKF ekf;
  ekf.SetDynamicsFunction(joint_dynamics);
  ekf.SetMeasurementFunction(meas_func_pos_clk);
  ekf.SetProcessNoiseFunction(proc_noise_func);
  ekf.SetAdaptiveProcessNoise(use_adaptive_proc);
  ekf.SetAdaptiveProcessNoiseCoeff(alpha_Q);
  ekf.SetOutlierThreshold(outlier_thresh);

  // Storage
  MatXd error_mat(4, Nt);
  VecXd N_meas(Nt);

  // Initilization
  VecX x_est = SampleMVN(joint_state.GetJointStateValue(), P0, 1, seed);
  ekf.Initialize(t0_tai, x_est, P0);
  VecXd est_err = ComputeEstimationErrors(moon_sat, &ekf);
  error_mat.col(0) = est_err;

  // Output
  std::string constellation_config = "_GPS";
  if (use_galileo) {
    constellation_config += "_GALILEO";
  }
  if (use_qzss) {
    constellation_config += "_QZSS";
  }

  std::string clock_str;
  switch (clk_model) {
    case ClockModel::kMiniRafs: clock_str = "MiniRafs"; break;
    case ClockModel::kMicrosemiCsac: clock_str = "Csac"; break;
    case ClockModel::kRafs: clock_str = "Rafs"; break;
    default: break;
  }

  std::string datafilename = "ExampleEKF" + constellation_config;
  auto output_path = std::filesystem::current_path() / "output" / datafilename / clock_str;

  /***********************************************
   * Main loop
   **********************************************/
  double t = 0;
  double t_tai = t0_tai;
  int t_idx = 0;
  int num_sat = 0;
  // tf = 50 * dt_sim;

  // Compute Estimation
  est_err = ComputeEstimationErrors(moon_sat, &ekf);  // pos, vel, clkb, clkd error

  PrintProgressHeader();
  PrintProgress(t, tf, est_err, ekf.GetCovariance(), num_sat, 0);

  while (t_tai < tf_tai) {
    t += dt_sim;
    t_tai += dt_sim;
    t_idx += 1;

    // moon_sat->Propagate(t_tai);
    // gnss_const.Propagate(t_tai);

    VecX z_true;
    auto measall = receiver->GetMeasurement(t_tai);
    auto meas_L1 = measall.ExtractSignal(signals);
    auto meas = meas_L1.ApplyIonoMask();
    num_sat = meas.GetTrackedSignalNum();
    N_meas(t_idx) = num_sat;

    bool with_noise = true;
    z_true = meas.GetGnssMeasurement(meas_types, with_noise, seed);
    if (use_range) z_true(0) += SampleRandNormal(0, sis_ure_std, seed);
    if (use_range_rate) z_true(1) += SampleRandNormal(0, sis_ure_rate_std, seed);

    ekf.Predict(t_tai);
    ekf.Update(z_true);

    est_err = ComputeEstimationErrors(moon_sat, &ekf);
    error_mat.col(t_idx) = est_err;

    int num_used_meas = ekf.GetMeasurementResidual().size();
    if (fmod(t, dt_print) < 1e-3) {
      PrintProgress(t, tf, est_err, ekf.GetCovariancePost(), num_sat, num_used_meas);
    }
  }

  PrintEstimationStatistics(N_meas, error_mat, double(1.0 / N_orbit));  // use last 30%
}
