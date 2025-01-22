#include <lupnt/lupnt.h>

using namespace lupnt;

void PrintProgressHeader() {
  std::cout << " " << std::endl;
  std::cout << std::string(80, '-') << std::endl;
  std::cout << "Time [min]  | Pos Err [m]         | Vel Err [mm/s]         | Clk Bias Err [m]      "
               "   | Num GPS Tracked"
            << std::endl;
  std::cout << std::string(80, '-') << std::endl;
}

VecXd ComputeEstimationErrors(const Ptr<Spacecraft> sat, EKF* ekf) {
  auto x_est = ekf->GetState();
  auto x_true = sat->GetStateVec();

  double r_err = 1000 * (x_true.segment(0, 3) - x_est.segment(0, 3)).norm().val();
  double v_err = 1e6 * (x_true.segment(3, 3) - x_est.segment(3, 3)).norm().val();
  double b_err = 1e3 * abs((x_true(6) - x_est(6)).val());
  double d_err = 1e3 * abs((x_true(7) - x_est(7)).val());

  VecXd x_err(4);
  x_err << r_err, v_err, b_err, d_err;
  return x_err;
}

void PrintProgress(Real t, Vec4 x_err, MatXd P, int num_sat, int num_used_meas) {
  std::cout.precision(5);

  int N_sigma = 3;
  Real sigma_r = 1e3 * sqrt(P(0, 0) + P(1, 1) + P(2, 2));  // [m]
  Real sigma_v = 1e6 * sqrt(P(3, 3) + P(4, 4) + P(5, 5));  // [mm/s]
  Real sigma_b = 1e3 * sqrt(P(6, 6));                      // [m]

  std::cout << std::left << std::setw(12) << t / 60 << " ";
  std::cout << std::left << std::setw(8) << x_err(0);
  std::cout << " (" << std::setw(7) << N_sigma * sigma_r << ")     ";
  std::cout << std::left << std::setw(10) << x_err(1);
  std::cout << " (" << std::setw(7) << N_sigma * sigma_v << ")     ";
  std::cout << std::left << std::setw(12) << x_err(2);
  std::cout << " (" << std::setw(7) << N_sigma * sigma_b << ")     ";
  std::cout << std::left << std::setw(4) << num_sat;
  std::cout << " (used meas:" << num_used_meas << ")";
  std::cout << std::endl;
};

/**
 * @brief Print Estimation Errors
 *
 * @param num_meas (n_time,)   Number of GPS measurements
 * @param error_mat (4, n_time)  Error Mat (Position, Velocity, Clock Bias,
 * Clock Drift)
 */
void PrintEstimationStatistics(VecXd num_meas, MatXd error_mat, double data_ratio = 1.0) {
  int n_time = num_meas.size();
  Vec4d rms, means, stds, p68, p95, p99;

  if (error_mat.rows() != 4) {
    std::cout << "Wrong Mat Size, Error Mat size must be (4 x timestep)" << std::endl;
    return;
  }

  // extract statistics range data
  int start_idx = (int)((1.0 - data_ratio) * n_time);
  int end_idx = n_time - 1;
  int n_range = end_idx - start_idx;

  VecXd num_meas_range(n_range);
  MatXd error_mat_range(4, n_range);

  num_meas_range = num_meas.segment(start_idx, n_range);
  error_mat_range = error_mat.block(0, start_idx, 4, n_range);

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

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "< Simulation Statistics (Last " << data_ratio * 100 << "%)>" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Statistics  | Position [m]  | Velocity [mm/s] | Clock Bias [ns] "
               "| Clk Drift [ns/s] "
            << std::endl;
  std::cout << std::string(80, '-') << std::endl;

  std::cout.precision(5);
  std::cout << "RMS         | " << std::left << std::setw(16) << rms(0) << "  " << std::left
            << std::setw(16) << rms(1) << "   " << std::left << std::setw(16) << rms(2) << std::left
            << std::setw(16) << rms(3) << std::endl;
  std::cout << "Mean+-Std   | " << std::left << means(0) << "+-" << std::left << stds(0) << "  "
            << std::left << means(1) << "+-" << std::left << stds(1) << "   " << std::left
            << means(2) << "+-" << std::left << stds(2) << "   " << std::left << means(3) << "+-"
            << std::left << stds(3) << std::endl;
  std::cout << "68%         | " << std::left << std::setw(16) << p68(0) << "  " << std::left
            << std::setw(16) << p68(1) << "   " << std::left << std::setw(16) << p68(2) << std::left
            << std::setw(16) << p68(3) << std::endl;
  std::cout << "95%         | " << std::left << std::setw(16) << p95(0) << "  " << std::left
            << std::setw(16) << p95(1) << "   " << std::left << std::setw(16) << p95(2) << std::left
            << std::setw(16) << p95(3) << std::endl;
  std::cout << "99%         | " << std::left << std::setw(16) << p99(0) << "  " << std::left
            << std::setw(16) << p99(1) << "   " << std::left << std::setw(16) << p99(2) << std::left
            << std::setw(16) << p99(3) << std::endl;
  std::cout << "  " << std::endl;
}

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
  double tf_tai = t0_tai + N_orbit * period;
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

    // ADD signal in space URE
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
  VecXd num_meas(Nt);

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
  Real t_tai = t0_tai;
  int time_index = 0;
  int num_sat = 0;
  // tf = 50 * dt_sim;

  // Compute Estimation
  est_err = ComputeEstimationErrors(moon_sat, &ekf);  // pos, vel, clkb, clkd error

  PrintProgressHeader();
  PrintProgress((t_tai - t0_tai), est_err, ekf.GetCovariance(), num_sat, 0);

  while (t_tai < tf_tai) {
    t_tai += dt_sim;
    time_index += 1;

    moon_sat->Propagate(t_tai);
    gnss_const.Propagate(t_tai);

    VecX z_true;
    auto measall = receiver->GetMeasurement(t_tai);
    auto meas_L1 = measall.ExtractSignal(signals);
    auto meas = meas_L1.ApplyIonoMask();
    num_sat = meas.GetTrackedSignalNum();
    num_meas(time_index) = num_sat;

    bool with_noise = true;
    z_true = meas.GetGnssMeasurement(meas_types, with_noise, seed);
    if (use_range) z_true(0) += SampleRandNormal(0, sis_ure_std, seed);
    if (use_range_rate) z_true(1) += SampleRandNormal(0, sis_ure_rate_std, seed);

    ekf.Predict(t_tai);
    ekf.Update(z_true);

    est_err = ComputeEstimationErrors(moon_sat, &ekf);
    error_mat.col(time_index) = est_err;

    int num_used_meas = ekf.GetMeasurementResidual().size();
    if (fmod((t_tai - t0_tai).val(), dt_print) < 1e-3) {
      PrintProgress((t_tai - t0_tai), est_err, ekf.GetCovariancePost(), num_sat, num_used_meas);
    }
  }

  PrintEstimationStatistics(num_meas, error_mat, double(1.0 / N_orbit));  // use last 30%
}
