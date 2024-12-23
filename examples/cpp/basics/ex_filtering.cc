#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {
  auto output_path = GetOutputPath("ex_filtering");
  cout << "Output path: " << output_path << endl;

  // Parameters
  bool use_autodiff = true;
  const Real Dt = 1;                             // Time step
  const int N_t = 100;                           // Number of steps
  Vec3 x{0, 0, 0};                               // State
  Mat3d Q = Vec3d(0.1, 0.1, 0.01).asDiagonal();  // Process noise
  Vec2 u{5, 0.1};                                // Control input
  Vec<N_t> tspan = arange(0, N_t * Dt, Dt);      // Time span

  // Stations
  const int N_sta = 4;
  Matd<N_sta, 2> pos{{-120, -120}, {-120, +120}, {+120, -120}, {+120, +120}};
  Matd<N_sta, N_sta> R = 10 * Matd<N_sta, N_sta>::Identity();

  // Measurement function
  auto g = [N_sta, pos](Vec3 x) -> Vec<N_sta> {
    Vec<N_sta> z;
    for (int i = 0; i < N_sta; i++) z(i) = (pos.row(i).transpose() - x.head(2)).norm();

    return z;
  };

  // Dynamics function
  auto f = [](Vec3 x, Vec2 u, Real Dt) -> Vec3 {
    return x + Dt * Vec3(u(0) * cos(x(2)), u(0) * sin(x(2)), u(1));
  };

  // Measurement model
  auto f_meas
      = [pos, R, g, use_autodiff](VecX x, MatXd *H = nullptr, MatXd *R_in = nullptr) -> VecX {
    VecX z;
    if (H == nullptr) {
      z = g(x);
    } else {
      if (use_autodiff) {
        VecX x_tmp = x.cast<double>();
        *H = jacobian(g, wrt(x_tmp), at(x_tmp), z);
      } else {
        Mat<N_sta, 3> H_tmp;
        H_tmp.col(0) = -(pos.col(0).array() - x(0)) / z.array();
        H_tmp.col(1) = -(pos.col(1).array() - x(1)) / z.array();
        H_tmp.col(2).setZero();
        *H = H_tmp.cast<double>();
      }
    }
    if (R_in != nullptr) *R_in = R;
    return z;
  };

  // Dynamics model
  auto f_dyn = [u, f, use_autodiff](VecX x, Real t_curr, Real t_end, MatXd *F) -> VecX {
    VecX xf;
    if (F == nullptr) {
      xf = f(x, u, t_end - t_curr);
    } else {
      if (use_autodiff) {
        VecX x_tmp = x.cast<double>();
        *F = jacobian(f, wrt(x_tmp), at(x_tmp, u, t_end - t_curr), xf);
      } else {
        xf = f(x, u, t_end - t_curr);
        Real Dt = t_end - t_curr;
        MatX F_tmp{{1, 0, -Dt * u[0] * sin(x[2])}, {0, 1, Dt * u[0] * cos(x[2])}, {0, 0, 1}};
        *F = F_tmp.cast<double>();
      }
    }
    return xf;
  };

  // Process noise
  auto f_noise = [Q](VecX x, Real t_curr, Real t_end) -> MatXd { return Q; };

  // Filter
  Vec3 x0 = SampleMVN(x, Q, 1);  // Initial estimate
  Mat3d P0 = Q;                  // Initial covariance estimate

  filters::EKF ekf;
  ekf.SetDynamicsFunction(f_dyn);
  ekf.SetMeasurementFunction(f_meas);
  ekf.SetProcessNoiseFunction(f_noise);
  ekf.Initialize(0, x0, Q);

  // History
  Mat<N_t, 3> x_hist;
  Mat<N_t, 3> x_pred_hist;
  Mat<N_t, 3> x_updt_hist;
  Mat<N_t, 3> P_hist;
  Mat<N_t, N_sta> z_hist;

  // Noise
  int seed = 1234;
  default_random_engine rng(seed);
  normal_distribution<double> noise(0, 1);
  // auto mvr = Eigen::Rand::makeMvNformalGen(mean, cov);

  // Propagate
  Real t = 0;
  for (int it = 0; it < N_t; it++) {
    t += Dt;
    // Vec2 w = gen.generate(rng);

    // Propagate
    x += Dt * Vec3(u(0) * cos(x(2)), u(0) * sin(x(2)), u(1));
    for (int i = 0; i < 3; i++) x(i) += sqrt(Q(i, i)) * noise(rng);

    // Measurements
    VecX z = g(x);
    for (int i = 0; i < N_sta; i++) z(i) += sqrt(R(i, i)) * noise(rng);

    ekf.Predict(t);
    ekf.Update(z);

    x_hist.row(it) = x;
    x_pred_hist.row(it) = ekf.GetPredictedStateEstimate();
    x_updt_hist.row(it) = ekf.GetUpdatedStateEstimate();
    P_hist.row(it) = ekf.GetPostCov().diagonal();
    z_hist.row(it) = z;
  }

  // Plot
  line_handle p;
  figure_handle fig;

  // Trajectory
  fig = figure(true);
  hold(true);
  Plot(x_hist.col(0), x_hist.col(1));
  Plot(x_updt_hist.col(0), x_updt_hist.col(1));
  xlabel("x [m]");
  ylabel("y [m]");
  grid(true);
  xlim({-120, 120});
  ylim({-50, 120});
  fig->save(output_path / "trajectory.png");

  // Measurement residual
  Mat<N_t, N_sta> res_prior, res_post;
  for (int it = 0; it < N_t; it++) {
    res_prior.row(it) = z_hist.row(it) - g(x_pred_hist.row(it)).transpose();
    res_post.row(it) = z_hist.row(it) - g(x_updt_hist.row(it)).transpose();
  }
  fig = figure(true);
  hold(true);
  for (int i = 0; i < N_sta; i++) Plot(tspan, res_prior.col(i), "bo");
  for (int i = 0; i < N_sta; i++) Plot(tspan, res_post.col(i), "rx");
  xlabel("t [s]");
  ylabel("Resiudal [m]");
  grid(true);
  fig->save(output_path / "residuals.png");

  // Errors
  double N_sig = 3;
  Mat<N_t, 3> error = x_hist - x_updt_hist;
  Mat<N_t, 3> error_std = sqrt(P_hist.array());
  error.col(2) *= DEG;
  error_std.col(2) *= DEG;

  fig = figure(true);
  vector<string> ylabels = {"x error [m]", "y error [m]", "theta error [deg]"};
  for (int i = 0; i < 3; i++) {
    fig->add_subplot(3, 1, i);
    hold(true);
    Plot(tspan, error.col(i));
    Plot(tspan, +N_sig * error_std.col(i));
    Plot(tspan, -N_sig * error_std.col(i));
    xlabel("t [s]");
    ylabel(ylabels[i]);
    grid(true);
  }
  fig->save(output_path / "errors.png");

  cout << "Done" << endl;
  return 0;
}
