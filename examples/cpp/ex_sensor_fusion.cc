#include <lupnt/lupnt.h>

#include <opencv2/opencv.hpp>

using namespace lupnt;
using std::cout;
using std::endl;
using std::vector;
using namespace matplot;

int main() {
  // Cache
  bool recompute = false;
  auto output_mpath = GetOutputPath("ex_sensor_fusion");
  auto cache_mpath = output_mpath / "cache.h5";
  H5Easy::File cache_file = GetCacheFile(cache_mpath, recompute);
  cout << "Cache path: " << cache_mpath << endl;

  // Time
  std::string t0_utc_string = "2028-06-25T12:00:00.0";
  Real t0_utc = Gregorian2Time(t0_utc_string);
  Real t0_tai = ConvertTime(t0_utc, Time::UTC, Time::TAI);
  const int Nsc = 2;

  // Classical orbital elements (a, e, i, Ω, ω, M) [km, -, rad, rad, rad, rad]
  Real a = 5740;                      // [km] Semi-Mjor axis
  Real e = 0.58;                      // [-] Eccentricity
  Real i = 54.856 * RAD;              // [rad] Inclination
  Real Omega = 0 * RAD;               // [rad] Right ascension of the ascending node
  Real w = 86.322 * RAD;              // [rad] Argument of periapsis
  Vec<Nsc> M = {0 * RAD, 180 * RAD};  // [rad] Mean anomaly

  // Time
  Real period = GetOrbitalPeriod(a, GM_MOON);  // [s] Orbital period
  Real dt = 10;                                // [s] Simulation time step
  Real dt_prop = 1;                            // [s] Propagation time step
  Real tf = 1 * SECS_DAY;                      // [s] Simulation final time
  VecX tspan = arange(0, tf + dt, dt);         // [s] Time sie first epoch
  VecX tspan_h = tspan / SECS_HOUR;            // [h] Time sie first epoch
  VecX tfs_tai = tspan.array() + t0_tai;       // [s] Epochs (TAI)
  int Nt = tspan.size();                       // [-] Number of time steps

  // Initial conditions
  MatX6 coe0_mop(Nsc, 6);
  for (int j = 0; j < Nsc; j++) coe0_mop.row(j) << a, e, i, Omega, w, M(j);
  MatX6 rv0_m2sc_mop = Classical2Cart(coe0_mop, GM_MOON);
  MatX6 rv0_m2sc_mci = ConvertFrame(t0_tai, rv0_m2sc_mop, Frame::MOON_OP, Frame::MOON_CI);

  // Dynamics
  NBodyDynamics dyn;
  dyn.AddBody(Body::Moon(100, 100));
  dyn.AddBody(Body::Earth());
  dyn.AddBody(Body::Sun());
  dyn.AddBody(Body::Mars());
  dyn.AddBody(Body::Venus());
  dyn.AddBody(Body::Jupiter());
  dyn.SetFrame(Frame::MOON_CI);
  dyn.SetTimeStep(dt_prop);
  dyn.SetPrintProgress(true);

  NBodyDynamics dyn_filter;
  dyn_filter.AddBody(Body::Moon(18, 18));
  dyn_filter.AddBody(Body::Earth());
  dyn_filter.SetFrame(Frame::MOON_CI);
  dyn_filter.SetTimeStep(dt_prop);

  // Propagation
  // rv_from2to_frame [km, km/s] (x, y, z, vx, vy, vz)
  vector<MatX6> rv_m2sc_mci(Nsc), rv_m2sc_mpa(Nsc);
#pragma omp parallel for
  for (int i = 0; i < Nsc; i++) {
    auto name = std::format("/rv_m2sc_mci_{}", i);
    auto func = [&]() { return dyn.Propagate(rv0_m2sc_mci.row(i), t0_tai, tfs_tai); };
    rv_m2sc_mci[i] = LoadOrRecompute<-1, 6, Real>(name, cache_file, recompute, func);
    rv_m2sc_mpa[i] = ConvertFrame(tfs_tai, rv_m2sc_mci[i], Frame::MOON_CI, Frame::MOON_PA);
  }

  // Earth and Sun
  MatX3 r_m2e_mci = GetBodyPos(tfs_tai, NaifId::MOON, NaifId::EARTH, Frame::MOON_CI);
  MatX3 r_m2s_mci = GetBodyPos(tfs_tai, NaifId::MOON, NaifId::SUN, Frame::MOON_CI);
  MatX3 r_e2s_eci = GetBodyPos(tfs_tai, NaifId::EARTH, NaifId::SUN, Frame::ECI);
  MatX3 r_m2e_mpa = ConvertFrame(tfs_tai, r_m2e_mci, Frame::MOON_CI, Frame::MOON_PA);
  MatX3 r_m2s_mpa = ConvertFrame(tfs_tai, r_m2s_mci, Frame::MOON_CI, Frame::MOON_PA);

  // Attitude
  vector<MatX3> e_sc2m(Nsc), e_sc2s(Nsc), e_sc2e(Nsc);
  vector<MatX3> ez_sc(Nsc), ey_sc(Nsc), ex_sc(Nsc);
  vector<MatX3> er_sc(Nsc), et_sc(Nsc), en_sc(Nsc);
  vector<VecX> sun_angle(Nsc);
  vector<vector<Mat3>> R_mci2sc(Nsc), R_mci2rtn(Nsc);
  vector<Mat3> R_mpa2ci(Nt);
  for (int i = 0; i < Nsc; i++) {
    // Moon, Sun, Earth
    e_sc2m[i] = -rv_m2sc_mci[i].leftCols(3).rowwise().normalized();
    e_sc2s[i] = (r_m2s_mci - rv_m2sc_mci[i].leftCols(3)).rowwise().normalized();
    e_sc2e[i] = (r_m2e_mci - rv_m2sc_mci[i].leftCols(3)).rowwise().normalized();

    // Yaw-Steering
    ez_sc[i].resize(Nt, 3);
    ey_sc[i].resize(Nt, 3);
    ex_sc[i].resize(Nt, 3);
    for (int t = 0; t < Nt; t++) {
      ez_sc[i].row(t) = e_sc2m[i].row(t);
      ey_sc[i].row(t) = ez_sc[i].row(t).cross(e_sc2s[i].row(t)).normalized();
      ex_sc[i].row(t) = ey_sc[i].row(t).cross(ez_sc[i].row(t)).normalized();
    }

    // Radial, tangential, normal (RTN)
    MatX3 r_i = rv_m2sc_mci[i].leftCols(3);
    MatX3 v_i = rv_m2sc_mci[i].rightCols(3);
    er_sc[i].resize(Nt, 3);
    et_sc[i].resize(Nt, 3);
    en_sc[i].resize(Nt, 3);
    for (int t = 0; t < Nt; t++) {
      er_sc[i].row(t) = r_i.row(t).normalized();
      en_sc[i].row(t) = r_i.row(t).cross(v_i.row(t)).normalized();
      et_sc[i].row(t) = en_sc[i].row(t).cross(er_sc[i].row(t)).normalized();
    }

    // Sun angle
    sun_angle[i] = e_sc2s[i] * e_sc2m[i].transpose();

    // Rotation Mtrices
    R_mci2sc[i].resize(Nt);
    R_mci2rtn[i].resize(Nt);
    for (int t = 0; t < Nt; t++) {
      R_mci2sc[i][t] << ex_sc[i].row(t), ey_sc[i].row(t), ez_sc[i].row(t);
      R_mci2rtn[i][t] << er_sc[i].row(t), et_sc[i].row(t), en_sc[i].row(t);

      auto [R_tmp, t_tmp] = GetFrameRotationTranslation(tfs_tai[t], Frame::MOON_CI, Frame::MOON_PA);
      R_mpa2ci[t] = R_tmp;
    }
  }

  // Spacecraft frame
  Vec3 cam_dir = {0, 0, 1};
  Vec3 cam_up = {1, 0, 0};
  Vec3 cam_right = {0, 1, 0};
  Vec3 ez_c = cam_dir;
  Vec3 ey_c = -cam_up;
  Vec3 ex_c = cam_right;
  Mat3 R_sc2ocv;
  R_sc2ocv << ex_c.transpose(), ey_c.transpose(), ez_c.transpose();

  // Plot
  auto fig = figure(true);
  title("Case 0: Moon orbit in CI frame");
  hold(true);
  PlotBody(NaifId::MOON);                                                                   // Moon
  Plot3(rv_m2sc_mci[0], "")->line_width(2).marker_indices({0});                             // Orbit
  PlotArrow3(Vec3::Zero(), r_m2e_mci.row(0).normalized() * 3 * R_MOON, "")->line_width(2);  // Earth
  PlotArrow3(Vec3::Zero(), r_m2s_mci.row(0).normalized() * 3 * R_MOON, "")->line_width(2);  // Sun
  for (int t = 0; t < int(period / dt); t += int(period / dt / 10))
    PlotFrame(rv_m2sc_mci[0].row(t).head(3), R_mci2sc[0][t] * 1 * R_MOON);  // Attitude
  legend({"Moon", "Orbit", "Sun", "Earth", "x_{sc}", "y_{sc}", "z_{sc}"});
  SetLim(12e3);
  fig->draw();
  fig->show();
}
