#include <lupnt/lupnt.h>

#include <chrono>

using namespace lupnt;
using namespace std::chrono;

int main() {
  // Time
  Real t0_utc = Gregorian2Time(2020, 1, 1, 12, 0, 0);
  Real t0_tai = ConvertTime(t0_utc, Time::UTC, Time::TAI);
  Real tf_tai = t0_tai + 90 * SECS_DAY;
  Real Dt = 10;
  int N_steps = static_cast<int>((tf_tai - t0_tai) / Dt);
  VecX tfs = arange(t0_tai, tf_tai + Dt, Dt);
  std::cout << N_steps << " steps" << std::endl;

  // Initial state
  Real a = 6541.4;      // [km] Semi-major axis
  Real e = 0.6000;      // [--] Eccentricity
  Real i = 56.2 * RAD;  // [deg] Inclination
  Real O = 0.00 * RAD;  // [deg] Right ascension of the ascending node
  Real w = 90.0 * RAD;  // [deg] Argument of perigee
  Real M = 0.00 * RAD;  // [deg] Mean anomaly
  Vec6 coe0_op(a, e, i, O, w, M);
  Vec6 rv0_op = Classical2Cart(coe0_op, GM_MOON);
  Vec6 rv0_ci = ConvertFrame(t0_tai, rv0_op, Frame::MOON_OP, Frame::MOON_CI);
  Vec6 rv0_gcrf = ConvertFrame(t0_tai, rv0_ci, Frame::MOON_CI, Frame::GCRF);

  // Integrator
  IntegratorParams params;
  params.max_iter = 10;
  params.abstol = 1e-8;
  params.reltol = 1e-8;

  // Dynamics
  NBodyDynamics dyn(IntegratorType::RK4);
  dyn.AddBody(Body::Moon(20, 20));
  dyn.AddBody(Body::Earth());
  dyn.AddBody(Body::Sun());
  dyn.SetFrame(Frame::GCRF);
  dyn.SetIntegratorParams(params);
  dyn.SetTimeStep(10);

  // Propagate
  auto start = high_resolution_clock::now();
  MatX6 rv = dyn.Propagate(rv0_gcrf, t0_tai, tfs, true);
  auto end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(end - start);

  std::cout << "Elapsed time: " << duration.count() / 1e6 << " s" << std::endl;
}
