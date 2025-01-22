#include <lupnt/lupnt.h>

#include <chrono>

using namespace lupnt;
using namespace std::chrono;

int main() {
  // Create dynamics
  // auto dyn= MakePtr<CartesianTwoBodyDynamics>(GM_MOON, IntegratorType::RK4);

  auto dyn = MakePtr<NBodyDynamics<Real>>(IntegratorType::RKF45);
  dyn->SetIntegratorParams(IntegratorParams(20, 1e-12, 1e-12));
  dyn->SetFrame(Frame::MOON_CI);
  dyn->AddBody(BodyT<Real>::Moon(10, 10));
  dyn->AddBody(BodyT<Real>::Earth());
  dyn->SetTimeStep(1.0);

  // Initial state
  Real a = 6541.4;
  Real e = 0.6;
  Real i = 65.5 * RAD;
  Real Omega = 0.0 * RAD;
  Real w = 90.0 * RAD;
  Real M = 0.0 * RAD;

  // Time
  Real et0_utc = Gregorian2Time(2025, 1, 1, 12, 0, 0).val();  // in UTC
  double et0 = UTC2TAI(et0_utc).val();                        // in TAI
  double dt = 30 * 60;                                        // Integration time step [s]

  // Initial state
  ClassicalOE coe_moon({a, e, i, Omega, w, M}, Frame::MOON_CI);
  CartesianOrbitState cart_op = Classical2Cart(coe_moon, GM_MOON);

  VecX x0 = cart_op.GetVec();

  dyn->SetTimeStep(1.0);

  MatXd Phi = MatXd::Zero(x0.size(), x0.size());
  VecX x_next = dyn->Propagate(x0, et0, et0 + dt, &Phi);

  // Compare with numerical Derivative
  VecX x_pert_plus = VecX::Zero(x0.size());
  VecX x_pert_minus = VecX::Zero(x0.size());
  VecX dx_pert = VecX::Zero(x0.size());
  MatXd J_num = MatXd::Zero(x0.size(), x0.size());
  MatXd F_dum_ = MatXd::Zero(x0.size(), x0.size());
  double eps = 1e-6;
  for (int i = 0; i < x0.size(); i++) {
    x_pert_plus = x0;
    x_pert_minus = x0;
    x_pert_plus(i) = x0(i) + eps;
    x_pert_minus(i) = x0(i) - eps;
    dx_pert = dyn->Propagate(x_pert_plus, et0, et0 + dt, &F_dum_)
              - dyn->Propagate(x_pert_minus, et0, et0 + dt, &F_dum_);
    J_num.col(i) = dx_pert.cast<double>() / (2 * eps);
  }

  std::cout << "<Dynamics Example>" << std::endl;
  std::cout << " Initial State: " << x0.transpose() << std::endl << std::endl;
  std::cout << " Propagated State: " << x_next.transpose() << std::endl << std::endl;

  std::cout << " Jacobian: " << std::endl << Phi << std::endl << std::endl;
  std::cout << " Numerical Jacobian: " << std::endl << J_num << std::endl << std::endl;

  // compute the ratio (F_ - J_num) / J_num
  MatXd diff = (Phi - J_num).array() / J_num.array();
  std::cout << " Difference in STM (ratio): " << std::endl << diff << std::endl << std::endl;

  return 0;
}
