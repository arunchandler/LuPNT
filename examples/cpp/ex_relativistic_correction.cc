/*
* ex_relativistic_correction.cc
* Tests relativistic clock correction for N-body dynamics
*/

#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {

    //initialize dynamics
    ClockOrbitDynamics clock_orbit;
    Frame frame = Frame::MOON_CI;
    clock_orbit.SetFrame_COD(frame);
    clock_orbit.AddBody_COD(Body::Moon());

    //start with initial orbital elements and time parameters
    Real t0 = 0.0;
    Real tf = 10.0;
    Real dt = 0.1;
    Real init_t_err = 0.0;
    int num_steps = static_cast<int>((tf - t0) / dt) + 1;
    clock_orbit.SetTimeStep(dt);

    Real a = R_MOON + 100.0;
    Real e = 0.0;
    Real i = 0.0;
    Real Omega = 0.0;
    Real omega = 0.0;
    Real M0 = 0.0;
    Vec6 elements = {a, e, i, Omega, omega, M0};
    Vec6 state0 = Classical2Cart(elements, GM_MOON);
    VecX init_state(7);
    init_state.head<6>() = state0;
    init_state[6] = init_t_err;

    //propagate state and timing error
    MatX state_history(num_steps, 7);
    state_history.row(0) = init_state;
    Real t = t0;
    for (int i = 1; i < num_steps; i++) {
        VecX x = state_history.row(i-1);
        state_history.row(i) = clock_orbit.Propagate(x, t, t + dt);
        t += dt;
    }

    //plot results
    figure();
    hold(on);
    VecX x_vals = state_history.col(0);
    VecX y_vals = state_history.col(1);
    VecX z_vals = state_history.col(2);
    Plot3(x_vals, y_vals, z_vals, "b", 0);
    Vec3 r_Moon = GetBodyPos(t, NaifId::MOON, frame);
    PlotBody(NaifId::MOON, r_Moon, 1);
    xlabel("X [km]");
    ylabel("Y [km]");
    zlabel("Z [km]");
    title("Relativistic Clock Correction for Lunar Orbit");
    show();

}
