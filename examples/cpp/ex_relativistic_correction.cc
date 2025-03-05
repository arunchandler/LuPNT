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

    //initialize ground station state at south pole of moon
    SurfaceStaticDynamics ground_asset(NaifId::MOON, frame);
    Vec6 rv_asset = {0.0, 0.0, -R_MOON, 0.0, 0.0, 0.0}; //position and velocity
    clock_orbit.SetAssetState_COD(rv_asset);

    //start with initial orbital elements and time parameters
    Real t0 = Gregorian2Time(2025, 11, 9, 0, 0, 0);
    Real tf = t0 + SECS_DAY;
    Real dt = 10.0;
    int num_steps = static_cast<int>((tf - t0) / dt) + 1;
    clock_orbit.SetTimeStep(dt);

    Real a = R_MOON + 100.0;
    Real e = 0.0;
    Real i = 90.0;
    Real Omega = 0.0;
    Real omega = 0.0;
    Real M0 = 0.0;
    Vec6 elements = {a, e, i, Omega, omega, M0};
    Vec6 state0 = Classical2Cart(elements, GM_MOON);
    VecX init_state(8);
    init_state.head<6>() = state0;
    init_state[6] = 0.0; //velocity time dilation
    init_state[7] = 0.0; //gravitational time dilation

    //propagate state and timing error
    MatX state_history(num_steps, 8);
    state_history.row(0) = init_state;

    Real t = t0;
    for (int i = 1; i < num_steps; i++) {
        VecX x = state_history.row(i-1);
        state_history.row(i) = clock_orbit.Propagate(x, t, t + dt);
        rv_asset = ground_asset.Propagate(rv_asset, t, t + dt);
        clock_orbit.SetAssetState_COD(rv_asset);
        t += dt;
    }

    //show results
    cout << "Time Correction from Velocity Time Dilation: " << state_history(num_steps-1, 6) - (tf-t0) << "s" << endl;
    cout << "Time Correction from Gravitational Time Dilation: " << state_history(num_steps-1, 7) - (tf-t0) << "s" << endl;

    bool plot = false;
    if (plot) {
        figure();
        hold(on);
        VecX x_vals = state_history.col(0);
        VecX y_vals = state_history.col(1);
        VecX z_vals = state_history.col(2);
        Plot3(x_vals, y_vals, z_vals, "b", 0);
        Vec3 r_Moon = GetBodyPos(t, NaifId::MOON, frame);
        PlotBody(NaifId::MOON, r_Moon, 0);
        xlabel("X [km]");
        ylabel("Y [km]");
        zlabel("Z [km]");
        title("Relativistic Clock Correction for LLO");
        SetLim(a, 0);
        show();
    }

}
