/*
* ex_lunar_relativistic.cc
* Tests relativistic clock correction for Lunar dynamics
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

    //initialize ground station state
    SurfaceStaticDynamics ground_observer(NaifId::MOON, frame);
    Vec6 rv_observer = {0.0, 0.0, -R_MOON, 0.0, 0.0, 0.0};
    clock_orbit.SetObserverState_COD(rv_observer);

    //start with initial orbital elements and time parameters
    Real t0 = Gregorian2Time(2025, 11, 9, 0, 0, 0);
    Real num_days = 30;
    Real tf = t0 + num_days * SECS_DAY;
    Real dt = 10.0;
    Real t_span = tf - t0;
    int num_steps = static_cast<int>(t_span / dt) + 1;
    clock_orbit.SetTimeStep(dt);

    //LLO
    // Real a = R_MOON + 100.0;
    // Real e = 0.0;
    // Real i = 90.0;
    // Real Omega = 0.0;
    // Real omega = 0.0;
    //ELFO - fill these out
    Real a = 6541.4;
    Real e = 0.6;
    Real i = 90.0;
    Real Omega = 0.0;
    Real omega = 90.0;
    //NRHO - fill these out
    // Real a = 0.0;
    // Real e = 0.0;
    // Real i = 0.0;
    // Real Omega = 0.0;
    // Real omega = 0.0;

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
    MatX observer_state_history(num_steps, 6);
    observer_state_history.row(0) = rv_observer;
    Real t = t0;
    for (int i = 1; i < num_steps; i++) {
        VecX x = state_history.row(i-1);
        VecX x_observer = observer_state_history.row(i-1);
        state_history.row(i) = clock_orbit.Propagate(x, t, t + dt);
        observer_state_history.row(i) = ground_observer.Propagate(x_observer, t, t + dt);
        clock_orbit.SetObserverState_COD(observer_state_history.row(i));
        t += dt;
    }

    //show results
    cout << "Time correction from velocity time dilation: " << (t_span - state_history(num_steps-1, 6))/num_days << "s/day" << endl;
    cout << "Time correction from gravitational time dilation: " << (t_span - state_history(num_steps-1, 7))/num_days << "s/day" << endl;

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
