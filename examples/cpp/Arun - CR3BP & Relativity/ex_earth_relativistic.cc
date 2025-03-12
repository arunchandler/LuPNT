/*
* ex_lunar_relativistic.cc
* Tests relativistic clock correction for Earth dynamics
*/

#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {

    //initialize dynamics
    ClockOrbitDynamics clock_orbit;
    Frame frame = Frame::ECI;
    clock_orbit.SetFrame_COD(frame);
    clock_orbit.AddBody_COD(Body::Earth());

    //initialize ground station state
    SurfaceStaticDynamics ground_observer(NaifId::EARTH, frame);
    Vec3 latlonalt_observer = {37.428230, -122.168861, 0.0};
    Vec3 xyz_observer = LatLonAlt2Cart(latlonalt_observer, R_EARTH, 0.0);
    Vec6 rv_observer = {xyz_observer[0], xyz_observer[1], xyz_observer[2], 0.0, 0.0, 0.0};
    clock_orbit.SetObserverState_COD(rv_observer);

    //start with initial orbital elements and time parameters
    Real t0 = Gregorian2Time(2025, 11, 9, 0, 0, 0);
    Real num_days = 30;
    Real tf = t0 + num_days * SECS_DAY;
    Real dt = 100.0;
    Real t_span = tf - t0;
    int num_steps = static_cast<int>(t_span / dt) + 1;
    clock_orbit.SetTimeStep(dt);

    //GPS
    // Real a = 26561.8;
    // Real e = 0.02;
    // Real i = 55.0;
    // Real Omega = 0.0;
    // Real omega = 0.0;
    //Galileo
    // Real a = 29994;
    // Real e = 0.02;
    // Real i = 56.0;
    // Real Omega = 0.0;
    // Real omega = 0.0;
    //GLONASS
    Real a = 25510;
    Real e = 0.02;
    Real i = 64.8;
    Real Omega = 0.0;
    Real omega = 0.0;

    Real M0 = 0.0;
    Vec6 elements = {a, e, i, Omega, omega, M0};
    Vec6 state0 = Classical2Cart(elements, GM_EARTH);
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
        Vec3 r_Earth = GetBodyPos(t, NaifId::EARTH, frame);
        PlotBody(NaifId::EARTH, r_Earth, 0);
        xlabel("X [km]");
        ylabel("Y [km]");
        zlabel("Z [km]");
        SetLim(a, 0);
        show();
    }
}
