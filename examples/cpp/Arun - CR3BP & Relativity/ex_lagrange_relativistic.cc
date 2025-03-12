/*
* ex_lagrange_relativistic.cc
* Tests relativistic corrections at Lagrange points
*/

#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main(){

    //initialize dynamics
    ClockOrbitDynamics clock_orbit;
    Frame frame = Frame::ECI;
    clock_orbit.SetFrame_COD(frame);
    clock_orbit.AddBody_COD(Body::Earth());
    clock_orbit.AddBody_COD(Body::Moon());
    Real mu = GM_MOON / (GM_EARTH + GM_MOON);

    //initialize ground station state
    SurfaceStaticDynamics ground_observer(NaifId::EARTH, frame);
    Vec3 latlonalt_observer = {37.428230, -122.168861, 0.0};
    Vec3 xyz_observer = LatLonAlt2Cart(latlonalt_observer, R_EARTH, 0.0);
    Vec6 rv_observer = {xyz_observer[0], xyz_observer[1], xyz_observer[2], 0.0, 0.0, 0.0};
    clock_orbit.SetObserverState_COD(rv_observer);

    //time parameters
    Real t0 = Gregorian2Time(2025, 11, 9, 0, 0, 0);
    Real num_days = 30;
    Real tf = t0 + num_days * SECS_DAY;
    Real dt = 100.0;
    Real t_span = tf - t0;
    int num_steps = static_cast<int>(t_span / dt) + 1;
    clock_orbit.SetTimeStep(dt);

    //Lagrange point states
    VecX x0_L1(8), x0_L2(8), x0_L3(8), x0_L4(8), x0_L5(8);
    Vec3 r_moon = GetBodyPos(t0, NaifId::MOON, frame);
    Vec3 omega(0.0, 0.0, OMEGA_EARTH_MOON);
    // L1
    Vec3 L1_pos = (0.8446819150716 + mu) * r_moon;
    Vec3 L1_vel = omega.cross(L1_pos);
    x0_L1 << L1_pos(0), L1_pos(1), L1_pos(2), //pos
          L1_vel(0), L1_vel(1), L1_vel(2), //vel
          0.0, 0.0; //time error

    // L2
    Vec3 L2_pos = (1.1495185441167 + mu) * r_moon;
    Vec3 L2_vel = omega.cross(L2_pos);
    x0_L2 << L2_pos(0), L2_pos(1), L2_pos(2), //pos
          L2_vel(0), L2_vel(1), L2_vel(2), //vel
          0.0, 0.0; //time error

    // L3
    Vec3 L3_pos = (-1.0044288297623 + mu) * r_moon;
    Vec3 L3_vel = omega.cross(L3_pos);
    x0_L3 << L3_pos(0), L3_pos(1), L3_pos(2), //pos
          L3_vel(0), L3_vel(1), L3_vel(2), //vel
          0.0, 0.0; //time error

    // L4
    Vec3 L4_pos = 0.5 * r_moon + sqrt(3) / 2 * Vec3(-r_moon(1), r_moon(0), 0.0);
    Vec3 L4_vel = omega.cross(L4_pos);
    x0_L4 << L4_pos(0), L4_pos(1), L4_pos(2), //pos
          L4_vel(0), L4_vel(1), L4_vel(2), //vel
          0.0, 0.0; //time error

    // L5
    Vec3 L5_pos = 0.5 * r_moon - sqrt(3) / 2 * Vec3(-r_moon(1), r_moon(0), 0.0);
    Vec3 L5_vel = omega.cross(L5_pos);
    x0_L5 << L5_pos(0), L5_pos(1), L5_pos(2), //pos
          L5_vel(0), L5_vel(1), L5_vel(2), //vel
          0.0, 0.0; //time error

    //propagate state and timing error
    MatX state_history_L1(num_steps, 8);
    MatX state_history_L2(num_steps, 8);
    MatX state_history_L3(num_steps, 8);
    MatX state_history_L4(num_steps, 8);
    MatX state_history_L5(num_steps, 8);
    
    state_history_L1.row(0) = x0_L1;
    state_history_L2.row(0) = x0_L2;
    state_history_L3.row(0) = x0_L3;
    state_history_L4.row(0) = x0_L4;
    state_history_L5.row(0) = x0_L5;

    MatX observer_state_history(num_steps, 6);
    observer_state_history.row(0) = rv_observer;
    Real t = t0;

    for (int i = 1; i < num_steps; i++) {

        VecX x_L1 = state_history_L1.row(i - 1);
        VecX x_L2 = state_history_L2.row(i - 1);
        VecX x_L3 = state_history_L3.row(i - 1);
        VecX x_L4 = state_history_L4.row(i - 1);
        VecX x_L5 = state_history_L5.row(i - 1);
        VecX x_observer = observer_state_history.row(i-1);
        VecX x_L1_temp = clock_orbit.Propagate(x_L1, t, t + dt);
        VecX x_L2_temp = clock_orbit.Propagate(x_L2, t, t + dt);
        VecX x_L3_temp = clock_orbit.Propagate(x_L3, t, t + dt);
        VecX x_L4_temp = clock_orbit.Propagate(x_L4, t, t + dt);
        VecX x_L5_temp = clock_orbit.Propagate(x_L5, t, t + dt);

        state_history_L1.row(i).head(6) = x_L1.head(6);
        state_history_L1.row(i).tail<2>() = x_L1_temp.tail(2);
        state_history_L2.row(i).head(6) = x_L2.head(6);
        state_history_L2.row(i).tail<2>() = x_L2_temp.tail(2);
        state_history_L3.row(i).head(6) = x_L3.head(6);
        state_history_L3.row(i).tail<2>() = x_L3_temp.tail(2);
        state_history_L4.row(i).head(6) = x_L4.head(6);
        state_history_L4.row(i).tail<2>() = x_L4_temp.tail(2);
        state_history_L5.row(i).head(6) = x_L5.head(6);
        state_history_L5.row(i).tail<2>() = x_L5_temp.tail(2);

        observer_state_history.row(i) = ground_observer.Propagate(x_observer, t, t + dt);
        clock_orbit.SetObserverState_COD(observer_state_history.row(i));
        t += dt;
    }

    //show results
    bool partial = false;
    if (partial) {
        cout << "L1 time correction from velocity time dilation: " << (t_span - state_history_L1(num_steps-1, 6))/num_days << "s/day" << endl;
        cout << "L1 time correction from gravitational time dilation: " << (t_span - state_history_L1(num_steps-1, 7))/num_days << "s/day" << endl;
        cout << "L2 time correction from velocity time dilation: " << (t_span - state_history_L2(num_steps-1, 6))/num_days << "s/day" << endl;
        cout << "L2 time correction from gravitational time dilation: " << (t_span - state_history_L2(num_steps-1, 7))/num_days << "s/day" << endl;
        cout << "L3 time correction from velocity time dilation: " << (t_span - state_history_L3(num_steps-1, 6))/num_days << "s/day" << endl;
        cout << "L3 time correction from gravitational time dilation: " << (t_span - state_history_L3(num_steps-1, 7))/num_days << "s/day" << endl;
        cout << "L4 time correction from velocity time dilation: " << (t_span - state_history_L4(num_steps-1, 6))/num_days << "s/day" << endl;
        cout << "L4 time correction from gravitational time dilation: " << (t_span - state_history_L4(num_steps-1, 7))/num_days << "s/day" << endl;
        cout << "L5 time correction from velocity time dilation: " << (t_span - state_history_L5(num_steps-1, 6))/num_days << "s/day" << endl;
        cout << "L5 time correction from gravitational time dilation: " << (t_span - state_history_L5(num_steps-1, 7))/num_days << "s/day" << endl;
    }
    cout << "L1 time correction from both effects: " << (t_span - state_history_L1(num_steps-1, 6) + t_span - state_history_L1(num_steps-1, 7))/num_days << "s/day" << endl;
    cout << "L2 time correction from both effects: " << (t_span - state_history_L2(num_steps-1, 6) + t_span - state_history_L2(num_steps-1, 7))/num_days << "s/day" << endl;
    cout << "L3 time correction from both effects: " << (t_span - state_history_L3(num_steps-1, 6) + t_span - state_history_L3(num_steps-1, 7))/num_days << "s/day" << endl;
    cout << "L4 time correction from both effects: " << (t_span - state_history_L4(num_steps-1, 6) + t_span - state_history_L4(num_steps-1, 7))/num_days << "s/day" << endl;
    cout << "L5 time correction from both effects: " << (t_span - state_history_L5(num_steps-1, 6) + t_span - state_history_L5(num_steps-1, 7))/num_days << "s/day" << endl;

    bool plot = false;
    if (plot) {
        figure();
        hold(on);
        VecX x_vals = state_history_L1.col(0);
        VecX y_vals = state_history_L1.col(1);
        VecX z_vals = state_history_L1.col(2);
        Plot3(x_vals, y_vals, z_vals, "b", 0);
        Vec3 r_Earth = GetBodyPos(t, NaifId::EARTH, frame);
        PlotBody(NaifId::EARTH, r_Earth, 0);
        xlabel("X [km]");
        ylabel("Y [km]");
        zlabel("Z [km]");
        show();
    }

}
