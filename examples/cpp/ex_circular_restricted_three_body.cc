/*
* ex_circular_restricted_three_body.cc
* Tests dynamics of CR3BP around the Earth-Moon Lagrange points
* Position is normalized by the Earth-Moon distance
* Velocity is normalized by the Earth-Moon distance and the Earth-Moon gravitational parameter
* Time is normalized by the Earth-Moon distance and the Earth-Moon gravitational parameter
*/

#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main() {

    IntegratorType integ = IntegratorType::RK4;
    Vec2 GMs = {Body::Earth().GM, Body::Moon().GM};
    CR3BPDynamics cr3bp(GMs, integ);
    Real mu;
    if (GMs[1] < GMs[0]) {
        mu = GMs[1] / (GMs[0] + GMs[1]);
    } else {
        mu = GMs[0] / (GMs[0] + GMs[1]);
    }
//     RelativisticClockDynamics time_error(integ);

    //initial state - uncomment one of the following initial conditions
    VecX x0(7);
    // Test
    // x0 << 1.51805, 0.0, 0.0, //pos
    //       0.0, -1.0, 0.0, //vel
    //       0.0; //time error

    // L1 Static

    // L2 Static

    // L3 Static

    // L4 Static

    // L5 Static

    // L1 Axial
    // x0 << 0.913013, 0.0, 0.0, //pos
    //       0.0, -0.313025, 0.434375, //vel
    //      0.0; //time error

    // L1 Vertical
    // x0 << 0.861879, 0.0, 0.0, //pos
    //       0.0, 0.085903, 0.432614, //vel
    //      0.0; //time error 

    // L2 Northern Butterfly
    // x0 << 1.038394, 0.0, 0.173741, //pos
    //       0.0, -0.078548, 0.0, //vel
    //       0.0; //time error

    // L3 Lyaupunov
    x0 << -0.463824, 0.0, 0.0, //pos
          0.0, -1.388737, 0.0, //vel
          0.0; //time error

    // L4 Short Period
    // x0 << 0.416475, 0.866025, 0.0, //pos
    //       -0.045831, 0.054473, 0.0, //vel
         //  0.0; //time error

    // L4 Axial
    // x0 << 0.570828, 0.522879, 0.0, //pos
    //       -0.015379, -0.156282, 1.038768, //vel
         //  0.0; //time error

    // L4 Vertical
    // x0 << 0.503828, 0.856388, 0.0, //pos
    //       0.059258, -0.034869, 0.366554, //vel
    //       0.0; //time error

    // L5 Long Period
    // x0 << 0.450191, -0.866025, 0.0, //pos
    //       0.023670, 0.014418, 0.0, //vel
    //       0.0; //time error

    //time parameters - normalized
    Real t0 = 0;
    Real tf = 30;
    Real dt = 0.1;
    Real t_span = tf - t0;
    int num_steps = static_cast<int>((tf - t0) / dt) + 1;

    MatX state_history(num_steps, 7);
    state_history.row(0) = x0;

    Real t = t0;
    VecX current_state = x0;

    // propagate state and tracking timing errors
    for (int i = 1; i < num_steps; i++) {

        VecX x = state_history.row(i-1);
        state_history.row(i).head(6) = cr3bp.Propagate(static_cast<Vec6>(x.head(6)), t, t + dt);
        state_history(i, 6) = 0.0;

        t += dt;
    }

//     cout << "Final Timing Error: " << state_history - t_span << endl;

    Vec3d earth_state(-mu.val(), 0.0, 0.0);      // Earth state in CR3BP
    Vec3d moon_state((1.0 - mu).val(), 0.0, 0.0);      // Moon state in CR3BP
    figure();
    hold(on);
    VecX x_vals = state_history.col(0);
    VecX y_vals = state_history.col(1);
    VecX z_vals = state_history.col(2);
    Plot3(x_vals,y_vals,z_vals, "b", 0);
    scatter3(std::vector<double>{earth_state(0)}, std::vector<double>{earth_state(1)}, std::vector<double>{earth_state(2)});
    scatter3(std::vector<double>{moon_state(0)}, std::vector<double>{moon_state(1)}, std::vector<double>{moon_state(2)});
    title("Circular Restricted Three-Body Problem");
    xlabel("X Position (Normalized)");
    ylabel("Y Position (Normalized)");
    zlabel("Z Position (Normalized)");
    show();

}
