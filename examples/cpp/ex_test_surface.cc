#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;
using namespace matplot;

int main(){

    Frame frame = Frame::MOON_CI;
    SurfaceStaticDynamics surface_dynamics(NaifId::MOON, frame);
    Vec6 init_surface_state = {0.0, 0.0, -R_MOON, 0.0, 0.0, 0.0};

    Real t0 = 0.0;
    Real tf = SECS_DAY;
    Real dt = 10.0;
    Real init_t_err = 0.0;
    int num_steps = static_cast<int>((tf - t0) / dt) + 1;

    MatX surface_state_history(num_steps, 6);
    surface_state_history.row(0) = init_surface_state;

    Real t = t0;
    for (int i = 1; i < num_steps; i++) {
        Vec6 x_surface = surface_state_history.row(i-1);
        surface_state_history.row(i) = surface_dynamics.Propagate(x_surface, t, t + dt);
        t += dt;
    }

    //show results
    cout << "Surface State: " << surface_state_history.row(num_steps-1) << endl;
}
