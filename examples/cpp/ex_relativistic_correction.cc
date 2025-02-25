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

    ClockOrbitDynamics clock_orbit(integ);

    clock_orbit.AddBody_COD(Body::Earth());

    VecX x0(7);
    x0 << 1.0, 0.0, 0.0, //pos
          0.0, 0.0, 0.0, //vel
          0.0; //initial time error

}
