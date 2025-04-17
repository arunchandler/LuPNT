#include <lupnt/lupnt.h>
#include <matplot/matplot.h>

using namespace lupnt;
using namespace std;

int main(){
    Real mjd_time = 51544.5;
    Real t_TT = MJD2Time(mjd_time);
    Real t_tai = TT2TAI(t_TT);
    cout << "t_tai: " << t_tai << endl;

    Mat3 R = GetFrameConversionMatrix(20122, Frame::ECI, Frame::ECEF);
    cout << "R: " << R << endl;

    Vec3 r_init = {6378.137, 0.0, 0.0};
    Vec3 r_final = R * r_init;
    cout << "r_final: " << r_final << endl;

}
