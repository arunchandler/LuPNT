#include "lupnt/core/plot.h"

#include <matplot/matplot.h>

#include "lupnt/core/constants.h"
#include "lupnt/numerics/math_utils.h"
#include "lupnt/physics/body.h"

namespace lupnt {
  matplot::line_handle Plot3(const VecX &x, const VecX &y, const VecX &z,
                             std::string_view line_spec, double scale) {
    scale = pow(10, scale);
    return matplot::plot3(ToDouble(x / scale), ToDouble(y / scale), ToDouble(z / scale), line_spec);
  }

  matplot::line_handle Plot3(const MatX &xyz, std::string_view line_spec, double scale) {
    VecX x = xyz.col(0), y = xyz.col(1), z = xyz.col(2);
    return Plot3(x, y, z, line_spec, scale);
  }

  matplot::line_handle Scatter3(const Vec3 &xyz, std::string_view line_spec, double scale) {
    Vec1 x(xyz(0)), y(xyz(1)), z(xyz(2));
    return Plot3(x, y, z, line_spec, scale);
  }

  matplot::line_handle PlotArrow3(const Vec3 &center, const Vec3 &dir, std::string_view line_spec,
                                  double scale) {
    Vec2 x(center(0), center(0) + dir(0));
    Vec2 y(center(1), center(1) + dir(1));
    Vec2 z(center(2), center(2) + dir(2));
    return Plot3(x, y, z, line_spec, scale);
  }

  matplot::line_handle Plot(const VecX &x, const VecX &y, std::string_view line_spec) {
    return matplot::plot(ToDouble(x), ToDouble(y), line_spec);
  }

  std::vector<matplot::line_handle> PlotFrame(const Vec3 &center, const Mat3 &R, double scale) {
    std::vector<std::string> line_specs = {"r-", "g-", "b-"};
    std::vector<matplot::line_handle> lines;
    for (int i = 0; i < 3; i++) lines.push_back(PlotArrow3(center, R.row(i), line_specs[i], scale));
    return lines;
  }

  matplot::surface_handle PlotBody(NaifId body, Vec3 r_body, double scale) {
    scale = pow(10, scale);
    Vec3d r_body_ = ToDouble(r_body / scale);
    double radius = GetBodyRadius(body) / scale;
    using namespace matplot;
    int n = 40;
    auto theta = linspace(0, pi, n);
    auto phi = linspace(0, 2 * pi, n);
    auto [T, P] = meshgrid(theta, phi);
    auto X = transform(T, P, [radius, r_body_](double theta, double phi) {
      return radius * sin(theta) * cos(phi) + r_body_(0);
    });
    auto Y = transform(T, P, [radius, r_body_](double theta, double phi) {
      return radius * sin(theta) * sin(phi) + r_body_(1);
    });
    auto Z = transform(T, P, [radius, r_body_](double theta, double phi) {
      (void)phi;
      return radius * cos(theta) + r_body_(2);
    });
    auto h = surf(X, Y, Z);
    h->edge_color("gray");
    xlabel("X [1e3 km]");
    ylabel("Y [1e3 km]");
    zlabel("Z [1e3 km]");
    // colormap(palette::gray());
    return h;
  }

  void SetLim(Real lim, double scale) {
    scale = pow(10, scale);
    double lim_ = lim.val() / scale;
    matplot::xlim({-lim_, lim_});
    matplot::ylim({-lim_, lim_});
    matplot::zlim({-lim_, lim_});
  }

}  // namespace lupnt
