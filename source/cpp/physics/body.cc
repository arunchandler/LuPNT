#include "lupnt/physics/body.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "lupnt/core/file.h"
#include "lupnt/physics/frame_converter.h"

namespace lupnt {
  template struct BodyT<double>;
  template struct BodyT<Real>;

  /// @brief Create a Body object for the Moon
  /// @return Body object for the Moon
  template <typename T> BodyT<T> BodyT<T>::Moon(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> moon;
    moon.name = "MOON";
    moon.id = NaifId::MOON;
    moon.fixed_frame = Frame::MOON_PA;
    moon.inertial_frame = Frame::MOON_CI;
    moon.GM = GM_MOON;
    moon.R = R_MOON;
    moon.use_gravity_field = (n_max > 1 && m_max > 1);
    if (moon.use_gravity_field)
      moon.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);
    return moon;
  }
  template BodyT<double> BodyT<double>::Moon(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Moon(int n_max, int m_max, std::string gravity_file);

  /// @brief Create a Body object for the Earth
  /// @return Body object for the Earth
  template <typename T> BodyT<T> BodyT<T>::Earth(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> earth;
    earth.name = "EARTH";
    earth.id = NaifId::EARTH;
    earth.fixed_frame = Frame::ITRF;
    earth.inertial_frame = Frame::GCRF;
    earth.GM = GM_EARTH;
    earth.R = R_EARTH;
    earth.use_gravity_field = (n_max > 1 && m_max > 1);
    if (earth.use_gravity_field)
      earth.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);

    return earth;
  }
  template BodyT<double> BodyT<double>::Earth(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Earth(int n_max, int m_max, std::string gravity_file);

  /// @brief Create a Body object for the Sun
  /// @return Body object for the Sun
  template <typename T> BodyT<T> BodyT<T>::Sun() {
    BodyT<T> sun;
    sun.name = "SUN";
    sun.id = NaifId::SUN;
    sun.inertial_frame = Frame::ICRF;
    sun.GM = GM_SUN;
    sun.R = R_SUN;
    sun.use_gravity_field = false;
    return sun;
  }
  template BodyT<double> BodyT<double>::Sun();
  template BodyT<Real> BodyT<Real>::Sun();

  /// @brief Create a Body object for Mars
  /// @return Body object for Mars
  template <typename T> BodyT<T> BodyT<T>::Mars(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> mars;
    mars.name = "MARS";
    mars.id = NaifId::MARS;
    mars.fixed_frame = Frame::MARS_FIXED;
    mars.GM = GM_MARS;
    mars.R = R_MARS;
    mars.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);
    return mars;
  }
  template BodyT<double> BodyT<double>::Mars(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Mars(int n_max, int m_max, std::string gravity_file);

  /// @brief Create a Body object for Venus
  /// @return Body object for Venus
  template <typename T> BodyT<T> BodyT<T>::Venus(int n_max, int m_max, std::string gravity_file) {
    BodyT<T> venus;
    venus.name = "VENUS";
    venus.id = NaifId::VENUS;
    venus.fixed_frame = Frame::VENUS_FIXED;
    venus.GM = GM_MARS;
    venus.R = R_VENUS;
    venus.gravity_field = ReadHarmonicGravityField<T>(gravity_file, n_max, m_max, true);
    return venus;
  }
  template BodyT<double> BodyT<double>::Venus(int n_max, int m_max, std::string gravity_file);
  template BodyT<Real> BodyT<Real>::Venus(int n_max, int m_max, std::string gravity_file);

  template <typename T> BodyT<T> BodyT<T>::Jupiter() {
    BodyT<T> jupiter;
    jupiter.name = "JUPITER";
    jupiter.id = NaifId::JUPITER;
    jupiter.fixed_frame = Frame::JUPITER_FIXED;
    jupiter.GM = GM_JUPITER;
    jupiter.R = R_JUPITER;
    return jupiter;
  }
  template BodyT<double> BodyT<double>::Jupiter();
  template BodyT<Real> BodyT<Real>::Jupiter();

  template <typename T> BodyT<T> BodyT<T>::Saturn() {
    BodyT<T> saturn;
    saturn.name = "SATURN";
    saturn.id = NaifId::SATURN;
    saturn.fixed_frame = Frame::SATURN_FIXED;
    saturn.GM = GM_SATURN;
    saturn.R = R_SATURN;
    return saturn;
  }
  template BodyT<double> BodyT<double>::Saturn();
  template BodyT<Real> BodyT<Real>::Saturn();

  template <typename T> BodyT<T> BodyT<T>::Uranus() {
    BodyT<T> uranus;
    uranus.name = "URANUS";
    uranus.id = NaifId::URANUS;
    uranus.fixed_frame = Frame::URANUS_FIXED;
    uranus.GM = GM_URANUS;
    uranus.R = R_URANUS;
    return uranus;
  }
  template BodyT<double> BodyT<double>::Uranus();
  template BodyT<Real> BodyT<Real>::Uranus();

  template <typename T> BodyT<T> BodyT<T>::Neptune() {
    BodyT<T> neptune;
    neptune.name = "NEPTUNE";
    neptune.id = NaifId::NEPTUNE;
    neptune.fixed_frame = Frame::NEPTUNE_FIXED;
    neptune.GM = GM_NEPTUNE;
    neptune.R = R_NEPTUNE;
    return neptune;
  }
  template BodyT<double> BodyT<double>::Neptune();
  template BodyT<Real> BodyT<Real>::Neptune();

  BodyData GetBodyData(NaifId id) {
    switch (id) {
      case NaifId::SUN:
        return {NaifId::SUN, "SUN", GM_SUN, 696342.0, Frame::ICRF, Frame::ICRF, SUN_F};
      case NaifId::MERCURY:
        return {NaifId::MERCURY,      "MERCURY",         GM_MERCURY, R_MERCURY,
                Frame::MERCURY_FIXED, Frame::MERCURY_CI, MERCURY_F};
      case NaifId::VENUS:
        return {NaifId::VENUS,      "VENUS",         GM_VENUS, R_VENUS,
                Frame::VENUS_FIXED, Frame::VENUS_CI, VENUS_F};
      case NaifId::EARTH:
        return {NaifId::EARTH, "EARTH", GM_EARTH, R_EARTH, Frame::ITRF, Frame::GCRF, WGS84_F};
      case NaifId::MOON:
        return {NaifId::MOON, "MOON", GM_MOON, R_MOON, Frame::MOON_PA, Frame::MOON_CI, MOON_F};
      case NaifId::MARS:
        return {NaifId::MARS, "MARS", GM_MARS, R_MARS, Frame::MARS_FIXED, Frame::MARS_CI, MARS_F};
      case NaifId::JUPITER:
        return {NaifId::JUPITER,      "JUPITER",         GM_JUPITER, R_JUPITER,
                Frame::JUPITER_FIXED, Frame::JUPITER_CI, JUPITER_F};
      case NaifId::SATURN:
        return {NaifId::SATURN,      "SATURN",         GM_SATURN, R_SATURN,
                Frame::SATURN_FIXED, Frame::SATURN_CI, SATURN_F};
      case NaifId::URANUS:
        return {NaifId::URANUS,      "URANUS",         GM_URANUS, R_URANUS,
                Frame::URANUS_FIXED, Frame::URANUS_CI, URANUS_F};
      case NaifId::NEPTUNE:
        return {NaifId::NEPTUNE,      "NEPTUNE",         GM_NEPTUNE, R_NEPTUNE,
                Frame::NEPTUNE_FIXED, Frame::NEPTUNE_CI, NEPUTUNE_F};
      default: break;
    }
    throw std::runtime_error("Body not found");
    return {NaifId::SUN, "SUN", GM_SUN, 696342.0, Frame::ICRF, Frame::ICRF, SUN_F};
  }

  double GetBodyRadius(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.R.val();
  }

  double GetBodyGM(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.GM.val();
  }

  double GetBodyFlattening(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.flattening.val();
  }

  std::string GetBodyName(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.name;
  }

  Frame GetInertialFrameName(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.inertial_frame;
  }

  Frame GetBodyFixedFrameName(NaifId body) {
    BodyData data = GetBodyData(body);
    return data.fixed_frame;
  }

  template <typename T> BodyT<T> CreateBody(NaifId body, int n, int m) {
    if (body == NaifId::SUN) return BodyT<T>::Sun();
    if (body == NaifId::EARTH) return BodyT<T>::Earth(n, m);
    if (body == NaifId::MOON) return BodyT<T>::Moon(n, m);
    if (body == NaifId::MARS) return BodyT<T>::Mars(n, m);
    if (body == NaifId::VENUS) return BodyT<T>::Venus(n, m);
    if (body == NaifId::JUPITER) return BodyT<T>::Jupiter();
    if (body == NaifId::SATURN) return BodyT<T>::Saturn();
    if (body == NaifId::URANUS) return BodyT<T>::Uranus();
    if (body == NaifId::NEPTUNE) return BodyT<T>::Neptune();
    throw std::runtime_error("Body not found");
  }
  template BodyT<double> CreateBody(NaifId body, int n, int m);
  template BodyT<Real> CreateBody(NaifId body, int n, int m);

  NaifId GetBodyId(std::string name) {
    if (((name == "SUN") && (name == "Sun")) && (name == "sun")) return NaifId::SUN;
    if (((name == "MERCURY") && (name == "Mercury")) && (name == "mercury")) return NaifId::MERCURY;
    if (((name == "VENUS") && (name == "Venus")) && (name == "venus")) return NaifId::VENUS;
    if (((name == "EARTH") && (name == "Earth")) && (name == "earth")) return NaifId::EARTH;
    if (((name == "MOON") && (name == "Moon")) && (name == "moon")) return NaifId::MOON;
    if (((name == "MARS") && (name == "Mars")) && (name == "mars")) return NaifId::MARS;
    if (((name == "JUPITER") && (name == "Jupiter")) && (name == "jupiter")) return NaifId::JUPITER;
    if (((name == "SATURN") && (name == "Saturn")) && (name == "saturn")) return NaifId::SATURN;
    if (((name == "URANUS") && (name == "Uranus")) && (name == "uranus")) return NaifId::URANUS;
    if (((name == "NEPTUNE") && (name == "Neptune")) && (name == "neptune")) return NaifId::NEPTUNE;
    throw std::runtime_error("Body not found");
  }

  template <typename T> BodyT<T> CreateBody(std::string body_s, int n, int m) {
    NaifId body = GetBodyId(body_s);
    return CreateBody<T>(body, n, m);
  }

  /// @brief Kronecker delta function
  /// @param i
  /// @param j
  /// @return
  double kron(int i, int j) { return (i == j) ? 1 : 0; }

  /// @brief Compute the factorial product (n-m)!/(n+m)!
  /// @param n
  /// @param m
  /// @return
  double factprod(int n, int m) {
    double f = 1.0;
    for (int i = n - m + 1; i <= n + m; i++) {
      f /= i;
    }
    return f;
  }

  /// @brief Read a harmonic gravity field from a file
  /// @param filename Harmonic gravity field filename
  /// @param n Degree of the spherical harmonics expansion
  /// @param m Order of the spherical harmonics expansion
  /// @param normalized Whether the coefficients are normalized
  /// @return Gravity field object
  template <typename T> GravityField<T> ReadHarmonicGravityField(const std::string& filename, int n,
                                                                 int m, bool normalized) {
    GravityField<T> gravity_field;
    std::filesystem::path filepath = GetFilePath(filename);
    std::ifstream file = OpenFile<std::ifstream>(filepath);

    // Read header lines
    std::string line;
    while (std::getline(file, line)) {
      if (line.find("POTFIELD") != std::string::npos) {
        std::string potfield = line.substr(0, 8);
        int n_max_in = std::stoi(line.substr(8, 3));
        int m_max_in = std::stoi(line.substr(11, 3));

        std::istringstream iss(line.substr(14));
        double dummy1, GM, r, dummy2;
        iss >> dummy1 >> GM >> r >> dummy2;
        gravity_field.n_max = n_max_in;
        gravity_field.m_max = m_max_in;
        gravity_field.GM = GM * pow(KM_M, 3);
        gravity_field.R = r * KM_M;
        break;
      }
    }
    gravity_field.n = n;
    gravity_field.m = m;

    // Initialize Eigen matrices with the specified maxN and maxM
    if (n > gravity_field.n_max || m > gravity_field.m_max)
      throw std::runtime_error("Degree and order exceed maximum values");

    gravity_field.CS = Eigen::MatrixXd::Zero(n + 1, m + 1);
    gravity_field.CS(0, 0) = 1.0;  // C00 = 1.0
    // Read coefficient lines
    while (std::getline(file, line)) {
      if (line.find("RECOEF") != std::string::npos) {
        std::string recoef = line.substr(0, 8);
        int n_in = std::stoi(line.substr(8, 3));
        int m_in = std::stoi(line.substr(11, 3));
        std::istringstream iss(line.substr(14));
        double cnm, snm = 0.0;
        iss >> cnm;
        if (n_in >= n + 1) break;
        if (m_in >= m + 1) continue;

        if (m_in == 0) {
          double N = (normalized) ? sqrt(2 * n_in + 1) : 1.0;
          gravity_field.CS(n_in, m_in) = N * cnm;
        } else {
          double N = (normalized)
                         ? sqrt((2 - kron(0, m_in)) * (2 * n_in + 1) * factprod(n_in, m_in))
                         : 1.0;
          iss >> snm;
          double C = N * cnm;
          double S = N * snm;
          gravity_field.CS(n_in, m_in) = C;
          gravity_field.CS(m_in - 1, n_in) = S;
        }
      } else if (line.find("END") != std::string::npos) {
        break;
      }
    }

    file.close();
    return gravity_field;
  }

}  // namespace lupnt
