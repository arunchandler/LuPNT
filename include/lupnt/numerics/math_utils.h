/**
 * @file math_utils.h
 * @author Stanford NAV LAB
 * @brief
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <random>
#include <tuple>
#include <utility>

#include "lupnt/core/constants.h"

namespace lupnt {

  template <typename Vec, std::size_t... Indices>
  auto unpackImpl(const Vec& vec, std::index_sequence<Indices...>) {
    return std::make_tuple(vec(Indices)...);
  }

  template <typename T, int Size> auto unpack(const Eigen::Matrix<T, Size, 1>& vec) {
    return unpackImpl(vec, std::make_index_sequence<Size>{});
  }

  template <typename T> VectorX<T> arange(T start, T stop, T step = 1);

  Real AngleBetweenVecs(const VecX& x, const VecX& y);

  double DegMinSec2DeciDeg(double degrees, double minutes, double seconds);

  Real Wrap2Pi(Real angle);
  Real Wrap2TwoPi(Real angle);
  VecX Wrap2Pi(VecX angle);
  VecX Wrap2TwoPi(VecX angle);

  Real Decimal2Decibel(Real x);
  VecX Decimal2Decibel(VecX x);
  MatX Decimal2Decibel(MatX x);

  Real Decibel2Decimal(Real x);
  VecX Decibel2Decimal(VecX x);
  MatX Decibel2Decimal(MatX x);

  Real Max(Real x, Real y);
  Real Min(Real x, Real y);
  double MaxD(double x, double y);
  double MinD(double x, double y);

  Real round(Real x, int n = 0);
  Real frac(Real x);
  Real ceil(Real x);
  Real floor(Real x);
  Real mod(Real x, Real y);

  Vec3 Degrees2DegMinSec(Real deg);
  Real DegMinSec2Degrees(Vec3 hms);

  Real sind(Real x);
  Real cosd(Real x);
  Real tand(Real x);

  Real safe_acos(Real x);
  Real safe_asin(Real x);

  /// @brief Compute the Bessel function of the first kind of order 0
  /// @param x
  /// @return
  /// @note https://en.wikipedia.org/wiki/Bessel_function
  template <typename T> T J0Bessel(T x) {
    // double J0 = 0.0;
    T y = 1.0;
    T sum = 1.0;
    for (int i = 1; i < 10; i++) {
      y = y * x * x / (4 * i * i);
      sum += y;
    }
    return sum;
  }

  /// @brief Compute the Bessel function of the first kind of order 1
  /// @tparam T
  /// @param x
  /// @return
  /// @note https://en.wikipedia.org/wiki/Bessel_function
  template <typename T> T J1Bessel(T x) {
    // double J1 = 0.0;
    T y = 1.0;
    T sum = 1.0;
    for (int i = 1; i < 10; i++) {
      y = y * x / (2 * i * (2 * i + 1));
      sum += y;
    }
    return sum;
  }

  /// @brief Compute the root mean square of a vector
  /// @param x Input vector
  /// @return Root mean square of the vector
  template <typename T> T::Scalar RootMeanSquare(const MatrixBase<T>& x) {
    return x.norm() / sqrt(x.size());
  }

  /// @brief Compute the pth percentile of a vector
  /// @param x Input vector
  /// @param p Percentile value
  template <typename T> T::Scalar Percentile(const MatrixBase<T>& x, double p) {
    auto evaluated = x.eval();
    std::vector<typename T::Scalar> data(evaluated.data(), evaluated.data() + evaluated.size());
    std::sort(data.begin(), data.end());
    size_t index = std::ceil(p * (data.size() - 1));
    return data[index];
  }

  /// @brief Compute the standard deviation of a vector
  /// @param x Input vector
  /// @return Standard deviation of the vector
  template <typename T> T::Scalar Std(const MatrixBase<T>& x) {
    return sqrt((x.array() - x.mean()).square().sum() / (x.size() - 1));
  }

  template <typename T> T erfc(T x) { return 1 - erf(x.val()); }
  template <typename T> T qfunc(T x) { return 0.5 * erfc(x / sqrt(2)); }

  MatX SampleMVN(const VecX& mean, const MatX& cov, int nn, int seed = 0);
  Real SampleRandNormal(Real mean, Real std, int seed = 0);

  template <typename T, int N1, int M1, int N2, int M2>
  Matrix<T, N1 + N2, M1 + N2> BlockDiagonal(const Matrix<T, N1, M1>& A,
                                            const Matrix<T, N2, M2>& B) {
    Matrix<T, N1 + N2, M1 + M2> C;
    C << A, Matrix<T, N1, M2>::Zero(), Matrix<T, N2, M1>::Zero(), B;
    return C;
  }

  template <typename T> Matrix<T, 3, 3> RotX(T angle);
  template <typename T> Matrix<T, 3, 3> RotY(T angle);
  template <typename T> Matrix<T, 3, 3> RotZ(T angle);
  template <typename T> Matrix<T, 3, 3> Skew(Vector<T, 3> x);

  VecXd ToDouble(const VecX& x);
  // MatXd ToDouble(const MatX& x);
  std::vector<double> ToDoubleVec(const VecX& x);
  std::vector<double> ToDoubleVec(const VecXd& x);
  std::vector<double> ToDoubleVec(const VecXi& x);

  Real RatioOfSectorToTriangleArea(Vec3 r1, Vec3 r2, Real tau);

  VecX arange(Real start, Real stop, Real step);

  VecXd SolveLinearEqSVD(const MatXd& A, const VecXd& b);  // Ax = b

  MatXd SolveLinearEqSVD(const MatXd& A, const MatXd& B);  // AX = B

  MatXd PseudoInverse(const MatXd& A);

}  // namespace lupnt
