/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_SYSTEMS_SINCOS_SEMIEXPLICIT_HH
#define TESTS_SYSTEMS_SINCOS_SEMIEXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;
using namespace Eigen;

template <typename Real = double>
class SinCosSemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using typename SemiExplicit<Real, 2, 0>::VectorF;
  using typename SemiExplicit<Real, 2, 0>::MatrixA;
  using typename SemiExplicit<Real, 2, 0>::TensorTA;
  using typename SemiExplicit<Real, 2, 0>::VectorB;
  using typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using typename SemiExplicit<Real, 2, 0>::VectorH;
  using typename SemiExplicit<Real, 2, 0>::MatrixJH;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  SinCosSemiExplicit() : SemiExplicit<Real, 2, 0>("SinCosSemiExplicit") {}

  ~SinCosSemiExplicit() {}

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(const VectorF & /*x*/, Real t) const override {
    VectorB b;
    b << std::cos(t), std::sin(t);
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJB::Zero();
  }

  static VectorF ics() {
    return VectorF::Zero();
  }

  static VectorF exact_solution(const Real t) {
    VectorF x;
    x << std::sin(t), Real(1.0) - std::cos(t);
    return x;
  }

  static MatrixX exact_solution(const VectorX &t) {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) = SinCosSemiExplicit::exact_solution(t(i));
    }
    return x;
  }
};

#endif  // TESTS_SYSTEMS_SINCOS_SEMIEXPLICIT_HH
