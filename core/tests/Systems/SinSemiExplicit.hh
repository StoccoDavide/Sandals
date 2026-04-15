/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef TESTS_SYSTEMS_SIN_SEMIEXPLICIT_HH
#define TESTS_SYSTEMS_SIN_SEMIEXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;
using namespace Eigen;

template <typename Real = double>
class SinSemiExplicit : public SemiExplicit<Real, 1, 0> {
 public:
  using VectorF  = typename SemiExplicit<Real, 1, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 1, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 1, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 1, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 1, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 1, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 1, 0>::MatrixJH;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX  = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

  SinSemiExplicit() : SemiExplicit<Real, 1, 0>("SinSemiExplicit") {}

  ~SinSemiExplicit() {}

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    return TA_x;
  }

  VectorB b(const VectorF & /*x*/, Real t) const override {
    VectorB b;
    b << std::cos(t);
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJB::Zero();
  }

  VectorH h(const VectorF & /*x*/, const Real /*t*/) const override {
    return VectorH::Zero();
  }

  MatrixJH Jh_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJH::Zero();
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }

  static VectorF ics() {
    return VectorF::Zero();
  }

  static VectorF exact_solution(const Real t) {
    VectorF x;
    x << std::sin(t);
    return x;
  }

  static MatrixX exact_solution(const VectorX &t) {
    MatrixX x(1, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) = SinSemiExplicit::exact_solution(t(i));
    }
    return x;
  }
};

#endif  // TESTS_SYSTEMS_SIN_SEMIEXPLICIT_HH
