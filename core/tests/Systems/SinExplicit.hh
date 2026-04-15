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

#ifndef TESTS_SYSTEMS_SIN_EXPLICIT_HH
#define TESTS_SYSTEMS_SIN_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;
using namespace Eigen;

template <typename Real = double>
class SinExplicit : public Explicit<Real, 1, 0> {
 public:
  using VectorF  = typename Explicit<Real, 1, 0>::VectorF;
  using MatrixJF = typename Explicit<Real, 1, 0>::MatrixJF;
  using VectorH  = typename Explicit<Real, 1, 0>::VectorH;
  using MatrixJH = typename Explicit<Real, 1, 0>::MatrixJH;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX  = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

  SinExplicit() : Explicit<Real, 1, 0>("SinExplicit") {}

  ~SinExplicit() {}

  VectorF f(const VectorF & /*x*/, Real t) const override {
    VectorF f;
    f << std::cos(t);
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJF::Zero();
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
      x.col(i) = SinExplicit::exact_solution(t(i));
    }
    return x;
  }
};

#endif  // TESTS_SYSTEMS_SIN_EXPLICIT_HH
