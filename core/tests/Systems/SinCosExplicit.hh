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

#ifndef TESTS_SYSTEMS_SINCOS_EXPLICIT_HH
#define TESTS_SYSTEMS_SINCOS_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;
using namespace Eigen;

template <typename Real = double>
class SinCosExplicit : public Explicit<Real, 2, 0> {
 public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  SinCosExplicit() : Explicit<Real, 2, 0>("SinCosExplicit") {}

  ~SinCosExplicit() {}

  VectorF f(const VectorF & /*x*/, Real t) const override {
    VectorF f;
    f << std::cos(t), std::sin(t);
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJF::Zero();
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
      x.col(i) = SinCosExplicit::exact_solution(t(i));
    }
    return x;
  }
};

#endif  // TESTS_SYSTEMS_SINCOS_EXPLICIT_HH
