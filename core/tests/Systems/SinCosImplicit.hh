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

#ifndef TESTS_SYSTEMS_SINCOS_IMPLICIT_HH
#define TESTS_SYSTEMS_SINCOS_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;
using namespace Eigen;

template <typename Real = double>
class SinCosImplicit : public Implicit<Real, 2, 0> {
 public:
  using VectorF  = typename Implicit<Real, 2, 0>::VectorF;
  using MatrixJF = typename Implicit<Real, 2, 0>::MatrixJF;
  using VectorH  = typename Implicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename Implicit<Real, 2, 0>::MatrixJH;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX  = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  SinCosImplicit() : Implicit<Real, 2, 0>("SinCosImplicit") {}

  ~SinCosImplicit() {}

  VectorF F(const VectorF & /*x*/,
            const VectorF &x_dot,
            Real t) const override {
    VectorF F;
    F << x_dot(0) - std::cos(t), x_dot(1) - std::sin(t);
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    return MatrixJF::Zero();
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
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
    x << std::sin(t), 1.0 - std::cos(t);
    return x;
  }

  static MatrixX exact_solution(const VectorX &t) {
    MatrixX x(2, t.size());
    for (int i = 0; i < t.size(); ++i) {
      x.col(i) = SinCosImplicit::exact_solution(t(i));
    }
    return x;
  }
};

#endif  // TESTS_SYSTEMS_SINCOS_IMPLICIT_HH
