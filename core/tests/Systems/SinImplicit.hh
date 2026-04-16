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

#ifndef TESTS_SYSTEMS_SIN_IMPLICIT_HH
#define TESTS_SYSTEMS_SIN_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;
using namespace Eigen;

template <typename Real = double>
class SinImplicit : public Implicit<Real, 1, 0> {
 public:
  using typename Implicit<Real, 1, 0>::VectorF;
  using typename Implicit<Real, 1, 0>::MatrixJF;
  using typename Implicit<Real, 1, 0>::VectorH;
  using typename Implicit<Real, 1, 0>::MatrixJH;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

  SinImplicit() : Implicit<Real, 1, 0>("SinImplicit") {}

  ~SinImplicit() {}

  VectorF F(const VectorF & /*x*/,
            const VectorF &x_dot,
            Real t) const override {
    VectorF F;
    F << x_dot(0) - std::cos(t);
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
      x.col(i) = SinImplicit::exact_solution(t(i));
    }
    return x;
  }
};

#endif  // TESTS_SYSTEMS_SIN_IMPLICIT_HH
