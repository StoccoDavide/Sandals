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

#ifndef TESTS_SYSTEMS_SIMPLE_PENDULUM_IMPLICIT_HH
#define TESTS_SYSTEMS_SIMPLE_PENDULUM_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

template <typename Real = double>
class SimplePendulumImplicit : public Implicit<Real, 2, 0> {
 public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;
  using VectorX = Eigen::Matrix<Real, 2, 1>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  Real m_l{1.0};   // Length of the pendulum (m)
  Real m_g{9.81};  // Gravity acceleration (m/s^2)

 public:
  SimplePendulumImplicit() : Implicit<Real, 2, 0>("SimplePendulumImplicit") {}

  ~SimplePendulumImplicit() {}

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) + this->m_g / this->m_l * std::sin(x(0));
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x;
    JF_x << 0.0, -1.0, this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }

  static VectorF ics() {
    return VectorF::Unit(0, 2);
  }
};

#endif  // TESTS_SYSTEMS_SIMPLE_PENDULUM_IMPLICIT_HH
