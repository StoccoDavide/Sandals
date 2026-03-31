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

#ifndef TESTS_SYSTEMS_SIMPLE_PENDULUM_EXPLICIT_HH
#define TESTS_SYSTEMS_SIMPLE_PENDULUM_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;

template <typename Real = double>
class SimplePendulumExplicit : public Explicit<Real, 2, 0> {
 public:
  using VectorF  = typename Explicit<Real, 2, 0>::VectorF;
  using MatrixJF = typename Explicit<Real, 2, 0>::MatrixJF;
  using VectorH  = typename Explicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename Explicit<Real, 2, 0>::MatrixJH;
  using VectorX  = Eigen::Matrix<Real, 2, 1>;
  using MatrixX  = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  Real m_l{1.0};   // length of the pendulum (m)
  Real m_g{9.81};  // gravity acceleration (m/s^2)

 public:
  SimplePendulumExplicit() : Explicit<Real, 2, 0>("SimplePendulumExplicit") {}

  ~SimplePendulumExplicit() {}

  VectorF f(const VectorF &x, const Real /*t*/) const override {
    VectorF f;
    f << x(1), -this->m_g / this->m_l * std::sin(x(0));
    return f;
  }

  MatrixJF Jf_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJF Jf_x;
    Jf_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return Jf_x;
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
    return VectorF::Unit(0, 2);
  }
};

#endif  // TESTS_SYSTEMS_SIMPLE_PENDULUM_EXPLICIT_HH
