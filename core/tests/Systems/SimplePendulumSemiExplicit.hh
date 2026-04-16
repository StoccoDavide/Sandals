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

#ifndef TESTS_SYSTEMS_SIMPLE_PENDULUM_SEMIEXPLICIT_HH
#define TESTS_SYSTEMS_SIMPLE_PENDULUM_SEMIEXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

template <typename Real = double>
class SimplePendulumSemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using typename SemiExplicit<double, 2, 0>::VectorF;
  using typename SemiExplicit<double, 2, 0>::MatrixA;
  using typename SemiExplicit<double, 2, 0>::TensorTA;
  using typename SemiExplicit<double, 2, 0>::VectorB;
  using typename SemiExplicit<double, 2, 0>::MatrixJB;
  using typename SemiExplicit<double, 2, 0>::VectorH;
  using typename SemiExplicit<double, 2, 0>::MatrixJH;
  using VectorX = Eigen::Matrix<Real, 2, 1>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  Real m_l{1.0};   // Length of the pendulum (m)
  Real m_g{9.81};  // Gravity acceleration (m/s^2)

 public:
  SimplePendulumSemiExplicit()
      : SemiExplicit<double, 2, 0>("SimplePendulumSemiExplicit") {}

  ~SimplePendulumSemiExplicit() {}

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(const VectorF &x, const Real /*t*/) const override {
    VectorB b;
    b << x(1), -this->m_g / this->m_l * std::sin(x(0));
    return b;
  }

  MatrixJB Jb_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJB Jb_x;
    Jb_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return Jb_x;
  }

  static VectorF ics() {
    return VectorF::Unit(0, 2);
  }
};

#endif  // TESTS_SYSTEMS_SIMPLE_PENDULUM_SEMIEXPLICIT_HH
