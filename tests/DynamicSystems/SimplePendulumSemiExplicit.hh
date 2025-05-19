/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef TESTS_SIMPLE_PENDULUM_SEMIEXPLICIT_HH
#define TESTS_SIMPLE_PENDULUM_SEMIEXPLICIT_HH

using namespace Sandals;

class SimplePendulumSemiExplicit : public SemiExplicit<double, 2, 0>
{
  Real    m_l{1.0};        // Length of the pendulum (m)
  Real    m_g{9.81};       // Gravity acceleration (m/s^2)
  VectorF m_ics{1.0, 0.0}; // Initial conditions

public:
  using VectorF  = typename SemiExplicit<double, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<double, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<double, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<double, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<double, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<double, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<double, 2, 0>::MatrixJH;

  SimplePendulumSemiExplicit() : SemiExplicit<double, 2, 0>("SimplePendulumSemiExplicit") {}

  ~SimplePendulumSemiExplicit() {}

  MatrixA A(VectorF const &/*x*/, Real /*t*/) const override
  {
    MatrixA A;
    A.setIdentity();
    return A;
  }

  TensorTA TA_x(VectorF const &/*x*/, Real /*t*/) const override
  {
    TensorTA TA_x(2);
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(VectorF const &x, Real /*t*/) const override
  {
    VectorB b;
    b << x(1), -this->m_g / this->m_l * std::sin(x(0));
    return b;
  }

  MatrixJB Jb_x(VectorF const &x, Real /*t*/) const override
  {
    MatrixJB Jb_x;
    Jb_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return Jb_x;
  }

  VectorH h(VectorF const &/*x*/, Real /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const &/*x*/, Real /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}

  VectorF const & ics() const {return this->m_ics;}

};

#endif // TESTS_SIMPLE_PENDULUM_SEMIEXPLICIT_HH
