/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef TESTS_SIMPLE_PENDULUM_EXPLICIT_HH
#define TESTS_SIMPLE_PENDULUM_EXPLICIT_HH

using namespace Sandals;

class SimplePendulumExplicit : public Explicit<2, 0>
{
  Real    m_l{1.0};        // length of the pendulum (m)
  Real    m_g{9.81};       // gravity acceleration (m/s^2)
  VectorF m_ics{1.0, 0.0}; // initial conditions

public:
  using VectorF  = typename Explicit<2, 0>::VectorF;
  using MatrixJF = typename Explicit<2, 0>::MatrixJF;
  using VectorH  = typename Explicit<2, 0>::VectorH;
  using MatrixJH = typename Explicit<2, 0>::MatrixJH;

  SimplePendulumExplicit() : Explicit<2, 0>("SimplePendulumExplicit") {}

  ~SimplePendulumExplicit() {}

  VectorF f(VectorF const &x, Real /*t*/) const override
  {
    VectorF f;
    f << x(1), -this->m_g / this->m_l * std::sin(x(0));
    return f;
  }

  Matrix2 Jf_x(VectorF const &x, Real /*t*/) const override
  {
    Matrix2 Jf_x;
    Jf_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return Jf_x;
  }

  VectorH h(VectorF const &/*x*/, Real /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const &/*x*/, Real /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}

  VectorF const & ics() const {return this->m_ics;}
};

#endif // TESTS_SIMPLE_PENDULUM_EXPLICIT_HH
