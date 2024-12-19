/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_IMPLICIT_SIMPLEPENDULUMIMPLICIT_HH

using namespace Sandals;

class SimplePendulumImplicit : public Implicit<2>
{
  Real    m_l{1.0};        // Length of the pendulum (m)
  Real    m_g{9.81};       // Gravity acceleration (m/s^2)
  Vector2 m_ics{1.0, 0.0}; // Initial conditions

public:
  SimplePendulumImplicit(void) : Implicit<2>("SimplePendulumImplicit") {}

  ~SimplePendulumImplicit(void) {}

  Vector2 F(Vector2 const &x, Vector2 const &x_dot, Real /*t*/) const override
  {
    Vector2 F;
    F << x_dot(0) - x(1), x_dot(1) + this->m_g / this->m_l * std::sin(x(0));
    return F;
  }

  Matrix2 JF_x(Vector2 const &x, Vector2 const &/*x_dot*/, Real /*t*/) const override
  {
    Matrix2 JF_x;
    JF_x << 0.0, -1.0, this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return JF_x;
  }

  Matrix2 JF_x_dot(Vector2 const &/*x*/, Vector2 const &/*x_dot*/, Real /*t*/) const override
  {return Matrix2::Identity();}

  bool in_domain(Vector2 const &/*x*/, Real /*t*/) const override {return true;}

  Vector2 const & ics(void) const {return this->m_ics;}
};

#endif // TESTS_IMPLICIT_SIMPLEPENDULUMIMPLICIT_HH
