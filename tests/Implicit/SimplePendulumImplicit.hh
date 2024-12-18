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
  real m_l{1.0};        // Length of the pendulum (m)
  real m_g{9.81};       // Gravity acceleration (m/s^2)
  vec2 m_ics{1.0, 0.0}; // Initial conditions

public:
  using ptr = std::shared_ptr<SimplePendulumImplicit>;

  SimplePendulumImplicit(void) : Implicit<2>("SimplePendulumImplicit") {}

  ~SimplePendulumImplicit(void) {}

  vec2 F(vec2 const &x, vec2 const &x_dot, real /*t*/) const override
  {
    vec2 F;
    F << x_dot(0) - x(1), x_dot(1) + this->m_g / this->m_l * std::sin(x(0));
    return F;
  }

  mat2 JF_x(vec2 const &x, vec2 const &/*x_dot*/, real /*t*/) const override
  {
    mat2 JF_x;
    JF_x << 0.0, -1.0, this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return JF_x;
  }

  mat2 JF_x_dot(vec2 const &/*x*/, vec2 const &/*x_dot*/, real /*t*/) const override
  {return mat2::Identity();}

  bool in_domain(vec2 const &/*x*/, real /*t*/) const override {return true;}

  vec2 const & ics(void) const {return this->m_ics;}
};

#endif // TESTS_IMPLICIT_SIMPLEPENDULUMIMPLICIT_HH
