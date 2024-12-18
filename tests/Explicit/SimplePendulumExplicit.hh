/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_EXPLICIT_SIMPLEPENDULUMIMPLICIT_HH

using namespace Sandals;

class SimplePendulumExplicit : public Explicit<2>
{
  real m_l{1.0};        // length of the pendulum (m)
  real m_g{9.81};       // gravity acceleration (m/s^2)
  vec2 m_x_0{1.0, 0.0}; // initial conditions

public:
  using ptr = std::shared_ptr<SimplePendulumExplicit>;

  SimplePendulumExplicit(void) : Explicit<2>("SimplePendulumExplicit") {}

  ~SimplePendulumExplicit(void) {}

  vec2 f(vec2 const &x, real /*t*/) const override
  {
    vec2 f;
    f << x(1), -this->m_g / this->m_l * std::sin(x(0));
    return f;
  }

  mat2 Jf_x(vec2 const &x, real /*t*/) const override
  {
    mat2 Jf_x;
    Jf_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return Jf_x;
  }

  bool in_domain(vec2 const &/*x*/, real /*t*/) const override {return true;}

  vec2 ics(void) const {return this->m_x_0;}
};

#endif // TESTS_EXPLICIT_SIMPLEPENDULUMIMPLICIT_HH
