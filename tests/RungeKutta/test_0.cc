/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"

using namespace Sandals;

template class Sandals::RK4<2>;

class SimplePendulum : public Implicit<2>
{
  real m_l{1.0};        // length of the pendulum (m)
  real m_g{9.81};       // gravity acceleration (m/s^2)
  vec2 m_x_0{0.0, 0.0}; // initial conditions

public:
  SimplePendulum(void) : Implicit<2>("SimplePendulum") {}

  ~SimplePendulum(void) {}

  vec2 ic(void) const {return this->m_x_0;}

  using Implicit<2>::F;
  vec2 F(vec2 const &x, vec2 const &x_dot, real /*t*/) const
  {
    vec2 F;
    F << x_dot(1) - x(2), x_dot(2) + this->m_g / this->m_l * std::sin(x(1));
    return F;
  }

  using Implicit<2>::JF_x;
  mat2 JF_x(vec2 const &x, vec2 const &/*x_dot*/, real /*t*/) const
  {
    mat2 JF_x;
    JF_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return JF_x;
  }

  using Implicit<2>::JF_x_dot;
  mat2 JF_x_dot(vec2 const &/*x*/, vec2 const &/*x_dot*/, real /*t*/) const
  {return mat2::Identity();}

  using Implicit<2>::in_domain;
  bool in_domain(vec2 const &/*x*/, real /*t*/) const {return true;}
};

int main(void) {
  SimplePendulum ode;
  RK4<2> rk;
  rk.ode(std::make_shared<const SimplePendulum>(ode));

  vecD t = Eigen::VectorXd::LinSpaced(100, 0.0, 10.0);
  RK4<2>::matND x;
  rk.solve(t, ode.ic(), x);

  std::cout << "Hi, I'm finding the sandal size for you!" << std::endl;
  return 0;
}
