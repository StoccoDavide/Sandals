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
#include <matplot/matplot.h>

using namespace Sandals;
using namespace matplot;

template class Sandals::RK4<2>;

class SimplePendulumImplicit : public Implicit<2>
{
  real m_l{1.0};        // length of the pendulum (m)
  real m_g{9.81};       // gravity acceleration (m/s^2)
  vec2 m_x_0{1.0, 0.0}; // initial conditions

public:
  using ptr = std::shared_ptr<SimplePendulumImplicit>;

  SimplePendulumImplicit(void) : Implicit<2>("SimplePendulumImplicit") {}

  ~SimplePendulumImplicit(void) {}

  vec2 ic(void) const {return this->m_x_0;}

  vec2 F(vec2 const &x, vec2 const &x_dot, real /*t*/) const override
  {
    vec2 F;
    F << x_dot(1) - x(2), x_dot(2) + this->m_g / this->m_l * std::sin(x(1));
    return F;
  }

  mat2 JF_x(vec2 const &x, vec2 const &/*x_dot*/, real /*t*/) const override
  {
    mat2 JF_x;
    JF_x << 0.0, 1.0, -this->m_g / this->m_l * std::cos(x(0)), 0.0;
    return JF_x;
  }

  mat2 JF_x_dot(vec2 const &/*x*/, vec2 const &/*x_dot*/, real /*t*/) const override
  {return mat2::Identity();}

  bool in_domain(vec2 const &/*x*/, real /*t*/) const override {return true;}
};

class SimplePendulumExplicit : public Explicit<2>
{
  real m_l{1.0};        // length of the pendulum (m)
  real m_g{9.81};       // gravity acceleration (m/s^2)
  vec2 m_x_0{1.0, 0.0}; // initial conditions

public:
  using ptr = std::shared_ptr<SimplePendulumExplicit>;

  SimplePendulumExplicit(void) : Explicit<2>("SimplePendulumExplicit") {}

  ~SimplePendulumExplicit(void) {}

  vec2 ic(void) const {return this->m_x_0;}

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
};

int main(void) {
  SimplePendulumExplicit::ptr system{std::make_shared<SimplePendulumExplicit>()};
  RK4<2> integrator(system);
  RK4<2>::time t = Eigen::VectorXd::LinSpaced(10, 0.0, 1.0);
  RK4<2>::solution x;
  integrator.solve(t, system->ic(), x);

  std::cout << "t = " << t.transpose() << std::endl;
  std::cout << "x = " << x << std::endl;

  std::vector<real> omega(x.row(0).data(), x.row(0).data() + x.row(0).size());
  std::vector<real> theta(x.row(1).data(), x.row(1).data() + x.row(1).size());

  auto p = plot(t, omega, t, theta);
  p[0]->line_width(2);
  p[1]->marker(line_spec::marker_style::asterisk);

  std::cout << "Hi, I'm finding the sandal size for you!" << std::endl;
  return 0;
}
