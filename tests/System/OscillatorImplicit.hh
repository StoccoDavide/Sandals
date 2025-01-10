/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_OSCILLATOR_EXPLICIT_HH

using namespace Sandals;

class OscillatorImplicit : public Implicit<2, 1>
{
  Real    m_m{1.0};        // Mass (kg)
  Real    m_k{1.0};        // Spring constant (N/m)
  Vector2 m_ics{1.0, 0.0}; // Initial conditions

public:
  using Matrix12 = Eigen::Matrix<Real, 1, 2>;

  OscillatorImplicit() : Implicit<2, 1>("OscillatorImplicit") {}

  ~OscillatorImplicit() {}

  Vector2 F(Vector2 const &x, Vector2 const &x_dot, Real /*t*/) const override
  {
    Vector2 F;
    F << x_dot(0)-x(1), x_dot(1)+this->m_k/this->m_m*x(0);
    return F;
  }

  Matrix2 JF_x(Vector2 const &/*x*/, Vector2 const &/*x_dot*/, Real /*t*/) const override {
    Matrix2 JF_x;
    JF_x << -1.0, 1.0, this->m_k/this->m_m, 0.0;
    return JF_x;
  }

  Matrix2 JF_x_dot(Vector2 const &/*x*/, Vector2 const &/*x_dot*/, Real /*t*/) const override
  {return Matrix2::Identity();}

  Real energy(Vector2 const &x) const {return this->m_m/2.0*x(1)*x(1) + this->m_k/2.0*x(0)*x(0);}

  Vector1 h(Vector2 const &x, Real /*t*/) const override {
    Vector1 h;
    h << this->energy(x) - this->energy(this->m_ics);
    return h;
  }

  Matrix12 Jh_x(Vector2 const &/*x*/, Real /*t*/) const override {
    Matrix12 Jh_x;
    Jh_x << this->m_k, 0.0;
    return Jh_x;
  }

  bool in_domain(Vector2 const &/*x*/, Real /*t*/) const override {return true;}

  Vector2 analytical_solution(Real t) const {
    Vector2 x;
    Real tmp{std::sqrt(this->m_k / this->m_m)};
    x << std::cos(tmp*t), -tmp*std::sin(tmp*t);
    return x;
  }

  Vector2 ics() const {return this->m_ics;}
};

#endif // TESTS_OSCILLATOR_EXPLICIT_HH
