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

#ifndef TESTS_OSCILLATOR_IMPLICIT_HH
#define TESTS_OSCILLATOR_IMPLICIT_HH

using namespace Sandals;

class OscillatorImplicit : public Implicit<2, 1>
{
  Real    m_m{1.0};        // Mass (kg)
  Real    m_k{1.0};        // Spring constant (N/m)
  VectorF m_ics{1.0, 0.0}; // Initial conditions

public:
  using VectorF  = typename Implicit<2, 1>::VectorF;
  using MatrixJF = typename Implicit<2, 1>::MatrixJF;
  using VectorH  = typename Implicit<2, 1>::VectorH;
  using MatrixJH = typename Implicit<2, 1>::MatrixJH;

  OscillatorImplicit() : Implicit<2, 1>("OscillatorImplicit") {}

  ~OscillatorImplicit() {}

  VectorF F(VectorF const &x, VectorF const &x_dot, Real /*t*/) const override
  {
    VectorF F;
    F << x_dot(0)-x(1), x_dot(1)+this->m_k/this->m_m*x(0);
    return F;
  }

  MatrixJF JF_x(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const override {
    MatrixJF JF_x;
    JF_x << -1.0, 1.0, this->m_k/this->m_m, 0.0;
    return JF_x;
  }

  MatrixJF JF_x_dot(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const override
  {return MatrixJF::Identity();}

  Real energy(VectorF const &x) const {return this->m_m/2.0*x(1)*x(1) + this->m_k/2.0*x(0)*x(0);}

  VectorH h(VectorF const &x, Real /*t*/) const override {
    VectorH h;
    h << this->energy(x) - this->energy(this->m_ics);
    return h;
  }

  MatrixJH Jh_x(VectorF const &x, Real /*t*/) const override {
    MatrixJH Jh_x;
    Jh_x << this->m_k*x(0), this->m_m*x(1);
    return Jh_x;
  }

  bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}

  VectorF analytical_solution(Real t) const {
    VectorF x;
    Real tmp{std::sqrt(this->m_k / this->m_m)};
    x << std::cos(tmp*t), -tmp*std::sin(tmp*t);
    return x;
  }

  MatrixX analytical_solution(VectorX const &t) const {
    MatrixX x(2, t.size());
    for (int i = 0; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }

  VectorF const & ics() const {return this->m_ics;}
};

#endif // TESTS_OSCILLATOR_IMPLICIT_HH
