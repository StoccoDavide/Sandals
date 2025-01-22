/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_OSCILLATOR_SEMIEXPLICIT_HH
#define TESTS_OSCILLATOR_SEMIEXPLICIT_HH

using namespace Sandals;

class OscillatorSemiExplicit : public SemiExplicit<2, 1>
{
  Real    m_m{1.0};        // Mass (kg)
  Real    m_k{1.0};        // Spring constant (N/m)
  VectorF m_ics{1.0, 0.0}; // Initial conditions

public:
  using VectorF  = typename SemiExplicit<2, 1>::VectorF;
  using MatrixA  = typename SemiExplicit<2, 1>::MatrixA;
  using TensorTA = typename SemiExplicit<2, 1>::TensorTA;
  using VectorB  = typename SemiExplicit<2, 1>::VectorB;
  using MatrixJB = typename SemiExplicit<2, 1>::MatrixJB;
  using VectorH  = typename SemiExplicit<2, 1>::VectorH;
  using MatrixJH = typename SemiExplicit<2, 1>::MatrixJH;

  OscillatorSemiExplicit() : SemiExplicit<2, 1>("OscillatorSemiExplicit") {}

  ~OscillatorSemiExplicit() {}

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
    VectorF b;
    b << x(1), -this->m_k/this->m_m*x(0);
    return b;
  }

  MatrixJB Jb_x(VectorF const &/*x*/, Real /*t*/) const override {
    MatrixJB Jb_x;
    Jb_x << 0.0, 1.0, -this->m_k/this->m_m, 0.0;
    return Jb_x;
  }

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

#endif // TESTS_OSCILLATOR_SEMIEXPLICIT_HH
