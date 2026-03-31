/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef TESTS_SYSTEMS_OSCILLATOR_EXPLICIT_HH
#define TESTS_SYSTEMS_OSCILLATOR_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;

template <typename Real = double>
class OscillatorExplicit : public Explicit<Real, 2, 1> {
 public:
  using VectorF  = typename Explicit<Real, 2, 1>::VectorF;
  using MatrixJF = typename Explicit<Real, 2, 1>::MatrixJF;
  using VectorH  = typename Explicit<Real, 2, 1>::VectorH;
  using MatrixJH = typename Explicit<Real, 2, 1>::MatrixJH;
  using VectorX  = Eigen::Matrix<Real, 2, 1>;
  using MatrixX  = Eigen::Matrix<Real, Eigen::Dynamic, 2>;

 private:
  Real m_m{1.0};  // Mass (kg)
  Real m_k{1.0};  // Spring constant (N/m)

 public:
  OscillatorExplicit() : Explicit<Real, 2, 1>("OscillatorExplicit") {}

  ~OscillatorExplicit() {}

  VectorF f(const VectorF &x, const Real /*t*/) const override {
    VectorF f;
    f << x(1), -this->m_k / this->m_m * x(0);
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real /*t*/) const override {
    MatrixJF Jf_x;
    Jf_x << 0.0, 1.0, -this->m_k / this->m_m, 0.0;
    return Jf_x;
  }

  Real energy(const VectorF &x) const {
    return this->m_m / 2.0 * x(1) * x(1) + this->m_k / 2.0 * x(0) * x(0);
  }

  VectorH h(const VectorF &x, const Real /*t*/) const override {
    VectorH h;
    h << this->energy(x) - this->energy(OscillatorExplicit::ics());
    return h;
  }

  MatrixJH Jh_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJH Jh_x;
    Jh_x << this->m_k * x(0), this->m_m * x(1);
    return Jh_x;
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }

  VectorF exact_solution(Real t) const {
    VectorF x;
    Real tmp{std::sqrt(this->m_k / this->m_m)};
    x << std::cos(tmp * t), -tmp * std::sin(tmp * t);
    return x;
  }

  MatrixX exact_solution(const VectorX &t) const {
    MatrixX x(2, t.size());
    for (int i = 0; i < t.size(); ++i) {
      x.col(i) = this->exact_solution(t(i));
    }
    return x;
  }

  static VectorF ics() {
    return VectorF::Unit(0, 2);
  }
};

#endif  // TESTS_SYSTEMS_OSCILLATOR_EXPLICIT_HH
