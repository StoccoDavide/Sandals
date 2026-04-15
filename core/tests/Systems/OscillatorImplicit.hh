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

#ifndef TESTS_SYSTEMS_OSCILLATOR_IMPLICIT_HH
#define TESTS_SYSTEMS_OSCILLATOR_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

template <typename Real = double>
class OscillatorImplicit : public Implicit<Real, 2, 1> {
 public:
  using VectorF  = typename Implicit<Real, 2, 1>::VectorF;
  using MatrixJF = typename Implicit<Real, 2, 1>::MatrixJF;
  using VectorH  = typename Implicit<Real, 2, 1>::VectorH;
  using MatrixJH = typename Implicit<Real, 2, 1>::MatrixJH;
  using VectorX  = Eigen::Matrix<Real, 2, 1>;
  using MatrixX  = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  Real m_m{1.0};  // Mass (kg)
  Real m_k{1.0};  // Spring constant (N/m)

 public:
  OscillatorImplicit() : Implicit<Real, 2, 1>("OscillatorImplicit") {}

  ~OscillatorImplicit() {}

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) + this->m_k / this->m_m * x(0);
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x;
    JF_x << -1.0, 1.0, this->m_k / this->m_m, 0.0;
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }

  Real energy(const VectorF &x) const {
    return this->m_m / 2.0 * x(1) * x(1) + this->m_k / 2.0 * x(0) * x(0);
  }

  VectorH h(const VectorF &x, const Real /*t*/) const override {
    VectorH h;
    h << this->energy(x) - this->energy(OscillatorImplicit::ics());
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
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) = OscillatorImplicit::exact_solution(t(i));
    }
    return x;
  }

  static VectorF ics() {
    return VectorF::Unit(0, 2);
  }
};

#endif  // TESTS_SYSTEMS_OSCILLATOR_IMPLICIT_HH
