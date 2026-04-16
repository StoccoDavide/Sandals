/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_SYSTEMS_OSCILLATOR_SEMIEXPLICIT_HH
#define TESTS_SYSTEMS_OSCILLATOR_SEMIEXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

template <typename Real = double>
class OscillatorSemiExplicit : public SemiExplicit<Real, 2, 1> {
 public:
  using typename SemiExplicit<Real, 2, 1>::VectorF;
  using typename SemiExplicit<Real, 2, 1>::MatrixA;
  using typename SemiExplicit<Real, 2, 1>::TensorTA;
  using typename SemiExplicit<Real, 2, 1>::VectorB;
  using typename SemiExplicit<Real, 2, 1>::MatrixJB;
  using typename SemiExplicit<Real, 2, 1>::VectorH;
  using typename SemiExplicit<Real, 2, 1>::MatrixJH;
  using typename SemiExplicit<Real, 2, 1>::TensorTH;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  Real m_m{1.0};  // Mass (kg)
  Real m_k{1.0};  // Spring constant (N/m)

 public:
  OscillatorSemiExplicit()
      : SemiExplicit<Real, 2, 1>("OscillatorSemiExplicit") {}

  ~OscillatorSemiExplicit() {}

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(const VectorF &x, const Real /*t*/) const override {
    VectorF b;
    b << x(1), -this->m_k / this->m_m * x(0);
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real /*t*/) const override {
    MatrixJB Jb_x;
    Jb_x << 0.0, 1.0, -this->m_k / this->m_m, 0.0;
    return Jb_x;
  }

  Real energy(const VectorF &x) const {
    return this->m_m / 2.0 * x(1) * x(1) + this->m_k / 2.0 * x(0) * x(0);
  }

  VectorH h(const VectorF &x, const Real /*t*/) const override {
    VectorH h;
    h << this->energy(x) - this->energy(OscillatorSemiExplicit::ics());
    return h;
  }

  MatrixJH Jh_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJH Jh_x;
    Jh_x << this->m_k * x(0), this->m_m * x(1);
    return Jh_x;
  }

  TensorTH Th_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTH Th_x;
    for (MatrixJH &m : Th_x) {
      m.setZero();
    }
    Th_x[0](0, 0) = this->m_k;
    Th_x[1](0, 1) = this->m_m;
    return Th_x;
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
      x.col(i) = OscillatorSemiExplicit::exact_solution(t(i));
    }
    return x;
  }

  static VectorF ics() {
    return VectorF::Unit(0, 2);
  }
};

#endif  // TESTS_SYSTEMS_OSCILLATOR_SEMIEXPLICIT_HH
