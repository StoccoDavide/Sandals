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

#ifndef TESTS_PROBLEMS_BVPT24_HH
#define TESTS_PROBLEMS_BVPT24_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT24Explicit : public Explicit<Real, 2, 0> {
 public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT24Explicit() : Explicit<Real, 2, 0>("BVPT24Explicit") {}
  ~BVPT24Explicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }
  Real lambda() const {
    return this->m_lambda;
  }

  VectorF f(const VectorF &x, const Real t) const override {
    const Real Ax{1.0 + t * t};
    const Real Apx{2.0 * t};
    const Real Ga{1.4};
    const Real y0{x(0)};
    const Real y1{x(1)};
    VectorF f;
    f << y1, (((1.0 + Ga) / 2.0 - this->m_lambda * Apx) * y0 * y1 - y1 / y0 -
              (Apx / Ax) * (1.0 - (Ga - 1.0) * y0 * y0 / 2.0)) /
                 (this->m_lambda * Ax * y0);
    return f;
  }

  MatrixJF Jf_x(const VectorF &x, const Real t) const override {
    const Real Ax{1.0 + t * t};
    const Real Apx{2.0 * t};
    const Real Ga{1.4};
    const Real y0{x(0)};
    const Real y1{x(1)};
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = (2.0 * y1 / (y0 * y0 * y0) + Apx / (Ax * y0 * y0) +
                  Apx * (Ga - 1.0) / (2.0 * Ax)) /
                 (this->m_lambda * Ax);
    Jf_x(1, 1) = (((1.0 + Ga) / 2.0 - this->m_lambda * Apx) * y0 - 1.0 / y0) /
                 (this->m_lambda * Ax * y0);
    return Jf_x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT24Implicit : public Implicit<Real, 2, 0> {
 public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT24Implicit() : Implicit<Real, 2, 0>("BVPT24Implicit") {}
  ~BVPT24Implicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }
  Real lambda() const {
    return this->m_lambda;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real t) const override {
    const Real Ax{1.0 + t * t};
    const Real Apx{2.0 * t};
    const Real Ga{1.4};
    const Real y0{x(0)};
    const Real y1{x(1)};
    const Real f1{(((1.0 + Ga) / 2.0 - this->m_lambda * Apx) * y0 * y1 -
                   y1 / y0 - (Apx / Ax) * (1.0 - (Ga - 1.0) * y0 * y0 / 2.0)) /
                  (this->m_lambda * Ax * y0)};
    VectorF F;
    F << x_dot(0) - y1, x_dot(1) - f1;
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real t) const override {
    const Real Ax{1.0 + t * t};
    const Real Apx{2.0 * t};
    const Real Ga{1.4};
    const Real y0{x(0)};
    const Real y1{x(1)};
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = -(2.0 * y1 / (y0 * y0 * y0) + Apx / (Ax * y0 * y0) +
                   Apx * (Ga - 1.0) / (2.0 * Ax)) /
                 (this->m_lambda * Ax);
    JF_x(1, 1) = -(((1.0 + Ga) / 2.0 - this->m_lambda * Apx) * y0 - 1.0 / y0) /
                 (this->m_lambda * Ax * y0);
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT24SemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using typename SemiExplicit<Real, 2, 0>::VectorF;
  using typename SemiExplicit<Real, 2, 0>::MatrixA;
  using typename SemiExplicit<Real, 2, 0>::TensorTA;
  using typename SemiExplicit<Real, 2, 0>::VectorB;
  using typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using typename SemiExplicit<Real, 2, 0>::VectorH;
  using typename SemiExplicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT24SemiExplicit() : SemiExplicit<Real, 2, 0>("BVPT24SemiExplicit") {}
  ~BVPT24SemiExplicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }
  Real lambda() const {
    return this->m_lambda;
  }

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(const VectorF &x, const Real t) const override {
    const Real Ax{1.0 + t * t};
    const Real Apx{2.0 * t};
    const Real Ga{1.4};
    const Real y0{x(0)};
    const Real y1{x(1)};
    VectorB b;
    b << y1, (((1.0 + Ga) / 2.0 - this->m_lambda * Apx) * y0 * y1 - y1 / y0 -
              (Apx / Ax) * (1.0 - (Ga - 1.0) * y0 * y0 / 2.0)) /
                 (this->m_lambda * Ax * y0);
    return b;
  }

  MatrixJB Jb_x(const VectorF &x, const Real t) const override {
    const Real Ax{1.0 + t * t};
    const Real Apx{2.0 * t};
    const Real Ga{1.4};
    const Real y0{x(0)};
    const Real y1{x(1)};
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = (2.0 * y1 / (y0 * y0 * y0) + Apx / (Ax * y0 * y0) +
                  Apx * (Ga - 1.0) / (2.0 * Ax)) /
                 (this->m_lambda * Ax);
    Jb_x(1, 1) = (((1.0 + Ga) / 2.0 - this->m_lambda * Apx) * y0 - 1.0 / y0) /
                 (this->m_lambda * Ax * y0);
    return Jb_x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class BVPT24Problem : public BVP<Real, 2, 0, Integrator> {
 public:
  using typename BVP<Real, 2, 0, Integrator>::VectorF;
  using typename BVP<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  BVPT24Problem()
      : BVP<Real, 2, 0, Integrator>("BVPT24Problem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}
  ~BVPT24Problem() {}

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 1.0;
  }

  void lambda(const Real lambda) {
    static_cast<System *>(this->integrator()->system())->lambda(lambda);
  }
  Real lambda() const {
    return static_cast<const System *>(this->integrator()->system())->lambda();
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b << x_ini(0) - 0.9129, x_end(0) - 0.375;
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini(MatrixJF::Zero());
    Jb_x_ini(0, 0) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end(MatrixJF::Zero());
    Jb_x_end(1, 0) = 1.0;
    return Jb_x_end;
  }

  static constexpr bool has_exact_solution() {
    return false;
  }

  Real exact_solution(const Real /*t*/) const {
    return 0.0;
  }

  VectorX exact_solution(const VectorX &t) const {
    return VectorX::Zero(t.size());
  }

  VectorF guess(const Real /*t*/) {
    return VectorF::Constant(0.5);
  }

  MatrixX guess(const VectorX &t) {
    return MatrixX::Constant(2, t.size(), 0.5);
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PROBLEMS_BVPT24_HH
