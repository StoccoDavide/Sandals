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

#ifndef TESTS_PROBLEMS_BVPT33_HH
#define TESTS_PROBLEMS_BVPT33_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT33Explicit : public Explicit<Real, 6, 0> {
 public:
  using typename Explicit<Real, 6, 0>::VectorF;
  using typename Explicit<Real, 6, 0>::MatrixJF;
  using typename Explicit<Real, 6, 0>::VectorH;
  using typename Explicit<Real, 6, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT33Explicit() : Explicit<Real, 6, 0>("BVPT33Explicit") {}
  ~BVPT33Explicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }
  Real lambda() const {
    return this->m_lambda;
  }

  VectorF f(const VectorF &x, const Real /*t*/) const override {
    VectorF f;
    f << x(1), (x(0) * x(3) - x(2) * x(1)) / this->m_lambda, x(3), x(4), x(5),
        (-x(2) * x(5) - x(0) * x(1)) / this->m_lambda;
    return f;
  }

  MatrixJF Jf_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = x(3) / this->m_lambda;
    Jf_x(1, 1) = -x(2) / this->m_lambda;
    Jf_x(1, 2) = -x(1) / this->m_lambda;
    Jf_x(1, 3) = x(0) / this->m_lambda;
    Jf_x(2, 3) = 1.0;
    Jf_x(3, 4) = 1.0;
    Jf_x(4, 5) = 1.0;
    Jf_x(5, 0) = -x(1) / this->m_lambda;
    Jf_x(5, 1) = -x(0) / this->m_lambda;
    Jf_x(5, 2) = -x(5) / this->m_lambda;
    Jf_x(5, 5) = -x(2) / this->m_lambda;
    return Jf_x;
  }

  VectorH h(const VectorF & /*x*/, const Real /*t*/) const override {
    return VectorH::Zero();
  }

  MatrixJH Jh_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJH::Zero();
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT33Implicit : public Implicit<Real, 6, 0> {
 public:
  using typename Implicit<Real, 6, 0>::VectorF;
  using typename Implicit<Real, 6, 0>::MatrixJF;
  using typename Implicit<Real, 6, 0>::VectorH;
  using typename Implicit<Real, 6, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT33Implicit() : Implicit<Real, 6, 0>("BVPT33Implicit") {}
  ~BVPT33Implicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }
  Real lambda() const {
    return this->m_lambda;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F << x_dot(0) - x(1),
        x_dot(1) - (x(0) * x(3) - x(2) * x(1)) / this->m_lambda,
        x_dot(2) - x(3), x_dot(3) - x(4), x_dot(4) - x(5),
        x_dot(5) - (-x(2) * x(5) - x(0) * x(1)) / this->m_lambda;
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = -x(3) / this->m_lambda;
    JF_x(1, 1) = x(2) / this->m_lambda;
    JF_x(1, 2) = x(1) / this->m_lambda;
    JF_x(1, 3) = -x(0) / this->m_lambda;
    JF_x(2, 3) = -1.0;
    JF_x(3, 4) = -1.0;
    JF_x(4, 5) = -1.0;
    JF_x(5, 0) = x(1) / this->m_lambda;
    JF_x(5, 1) = x(0) / this->m_lambda;
    JF_x(5, 2) = x(5) / this->m_lambda;
    JF_x(5, 5) = x(2) / this->m_lambda;
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }

  VectorH h(const VectorF & /*x*/, const Real /*t*/) const override {
    return VectorH::Zero();
  }

  MatrixJH Jh_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJH::Zero();
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT33SemiExplicit : public SemiExplicit<Real, 6, 0> {
 public:
  using VectorF  = typename SemiExplicit<Real, 6, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 6, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 6, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 6, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 6, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 6, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 6, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT33SemiExplicit() : SemiExplicit<Real, 6, 0>("BVPT33SemiExplicit") {}
  ~BVPT33SemiExplicit() {}

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
    TA_x[2].setZero();
    TA_x[3].setZero();
    TA_x[4].setZero();
    TA_x[5].setZero();
    return TA_x;
  }

  VectorB b(const VectorF &x, const Real /*t*/) const override {
    VectorB b;
    b << x(1), (x(0) * x(3) - x(2) * x(1)) / this->m_lambda, x(3), x(4), x(5),
        (-x(2) * x(5) - x(0) * x(1)) / this->m_lambda;
    return b;
  }

  MatrixJB Jb_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = x(3) / this->m_lambda;
    Jb_x(1, 1) = -x(2) / this->m_lambda;
    Jb_x(1, 2) = -x(1) / this->m_lambda;
    Jb_x(1, 3) = x(0) / this->m_lambda;
    Jb_x(2, 3) = 1.0;
    Jb_x(3, 4) = 1.0;
    Jb_x(4, 5) = 1.0;
    Jb_x(5, 0) = -x(1) / this->m_lambda;
    Jb_x(5, 1) = -x(0) / this->m_lambda;
    Jb_x(5, 2) = -x(5) / this->m_lambda;
    Jb_x(5, 5) = -x(2) / this->m_lambda;
    return Jb_x;
  }

  VectorH h(const VectorF & /*x*/, const Real /*t*/) const override {
    return VectorH::Zero();
  }

  MatrixJH Jh_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJH::Zero();
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class BVPT33Problem : public BVP<Real, 6, 0, Integrator> {
 public:
  using typename BVP<Real, 6, 0, Integrator>::VectorF;
  using typename BVP<Real, 6, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 6, Eigen::Dynamic>;

  BVPT33Problem()
      : BVP<Real, 6, 0, Integrator>("BVPT33Problem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}
  ~BVPT33Problem() {}

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
    b << x_ini(0) + 1.0, x_ini(2), x_ini(3), x_end(0) - 1.0, x_end(2), x_end(3);
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini(MatrixJF::Zero());
    Jb_x_ini(0, 0) = 1.0;
    Jb_x_ini(1, 2) = 1.0;
    Jb_x_ini(2, 3) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end(MatrixJF::Zero());
    Jb_x_end(3, 0) = 1.0;
    Jb_x_end(4, 2) = 1.0;
    Jb_x_end(5, 3) = 1.0;
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
    return VectorF::Zero();
  }

  MatrixX guess(const VectorX &t) {
    return MatrixX::Zero(6, t.size());
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PROBLEMS_BVPT33_HH
