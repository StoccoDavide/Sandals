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

#ifndef TESTS_PROBLEMS_BVPT9_HH
#define TESTS_PROBLEMS_BVPT9_HH

#include "Sandals/System/BoundaryValueProblem.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT9Explicit : public Explicit<Real, 2, 0> {
 public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT9Explicit() : Explicit<Real, 2, 0>("BVPT9Explicit") {}

  ~BVPT9Explicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }

  Real lambda() const {
    return this->m_lambda;
  }

  VectorF f(const VectorF &x, const Real t) const override {
    const Real den{this->m_lambda + t * t};
    VectorF f;
    f << x(1), -(4.0 * t * x(1) + 2.0 * x(0)) / den;
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real t) const override {
    const Real den{this->m_lambda + t * t};
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = -2.0 / den;
    Jf_x(1, 1) = -4.0 * t / den;
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
class BVPT9Implicit : public Implicit<Real, 2, 0> {
 public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT9Implicit() : Implicit<Real, 2, 0>("BVPT9Implicit") {}

  ~BVPT9Implicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }

  Real lambda() const {
    return this->m_lambda;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real t) const override {
    const Real den{this->m_lambda + t * t};
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) + (4.0 * t * x(1) + 2.0 * x(0)) / den;
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real t) const override {
    const Real den{this->m_lambda + t * t};
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = 2.0 / den;
    JF_x(1, 1) = 4.0 * t / den;
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
class BVPT9SemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_lambda{1.0e-3};

 public:
  BVPT9SemiExplicit() : SemiExplicit<Real, 2, 0>("BVPT9SemiExplicit") {}

  ~BVPT9SemiExplicit() {}

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
    const Real den{this->m_lambda + t * t};
    VectorB b;
    b << x(1), -(4.0 * t * x(1) + 2.0 * x(0)) / den;
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real t) const override {
    const Real den{this->m_lambda + t * t};
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = -2.0 / den;
    Jb_x(1, 1) = -4.0 * t / den;
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
class BVPT9Problem : public BoundaryValueProblem<Real, 2, 0, Integrator> {
 public:
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::SystemPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::VectorF;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  BVPT9Problem()
      : BoundaryValueProblem<Real, 2, 0, Integrator>(
            "BVPT9Problem",
            std::make_unique<System>(),
            std::make_unique<Integrator>()) {}

  ~BVPT9Problem() {}

  static Real time_start() {
    return -1.0;
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
    const Real eta{1.0 / (1.0 + this->lambda())};
    VectorF b;
    b << x_ini(0) - eta, x_end(0) - eta;
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
    return true;
  }

  Real exact_solution(const Real t) const {
    return 1.0 / (t * t + this->lambda());
  }

  VectorX exact_solution(const VectorX &t) const {
    VectorX x(t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x(i) = BVPT9Problem<Real, System, Integrator>::exact_solution(t(i));
    }
    return x;
  }

  VectorF guess(const Real /*t*/) {
    return VectorF::Zero();
  }

  MatrixX guess(const VectorX &t) {
    return MatrixX::Zero(2, t.size());
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PROBLEMS_BVPT9_HH
