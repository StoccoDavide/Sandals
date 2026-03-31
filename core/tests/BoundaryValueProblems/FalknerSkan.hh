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

#ifndef TESTS_PROBLEMS_FALKNERSKAN_HH
#define TESTS_PROBLEMS_FALKNERSKAN_HH

#include "Sandals/System/BoundaryValueProblem.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class FalknerSkanExplicit : public Explicit<Real, 3, 0> {
 public:
  using typename Explicit<Real, 3, 0>::VectorF;
  using typename Explicit<Real, 3, 0>::MatrixJF;
  using typename Explicit<Real, 3, 0>::VectorH;
  using typename Explicit<Real, 3, 0>::MatrixJH;

 private:
  Real m_beta{0.5};

 public:
  FalknerSkanExplicit() : Explicit<Real, 3, 0>("FalknerSkanExplicit") {}

  ~FalknerSkanExplicit() {}

  void beta(const Real beta) {
    this->m_beta = beta;
  }

  Real beta() const {
    return this->m_beta;
  }

  VectorF f(const VectorF &x, const Real /*t*/) const override {
    VectorF f;
    f << x(1), x(2), -x(0) * x(2) - this->m_beta * (1.0 - x(1) * x(1));
    return f;
  }

  MatrixJF Jf_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJF Jf_x;
    Jf_x << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, -x(2), 2.0 * this->m_beta * x(1),
        -x(0);
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
class FalknerSkanImplicit : public Implicit<Real, 3, 0> {
 public:
  using typename Implicit<Real, 3, 0>::VectorF;
  using typename Implicit<Real, 3, 0>::MatrixJF;
  using typename Implicit<Real, 3, 0>::VectorH;
  using typename Implicit<Real, 3, 0>::MatrixJH;

 private:
  Real m_beta{0.5};

 public:
  FalknerSkanImplicit() : Implicit<Real, 3, 0>("FalknerSkanImplicit") {}

  ~FalknerSkanImplicit() {}

  void beta(const Real beta) {
    this->m_beta = beta;
  }

  Real beta() const {
    return this->m_beta;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) - x(2),
        x_dot(2) + x(0) * x(2) + this->m_beta * (1.0 - x(1) * x(1));
    return F;
  }
  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x << 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, x(2), -2.0 * this->m_beta * x(1),
        x(0);
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
class FalknerSkanSemiExplicit : public SemiExplicit<Real, 3, 0> {
 public:
  using VectorF  = typename SemiExplicit<Real, 3, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 3, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 3, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 3, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 3, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 3, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 3, 0>::MatrixJH;

 private:
  Real m_beta{0.5};

 public:
  FalknerSkanSemiExplicit()
      : SemiExplicit<Real, 3, 0>("FalknerSkanSemiExplicit") {}

  ~FalknerSkanSemiExplicit() {}

  void beta(const Real beta) {
    this->m_beta = beta;
  }

  Real beta() const {
    return this->m_beta;
  }

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    TA_x[2].setZero();
    return TA_x;
  }

  VectorB b(const VectorF &x, const Real /*t*/) const override {
    VectorB b;
    b << x(1), x(2), -x(0) * x(2) - this->m_beta * (1.0 - x(1) * x(1));
    return b;
  }

  MatrixJB Jb_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJB Jb_x;
    Jb_x << 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, -x(2), 2.0 * this->m_beta * x(1),
        -x(0);
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
class FalknerSkanProblem : public BoundaryValueProblem<Real, 3, 0, Integrator> {
 public:
  using typename BoundaryValueProblem<Real, 3, 0, Integrator>::SystemPtr;
  using typename BoundaryValueProblem<Real, 3, 0, Integrator>::IntegratorPtr;
  using typename BoundaryValueProblem<Real, 3, 0, Integrator>::SolutionPtr;
  using typename BoundaryValueProblem<Real, 3, 0, Integrator>::VectorF;
  using typename BoundaryValueProblem<Real, 3, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 3, Eigen::Dynamic>;

  FalknerSkanProblem()
      : BoundaryValueProblem<Real, 3, 0, Integrator>(
            "FalknerSkanProblem",
            std::make_unique<System>(),
            std::make_unique<Integrator>()) {}

  ~FalknerSkanProblem() {}

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 3.0;
  }

  void beta(const Real beta) {
    static_cast<System *>(this->integrator()->system())->beta(beta);
  }

  Real beta() const {
    return static_cast<const System *>(this->integrator()->system())->beta();
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b << x_ini(0), x_ini(1), x_end(1) - 1.0;
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini(MatrixJF::Zero());
    Jb_x_ini(0, 0) = 1.0;
    Jb_x_ini(1, 1) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end(MatrixJF::Zero());
    Jb_x_end(2, 1) = 1.0;
    return Jb_x_end;
  }

  VectorF exact_solution(const Real /*t*/) const {
    return VectorF::Zero();
  }

  MatrixX exact_solution(const VectorX &t) const {
    return MatrixX::Zero(3, t.size());
  }

  static VectorF guess(const Real /*t*/) {
    return VectorF::Ones();
  }

  static MatrixX guess(const VectorX &t) {
    return MatrixX::Ones(3, t.size());
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PROBLEMS_FALKNERSKAN_HH
