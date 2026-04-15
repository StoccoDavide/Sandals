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

#ifndef TESTS_PROBLEMS_SHAMPINE2_HH
#define TESTS_PROBLEMS_SHAMPINE2_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class Shampine2Explicit : public Explicit<Real, 2, 0> {
 public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_p{1.0e-5};

 public:
  Shampine2Explicit() : Explicit<Real, 2, 0>("Shampine2Explicit") {}

  ~Shampine2Explicit() {}

  void p(const Real p) {
    this->m_p = p;
  }

  Real p() const {
    return this->m_p;
  }

  VectorF f(const VectorF &x, const Real t) const override {
    VectorF f;
    f << x(1), -3.0 * this->m_p * x(0) / std::pow(this->m_p + t * t, 2.0);
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real t) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = -3.0 * this->m_p / std::pow(this->m_p + t * t, 2.0);
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
class Shampine2Implicit : public Implicit<Real, 2, 0> {
 public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_p{1.0e-5};

 public:
  Shampine2Implicit() : Implicit<Real, 2, 0>("Shampine2Implicit") {}

  ~Shampine2Implicit() {}

  void p(const Real p) {
    this->m_p = p;
  }

  Real p() const {
    return this->m_p;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real t) const override {
    VectorF F;
    F << x_dot(0) - x(1),
        x_dot(1) + 3.0 * this->m_p * x(0) / std::pow(this->m_p + t * t, 2.0);
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real t) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = 3.0 * this->m_p / std::pow(this->m_p + t * t, 2.0);
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
class Shampine2SemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;

 private:
  Real m_p{1.0e-5};

 public:
  Shampine2SemiExplicit() : SemiExplicit<Real, 2, 0>("Shampine2SemiExplicit") {}

  ~Shampine2SemiExplicit() {}

  void p(const Real p) {
    this->m_p = p;
  }

  Real p() const {
    return this->m_p;
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
    VectorB b;
    b << x(1), -3.0 * this->m_p * x(0) / std::pow(this->m_p + t * t, 2.0);
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real t) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = -3.0 * this->m_p / std::pow(this->m_p + t * t, 2.0);
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
class Shampine2Problem : public BVP<Real, 2, 0, Integrator> {
 public:
  using typename BVP<Real, 2, 0, Integrator>::VectorF;
  using typename BVP<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  Shampine2Problem()
      : BVP<Real, 2, 0, Integrator>("Shampine2Problem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}

  ~Shampine2Problem() {}

  static Real time_start() {
    return -0.1;
  }

  static Real time_end() {
    return 0.1;
  }

  void p(const Real p) {
    static_cast<System *>(this->integrator()->system())->p(p);
  }

  Real p() const {
    return static_cast<const System *>(this->integrator()->system())->p();
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b << x_ini(0) + 0.1 / std::sqrt(this->p() + 0.01),
        x_end(0) - 0.1 / std::sqrt(this->p() + 0.01);
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

  VectorF exact_solution(const Real t) const {
    VectorF x;
    x << t / std::sqrt(this->p() + t * t),
        this->p() / std::pow(this->p() + t * t, 1.5);
    return x;
  }

  MatrixX exact_solution(const VectorX &t) const {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) =
          Shampine2Problem<Real, System, Integrator>::exact_solution(t(i));
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

#endif  // TESTS_PROBLEMS_SHAMPINE2_HH
