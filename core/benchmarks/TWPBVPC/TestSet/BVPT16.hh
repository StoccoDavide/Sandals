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

#ifndef TESTS_PROBLEMS_BVPT16_HH
#define TESTS_PROBLEMS_BVPT16_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT16Explicit : public Explicit<Real, 2, 0> {
 public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};
  Real m_lambda{1.0e-3};

 public:
  BVPT16Explicit() : Explicit<Real, 2, 0>("BVPT16Explicit") {}
  ~BVPT16Explicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }
  Real lambda() const {
    return this->m_lambda;
  }

  VectorF f(const VectorF &x, const Real /*t*/) const override {
    VectorF f;
    f << x(1), -x(0) * pi * pi / (4.0 * this->m_lambda);
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = -pi * pi / (4.0 * this->m_lambda);
    return Jf_x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT16Implicit : public Implicit<Real, 2, 0> {
 public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};
  Real m_lambda{1.0e-3};

 public:
  BVPT16Implicit() : Implicit<Real, 2, 0>("BVPT16Implicit") {}
  ~BVPT16Implicit() {}

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
    F << x_dot(0) - x(1), x_dot(1) + x(0) * pi * pi / (4.0 * this->m_lambda);
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = pi * pi / (4.0 * this->m_lambda);
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
class BVPT16SemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using typename SemiExplicit<Real, 2, 0>::VectorF;
  using typename SemiExplicit<Real, 2, 0>::MatrixA;
  using typename SemiExplicit<Real, 2, 0>::TensorTA;
  using typename SemiExplicit<Real, 2, 0>::VectorB;
  using typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using typename SemiExplicit<Real, 2, 0>::VectorH;
  using typename SemiExplicit<Real, 2, 0>::MatrixJH;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};
  Real m_lambda{1.0e-3};

 public:
  BVPT16SemiExplicit() : SemiExplicit<Real, 2, 0>("BVPT16SemiExplicit") {}
  ~BVPT16SemiExplicit() {}

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

  VectorB b(const VectorF &x, const Real /*t*/) const override {
    VectorB b;
    b << x(1), -x(0) * pi * pi / (4.0 * this->m_lambda);
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = -pi * pi / (4.0 * this->m_lambda);
    return Jb_x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class BVPT16Problem : public BVP<Real, 2, 0, Integrator> {
 public:
  using typename BVP<Real, 2, 0, Integrator>::VectorF;
  using typename BVP<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};

 public:
  BVPT16Problem()
      : BVP<Real, 2, 0, Integrator>("BVPT16Problem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}
  ~BVPT16Problem() {}

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
    b << x_ini(0), x_end(0) - std::sin(pi / (2.0 * std::sqrt(this->lambda())));
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
    return std::sin(pi * t / (2.0 * std::sqrt(this->lambda())));
  }

  VectorX exact_solution(const VectorX &t) const {
    VectorX x(t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x(i) = BVPT16Problem<Real, System, Integrator>::exact_solution(t(i));
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

#endif  // TESTS_PROBLEMS_BVPT16_HH
