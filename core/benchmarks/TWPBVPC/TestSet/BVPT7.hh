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

#ifndef TESTS_PROBLEMS_BVPT7_HH
#define TESTS_PROBLEMS_BVPT7_HH

#include "Sandals/System/BoundaryValueProblem.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/**
  The problem is

  $$
  \lambda z^{\prime \prime} = -t z^{\prime} + z - (1 + \lambda \pi^2) \cos(\pi
  t)
  - \pi t \sin(\pi t),
  $$

  with

  $$
  z(-1) = -1, \quad z(1) = 1, \quad t \in [-1, 1].
  $$

  Exact solution:

  $$
  z(t) = \cos(\pi t) + t + \frac{t
  \operatorname{erf}\!\left(\frac{t}{\sqrt{2\lambda}}\right)
  + \sqrt{\frac{2\lambda}{\pi}} e^{-t^2/(2\lambda)}}
  {\operatorname{erf}\!\left(\frac{1}{\sqrt{2\lambda}}\right)
  + \sqrt{\frac{2\lambda}{\pi}} e^{-1/(2\lambda)}}.
  $$
*/

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BVPT7Explicit : public Explicit<Real, 2, 0> {
 public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};
  Real m_lambda{1.0e-3};

 public:
  BVPT7Explicit() : Explicit<Real, 2, 0>("BVPT7Explicit") {}

  ~BVPT7Explicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }

  Real lambda() const {
    return this->m_lambda;
  }

  VectorF f(const VectorF &x, const Real t) const override {
    const Real tmp{static_cast<Real>(EIGEN_PI) * t};
    VectorF f;
    f << x(1),
        (-t * x(1) + x(0) - (1.0 + this->m_lambda * pi * pi) * std::cos(tmp) -
         pi * t * std::sin(tmp)) /
            this->m_lambda;
    return f;
  }

  MatrixJF Jf_x(const VectorF & /*x*/, const Real t) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = 1.0 / this->m_lambda;
    Jf_x(1, 1) = -t / this->m_lambda;
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
class BVPT7Implicit : public Implicit<Real, 2, 0> {
 public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};
  Real m_lambda{1.0e-3};

 public:
  BVPT7Implicit() : Implicit<Real, 2, 0>("BVPT7Implicit") {}

  ~BVPT7Implicit() {}

  void lambda(const Real lambda) {
    this->m_lambda = lambda;
  }

  Real lambda() const {
    return this->m_lambda;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real t) const override {
    const Real tmp{static_cast<Real>(EIGEN_PI) * t};
    VectorF F;
    F << x_dot(0) - x(1),
        x_dot(1) - (-t * x(1) + x(0) -
                    (1.0 + this->m_lambda * pi * pi) * std::cos(tmp) -
                    pi * t * std::sin(tmp)) /
                       this->m_lambda;
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real t) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = -1.0 / this->m_lambda;
    JF_x(1, 1) = t / this->m_lambda;
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
class BVPT7SemiExplicit : public SemiExplicit<Real, 2, 0> {
 public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};
  Real m_lambda{1.0e-3};

 public:
  BVPT7SemiExplicit() : SemiExplicit<Real, 2, 0>("BVPT7SemiExplicit") {}

  ~BVPT7SemiExplicit() {}

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
    const Real tmp{pi * t};
    VectorB b;
    b << x(1),
        (-t * x(1) + x(0) - (1.0 + this->m_lambda * pi * pi) * std::cos(tmp) -
         pi * t * std::sin(tmp)) /
            this->m_lambda;
    return b;
  }

  MatrixJB Jb_x(const VectorF & /*x*/, const Real t) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = 1.0 / this->m_lambda;
    Jb_x(1, 1) = -t / this->m_lambda;
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
class BVPT7Problem : public BoundaryValueProblem<Real, 2, 0, Integrator> {
 public:
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::SystemPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::VectorF;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

 private:
  static constexpr Real pi{static_cast<Real>(EIGEN_PI)};

 public:
  BVPT7Problem()
      : BoundaryValueProblem<Real, 2, 0, Integrator>(
            "BVPT7Problem",
            std::make_unique<System>(),
            std::make_unique<Integrator>()) {}

  ~BVPT7Problem() {}

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
    VectorF b;
    b << x_ini(0) + 1.0, x_end(0) - 1.0;
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
    const Real lambda{this->lambda()};
    const Real A{t * t / (2.0 * lambda)};
    const Real B{1.0 / (2.0 * lambda)};
    const Real C{std::exp(-A)};
    const Real D{std::exp(-B)};
    const Real S{std::sqrt(2.0 * lambda)};
    return std::cos(pi * t) + t +
           (t * std::erf(t / S) + std::sqrt(2.0 * lambda / pi) * C) /
               (std::erf(1.0 / S) + std::sqrt(2.0 * lambda / pi) * D);
  }

  VectorX exact_solution(const VectorX &t) const {
    VectorX x(t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x(i) = BVPT7Problem<Real, System, Integrator>::exact_solution(t(i));
    }
    return x;
  }

  static VectorF guess(const Real /*t*/) {
    return VectorF::Zero();
  }

  static MatrixX guess(const VectorX &t) {
    return MatrixX::Zero(2, t.size());
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PROBLEMS_BVPT7_HH
