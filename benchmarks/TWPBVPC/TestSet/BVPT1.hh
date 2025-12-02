/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef TESTS_PROBLEMS_BVPT1_HH
#define TESTS_PROBLEMS_BVPT1_HH

#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"
#include "Sandals/BoundaryValueProblem.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/**
  The problem is

  $$
  \lambda z^{\prime \prime}=z, \quad z(0)=1, \quad z(1)=0,
  $$

  with

  $$
  z \in \mathbb{R}, \quad t \in[0,1] .
  $$


  We write this problem in first order form by defining $y_1=z$, and $y_2=z^{\prime}$, yielding a system of differential equations of the form

  $$
  \binom{y_1}{y_2}^{\prime}=\binom{y_2}{\frac{1}{\lambda} f\left(y_1\right)},
  $$

  where

  $$
  f(z)=z,
  $$

  and

  $$
  \left(y_1, y_2\right)^T \in \mathbb{R}^2, \quad t \in[0,1] .
  $$


  The boundary conditions are obtained from

  $$
  \left(\begin{array}{ll}
  1 & 0 \\
  0 & 0
  \end{array}\right)\binom{y_1(0)}{y_2(0)}+\left(\begin{array}{ll}
  0 & 0 \\
  1 & 0
  \end{array}\right)\binom{y_1(1)}{y_2(1)}=\binom{1}{0} .
  $$


  Exact solution

  $$
  z(t)=(\exp (-t / \sqrt{\lambda})-\exp ((t-2) / \sqrt{\lambda})) /(1-\exp (-2 / \sqrt{\lambda})) .
  $$
*/

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class BVPT1Explicit : public Explicit<Real, 2, 0>
{
public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

private:
  Real m_lambda{1.0e-3};

public:
  BVPT1Explicit() : Explicit<Real, 2, 0>("BVPT1Explicit") {}

  ~BVPT1Explicit() {}

  void lambda(Real const lambda) {this->m_lambda = lambda;}

  Real lambda() const {return this->m_lambda;}

  VectorF f(VectorF const & x, Real const /*t*/) const override
  {
    VectorF f;
    f << x(1), x(0)/this->m_lambda;
    return f;
  }

  MatrixJF Jf_x(VectorF const & /*x*/, Real const /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = 1.0/this->m_lambda;
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


template<typename Real = double>
class BVPT1Implicit : public Implicit<Real, 2, 0>
{
public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

private:
  Real m_lambda{1.0e-3};

public:
  BVPT1Implicit() : Implicit<Real, 2, 0>("BVPT1Implicit") {}

  ~BVPT1Implicit() {}

  void lambda(Real const lambda) {this->m_lambda = lambda;}

  Real lambda() const {return this->m_lambda;}

  VectorF F(VectorF const & x, VectorF const & x_dot, Real const /*t*/) const override
  {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) - x(0)/this->m_lambda;
    return F;
  }

  MatrixJF JF_x(VectorF const & /*x*/, VectorF const & /*x_dot*/, Real const /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = -1.0/this->m_lambda;
    return JF_x;
  }

  MatrixJF JF_x_dot(VectorF const & /*x*/, VectorF const & /*x_dot*/, Real const /*t*/) const override
  {return MatrixJF::Identity();}

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class BVPT1SemiExplicit : public SemiExplicit<Real, 2, 0>
{
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
  BVPT1SemiExplicit() : SemiExplicit<Real, 2, 0>("BVPT1SemiExplicit") {}

  ~BVPT1SemiExplicit() {}

  void lambda(Real const lambda) {this->m_lambda = lambda;}

  Real lambda() const {return this->m_lambda;}

  MatrixA A(VectorF const & /*x*/, Real const /*t*/) const override
  {return MatrixA::Identity();}

  TensorTA TA_x(VectorF const & /*x*/, Real const /*t*/) const override
  {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(VectorF const & x, Real const /*t*/) const override
  {
    VectorB b;
    b << x(1), x(0)/this->m_lambda;
    return b;
  }

  MatrixJB Jb_x(VectorF const & /*x*/, Real const /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = 1.0/this->m_lambda;
    return Jb_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real, typename System, typename Integrator>
class BVPT1Problem : public BoundaryValueProblem<Real, 2, 0, Integrator>
{
public:
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::SystemPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::VectorF;
  using typename BoundaryValueProblem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  BVPT1Problem()
    : BoundaryValueProblem<Real, 2, 0, Integrator>("BVPT1Problem", std::make_unique<System>(),
      std::make_unique<Integrator>()) {}

  ~BVPT1Problem() {}

  static Real time_start() {return 0.0;}

  static Real time_end() {return 1.0;}

  void lambda(Real const lambda) {static_cast<System*>(this->integrator()->system())->lambda(lambda);}

  Real lambda() const {return static_cast<const System*>(this->integrator()->system())->lambda();}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b << x_ini(0) - 1.0, x_end(0);
    return b;
  }

  MatrixJF Jb_x_ini(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
  {
    MatrixJF Jb_x_ini(MatrixJF::Zero());
    Jb_x_ini(0, 0) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
  {
    MatrixJF Jb_x_end(MatrixJF::Zero());
    Jb_x_end(1, 0) = 1.0;
    return Jb_x_end;
  }

  Real exact_solution(Real const t) const {
    Real sqrt_lambda{std::sqrt(this->lambda())};
    return (std::exp(-t/sqrt_lambda) - std::exp((t - 2.0)/sqrt_lambda)) /
      (1.0 - std::exp(-2/sqrt_lambda));
  }

  VectorX exact_solution(VectorX const & t) const {
    VectorX x(t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x(i) = BVPT1Problem<Real, System, Integrator>::exact_solution(t(i));
    }
    return x;
  }

  static VectorF guess(Real const /*t*/) {return VectorF::Zero();}

  static MatrixX guess(VectorX const & t) {return MatrixX::Zero(2, t.size());}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // TESTS_PROBLEMS_BVPT1_HH
