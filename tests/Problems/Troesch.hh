/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PROBLEMS_TROESCH0_HH
#define TESTS_PROBLEMS_TROESCH0_HH

#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"
#include "Sandals/Problem.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class TroeschExplicit : public Explicit<Real, 2, 0>
{
public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

private:
  Real m_lambda{1.0e-2};

public:
  TroeschExplicit() : Explicit<Real, 2, 0>("TroeschExplicit") {}

  ~TroeschExplicit() {}

  void lambda(Real const lambda) {this->m_lambda = lambda;}

  Real lambda() const {return this->m_lambda;}

  VectorF f(VectorF const & x, Real const /*t*/) const override
  {
    VectorF f;
    f << x(1), this->m_lambda*std::sinh(this->m_lambda*x(0));
    return f;
  }

  MatrixJF Jf_x(VectorF const & x, Real const /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = this->m_lambda*this->m_lambda*std::cosh(this->m_lambda*x(0));
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


template<typename Real = double>
class TroeschImplicit : public Implicit<Real, 2, 0>
{
public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

private:
  Real m_lambda{1.0e-2};

public:
  TroeschImplicit() : Implicit<Real, 2, 0>("TroeschImplicit") {}

  ~TroeschImplicit() {}

  void lambda(Real const lambda) {this->m_lambda = lambda;}

  Real lambda() const {return this->m_lambda;}

  VectorF F(VectorF const & x, VectorF const & x_dot, Real const /*t*/) const override
  {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) - this->m_lambda*std::sinh(this->m_lambda*x(0));
    return F;
  }

  MatrixJF JF_x(VectorF const & x, VectorF const & /*x_dot*/, Real const /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = -this->m_lambda*this->m_lambda*std::cosh(this->m_lambda*x(0));
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
class TroeschSemiExplicit : public SemiExplicit<Real, 2, 0>
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
  Real m_lambda{1.0e-2};

public:
  TroeschSemiExplicit() : SemiExplicit<Real, 2, 0>("TroeschSemiExplicit") {}

  ~TroeschSemiExplicit() {}

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
    VectorF b;
    b << x(1), this->m_lambda*std::sinh(this->m_lambda*x(0));
    return b;
  }

  MatrixJB Jb_x(VectorF const & x, Real const /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = this->m_lambda*this->m_lambda*std::cosh(this->m_lambda*x(0));
    return Jb_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real, typename System, typename Integrator>
class TroeschProblem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  TroeschProblem()
    : Problem<Real, 2, 0, Integrator>("TroeschProblem", std::make_unique<System>(),
      std::make_unique<Integrator>()) {}

  ~TroeschProblem() {}

  static Real time_start() {return 0.0;}

  static Real time_end() {return 1.0;}

  void lambda(Real const lambda) {static_cast<System*>(this->integrator()->system())->lambda(lambda);}

  Real lambda() const {return static_cast<const System*>(this->integrator()->system())->lambda();}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b << x_ini(0), x_end(0) - 1.0;
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

  VectorF analytical_solution(Real const t) const {
    VectorF x;
    Real lambda{this->lambda()};
    Real c{1.0 - std::sinh(lambda)/lambda};
    x <<
      std::sinh(t*lambda)/lambda + t*c,
      std::cosh(t*lambda) + c;
    return x;
  }

  MatrixX analytical_solution(VectorX const & t) const {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) = TroeschProblem<Real, System, Integrator>::analytical_solution(t(i));
    }
    return x;
  }

  static VectorF guess(Real const /*t*/) {return VectorF::Zero();}

  static MatrixX guess(VectorX const & t) {return MatrixX::Zero(2, t.size());}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // TESTS_PROBLEMS_TROESCH0_HH
