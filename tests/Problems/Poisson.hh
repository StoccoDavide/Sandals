/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PROBLEMS_POISSON_HH
#define TESTS_PROBLEMS_POISSON_HH

#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"
#include "Sandals/Problem.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class PoissonExplicit : public Explicit<Real, 2, 0>
{
public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

  PoissonExplicit() : Explicit<Real, 2, 0>("PoissonExplicit") {}

  ~PoissonExplicit() {}

  VectorF f(VectorF const & x, Real const /*t*/) const override
  {
    VectorF f;
    f << x(1), -5.0*std::sin(10.0*M_PI*x(0)) - std::sin(2.0*M_PI*x(0));
    return f;
  }

  MatrixJF Jf_x(VectorF const & x, Real const /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = -50.0*M_PI*std::cos(10.0*M_PI*x(0)) - 2.0*M_PI*std::cos(2.0*M_PI*x(0));
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class PoissonImplicit : public Implicit<Real, 2, 0>
{
public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

  PoissonImplicit() : Implicit<Real, 2, 0>("PoissonImplicit") {}

  ~PoissonImplicit() {}

  VectorF F(VectorF const & x, VectorF const & x_dot, Real const /*t*/) const override
  {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) + 5.0*std::sin(10.0*M_PI*x(0)) + std::sin(2.0*M_PI*x(0));
    return F;
  }

  MatrixJF JF_x(VectorF const & x, VectorF const & /*x_dot*/, Real const /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    JF_x(1, 0) = 50.0*M_PI*std::cos(10.0*M_PI*x(0)) + 2.0*M_PI*std::cos(2.0*M_PI*x(0));
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
class PoissonSemiExplicit : public SemiExplicit<Real, 2, 0>
{
public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;

  PoissonSemiExplicit() : SemiExplicit<Real, 2, 0>("PoissonSemiExplicit") {}

  ~PoissonSemiExplicit() {}

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
    b << x(1), -5.0*std::sin(10.0*M_PI*x(0)) - std::sin(2.0*M_PI*x(0));
    return b;
  }

  MatrixJB Jb_x(VectorF const & x, Real const /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    Jb_x(1, 0) = -50.0*M_PI*std::cos(10.0*M_PI*x(0)) - 2.0*M_PI*std::cos(2.0*M_PI*x(0));
    return Jb_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real, typename System, typename Integrator>
class PoissonProblem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  PoissonProblem()
    : Problem<Real, 2, 0, Integrator>("PoissonProblem", std::make_unique<System>(), std::make_unique<Integrator>()) {}

  ~PoissonProblem() {}

  static Real time_start() {return 0.0;}

  static Real time_end() {return 1.0;}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b(0) = x_ini(0);
    b(1) = x_end(0);
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

  static VectorF analytical_solution(Real const t) {
    VectorF x;
    x <<
      std::sin(2.0*M_PI*t)/(4.0*M_PI*M_PI) + std::sin(10.0*M_PI*t)/(20.0*M_PI*M_PI),
      std::cos(2.0*M_PI*t)/(2.0*M_PI) + std::cos(10.0*M_PI*t)/(2.0*M_PI);
    return x;
  }

  static MatrixX analytical_solution(VectorX const & t) {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) = PoissonProblem<Real, System, Integrator>::analytical_solution(t(i));
    }
    return x;
  }

  static VectorF guess(Real const /*t*/) {return VectorF::Ones();}

  static MatrixX guess(VectorX const & t) {return MatrixX::Ones(2, t.size());}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // TESTS_PROBLEMS_POISSON_HH
