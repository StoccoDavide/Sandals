/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PROBLEMS_HYPERSENSITIVE1_HH
#define TESTS_PROBLEMS_HYPERSENSITIVE1_HH

#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"
#include "Sandals/Problem.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class Hypersensitive1Explicit : public Explicit<Real, 2, 0>
{
public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

  Hypersensitive1Explicit() : Explicit<Real, 2, 0>("Hypersensitive1Explicit") {}

  ~Hypersensitive1Explicit() {}

  VectorF f(VectorF const & x, Real const /*t*/) const override
  {
    VectorF f;
    f << x(0)*x(0)*x(0) - 0.5*x(1), -2.0*x(0) + 3.0*x(0)*x(1)*x(1);
    return f;
  }

  MatrixJF Jf_x(VectorF const & x, Real const /*t*/) const override {
    MatrixJF Jf_x;
    Jf_x <<
      3.0*x(0)*x(0),        -0.5,
      -2.0 + 3.0*x(1)*x(1), 6.0*x(0)*x(1);
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real = double>
class Hypersensitive1Implicit : public Implicit<Real, 2, 0>
{
public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

  Hypersensitive1Implicit() : Implicit<Real, 2, 0>("Hypersensitive1Implicit") {}

  ~Hypersensitive1Implicit() {}

  VectorF F(VectorF const & x, VectorF const & x_dot, Real const /*t*/) const override
  {
    VectorF F;
    F << x_dot(0) - x(0)*x(0)*x(0) + 0.5*x(1), x_dot(1) + 2.0*x(0) - 3.0*x(0)*x(1)*x(1);
    return F;
  }

  MatrixJF JF_x(VectorF const & x, VectorF const & /*x_dot*/, Real const /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x <<
      -3.0*x(0)*x(0),      0.5,
      2.0 - 3.0*x(1)*x(1), -6.0*x(0)*x(1);
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
class Hypersensitive1SemiExplicit : public SemiExplicit<Real, 2, 0>
{
public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;

  Hypersensitive1SemiExplicit() : SemiExplicit<Real, 2, 0>("Hypersensitive1SemiExplicit") {}

  ~Hypersensitive1SemiExplicit() {}

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
    b << x(0)*x(0)*x(0) - 0.5*x(1), -2.0*x(0) + 3.0*x(0)*x(1)*x(1);
    return b;
  }

  MatrixJB Jb_x(VectorF const & x, Real const /*t*/) const override {
    MatrixJB Jb_x;
    Jb_x <<
      3.0*x(0)*x(0),        -0.5,
      -2.0 + 3.0*x(1)*x(1), 6.0*x(0)*x(1);
    return Jb_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename Real, typename System, typename Integrator>
class Hypersensitive1Problem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  Hypersensitive1Problem()
    : Problem<Real, 2, 0, Integrator>("Hypersensitive1Problem", std::make_unique<System>(),
      std::make_unique<Integrator>()) {}

  ~Hypersensitive1Problem() {}

  static Real time_start() {return 0.0;}

  static Real time_end() {return 1.0;}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b <<
      x_ini(0) - 1.0,
      x_end(0) - 1.5;
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

  VectorF analytical_solution(Real const /*t*/) const {return VectorF::Zero();}

  MatrixX analytical_solution(VectorX const & t) const {return MatrixX::Zero(2, t.size());}

  static VectorF guess(Real const /*t*/) {return VectorF::Ones();}

  static MatrixX guess(VectorX const & t) {return MatrixX::Ones(2, t.size());}

};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // TESTS_PROBLEMS_HYPERSENSITIVE1_HH
