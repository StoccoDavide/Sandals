/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PROBLEMS_BASIC_IMPLICIT_HH
#define TESTS_PROBLEMS_BASIC_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

template<typename Real = double>
class BasicImplicit : public Implicit<Real, 2, 0>
{
public:
  using typename Implicit<Real, 2, 0>::VectorF;
  using typename Implicit<Real, 2, 0>::MatrixJF;
  using typename Implicit<Real, 2, 0>::VectorH;
  using typename Implicit<Real, 2, 0>::MatrixJH;

  BasicImplicit() : Implicit<Real, 2, 0>("BasicImplicit") {}

  ~BasicImplicit() {}

  VectorF F(VectorF const & x,  VectorF const &x_dot, Real const /*t*/) const override
  {
    VectorF F;
    F << x_dot(0) - x(1), x_dot(1) - 1.0;
    return F;
  }

  MatrixJF JF_x(VectorF const & /*x*/,VectorF const & /*x_dot*/, Real const /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 1) = -1.0;
    return JF_x;
  }

  MatrixJF JF_x_dot(VectorF const & /*x*/, VectorF const & /*x_dot*/, Real const /*t*/) const override
  {return MatrixJF::Identity();}

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

template<typename Real, typename Integrator>
class BasicImplicitProblem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  BasicImplicitProblem(IntegratorPtr rk)
    : Problem<Real, 2, 0, Integrator>("BasicImplicitProblem", std::make_shared<BasicImplicit<Real>>(), rk)
  {
    rk->system(this->system());
  }

  ~BasicImplicitProblem() {}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b(0) = x_ini(0);
    b(1) = x_end(1);
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
    Jb_x_end(1, 1) = 1.0;
    return Jb_x_end;
  }

  VectorF analytical_solution(Real const t) const {
    VectorF x;
    constexpr Real c0{0.0};
    constexpr Real c1{-1.0};
    x << 0.5*t*t + c1*t + c0, t + c1;
    return x;
  }

  MatrixX analytical_solution(VectorX const & t) const {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }
};

template<typename Integrator>
BasicImplicitProblem(std::shared_ptr<Integrator>)
    -> BasicImplicitProblem<typename Integrator::real_type, Integrator>;

#endif // TESTS_PROBLEMS_BASIC_IMPLICIT_HH
