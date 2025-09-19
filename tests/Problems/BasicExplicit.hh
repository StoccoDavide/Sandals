/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PROBLEMS_BASIC_EXPLICIT_HH
#define TESTS_PROBLEMS_BASIC_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;

template<typename Real = double>
class BasicExplicit : public Explicit<Real, 2, 0>
{
public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

  BasicExplicit() : Explicit<Real, 2, 0>("BasicExplicit") {}

  ~BasicExplicit() {}

  VectorF f(VectorF const & x, Real const /*t*/) const override
  {
    VectorF f;
    f << x(1), 1.0;
    return f;
  }

  MatrixJF Jf_x(VectorF const & /*x*/, Real const /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

template<typename Real, typename Integrator>
class BasicExplicitProblem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  BasicExplicitProblem(IntegratorPtr rk)
    : Problem<Real, 2, 0, Integrator>("BasicExplicitProblem", std::make_shared<BasicExplicit<Real>>(), rk)
  {
    rk->system(this->system());
  }

  ~BasicExplicitProblem() {}

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

  VectorF ics() const {
    VectorF ics;
    ics << 0.0, 0.5;
    return ics;
  }

  VectorF analytical_solution(Real const t) const {
    VectorF x;
    x << 0.5*t*t - t, t - 1.0;
    return x;
  }

  MatrixX analytical_solution(VectorX const & t) const {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }
};

template<typename Integrator>
BasicExplicitProblem(std::shared_ptr<Integrator>)
    -> BasicExplicitProblem<typename Integrator::real_type, Integrator>;

#endif // TESTS_PROBLEMS_BASIC_EXPLICIT_HH
