/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef D8015E73_D378_4028_AA6F_9784F643F64F
#define D8015E73_D378_4028_AA6F_9784F643F64F

#ifndef TESTS_PROBLEMS_BASIC_SEMIEXPLICIT_HH
#define TESTS_PROBLEMS_BASIC_SEMIEXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

template<typename Real = double>
class BasicSemiExplicit : public SemiExplicit<Real, 2, 0>
{
public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;

  BasicSemiExplicit() : SemiExplicit<Real, 2, 0>("BasicSemiExplicit") {}

  ~BasicSemiExplicit() {}

  MatrixA A(VectorF const & /*x*/, Real const /*t*/) const override
  {
    MatrixA A;
    A.setIdentity();
    return A;
  }

  TensorTA TA_x(VectorF const & /*x*/, Real const /*t*/) const override
  {
    TensorTA TA_x(2);
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(VectorF const &x, Real const /*t*/) const override
  {
    VectorF b;
    b <<  x(1), 1.0;
    return b;
  }

  MatrixJB Jb_x(VectorF const & /*x*/, Real const /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 1) = 1.0;
    return Jb_x;
  }

  //virtual VectorH h(VectorF const & x, Real const t) const = 0;
  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

template<typename Real, typename Integrator>
class BasicSemiExplicitProblem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  BasicSemiExplicitProblem(IntegratorPtr rk)
    : Problem<Real, 2, 0, Integrator>("BasicSemiExplicitProblem", std::make_shared<BasicSemiExplicit<Real>>(), rk)
  {
    rk->system(this->system());
  }

  ~BasicSemiExplicitProblem() {}

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
BasicSemiExplicitProblem(std::shared_ptr<Integrator>)
    -> BasicSemiExplicitProblem<typename Integrator::real_type, Integrator>;

#endif // TESTS_PROBLEMS_BASIC_SEMIEXPLICIT_HH


#endif /* D8015E73_D378_4028_AA6F_9784F643F64F */
