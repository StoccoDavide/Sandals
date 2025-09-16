/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_MAZZIA_EXPLICIT_HH
#define TESTS_MAZZIA_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;

template<typename Real = double>
class MazziaExplicit : public Explicit<Real, 4, 0>
{
public:
  using VectorF  = typename Explicit<Real, 4, 0>::VectorF;
  using MatrixJF = typename Explicit<Real, 4, 0>::MatrixJF;
  using VectorH  = typename Explicit<Real, 4, 0>::VectorH;
  using MatrixJH = typename Explicit<Real, 4, 0>::MatrixJH;
  using VectorX  = Eigen::Matrix<Real, 4, 0>;
  using MatrixX  = Eigen::Matrix<Real, Eigen::Dynamic, 2>;

  MazziaExplicit() : Explicit<Real, 4, 0>("MazziaExplicit") {}

  ~MazziaExplicit() {}

  VectorF f(VectorF const & x, Real /*t*/) const override
  {
    VectorF f;
    f <<
       x(1),
      -x(0) - x(0)*x(0)*x(0) - x(3),
      -x(0) + x(3)*(1.0 - 3.0*x(0)*x(0)),
      -x(1) - x(3);
    return f;
  }

  MatrixJF Jf_x(VectorF const & x, Real /*t*/) const override {
    MatrixJF Jf_x;
    Jf_x <<
      0.0, 1.0, 0.0, 0.0,
      -3.0*x(0)*x(0) - 1.0, 0.0, 1.0, 0.0,
      -1.0 + x(3)*(1.0 - 6.0*x(0)), 0.0, 0.0, x(3)*(1.0 - 3.0*x(0)*x(0)),
      0.0, -1.0, 0.0, -1.0;
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real /*t*/) const override {return true;}

};

template<typename Real = double, Integer S>
class MazziaExplicitProblem : public Problem<Real, 4, 0, S>
{
public:
  using typename Problem<Real, 4, 0, S>::SystemPtr;
  using typename Problem<Real, 4, 0, S>::RungeKuttaPtr;
  using typename Problem<Real, 4, 0, S>::SolutionPtr;
  using typename Problem<Real, 4, 0, S>::VectorX;
  using typename Problem<Real, 4, 0, S>::VectorF;
  using typename Problem<Real, 4, 0, S>::MatrixJF;

  MazziaExplicitProblem(RungeKuttaPtr rk)
    : Problem<Real, 4, 0, S>("MazziaExplicitProblem", std::make_shared<MazziaExplicit<Real>>(), rk)
  {
    rk->system(this->system());
  }

  ~MazziaExplicitProblem() {}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b(0) = x_ini(0) - 1.0;
    b(1) = x_ini(1) - 0.0;
    b(2) = x_end(0) - 0.75;
    b(3) = x_end(1) + 0.0;
    return b;
  }

  MatrixJF Jb_x_ini(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
  {
    MatrixJF Jb_x_ini = MatrixJF::Zero();
    Jb_x_ini(0, 0) = 1.0;
    Jb_x_ini(1, 1) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
  {
    MatrixJF Jb_x_end = MatrixJF::Zero();
    Jb_x_end(2, 0) = 1.0;
    Jb_x_end(3, 1) = 1.0;
    return Jb_x_end;
  }
};

#endif // TESTS_MAZZIA_EXPLICIT_HH
