/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef TESTS_GERDTS_ACADEMIC_HH
#define TESTS_GERDTS_ACADEMIC_HH

#include "DAEProblemsAD.hh"

// Academic examples: small DAE optimal control problems with exact or extrapolated optima, from M.
// Gerdts, Optimal Control of ODEs and DAEs, 2nd edition, De Gruyter, 2023 (the example numbers
// refer to the book). Conventions as in DAEProblemsAD.hh: the node values of the algebraic states
// at t = 0 are gauges, fixed by a boundary condition.
namespace TestOCP
{

  /*\
   |  Cost on the algebraic state (Example 3.1.13, from Backes): the necessary conditions written with
   |  the "naive" Hamilton function (cost + adjoints times the DAE) are not solvable, but a direct
   |  method does not care:
   |    min  1/2 x1(tf)^2 + 1/2 int_0^tf (y^2 + u^2) dt
   |    s.t. x1' = u,  x2' = -y + u,  0 = x2       (index 2, x1(0) = a, x2(0) = 0)
   |  Index 1: 0 = u - y with the invariant h = x2.
   |  Optimum: u = y = -a/(2 + tf), x1 = a (1 - t/(2 + tf)), J* = a^2/(2 + tf). The control is
   |  constant, hence the piecewise constant controls contain the exact solution.
  \*/

  template <Integer Index>
  class AlgebraicCost : public OptimalControlProblemAD<AlgebraicCost<Index>, 3, 1, Index == 2 ? 0 : 1, 3>
  {
    using Base = OptimalControlProblemAD<AlgebraicCost<Index>, 3, 1, Index == 2 ? 0 : 1, 3>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;

    Real a{1.0};

    AlgebraicCost(Real const t_end) : Base("AlgebraicCostIndex" + std::to_string(Index), 0.0, t_end) {}

    Real optimal_cost() const {return a*a/(2.0 + this->t_end());}
    Real optimal_control() const {return -a/(2.0 + this->t_end());}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = x_dot(0) - u(0);
      F(1) = x_dot(1) + x(2) - u(0);
      F(2) = Index == 2 ? x(1) : u(0) - x(2);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Index == 1) {h(0) = x(1);}
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const &) const {return 0.5*(x(2)*x(2) + u(0)*u(0));}

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return 0.5*x_end(0)*x_end(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {return BoundaryVector<T>(x_ini(0) - a, x_ini(1), x_ini(2) /* gauge */);}
  };

  /*\
   |  Bang-singular-bang control of an index-two DAE (Example 7.1.13):
   |    min  1/2 int_0^3 x2^2 dt
   |    s.t. x1' = u,  x2' = y,  0 = x1 - x2,  -1 <= u <= 1     (index 2)
   |         x1(0) = x1(3) = 1,  x2(0) = 1
   |  Index 1: 0 = u - y with the invariant h = x1 - x2.
   |  Optimum: u = y = -1, 0, 1 on [0, 1), [1, 2), [2, 3], x1 = x2 = 1 - t, 0, t - 2, J* = 1/3. With
   |  a mesh multiple of three intervals, the piecewise constant controls contain the exact solution.
  \*/

  template <Integer Index>
  class BangSingularBang : public OptimalControlProblemAD<BangSingularBang<Index>, 3, 1, Index == 2 ? 0 : 1, 4>
  {
    using Base = OptimalControlProblemAD<BangSingularBang<Index>, 3, 1, Index == 2 ? 0 : 1, 4>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    using VectorU = typename Base::VectorU;

    BangSingularBang() : Base("BangSingularBangIndex" + std::to_string(Index), 0.0, 3.0)
    {this->u_bounds(VectorU(-1.0), VectorU(1.0));}

    static Real optimal_control(Real const t) {return t < 1.0 ? -1.0 : (t < 2.0 ? 0.0 : 1.0);}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = x_dot(0) - u(0);
      F(1) = x_dot(1) - x(2);
      F(2) = Index == 2 ? x(0) - x(1) : u(0) - x(2);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Index == 1) {h(0) = x(0) - x(1);}
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const &, T const &) const {return 0.5*x(1)*x(1);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {return BoundaryVector<T>(x_ini(0) - 1.0, x_ini(1) - 1.0, x_ini(2) /* gauge */, x_end(0) - 1.0);}
  };
  /*\
   |  Smoothed control in the algebraic constraint (Examples 1.1.17 and 5.7.4, from O. Kostyukova):
   |    min  int_0^3 (2 x1 + 1/2 y^2) dt
   |    s.t. x1' = x2,  x2' = y + u1,  x3' = x4,  x4' = u2,  0 = x2 - x4     (index 2)
   |         x(0) = (1, 2/3, 1, 2/3),  -1 <= u1 <= 1,  0 <= x3 <= 200
   |  Index 1: 0 = u2 - y - u1 with the invariant h = x2 - x4.
   |  Since x3 = x1, the problem is a double integrator x1'' = v driven to the state constraint x1 >= 0,
   |  where the acceleration v = y + u1 costs 1/2 y^2 = 1/2 max(|v| - 1, 0)^2: the constraint becomes
   |  active (x1 = x2 = 0) and the controls switch to y = u1 = 0 there.
  \*/

  template <Integer Index>
  class Kostyukova : public OptimalControlProblemAD<Kostyukova<Index>, 5, 2, Index == 2 ? 0 : 1, 5>
  {
    using Base = OptimalControlProblemAD<Kostyukova<Index>, 5, 2, Index == 2 ? 0 : 1, 5>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    using VectorF = typename Base::VectorF;
    using VectorU = typename Base::VectorU;

    Kostyukova() : Base("KostyukovaIndex" + std::to_string(Index), 0.0, 3.0)
    {
      VectorF lower(VectorF::Constant(-INF)), upper(VectorF::Constant(INF));
      lower(2) = 0.0;
      upper(2) = 200.0;
      this->x_bounds(lower, upper);
      this->u_bounds(VectorU(-1.0, -INF), VectorU(1.0, INF));
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = x_dot(0) - x(1);
      F(1) = x_dot(1) - x(4) - u(0);
      F(2) = x_dot(2) - x(3);
      F(3) = x_dot(3) - u(1);
      F(4) = Index == 2 ? x(1) - x(3) : u(1) - x(4) - u(0);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Index == 1) {h(0) = x(1) - x(3);}
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const &, T const &) const {return 2.0*x(0) + 0.5*x(4)*x(4);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - 1.0, x_ini(1) - 2.0/3.0, x_ini(2) - 1.0, x_ini(3) - 2.0/3.0, x_ini(4) /* gauge */;
      return b;
    }
  };

  /*\
   |  Degenerate problem (Example 7.1.6): every feasible control is optimal.
   |    min  int_0^1 (x^2 + alpha (y - u)^2) dt
   |    s.t. x' = y - u,  0 = x,  -1 <= u <= 1          (index 2, x(0) = 0)
   |  Index 1: 0 = y - u with the invariant h = x.
   |  On the feasible set y = u and x = 0, so J = 0 for every control: the reduced Hessian vanishes in
   |  the control directions, and for alpha < 0 (the case of the book) the Hessian of the Lagrangian
   |  is indefinite. The book uses it to show that the Hamilton function must be minimized jointly in
   |  (u, y) on the set y = u.
  \*/

  template <Integer Index>
  class Degenerate : public OptimalControlProblemAD<Degenerate<Index>, 2, 1, Index == 2 ? 0 : 1, 2>
  {
    using Base = OptimalControlProblemAD<Degenerate<Index>, 2, 1, Index == 2 ? 0 : 1, 2>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    using VectorU = typename Base::VectorU;

    Real alpha{-1.0};

    Degenerate() : Base("DegenerateIndex" + std::to_string(Index), 0.0, 1.0) {this->u_bounds(VectorU(-1.0), VectorU(1.0));}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = x_dot(0) - x(1) + u(0);
      F(1) = Index == 2 ? x(0) : x(1) - u(0);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Index == 1) {h(0) = x(0);}
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const &) const
    {return x(0)*x(0) + alpha*(x(1) - u(0))*(x(1) - u(0));}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {return BoundaryVector<T>(x_ini(0), x_ini(1) /* gauge */);}
  };

  /*\
   |  Kostyukova's problem as first stated (Example 1.1.17), with the control u2 in the algebraic
   |  constraint (no smoothing):
   |    min  int_0^3 (2 x1 + 1/2 y^2) dt
   |    s.t. x1' = x2,  x2' = y + u1,  x3' = u2,  0 = x2 - u2       (index 2 at least)
   |         x(0) = (1, 2/3, 1),  -1 <= u1 <= 1,  0 <= x3 <= 200
   |  The transcription with piecewise constant controls is inconsistent: on each interval the Runge-
   |  Kutta stages of x2 must all equal the control, so x2 changes by stage derivatives proportional to
   |  A^-1 1 (A the Runge-Kutta matrix); the constant u1 absorbs their mean (the rate of change of x2),
   |  but not their spread, which remains in y and costs ~ x2'^2 independently of the step size. The
   |  discrete optima converge (at first order) to about 8.963 instead of J* = 2.8 of the equivalent
   |  smoothed problem (5.7.4): this is why the book smooths u2.
  \*/

  class KostyukovaUnsmoothed : public OptimalControlProblemAD<KostyukovaUnsmoothed, 4, 2, 0, 4>
  {
  public:
    KostyukovaUnsmoothed() : OptimalControlProblemAD("KostyukovaUnsmoothed", 0.0, 3.0)
    {
      VectorF lower(VectorF::Constant(-INF)), upper(VectorF::Constant(INF));
      lower(2) = 0.0;
      upper(2) = 200.0;
      this->x_bounds(lower, upper);
      this->u_bounds(VectorU(-1.0, -INF), VectorU(1.0, INF));
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = x_dot(0) - x(1);
      F(1) = x_dot(1) - x(3) - u(0);
      F(2) = x_dot(2) - u(1);
      F(3) = x(1) - u(1);
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const &, T const &) const {return 2.0*x(0) + 0.5*x(3)*x(3);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - 1.0, x_ini(1) - 2.0/3.0, x_ini(2) - 1.0, x_ini(3) /* gauge */;
      return b;
    }
  };

} // namespace TestOCP

#endif // TESTS_GERDTS_ACADEMIC_HH
