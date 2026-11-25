/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_GERDTS_PROBLEMS_AD_HH
#define TESTS_GERDTS_PROBLEMS_AD_HH

#include "DAEProblemsAD.hh"

// Optimal control problems of DAEs from M. Gerdts, Optimal Control of ODEs and DAEs, 2nd edition,
// De Gruyter, 2023 (the example numbers refer to the book). Each DAE is written as in the book and,
// where meaningful, index-reduced with the lost constraints as invariants (same conventions as in
// DAEProblemsAD.hh: the node values of the algebraic states at t = 0 are gauges, fixed by a boundary
// condition).
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
   |  Pendulum in the Gear-Gupta-Leimkuhler stabilized formulation (Examples 1.1.19, 4.6.5, 6.1.9):
   |    min  int_0^tf u^2 dt
   |    s.t. x1' = x3 - 2 x1 y2,       x2' = x4 - 2 x2 y2,
   |         x3' = -2 x1 y1 + u x2,    x4' = -g - 2 x2 y1 - u x1,
   |         0 = x1^2 + x2^2 - 1,      0 = x1 x3 + x2 x4          (index 2, y2 = 0 on the solutions)
   |         x(0) = (x1_ini, x2_ini, 0, 0),  x1(tf) = x3(tf) = 0,  x2(tf) <= 0 (bottom)
   |  The same problem as GerdtsPendulum<Index> (Example 1.1.9), with both the position and the
   |  velocity constraints kept in the DAE through the additional multiplier y2.
  \*/

  class GerdtsPendulumGGL : public OptimalControlProblemAD<GerdtsPendulumGGL, 6, 1, 0, 9>
  {
  public:
    Real g{9.81};
    Real x1_ini{1.0};
    Real x2_ini{0.0};

    GerdtsPendulumGGL(Real const t_end = 3.0) : OptimalControlProblemAD("GerdtsPendulumGGL", 0.0, t_end) {}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = x_dot(0) - x(2) + 2.0*x(0)*x(5);
      F(1) = x_dot(1) - x(3) + 2.0*x(1)*x(5);
      F(2) = x_dot(2) + 2.0*x(0)*x(4) - u(0)*x(1);
      F(3) = x_dot(3) + g + 2.0*x(1)*x(4) + u(0)*x(0);
      F(4) = x(0)*x(0) + x(1)*x(1) - 1.0;
      F(5) = x(0)*x(2) + x(1)*x(3);
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - x1_ini, x_ini(1) - x2_ini, x_ini(2), x_ini(3), x_ini(4) /* gauge */, x_ini(5) /* gauge */,
           x_end(0), x_end(2), x_end(1) /* <= 0 */;
      return b;
    }

    VectorX boundary_conditions_lower() const override
    {
      VectorX lower(VectorX::Zero(9));
      lower(8) = -INF;
      return lower;
    }
  };

  /*\
   |  Ski driver (Example 5.7.3): a skier of mass m on the plane z = x tan(psi) steers with the angle
   |  delta (its rate u is a control) and pushes with the force F, and reaches the final position in
   |  minimum time. States (x, y, z, vx, vy, vz, delta, lambda, T), with the time scaled by T:
   |    x' = T vx,  y' = T vy,  z' = T vz,  delta' = T u,  T' = 0,
   |    m vx' = T (tan(psi) lambda + F cos(delta) cos(psi)),
   |    m vy' = T F sin(delta),
   |    m vz' = T (-m g - lambda + F cos(delta) sin(psi)),
   |    0 = z - x tan(psi)                                              (index 3)
   |    x(0) = y(0) = z(0) = delta(0) = 0,  v(0) = (-5, 5, -5 tan(psi)),
   |    x(tf) = -100,  y(tf) = -20,  vy(tf) = 0,  |u| <= 10,  |F| <= 200
   |  Index reduction (the lost constraints become the invariants):
   |    index 2: 0 = vz - vx tan(psi),    h = (z - x tan(psi))
   |    index 1: 0 = vz' - vx' tan(psi),  h = (z - x tan(psi), vz - vx tan(psi))
   |  The normal force is constant, lambda = -m g cos(psi)^2. The book reports tf = 8.426448607129171
   |  (m = 80, psi = 10 deg, 401 grid points), with two singular arcs of the steering rate u.
  \*/

  template <Integer Index>
  class SkiDriver : public OptimalControlProblemAD<SkiDriver<Index>, 9, 2, 3 - Index, 11>
  {
    using Base = OptimalControlProblemAD<SkiDriver<Index>, 9, 2, 3 - Index, 11>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    using VectorF = typename Base::VectorF;
    using VectorU = typename Base::VectorU;

    static constexpr Real m{80.0};
    static constexpr Real g{9.81};
    static constexpr Real REFERENCE_TIME{8.426448607129171};
    Real const tan_psi{std::tan(10.0*M_PI/180.0)};
    Real const cos_psi{std::cos(10.0*M_PI/180.0)};
    Real const sin_psi{std::sin(10.0*M_PI/180.0)};

    SkiDriver() : Base("SkiDriverIndex" + std::to_string(Index), 0.0, 1.0)
    {
      VectorF lower(VectorF::Constant(-INF)), upper(VectorF::Constant(INF));
      lower(8) = 1.0; // the duration
      this->x_bounds(lower, upper);
      this->u_bounds(VectorU(-10.0, -200.0), VectorU(10.0, 200.0));
    }

    Real normal_force() const {return -m*g*cos_psi*cos_psi;}

    // Consistent initial states, with the duration guess T
    VectorF initial_states(Real const T) const
    {
      VectorF x;
      x << 0.0, 0.0, 0.0, -5.0, 5.0, -5.0*tan_psi, 0.0, this->normal_force(), T;
      return x;
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      T const lambda{x(7)}, duration{x(8)}, steering{u(0)}, force{u(1)};
      StateVector<T> F;
      F(0) = x_dot(0) - duration*x(3);
      F(1) = x_dot(1) - duration*x(4);
      F(2) = x_dot(2) - duration*x(5);
      F(3) = m*x_dot(3) - duration*(tan_psi*lambda + force*cos(x(6))*cos_psi);
      F(4) = m*x_dot(4) - duration*force*sin(x(6));
      F(5) = m*x_dot(5) - duration*(-m*g - lambda + force*cos(x(6))*sin_psi);
      F(6) = x_dot(6) - duration*steering;
      if constexpr (Index == 3) {F(7) = x(2) - x(0)*tan_psi;}
      if constexpr (Index == 2) {F(7) = x(5) - x(3)*tan_psi;}
      if constexpr (Index == 1) {F(7) = x_dot(5) - x_dot(3)*tan_psi;}
      F(8) = x_dot(8);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Index < 3) {h(0) = x(2) - x(0)*tan_psi;}
      if constexpr (Index < 2) {h(1) = x(5) - x(3)*tan_psi;}
      return h;
    }

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return x_end(8);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      VectorF const start{this->initial_states(0.0)};
      BoundaryVector<T> b;
      for (Integer i{0}; i < 8; ++i) {b(i) = x_ini(i) - start(i);} // lambda(0) is a gauge
      b(8)  = x_end(0) + 100.0;
      b(9)  = x_end(1) + 20.0;
      b(10) = x_end(4);
      return b;
    }
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
   |  Semi-batch chemical reactor (Example 8.2.8): reactions A + B -> C and B + C -> D, fed with B at
   |  the rate FB and cooled with the power Q, maximizing the product C:
   |    min  -MC(tf) + 1e-2 int_0^20 (FB^2 + Q^2) dt
   |    s.t. MA' = -V r1,  MB' = FB - V (r1 + r2),  MC' = V (r1 - r2),  MD' = V r2,
   |         H' = 20 FB - Q + V (r1 + 75 r2),
   |         0 = H - sum_i Mi (alpha_i (TR - Tref) + beta_i/2 (TR^2 - Tref^2))     (index 1)
   |  with the reactor temperature TR algebraic, the volume V = sum_i Mi/rho_i, the concentrations
   |  Ci = Mi/V, and the rates r1 = A1 exp(-E1/TR) CA CB, r2 = A2 exp(-E2/TR) CB CC.
  \*/

  class ChemicalReactor : public OptimalControlProblemAD<ChemicalReactor, 6, 2, 0, 6>
  {
  public:
    static constexpr Real T_REF{298.0};
    static constexpr Real ALPHA[4]{0.1723, 0.2, 0.16, 0.155};
    static constexpr Real BETA[4]{0.000473, 0.0005, 0.00055, 0.000323};
    static constexpr Real RHO[4]{11250.0, 16000.0, 10400.0, 10000.0};
    static constexpr Real A1{0.8}, A2{0.002}, E1{3000.0}, E2{2400.0};

    ChemicalReactor() : OptimalControlProblemAD("ChemicalReactor", 0.0, 20.0)
    {
      VectorF nominal;
      nominal << 1.0e4, 10.0, 1.0e3, 1.0, 1.0e5, 100.0; // masses, energy, temperature
      this->nominal(nominal, VectorU(10.0, 0.01));
    }

    // The temperature TR from the energy balance, a quadratic equation (positive root)
    static Real temperature(Eigen::Vector<Real, 4> const & mass, Real const H)
    {
      Real a{0.0}, b{0.0};
      for (Integer i{0}; i < 4; ++i) {a += mass(i)*BETA[i]/2.0; b += mass(i)*ALPHA[i];}
      Real const c{-H - a*T_REF*T_REF - b*T_REF};
      return (-b + std::sqrt(b*b - 4.0*a*c))/(2.0*a);
    }

    // Consistent initial states
    static VectorF initial_states()
    {
      VectorF x;
      x << 9000.0, 0.0, 0.0, 0.0, 152509.97, 0.0;
      x(5) = temperature(x.head<4>(), x(4));
      return x;
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      T const TR{x(5)}, FB{u(0)}, Q{u(1)};
      T V{0.0}, energy{0.0};
      for (Integer i{0}; i < 4; ++i) {
        V += x(i)/RHO[i];
        energy += x(i)*(ALPHA[i]*(TR - T_REF) + BETA[i]/2.0*(TR*TR - T_REF*T_REF));
      }
      T const CA{x(0)/V}, CB{x(1)/V}, CC{x(2)/V};
      T const r1{A1*exp(-E1/TR)*CA*CB}, r2{A2*exp(-E2/TR)*CB*CC};
      StateVector<T> F;
      F(0) = x_dot(0) + V*r1;
      F(1) = x_dot(1) - FB + V*(r1 + r2);
      F(2) = x_dot(2) - V*(r1 - r2);
      F(3) = x_dot(3) - V*r2;
      F(4) = x_dot(4) - 20.0*FB + Q - V*(r1 + 75.0*r2);
      F(5) = x(4) - energy;
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 1.0e-2*(u(0)*u(0) + u(1)*u(1));}

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return -x_end(2);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {
      VectorF const start{initial_states()};
      BoundaryVector<T> b;
      for (Integer i{0}; i < 6; ++i) {b(i) = x_ini(i) - start(i);} // TR(0) is a gauge
      return b;
    }
  };

  /*\
   |  Emergency landing of a hypersonic flight system (Examples 1.1.10, 5.4.6, 6.1.10): after an engine
   |  failure above Bayreuth, glide to the altitude of 500 m with the maximum range, controlled by the
   |  lift coefficient CL and the bank angle mu. States (v, gamma, chi, h, Lambda, Theta, T) (velocity,
   |  inclination, azimuth, altitude, latitude, longitude, and the duration T scaling the time):
   |    min  -((Lambda(tf) - Lambda(0))/Lambda(0))^2 - ((Theta(tf) - Theta(0))/Theta(0))^2
   |    s.t. the equations of motion over the rotating Earth (EmergencyLandingModel::dynamics),
   |         0.01 <= CL <= 0.18326,  |mu| <= pi/2,  dynamic pressure q(v, h) <= 60000,
   |         h(0) = 33900 (or h(0) >= 33900 with free_altitude),  h(tf) = 500.
   |  The book reports tf = 726.57 and the cost -0.765 (Example 6.1.10, 201 grid points, classic
   |  Runge-Kutta, piecewise constant controls), with the dynamic pressure constraint active at the
   |  normalized time 0.19.
  \*/

  struct EmergencyLandingModel
  {
    static constexpr Real AREA{305.0}, R0{6.371e6}, CD0{0.017}, K{2.0}, RHO0{1.249512}, BETA{1.0/6900.0};
    static constexpr Real G0{9.80665}, OMEGA{7.27e-5}, MASS{115000.0}, Q_MAX{60000.0};
    static constexpr Real CL_MIN{0.01}, CL_MAX{0.18326};
    static constexpr Real V0{2150.5452900}, GAMMA0{0.1520181770}, CHI0{2.2689279889}, H0{33900.0};
    static constexpr Real LAMBDA0{0.8651597102}, THETA0{0.1980948701}, H_END{500.0};

    // Dynamic pressure
    template <typename T>
    static T pressure(T const & v, T const & h) {return 0.5*RHO0*exp(-BETA*h)*v*v;}

    // Time derivatives of (v, gamma, chi, h, Lambda, Theta)
    template <typename T>
    static Eigen::Vector<T, 6> dynamics(Eigen::Vector<T, 6> const & x, T const & CL, T const & mu)
    {
      T const v{x(0)}, gamma{x(1)}, chi{x(2)}, h{x(3)}, Lambda{x(4)};
      T const R{R0 + h}, g{G0*(R0/R)*(R0/R)}, q{pressure(v, h)};
      T const L{q*AREA*CL}, D{q*AREA*(CD0 + K*CL*CL)};
      T const w2{OMEGA*OMEGA*cos(Lambda)};
      Eigen::Vector<T, 6> f;
      f(0) = -D/MASS - g*sin(gamma) + w2*(sin(gamma)*cos(Lambda) - cos(gamma)*sin(chi)*sin(Lambda))*R;
      f(1) = L*cos(mu)/(MASS*v) - (g/v - v/R)*cos(gamma) + 2.0*OMEGA*cos(chi)*cos(Lambda)
           + w2*(sin(gamma)*sin(chi)*sin(Lambda) + cos(gamma)*cos(Lambda))*R/v;
      f(2) = L*sin(mu)/(MASS*v*cos(gamma)) - v/R*cos(gamma)*cos(chi)*tan(Lambda)
           + 2.0*OMEGA*(sin(chi)*cos(Lambda)*tan(gamma) - sin(Lambda)) - w2*sin(Lambda)*cos(chi)*R/(v*cos(gamma));
      f(3) = v*sin(gamma);
      f(4) = v/R*cos(gamma)*sin(chi);
      f(5) = v/R*cos(gamma)*cos(chi)/cos(Lambda);
      return f;
    }

    // The relative range (to be maximized)
    template <typename T>
    static T range(T const & Lambda, T const & Theta)
    {
      T const dLambda{(Lambda - LAMBDA0)/LAMBDA0}, dTheta{(Theta - THETA0)/THETA0};
      return dLambda*dLambda + dTheta*dTheta;
    }
  };

  class EmergencyLanding : public OptimalControlProblemAD<EmergencyLanding, 7, 2, 0, 7, 1>, public EmergencyLandingModel
  {
  public:
    static constexpr Real REFERENCE_TIME{726.57}, REFERENCE_COST{-0.765};
    bool free_altitude{false}; // h(0) >= 33900 instead of h(0) = 33900

    EmergencyLanding() : OptimalControlProblemAD("EmergencyLanding", 0.0, 1.0)
    {
      VectorF lower(VectorF::Constant(-INF)), upper(VectorF::Constant(INF));
      lower(6) = 100.0; // the duration
      this->x_bounds(lower, upper);
      this->u_bounds(VectorU(CL_MIN, -M_PI/2.0), VectorU(CL_MAX, M_PI/2.0));
      VectorF nominal;
      nominal << 1000.0, 1.0, 1.0, 10000.0, 1.0, 1.0, 100.0; // m/s, rad, rad, m, rad, rad, s
      this->nominal(nominal, VectorU(0.1, 1.0));
    }

    static VectorF initial_states(Real const T)
    {
      VectorF x;
      x << V0, GAMMA0, CHI0, H0, LAMBDA0, THETA0, T;
      return x;
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      Eigen::Vector<T, 6> const f{dynamics<T>(x.template head<6>(), u(0), u(1))};
      StateVector<T> F;
      F << x_dot.template head<6>() - x(6)*f, x_dot(6);
      return F;
    }

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return -range(x_end(4), x_end(5));}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - V0, x_ini(1) - GAMMA0, x_ini(2) - CHI0, x_ini(3) - H0, x_ini(4) - LAMBDA0,
           x_ini(5) - THETA0, x_end(3) - H_END;
      return b;
    }

    VectorX boundary_conditions_upper() const override
    {
      VectorX upper(VectorX::Zero(7));
      if (free_altitude) {upper(3) = INF;}
      return upper;
    }

    template <typename T>
    PathVector<T> path(StateVector<T> const & x, ControlVector<T> const &, T const &) const
    {return PathVector<T>(pressure(x(0), x(3))/Q_MAX - 1.0);}
  };

  /*\
   |  Emergency landing on a prescribed path (Example 1.1.10): the ground track is the straight line
   |  Lambda + Theta = Lambda(0) + Theta(0), and the bank angle mu becomes an algebraic state, of index 3.
   |  States (v, gamma, chi, h, Lambda, Theta, mu, T), from the altitude 33900 with the initial azimuth
   |  along the path (tan(chi) = -1/cos(Lambda)), with the dynamic pressure constraint. Index reduction (the lost constraints become the invariants):
   |    index 3: 0 = Lambda + Theta - Lambda(0) - Theta(0)
   |    index 2: 0 = sin(chi) cos(Lambda) + cos(chi)               (Lambda' + Theta' = 0, times R cos(Lambda)/(v cos(gamma)))
   |    index 1: 0 = (cos(chi) cos(Lambda) - sin(chi)) chi' - sin(chi) sin(Lambda) Lambda'
   |  The book reports the consistent initial values chi(0) = -0.9955 and mu(0) = 0.3695.
  \*/

  template <Integer Index>
  class EmergencyLandingPath
    : public OptimalControlProblemAD<EmergencyLandingPath<Index>, 8, 1, 3 - Index, 8, 1>, public EmergencyLandingModel
  {
    using Base = OptimalControlProblemAD<EmergencyLandingPath<Index>, 8, 1, 3 - Index, 8, 1>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    template <typename T> using PathVector      = typename Base::template PathVector<T>;
    using VectorF = typename Base::VectorF;
    using VectorU = typename Base::VectorU;

    static constexpr Real REFERENCE_CHI0{-0.9955}, REFERENCE_MU0{0.3695};

    EmergencyLandingPath() : Base("EmergencyLandingPathIndex" + std::to_string(Index), 0.0, 1.0)
    {
      VectorF lower(VectorF::Constant(-INF)), upper(VectorF::Constant(INF));
      lower(7) = 100.0; // the duration
      this->x_bounds(lower, upper);
      this->u_bounds(VectorU(CL_MIN), VectorU(CL_MAX));
      VectorF nominal;
      nominal << 1000.0, 1.0, 1.0, 10000.0, 1.0, 1.0, 1.0, 100.0; // m/s, rad, rad, m, rad, rad, rad, s
      this->nominal(nominal, VectorU(0.1));
    }

    // The azimuth along the path: tan(chi) = -1/cos(Lambda), heading south-east
    static Real path_azimuth(Real const Lambda) {return std::atan(-1.0/std::cos(Lambda));}

    static VectorF initial_states(Real const T)
    {
      VectorF x;
      x << V0, GAMMA0, path_azimuth(LAMBDA0), H0, LAMBDA0, THETA0, 0.0, T;
      return x;
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      Eigen::Vector<T, 6> const f{dynamics<T>(x.template head<6>(), u(0), x(6))};
      T const chi{x(2)}, Lambda{x(4)};
      StateVector<T> F;
      F.template head<6>() = x_dot.template head<6>() - x(7)*f;
      if constexpr (Index == 3) {F(6) = Lambda + x(5) - LAMBDA0 - THETA0;}
      if constexpr (Index == 2) {F(6) = sin(chi)*cos(Lambda) + cos(chi);}
      if constexpr (Index == 1) {F(6) = (cos(chi)*cos(Lambda) - sin(chi))*x_dot(2) - sin(chi)*sin(Lambda)*x_dot(4);}
      F(7) = x_dot(7);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Index < 3) {h(0) = x(4) + x(5) - LAMBDA0 - THETA0;}
      if constexpr (Index < 2) {h(1) = sin(x(2))*cos(x(4)) + cos(x(2));}
      return h;
    }

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return -range(x_end(4), x_end(5));}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - V0, x_ini(1) - GAMMA0, x_ini(2) - path_azimuth(LAMBDA0), x_ini(3) - H0, x_ini(4) - LAMBDA0,
           x_ini(5) - THETA0, x_ini(6) /* gauge */, x_end(3) - H_END;
      return b;
    }

    template <typename T>
    PathVector<T> path(StateVector<T> const & x, ControlVector<T> const &, T const &) const
    {return PathVector<T>(pressure(x(0), x(3))/Q_MAX - 1.0);}
  };

} // namespace TestOCP

#endif // TESTS_GERDTS_PROBLEMS_AD_HH
