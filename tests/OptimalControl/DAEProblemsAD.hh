/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_DAE_PROBLEMS_AD_HH
#define TESTS_DAE_PROBLEMS_AD_HH

#include "Sandals/OptimalControlAD.hh"

// Optimal control problems of DAEs, index-reduced or not, written once as templates (the derivatives
// are computed by CppAD). Notes common to all the problems:
//  - the node value of an algebraic state only shifts the stage variables of the Runge-Kutta method,
//    so its initial value is a gauge, fixed by a boundary condition (to its consistent value);
//  - with invariants, the boundary conditions act on the incoming state of the first node, which is
//    then corrected onto the invariants, while at the final node the invariants already hold: the
//    final conditions only fix the components transversal to them (otherwise they would be redundant).
namespace TestOCP
{
  using namespace Sandals;
  using Real = double;

  inline constexpr Real INF{std::numeric_limits<Real>::infinity()};

  /*\
   |  Minimum energy problem (Example 1.0.1 in Gerdts, Optimal Control of ODEs and DAEs, 2023)
   |    min  1/2 int_0^1 u^2 dt
   |    s.t. y1' = u - y4,  y2' = u,  y3' = -y2,  0 = y1 + y3   (index 2, y4 algebraic)
   |         y1(0) = y1(1) = 0,  y2(0) = -y2(1) = 1,  y3(0) = 0
   |  Optimum: u = -2, J* = 2, y1 = t - t^2, y2 = 1 - 2t, y3 = t^2 - t, y4 = 2t - 3.
  \*/

  class MinimumEnergyIndex2 : public OptimalControlProblemAD<MinimumEnergyIndex2, 4, 1, 0, 6>
  {
  public:
    MinimumEnergyIndex2() : OptimalControlProblemAD("MinimumEnergyIndex2", 0.0, 1.0) {}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & y, StateVector<T> const & y_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = y_dot(0) - u(0) + y(3);
      F(1) = y_dot(1) - u(0);
      F(2) = y_dot(2) + y(1);
      F(3) = y(0) + y(2);
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & y_ini, StateVector<T> const & y_end) const
    {
      BoundaryVector<T> b;
      b << y_ini(0), y_ini(1) - 1.0, y_ini(2), y_ini(3) + 3.0 /* gauge */, y_end(0), y_end(1) + 1.0;
      return b;
    }
  };

  // Index 1: the constraint differentiated once, 0 = u - y2 - y4, with the invariant h = y1 + y3.
  class MinimumEnergyIndex1 : public OptimalControlProblemAD<MinimumEnergyIndex1, 4, 1, 1, 6>
  {
  public:
    MinimumEnergyIndex1() : OptimalControlProblemAD("MinimumEnergyIndex1", 0.0, 1.0) {}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & y, StateVector<T> const & y_dot, ControlVector<T> const & u, T const &) const
    {
      StateVector<T> F;
      F(0) = y_dot(0) - u(0) + y(3);
      F(1) = y_dot(1) - u(0);
      F(2) = y_dot(2) + y(1);
      F(3) = u(0) - y(1) - y(3);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & y, T const &) const {return InvariantVector<T>(y(0) + y(2));}

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & y_ini, StateVector<T> const & y_end) const
    {
      BoundaryVector<T> b;
      b << y_ini(0), y_ini(1) - 1.0, y_ini(2), y_ini(3) + 3.0 /* gauge */, y_end(0), y_end(1) + 1.0;
      return b;
    }
  };

  /*\
   |  Pendulum (Example 1.1.9 in Gerdts, Optimal Control of ODEs and DAEs, 2023), control u acting
   |  tangentially on the bob, from horizontal at rest to vertical at rest (bottom: y2(3) <= 0):
   |    min  int_0^3 u^2 dt
   |    s.t. y1' = y3,  y2' = y4,  y3' = -2 y1 y5 + u y2,  y4' = -g - 2 y2 y5 - u y1,
   |         0 = y1^2 + y2^2 - 1                                   (index 3, y5 algebraic)
   |         y1(0) = 1,  y2(0) = y3(0) = y4(0) = 0,  y1(3) = y3(3) = 0
   |  Index reduction (the lost constraints become the invariants):
   |    index 2: 0 = y1 y3 + y2 y4,                       h = (y1^2 + y2^2 - 1)
   |    index 1: 0 = y3^2 + y4^2 + y1 y3' + y2 y4',       h = (y1^2 + y2^2 - 1, y1 y3 + y2 y4)
   |    index 0: y5 = (y3^2 + y4^2 - g y2)/(2 (y1^2 + y2^2)) eliminated, same invariants
  \*/

  template <Integer Index>
  struct GerdtsPendulumSizes
  {
    static constexpr Integer N{Index == 0 ? 4 : 5};                    // states (y5 only if Index > 0)
    static constexpr Integer M{Index == 3 ? 0 : (Index == 2 ? 1 : 2)}; // invariants
    static constexpr Integer NB{N + 3};                                // initial states, y1(3), y3(3), y2(3)
  };

  template <Integer Index>
  class GerdtsPendulum : public OptimalControlProblemAD<GerdtsPendulum<Index>, GerdtsPendulumSizes<Index>::N, 1,
    GerdtsPendulumSizes<Index>::M, GerdtsPendulumSizes<Index>::NB>
  {
    using Sizes = GerdtsPendulumSizes<Index>;
    using Base  = OptimalControlProblemAD<GerdtsPendulum<Index>, Sizes::N, 1, Sizes::M, Sizes::NB>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    using VectorX = typename Base::VectorX;

    Real g{9.81};
    Real y1_ini{1.0}; // initial position (y1_ini^2 + y2_ini^2 != 1 for inconsistent initial data)
    Real y2_ini{0.0};

    GerdtsPendulum(Real const t_end = 3.0) : Base("GerdtsPendulumIndex" + std::to_string(Index), 0.0, t_end) {}

    // The multiplier y5: a state for Index > 0, eliminated for Index = 0
    template <typename T>
    T multiplier(StateVector<T> const & y) const
    {
      if constexpr (Index == 0) {
        return (y(2)*y(2) + y(3)*y(3) - g*y(1)) / (2.0*(y(0)*y(0) + y(1)*y(1)));
      } else {
        return y(4);
      }
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & y, StateVector<T> const & y_dot, ControlVector<T> const & u, T const &) const
    {
      T const lambda{this->multiplier(y)};
      StateVector<T> F;
      F(0) = y_dot(0) - y(2);
      F(1) = y_dot(1) - y(3);
      F(2) = y_dot(2) + 2.0*y(0)*lambda - u(0)*y(1);
      F(3) = y_dot(3) + g + 2.0*y(1)*lambda + u(0)*y(0);
      if constexpr (Index == 3) {F(4) = y(0)*y(0) + y(1)*y(1) - 1.0;}
      if constexpr (Index == 2) {F(4) = y(0)*y(2) + y(1)*y(3);}
      if constexpr (Index == 1) {F(4) = y(2)*y(2) + y(3)*y(3) + y(0)*y_dot(2) + y(1)*y_dot(3);}
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & y, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (Sizes::M > 0) {h(0) = y(0)*y(0) + y(1)*y(1) - 1.0;}
      if constexpr (Sizes::M > 1) {h(1) = y(0)*y(2) + y(1)*y(3);}
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & y_ini, StateVector<T> const & y_end) const
    {
      BoundaryVector<T> b;
      b(0) = y_ini(0) - y1_ini;
      b(1) = y_ini(1) - y2_ini;
      for (Integer i{2}; i < Sizes::N; ++i) {b(i) = y_ini(i);} // at rest, and y5(0) = 0 (gauge)
      b(Sizes::N)     = y_end(0);
      b(Sizes::N + 1) = y_end(2);
      b(Sizes::N + 2) = y_end(1);                               // y2(3) <= 0 (bottom)
      return b;
    }

    VectorX boundary_conditions_lower() const override
    {
      VectorX lower(VectorX::Zero(Sizes::NB));
      lower(Sizes::N + 2) = -INF;
      return lower;
    }
  };

  /*\
   |  Rest-to-rest attitude slew of a rigid body (quaternion q, scalar first, angular velocity w):
   |    min  1/2 int_0^T |tau|^2 dt
   |    s.t. q' = 1/2 q (x) (0, w),  J w' = tau - w x J w,  invariant |q|^2 - 1 = 0
   |  From the identity to the rotation of an angle theta about the principal axis z. The rotation
   |  about z is a stationary point, with theta'' = 6 theta/T^2 (1 - 2t/T) and J* = 6 Jz^2 theta^2/T^3.
  \*/

  class AttitudeSlew : public OptimalControlProblemAD<AttitudeSlew, 7, 3, 1, 14>
  {
  public:
    Eigen::Vector3d J{1.0, 2.0, 3.0}; // principal moments of inertia
    Real theta;                       // rotation angle about z
    Real final_time;                  // duration of the slew

    AttitudeSlew(Real const t_theta = 1.0, Real const t_final_time = 2.0)
      : OptimalControlProblemAD("AttitudeSlew", 0.0, t_final_time), theta(t_theta), final_time(t_final_time) {}

    Real optimal_cost() const {return 6.0*J(2)*J(2)*theta*theta/std::pow(final_time, 3);}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & tau, T const &) const
    {
      Eigen::Vector<T, 4> const q{x.template head<4>()};
      Eigen::Vector<T, 3> const w{x.template tail<3>()};
      Eigen::Vector<T, 3> const Jw{J(0)*w(0), J(1)*w(1), J(2)*w(2)};

      Eigen::Vector<T, 4> q_dot;
      q_dot(0) = -0.5*(q(1)*w(0) + q(2)*w(1) + q(3)*w(2));
      q_dot(1) =  0.5*(q(0)*w(0) + q(2)*w(2) - q(3)*w(1));
      q_dot(2) =  0.5*(q(0)*w(1) + q(3)*w(0) - q(1)*w(2));
      q_dot(3) =  0.5*(q(0)*w(2) + q(1)*w(1) - q(2)*w(0));
      Eigen::Vector<T, 3> const gyroscopic{w.cross(Jw)};

      StateVector<T> F;
      F.template head<4>() = x_dot.template head<4>() - q_dot;
      for (Integer i{0}; i < 3; ++i) {F(4 + i) = J(i)*x_dot(4 + i) - tau(i) + gyroscopic(i);}
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {return InvariantVector<T>(x.template head<4>().squaredNorm() - 1.0);}

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & tau, T const &) const {return 0.5*tau.squaredNorm();}

    // Initial: identity at rest (all the states). Final: vector part of the quaternion and angular
    // velocity (the scalar part follows from the invariant, its sign from q0 >= 0).
    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b(0) = x_ini(0) - 1.0;
      for (Integer i{1}; i < 7; ++i) {b(i) = x_ini(i);}
      b(7)  = x_end(1);
      b(8)  = x_end(2);
      b(9)  = x_end(3) - sin(0.5*theta);
      b(10) = x_end(4);
      b(11) = x_end(5);
      b(12) = x_end(6);
      b(13) = x_end(0); // q0 >= 0
      return b;
    }

    VectorX boundary_conditions_upper() const override
    {
      VectorX upper(VectorX::Zero(14));
      upper(13) = INF;
      return upper;
    }
  };

  /*\
   |  Kepler problem (mu = 1) from the periapsis at distance 1, with eccentricity e, over a whole number
   |  of periods, with the energy (non-quadratic) and the angular momentum as invariants. The control
   |  does not act on the dynamics (the optimal control is zero): the problem tests the drift of the
   |  invariants over a long horizon. After every period the exact solution is back at the periapsis.
   |  On a circular orbit (e = 0) the two invariants have parallel gradients (the circular orbit has the
   |  minimum energy for its angular momentum): they are dependent there.
  \*/

  template <Integer M>
  class KeplerOrbit : public OptimalControlProblemAD<KeplerOrbit<M>, 4, 1, M, 4>
  {
    using Base = OptimalControlProblemAD<KeplerOrbit<M>, 4, 1, M, 4>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;

    Real e; // eccentricity

    KeplerOrbit(Real const t_e, Integer const periods) : Base("KeplerOrbit", 0.0, periods*period(t_e)), e(t_e) {}

    static Real period(Real const e) {return 2.0*M_PI*std::pow(1.0/(1.0 - e), 1.5);}

    Eigen::Vector4d periapsis() const {return Eigen::Vector4d(1.0, 0.0, 0.0, std::sqrt(1.0 + e));}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const &, T const &) const
    {
      T const r3{pow(x(0)*x(0) + x(1)*x(1), 1.5)};
      StateVector<T> F;
      F(0) = x_dot(0) - x(2);
      F(1) = x_dot(1) - x(3);
      F(2) = x_dot(2) + x(0)/r3;
      F(3) = x_dot(3) + x(1)/r3;
      return F;
    }

    // Energy and angular momentum, minus their values on the orbit
    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> h;
      if constexpr (M > 0) {h(0) = 0.5*(x(2)*x(2) + x(3)*x(3)) - 1.0/sqrt(x(0)*x(0) + x(1)*x(1)) - 0.5*(e - 1.0);}
      if constexpr (M > 1) {h(1) = x(0)*x(3) - x(1)*x(2) - sqrt(1.0 + e);}
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - 1.0, x_ini(1), x_ini(2), x_ini(3) - sqrt(1.0 + e);
      return b;
    }
  };

  /*\
   |  Double pendulum driven by a horizontal force u on the first mass, from hanging at rest to both
   |  links rotated by the angle a at rest:  min 1/2 int_0^T u^2 dt.
   |  - Cartesian coordinates, index 1 (acceleration-level constraints), with the position and velocity
   |    constraints as invariants. States (x1, y1, x2, y2, vx1, vy1, vx2, vy2, l1, l2).
   |  - Minimal coordinates (angles from the downward vertical), as reference.
  \*/

  struct DoublePendulumParameters
  {
    Real m1{1.0}, m2{1.0}, L1{1.0}, L2{1.0}, g{9.81};
    Real a{0.3}; // final angle of both links
    Real T{2.0}; // final time
  };

  class DoublePendulumCartesian
    : public OptimalControlProblemAD<DoublePendulumCartesian, 10, 1, 4, 14>, public DoublePendulumParameters
  {
  public:
    enum : Integer {X1 = 0, Y1, X2, Y2, VX1, VY1, VX2, VY2, L1_, L2_};

    DoublePendulumCartesian() : OptimalControlProblemAD("DoublePendulumCartesian", 0.0, DoublePendulumParameters().T) {}

    // Position constraints g = (|p1|^2 - L1^2, |p2 - p1|^2 - L2^2), and their time derivatives
    template <typename T>
    StateVector<T> dae(StateVector<T> const & s, StateVector<T> const & s_dot, ControlVector<T> const & u, T const &) const
    {
      T const dx{s(X2) - s(X1)}, dy{s(Y2) - s(Y1)};
      T const dvx{s(VX2) - s(VX1)}, dvy{s(VY2) - s(VY1)};
      T const dax{s_dot(VX2) - s_dot(VX1)}, day{s_dot(VY2) - s_dot(VY1)};

      StateVector<T> F;
      for (Integer i{0}; i < 4; ++i) {F(i) = s_dot(i) - s(VX1 + i);}
      // m p'' = -G^T lambda + forces, with G the Jacobian of the position constraints
      F(4) = m1*s_dot(VX1) + 2.0*s(X1)*s(L1_) - 2.0*dx*s(L2_) - u(0);
      F(5) = m1*s_dot(VY1) + 2.0*s(Y1)*s(L1_) - 2.0*dy*s(L2_) + m1*g;
      F(6) = m2*s_dot(VX2) + 2.0*dx*s(L2_);
      F(7) = m2*s_dot(VY2) + 2.0*dy*s(L2_) + m2*g;
      // Acceleration-level constraints (index 1)
      F(8) = s(VX1)*s(VX1) + s(VY1)*s(VY1) + s(X1)*s_dot(VX1) + s(Y1)*s_dot(VY1);
      F(9) = dvx*dvx + dvy*dvy + dx*dax + dy*day;
      return F;
    }

    // Position and velocity constraints
    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & s, T const &) const
    {
      T const dx{s(X2) - s(X1)}, dy{s(Y2) - s(Y1)};
      InvariantVector<T> h;
      h(0) = s(X1)*s(X1) + s(Y1)*s(Y1) - L1*L1;
      h(1) = dx*dx + dy*dy - L2*L2;
      h(2) = s(X1)*s(VX1) + s(Y1)*s(VY1);
      h(3) = dx*(s(VX2) - s(VX1)) + dy*(s(VY2) - s(VY1));
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    // Initial: hanging at rest, with the static multipliers (gauges). Final: horizontal positions and
    // velocities (transversal to the invariants).
    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & s_ini, StateVector<T> const & s_end) const
    {
      BoundaryVector<T> b;
      b(0) = s_ini(X1);
      b(1) = s_ini(Y1) + L1;
      b(2) = s_ini(X2);
      b(3) = s_ini(Y2) + L1 + L2;
      for (Integer i{0}; i < 4; ++i) {b(4 + i) = s_ini(VX1 + i);}
      b(8) = s_ini(L1_) - (m1 + m2)*g/(2.0*L1);
      b(9) = s_ini(L2_) - m2*g/(2.0*L2);
      b(10) = s_end(X1) - L1*sin(a);
      b(11) = s_end(X2) - (L1 + L2)*sin(a);
      b(12) = s_end(VX1);
      b(13) = s_end(VX2);
      return b;
    }
  };

  class DoublePendulumAngles
    : public OptimalControlProblemAD<DoublePendulumAngles, 4, 1, 0, 8>, public DoublePendulumParameters
  {
  public:
    DoublePendulumAngles() : OptimalControlProblemAD("DoublePendulumAngles", 0.0, DoublePendulumParameters().T) {}

    // States (theta1, theta2, omega1, omega2): mass matrix times accelerations equals the forces
    template <typename T>
    StateVector<T> dae(StateVector<T> const & q, StateVector<T> const & q_dot, ControlVector<T> const & u, T const &) const
    {
      T const c{cos(q(0) - q(1))}, s{sin(q(0) - q(1))};
      StateVector<T> F;
      F(0) = q_dot(0) - q(2);
      F(1) = q_dot(1) - q(3);
      F(2) = (m1 + m2)*L1*L1*q_dot(2) + m2*L1*L2*c*q_dot(3)
           + m2*L1*L2*s*q(3)*q(3) + (m1 + m2)*g*L1*sin(q(0)) - u(0)*L1*cos(q(0));
      F(3) = m2*L1*L2*c*q_dot(2) + m2*L2*L2*q_dot(3)
           - m2*L1*L2*s*q(2)*q(2) + m2*g*L2*sin(q(1));
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & q_ini, StateVector<T> const & q_end) const
    {
      BoundaryVector<T> b;
      b << q_ini(0), q_ini(1), q_ini(2), q_ini(3), q_end(0) - a, q_end(1) - a, q_end(2), q_end(3);
      return b;
    }
  };

  /*\
   |  Stiff tracking problem: the state follows the control with the fast rate 50,
   |    min  int_0^1 ((x - 1)^2 + 0.1 u^2) dt,  s.t. x' = -50 (x - u),  x(0) = 0.
   |  Explicit methods need small substeps (|h lambda| < 2.8 for RK4).
  \*/

  class StiffTracking : public OptimalControlProblemAD<StiffTracking, 1, 1, 0, 1>
  {
  public:
    StiffTracking() : OptimalControlProblemAD("StiffTracking", 0.0, 1.0) {}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {return StateVector<T>(x_dot(0) + 50.0*(x(0) - u(0)));}

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const &) const
    {return (x(0) - 1.0)*(x(0) - 1.0) + 0.1*u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {return BoundaryVector<T>(x_ini(0));}
  };

} // namespace TestOCP

#endif // TESTS_DAE_PROBLEMS_AD_HH
