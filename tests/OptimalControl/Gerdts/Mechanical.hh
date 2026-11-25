/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef TESTS_GERDTS_MECHANICAL_HH
#define TESTS_GERDTS_MECHANICAL_HH

#include "DAEProblemsAD.hh"

// Mechanical examples: the pendulum in GGL form and the ski driver (index-3 DAEs of motion), from
// M. Gerdts, Optimal Control of ODEs and DAEs, 2nd edition, De Gruyter, 2023 (the example numbers
// refer to the book). Conventions as in DAEProblemsAD.hh: the node values of the algebraic states
// at t = 0 are gauges, fixed by a boundary condition.
namespace TestOCP
{

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
   |  Semi-active chassis control of a half-car (Example 5.7.6): the chassis (height zc, pitch phi) on
   |  the front and rear wheels (heights zf, zr), driving at v = 20 m/s over a pothole of depth 4 cm,
   |  with the damping of the front and rear suspensions as controls uf, ur in [0, 1]. States
   |  (zc, zf, zr, phi, vc, vf, vr, omega, lambda_f, lambda_r), with the wheel contact forces lambda:
   |    mc vc' = -mc g + Ff + Fr,     mf vf' = -mf g - Ff - lambda_f,   mr vr' = -mr g - Fr - lambda_r,
   |    Jc omega' = cos(phi) (lr Fr - lf Ff),
   |    0 = zf - xi(t),  0 = zr - xi(t - l cos(phi)/v)                  (index 3)
   |  with the suspension forces
   |    Ff = cf0 - cf (zc - lf sin(phi) - zf) - uf df (vc - lf omega cos(phi) - vf),
   |    Fr = cr0 - cr (zc + lr sin(phi) - zr) - ur dr (vc + lr omega cos(phi) - vr),
   |  the road xi(s) = 0.02 (cos(2 pi (s - 0.3)/0.2) - 1) on [0.3, 0.5] (zero elsewhere, once
   |  continuously differentiable, so the contact forces jump), and the rear wheel delayed by
   |  l cos(phi)/v. Minimize the comfort, handling and regularization criteria
   |    a1 int (zc''^2 + phi''^2) + a2 int (lambda_f^2 + lambda_r^2) + a3 int ((uf - 1/2)^2 + (ur - 1/2)^2)
   |  on [0, 3], from rest at the equilibrium. The spring offsets cf0, cr0 are not given in the book:
   |  they carry the static load (cf0 = lr mc g/l, cr0 = lf mc g/l), so the equilibrium is zc = phi = 0.
   |  Index reduction (the lost constraints become the invariants, which depend on time):
   |    index 2: 0 = vf - xi'(t),   0 = vr - d/dt xi(t - l cos(phi)/v)
   |    index 1: the accelerations,  with the positions and velocities constraints as invariants.
   |  The book shows the solution in figures only (BDF with GGL stabilization, 101 grid points).
  \*/

  template <Integer Index>
  class HalfCar : public OptimalControlProblemAD<HalfCar<Index>, 10, 2, 2*(3 - Index), 10>
  {
    using Base = OptimalControlProblemAD<HalfCar<Index>, 10, 2, 2*(3 - Index), 10>;

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;
    using VectorF = typename Base::VectorF;
    using VectorU = typename Base::VectorU;

    // Weights and parameters of the book
    static constexpr Real A1{1.0}, A2{1.0e-9}, A3{1.0e-2};
    static constexpr Real LF{1.1}, LR{1.4}, L{LF + LR}, MC{1300.0}, MF{40.0}, MR{40.0}, JC{1659.97};
    static constexpr Real CF{25000.0}, CR{25000.0}, DF{2500.0}, DR{2500.0}, V{20.0}, G{9.81};
    static constexpr Real CF0{LR*MC*G/L}, CR0{LF*MC*G/L}; // spring offsets carrying the static load

    HalfCar() : Base("HalfCarIndex" + std::to_string(Index), 0.0, 3.0)
    {
      this->u_bounds(VectorU(0.0, 0.0), VectorU(1.0, 1.0));
      VectorF nominal;
      nominal << 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 1.0e4, 1.0e4; // m, rad, m/s, rad/s, N
      this->nominal(nominal, VectorU(1.0, 1.0));
    }

    // Equilibrium at rest on the flat road
    static VectorF initial_states()
    {
      VectorF x(VectorF::Zero());
      x(8) = -(MF*G + CF0);
      x(9) = -(MR*G + CR0);
      return x;
    }

    // The road profile under the front wheel and its first two derivatives at the position s (in time)
    template <typename T>
    static Eigen::Vector<T, 3> road(T const & s)
    {
      Real const depth{0.02}, w{2.0*M_PI/0.2};
      T const a{w*(s - 0.3)};
      Eigen::Vector<T, 3> bump(depth*(cos(a) - 1.0), -depth*w*sin(a), -depth*w*w*cos(a));
      Eigen::Vector<T, 3> xi;
      for (Integer i{0}; i < 3; ++i) {
        xi(i) = CppAD::CondExpLt(s, T(0.3), T(0.0), CppAD::CondExpGt(s, T(0.5), T(0.0), bump(i)));
      }
      return xi;
    }

    // The suspension forces (front, rear)
    template <typename T>
    static Eigen::Vector<T, 2> forces(StateVector<T> const & x, ControlVector<T> const & u)
    {
      T const zc{x(0)}, phi{x(3)}, vc{x(4)}, omega{x(7)};
      return Eigen::Vector<T, 2>(
        CF0 - CF*(zc - LF*sin(phi) - x(1)) - u(0)*DF*(vc - LF*omega*cos(phi) - x(5)),
        CR0 - CR*(zc + LR*sin(phi) - x(2)) - u(1)*DR*(vc + LR*omega*cos(phi) - x(6)));
    }

    // The chassis accelerations (zc'', phi'')
    template <typename T>
    static Eigen::Vector<T, 2> chassis_accelerations(StateVector<T> const & x, ControlVector<T> const & u)
    {
      Eigen::Vector<T, 2> const F{forces(x, u)};
      return Eigen::Vector<T, 2>((F(0) + F(1))/MC - G, cos(x(3))*(LR*F(1) - LF*F(0))/JC);
    }

    // The constraints of the wheels on the road: positions, velocities (and accelerations, given phi'')
    template <typename T>
    static Eigen::Vector<T, 6> contact(StateVector<T> const & x, T const & t, T const & phi_dd)
    {
      T const phi{x(3)}, omega{x(7)};
      T const s{t - L*cos(phi)/V}, s_dot{1.0 + L*sin(phi)*omega/V};
      Eigen::Vector<T, 3> const front{road(t)}, rear{road(s)};
      Eigen::Vector<T, 6> c;
      c << x(1) - front(0), x(2) - rear(0),
           x(5) - front(1), x(6) - rear(1)*s_dot,
           -front(2), -(rear(2)*s_dot*s_dot + rear(1)*L/V*(cos(phi)*omega*omega + sin(phi)*phi_dd));
      return c;
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const & t) const
    {
      Eigen::Vector<T, 2> const F{forces(x, u)};
      StateVector<T> R;
      R(0) = x_dot(0) - x(4);
      R(1) = x_dot(1) - x(5);
      R(2) = x_dot(2) - x(6);
      R(3) = x_dot(3) - x(7);
      R(4) = MC*x_dot(4) + MC*G - F(0) - F(1);
      R(5) = MF*x_dot(5) + MF*G + F(0) + x(8);
      R(6) = MR*x_dot(6) + MR*G + F(1) + x(9);
      R(7) = JC*x_dot(7) - cos(x(3))*(LR*F(1) - LF*F(0));
      Eigen::Vector<T, 6> const c{contact(x, t, x_dot(7))};
      if constexpr (Index == 3) {R(8) = c(0); R(9) = c(1);}
      if constexpr (Index == 2) {R(8) = c(2); R(9) = c(3);}
      if constexpr (Index == 1) {R(8) = x_dot(5) + c(4); R(9) = x_dot(6) + c(5);}
      return R;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const & t) const
    {
      InvariantVector<T> h;
      if constexpr (Index < 3) {
        Eigen::Vector<T, 6> const c{contact(x, t, T(0.0))};
        h.template head<2>() = c.template head<2>();
        if constexpr (Index < 2) {h.template tail<2>() = c.template segment<2>(2);}
      }
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const &) const
    {
      Eigen::Vector<T, 2> const a{chassis_accelerations(x, u)};
      return A1*a.squaredNorm() + A2*(x(8)*x(8) + x(9)*x(9))
           + A3*((u(0) - 0.5)*(u(0) - 0.5) + (u(1) - 0.5)*(u(1) - 0.5));
    }

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {
      VectorF const start{initial_states()};
      BoundaryVector<T> b;
      for (Integer i{0}; i < 10; ++i) {b(i) = x_ini(i) - start(i);} // the contact forces are gauges
      return b;
    }
  };

} // namespace TestOCP

#endif // TESTS_GERDTS_MECHANICAL_HH
