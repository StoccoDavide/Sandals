/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_HARD_PROBLEMS_AD_HH
#define TESTS_HARD_PROBLEMS_AD_HH

#include "Sandals/OptimalControlAD.hh"

// Harder optimal control problems, with known optima or independent references.
//
// Free final time: the duration T is a constant state (T' = 0), and the dynamics is written in the
// scaled time tau in [0, 1], where d/dtau = T d/dt.
namespace TestOCP
{
  using namespace Sandals;
  using Real = double;

  /*\
   |  Goddard rocket (COPS benchmark, normalized units): maximize the final altitude with free final
   |  time, thrust 0 <= u <= 3.5, and the final mass 0.6:
   |    h' = v,  v' = (u - D(h, v))/m - 1/h^2,  m' = -u/c,
   |    D = 310 v^2 exp(-500 (h - 1)),  c = 0.5,  h(0) = 1, v(0) = 0, m(0) = 1, m(T) = 0.6.
   |  The optimal thrust is bang-singular-bang (the singular arc balances thrust and drag). Optimum:
   |  h(T) = 1.01283 (Dolan, More, Munson, COPS 3.0).
   |  States (h, v, m, T).
  \*/

  class GoddardRocket : public OptimalControlProblemAD<GoddardRocket, 4, 1, 0, 4>
  {
  public:
    static constexpr Real U_MAX{3.5}, C{0.5}, DRAG{310.0}, SCALE_HEIGHT{500.0}, FINAL_MASS{0.6};
    static constexpr Real OPTIMAL_ALTITUDE{1.01283};

    GoddardRocket() : OptimalControlProblemAD("GoddardRocket", 0.0, 1.0)
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      this->x_bounds(VectorF(1.0, 0.0, FINAL_MASS, 0.01), VectorF(inf, inf, 1.0, inf));
      this->u_bounds(VectorU(0.0), VectorU(U_MAX));
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      T const h{x(0)}, v{x(1)}, m{x(2)}, duration{x(3)};
      T const drag{DRAG*v*v*exp(-SCALE_HEIGHT*(h - 1.0))};
      StateVector<T> F;
      F(0) = x_dot(0) - duration*v;
      F(1) = x_dot(1) - duration*((u(0) - drag)/m - 1.0/(h*h));
      F(2) = x_dot(2) + duration*u(0)/C;
      F(3) = x_dot(3);
      return F;
    }

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return -x_end(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - 1.0, x_ini(1), x_ini(2) - 1.0, x_end(2) - FINAL_MASS;
      return b;
    }
  };

  /*\
   |  Low-thrust orbit raising (Bryson and Ho, Applied Optimal Control, 1975): maximize the radius of
   |  the final circular orbit reached in the time 3.32, with the constant thrust 0.1405, the mass
   |  flow 0.0749, and the thrust direction as control (mu = 1, initial circular orbit of radius 1).
   |  Optimum: r(T) = 1.525 (reported to three digits; the Richardson extrapolation of the two
   |  formulations below gives 1.52528). Two formulations of the same problem:
   |  - polar coordinates (r, theta, v_r, v_theta), thrust angle phi from the local horizontal;
   |  - Cartesian coordinates (x, y, vx, vy), thrust angle psi from the x axis.
  \*/

  struct OrbitRaisingData
  {
    static constexpr Real THRUST{0.1405}, MASS_FLOW{0.0749}, DURATION{3.32}, OPTIMAL_RADIUS{1.525};
    template <typename T> static T acceleration(T const & t) {return THRUST/(1.0 - MASS_FLOW*t);}
  };

  class OrbitRaisingPolar : public OptimalControlProblemAD<OrbitRaisingPolar, 4, 1, 0, 6>, public OrbitRaisingData
  {
  public:
    OrbitRaisingPolar() : OptimalControlProblemAD("OrbitRaisingPolar", 0.0, DURATION) {}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const & t) const
    {
      T const r{x(0)}, v_r{x(2)}, v_t{x(3)}, a{acceleration(t)};
      StateVector<T> F;
      F(0) = x_dot(0) - v_r;
      F(1) = x_dot(1) - v_t/r;
      F(2) = x_dot(2) - (v_t*v_t/r - 1.0/(r*r) + a*sin(u(0)));
      F(3) = x_dot(3) - (-v_r*v_t/r + a*cos(u(0)));
      return F;
    }

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return -x_end(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b << x_ini(0) - 1.0, x_ini(1), x_ini(2), x_ini(3) - 1.0, x_end(2), x_end(3) - 1.0/sqrt(x_end(0));
      return b;
    }
  };

  class OrbitRaisingCartesian : public OptimalControlProblemAD<OrbitRaisingCartesian, 4, 1, 0, 6>, public OrbitRaisingData
  {
  public:
    OrbitRaisingCartesian() : OptimalControlProblemAD("OrbitRaisingCartesian", 0.0, DURATION) {}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const & t) const
    {
      T const r3{pow(x(0)*x(0) + x(1)*x(1), 1.5)}, a{acceleration(t)};
      StateVector<T> F;
      F(0) = x_dot(0) - x(2);
      F(1) = x_dot(1) - x(3);
      F(2) = x_dot(2) - (-x(0)/r3 + a*cos(u(0)));
      F(3) = x_dot(3) - (-x(1)/r3 + a*sin(u(0)));
      return F;
    }

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const
    {return -sqrt(x_end(0)*x_end(0) + x_end(1)*x_end(1));}

    // Final circular orbit: zero radial velocity and |v|^2 = 1/r
    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      T const r{sqrt(x_end(0)*x_end(0) + x_end(1)*x_end(1))};
      BoundaryVector<T> b;
      b << x_ini(0) - 1.0, x_ini(1), x_ini(2), x_ini(3) - 1.0,
           x_end(0)*x_end(2) + x_end(1)*x_end(3), x_end(2)*x_end(2) + x_end(3)*x_end(3) - 1.0/r;
      return b;
    }
  };

  /*\
   |  Cart-pole swing-up: a cart of mass M on a track |p| <= p_max, driven by the force |u| <= u_max,
   |  carries a pendulum of mass m and length l, from hanging at rest to upright at rest with the cart
   |  back at the origin, minimizing 1/2 int u^2 dt.
   |  - Minimal coordinates (p, theta, p', theta'), theta from the downward vertical:
   |      (M + m) p'' + m l (theta'' cos(theta) - theta'^2 sin(theta)) = u
   |      m l p'' cos(theta) + m l^2 theta'' + m g l sin(theta) = 0
   |  - Cartesian coordinates, index 1 with the position and velocity constraints as invariants: states
   |    (p, x, y, p', vx, vy, lambda), with the rod constraint (x - p)^2 + y^2 = l^2.
  \*/

  struct CartPoleParameters
  {
    static constexpr Real CART{1.0}, POLE{0.3}, LENGTH{0.5}, GRAVITY{9.81};
    static constexpr Real DURATION{2.0}, TRACK{0.4}, FORCE{20.0};
  };

  class CartPoleAngles : public OptimalControlProblemAD<CartPoleAngles, 4, 1, 0, 8>, public CartPoleParameters
  {
  public:
    CartPoleAngles() : OptimalControlProblemAD("CartPoleAngles", 0.0, DURATION)
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      this->x_bounds(VectorF(-TRACK, -inf, -inf, -inf), VectorF(TRACK, inf, inf, inf));
      this->u_bounds(VectorU(-FORCE), VectorU(FORCE));
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & q, StateVector<T> const & q_dot, ControlVector<T> const & u, T const &) const
    {
      T const c{cos(q(1))}, s{sin(q(1))};
      StateVector<T> F;
      F(0) = q_dot(0) - q(2);
      F(1) = q_dot(1) - q(3);
      F(2) = (CART + POLE)*q_dot(2) + POLE*LENGTH*(q_dot(3)*c - q(3)*q(3)*s) - u(0);
      F(3) = POLE*LENGTH*q_dot(2)*c + POLE*LENGTH*LENGTH*q_dot(3) + POLE*GRAVITY*LENGTH*s;
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & q_ini, StateVector<T> const & q_end) const
    {
      BoundaryVector<T> b;
      b << q_ini(0), q_ini(1), q_ini(2), q_ini(3), q_end(0), q_end(1) - M_PI, q_end(2), q_end(3);
      return b;
    }
  };

  class CartPoleCartesian : public OptimalControlProblemAD<CartPoleCartesian, 7, 1, 2, 12>, public CartPoleParameters
  {
  public:
    enum : Integer {P = 0, X, Y, VP, VX, VY, LAMBDA};

    CartPoleCartesian() : OptimalControlProblemAD("CartPoleCartesian", 0.0, DURATION)
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      VectorF lower(VectorF::Constant(-inf)), upper(VectorF::Constant(inf));
      lower(P) = -TRACK;
      upper(P) = TRACK;
      this->x_bounds(lower, upper);
      this->u_bounds(VectorU(-FORCE), VectorU(FORCE));
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & s, StateVector<T> const & s_dot, ControlVector<T> const & u, T const &) const
    {
      T const dx{s(X) - s(P)}, dvx{s(VX) - s(VP)}, dax{s_dot(VX) - s_dot(VP)};
      StateVector<T> F;
      F(0) = s_dot(P) - s(VP);
      F(1) = s_dot(X) - s(VX);
      F(2) = s_dot(Y) - s(VY);
      F(3) = CART*s_dot(VP) - 2.0*dx*s(LAMBDA) - u(0);
      F(4) = POLE*s_dot(VX) + 2.0*dx*s(LAMBDA);
      F(5) = POLE*s_dot(VY) + 2.0*s(Y)*s(LAMBDA) + POLE*GRAVITY;
      F(6) = dvx*dvx + s(VY)*s(VY) + dx*dax + s(Y)*s_dot(VY); // rod constraint, index 1
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & s, T const &) const
    {
      T const dx{s(X) - s(P)};
      InvariantVector<T> h;
      h(0) = dx*dx + s(Y)*s(Y) - LENGTH*LENGTH;
      h(1) = dx*(s(VX) - s(VP)) + s(Y)*s(VY);
      return h;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 0.5*u(0)*u(0);}

    // Initial: hanging at rest with the static multiplier (gauge). Final: cart at the origin at rest,
    // pole above the cart at rest (the transversal components), upright (y >= 0).
    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & s_ini, StateVector<T> const & s_end) const
    {
      BoundaryVector<T> b;
      b << s_ini(P), s_ini(X), s_ini(Y) + LENGTH, s_ini(VP), s_ini(VX), s_ini(VY), s_ini(LAMBDA) - POLE*GRAVITY/(2.0*LENGTH),
           s_end(P), s_end(VP), s_end(X) - s_end(P), s_end(VX) - s_end(VP), s_end(Y);
      return b;
    }

    VectorX boundary_conditions_upper() const override
    {
      VectorX upper(VectorX::Zero(12));
      upper(11) = std::numeric_limits<Real>::infinity();
      return upper;
    }
  };

  /*\
   |  Hypersensitive problem (Rao and Mease, 2000):
   |    min 1/2 int_0^T (x^2 + u^2) dt,  x' = -x^3 + u,  x(0) = 1.5,  x(T) = 1.
   |  For a long horizon the solution has two boundary layers and x ~ 0 in between, and the optimal cost
   |  tends to V(1.5) + W(1), the costs of the two layers, from the stationary Hamilton-Jacobi equation:
   |    V(a) = (a^2 sqrt(a^4 + 1) + asinh(a^2))/4 - a^4/4,  W(b) = (b^2 sqrt(b^4 + 1) + asinh(b^2))/4 + b^4/4.
  \*/

  class Hypersensitive : public OptimalControlProblemAD<Hypersensitive, 1, 1, 0, 2>
  {
  public:
    Hypersensitive(Real const T) : OptimalControlProblemAD("Hypersensitive", 0.0, T) {}

    static Real layer(Real const a) {return (a*a*std::sqrt(a*a*a*a + 1.0) + std::asinh(a*a))/4.0;}
    static Real asymptotic_cost() {return layer(1.5) - std::pow(1.5, 4)/4.0 + layer(1.0) + 1.0/4.0;}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {return StateVector<T>(x_dot(0) + x(0)*x(0)*x(0) - u(0));}

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const &) const {return 0.5*(x(0)*x(0) + u(0)*u(0));}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {return BoundaryVector<T>(x_ini(0) - 1.5, x_end(0) - 1.0);}
  };

  /*\
   |  Brachistochrone with the energy invariant and free final time: a bead slides without friction
   |  from rest at the origin to (1, -1) in minimum time, the control being the direction of motion.
   |  States (x, y, v, T), with the energy conservation v^2 + 2 g y = 0 as invariant:
   |    x' = T v sin(theta),  y' = -T v cos(theta),  v' = T g cos(theta).
   |  Optimum: the cycloid x = a (phi - sin(phi)), y = -a (1 - cos(phi)), with phi_f solving
   |  phi - sin(phi) = 1 - cos(phi), a = 1/(1 - cos(phi_f)), and the time phi_f sqrt(a/g).
  \*/

  class Brachistochrone : public OptimalControlProblemAD<Brachistochrone, 4, 1, 1, 5>
  {
  public:
    static constexpr Real G{9.81};

    Brachistochrone() : OptimalControlProblemAD("Brachistochrone", 0.0, 1.0)
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      this->x_bounds(VectorF(-inf, -inf, -inf, 0.01), VectorF(inf, inf, inf, inf));
    }

    static Real optimal_time()
    {
      Real phi{2.4}; // Newton on phi - sin(phi) - (1 - cos(phi)) = 0
      for (Integer i{0}; i < 50; ++i) {phi -= (phi - std::sin(phi) - 1.0 + std::cos(phi)) / (1.0 - std::cos(phi) - std::sin(phi));}
      Real const a{1.0/(1.0 - std::cos(phi))};
      return phi*std::sqrt(a/G);
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      T const v{x(2)}, duration{x(3)};
      StateVector<T> F;
      F(0) = x_dot(0) - duration*v*sin(u(0));
      F(1) = x_dot(1) + duration*v*cos(u(0));
      F(2) = x_dot(2) - duration*G*cos(u(0));
      F(3) = x_dot(3);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const {return InvariantVector<T>(x(2)*x(2) + 2.0*G*x(1));}

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return x_end(3);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
    {
      BoundaryVector<T> b;
      b << x_ini(0), x_ini(1), x_ini(2), x_end(0) - 1.0, x_end(1) + 1.0;
      return b;
    }
  };

} // namespace TestOCP

#endif // TESTS_HARD_PROBLEMS_AD_HH
