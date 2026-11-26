/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef TESTS_GERDTS_FLIGHT_HH
#define TESTS_GERDTS_FLIGHT_HH

#include "DAEProblemsAD.hh"

// Flight example: emergency landing of a hypersonic flight system (ODE, and index-3 DAE on a
// prescribed path), from M. Gerdts, Optimal Control of ODEs and DAEs, 2nd edition, De Gruyter, 2023
// (the example numbers refer to the book). Conventions as in DAEProblemsAD.hh: the node values of
// the algebraic states at t = 0 are gauges, fixed by a boundary condition.
namespace TestOCP
{

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
    bool pressure_limit{true}; // the dynamic pressure constraint (Figure 1.5 compares without it)

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

    VectorX path_constraints_upper() const override {return VectorX::Constant(1, pressure_limit ? 0.0 : INF);}
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

#endif // TESTS_GERDTS_FLIGHT_HH
