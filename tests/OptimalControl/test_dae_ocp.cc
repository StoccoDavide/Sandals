/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimal control of DAEs (index-reduced or not) by direct multiple shooting with IPOPT, validated
// against analytical optima, equivalent formulations, exact solutions, and the collocation.

#include <gtest/gtest.h>

#include "Sandals/Ipopt.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/GaussLegendre4.hh"
#include "Sandals/RungeKutta/RK4.hh"

#include "DAEProblemsAD.hh"
#include "PendulumProblems.hh"
#include "OptimalControlProblems.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using Vector1 = Eigen::Vector<Real, 1>;

// Quiet IPOPT with the default tolerance (1e-12); SANDALS_IPOPT_DEBUG=1 prints the iterations.
IpoptOptions options()
{
  IpoptOptions opts;
  opts.integers["print_level"] = std::getenv("SANDALS_IPOPT_DEBUG") ? 5 : 0;
  opts.strings["sb"] = "yes";
  opts.integers["max_iter"] = 500;
  return opts;
}

// A problem transcribed by multiple shooting and solved by IPOPT, starting from a forward simulation
// with zero controls (or from a states guess function), optionally with a continuation on the control
// bounds.
template <typename TableauType, Integer N, Integer U, Integer M, typename Problem>
struct ShootingSolve
{
  static constexpr Integer S{TableauType::Vector::RowsAtCompileTime};

  std::shared_ptr<Problem> problem;
  MultipleShooting<Real, S, N, U, M> nlp;
  bool succeeded{true};
  Real J{0.0};

  using FunctionX = std::function<Eigen::Vector<Real, N>(Real)>;

  ShootingSolve(std::shared_ptr<Problem> p, Integer const intervals, Integer const substeps,
    Eigen::Vector<Real, N> const & x_ini, Real const epsilon = 0.0, std::vector<Real> const & u_max_sequence = {},
    FunctionX const & x_guess = nullptr)
    : problem(p), nlp(p, TableauType(), intervals, substeps, epsilon)
  {
    auto const zero_controls = [] (Real) {return Eigen::Vector<Real, U>::Zero().eval();};
    if (x_guess) {
      nlp.guess(x_guess, zero_controls);
    } else if (!nlp.simulate_guess(x_ini, zero_controls)) {
      succeeded = false;
      return;
    }
    if (u_max_sequence.empty()) {
      succeeded = ipopt_solve(nlp, options(), J) == Ipopt::Solve_Succeeded;
    }
    for (Real const u_max : u_max_sequence) {
      p->u_bounds(Eigen::Vector<Real, U>::Constant(-u_max), Eigen::Vector<Real, U>::Constant(u_max));
      succeeded = succeeded && ipopt_solve(nlp, options(), J) == Ipopt::Solve_Succeeded;
    }
  }

  // Largest violation of the invariants on the nodes
  Real invariants_violation() const
  {
    Real violation{0.0};
    for (Integer k{0}; k <= nlp.intervals(); ++k) {
      violation = std::max(violation, problem->system()->h(nlp.x(nlp.guess(), k), nlp.t_mesh()(k)).template lpNorm<Eigen::Infinity>());
    }
    return violation;
  }

  // Largest projection multiplier (the size of the least-squares correction)
  Real largest_correction() const
  {
    Real largest{0.0};
    for (Integer k{0}; k <= nlp.intervals(); ++k) {
      if constexpr (M > 0) {largest = std::max(largest, nlp.mu(nlp.guess(), k).template lpNorm<Eigen::Infinity>());}
    }
    return largest;
  }
};

template <typename Problem, Integer N, Integer U, Integer M> using Radau = ShootingSolve<RadauIIA5Tableau<Real>, N, U, M, Problem>;

Real const UNBOUNDED{std::numeric_limits<Real>::infinity()};

/*\
 |  Minimum energy problem (Gerdts 1.0.1): index 2, and index 1 with the invariant
\*/

TEST(DAEOptimalControl, MinimumEnergy) {
  Eigen::Vector4d const y_ini(0.0, 1.0, 0.0, -3.0);

  Radau<MinimumEnergyIndex2, 4, 1, 0> index2(std::make_shared<MinimumEnergyIndex2>(), 20, 1, y_ini);
  ASSERT_TRUE(index2.succeeded);
  EXPECT_NEAR(index2.J, 2.0, 1.0e-10);
  EXPECT_LT((index2.nlp.u_mesh(index2.nlp.guess()).array() + 2.0).abs().maxCoeff(), 1.0e-9);

  Radau<MinimumEnergyIndex1, 4, 1, 1> index1(std::make_shared<MinimumEnergyIndex1>(), 20, 1, y_ini);
  ASSERT_TRUE(index1.succeeded);
  EXPECT_NEAR(index1.J, 2.0, 1.0e-10);
  EXPECT_LT(index1.invariants_violation(), 1.0e-12);

  // Collocation of the index-2 formulation
  Collocation<Real, 3, 4, 1> collocation(std::make_shared<MinimumEnergyIndex2>(), RadauIIA5Tableau<Real>(), 20);
  Real J_collocation;
  ASSERT_EQ(ipopt_solve(collocation, options(), J_collocation), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_collocation, 2.0, 1.0e-10);
}

/*\
 |  Pendulum (Gerdts 1.1.9): index 3, 2, 1, 0 agree with each other and with the collocation
\*/

Eigen::Vector<Real, 5> const GERDTS_START{(Eigen::Vector<Real, 5>() << 1.0, 0.0, 0.0, 0.0, 0.0).finished()};

TEST(DAEOptimalControl, GerdtsPendulumFormulationsAgree) {
  Integer const n{30}, substeps{2};
  Radau<GerdtsPendulum<3>, 5, 1, 0> index3(std::make_shared<GerdtsPendulum<3>>(), n, substeps, GERDTS_START);
  Radau<GerdtsPendulum<2>, 5, 1, 1> index2(std::make_shared<GerdtsPendulum<2>>(), n, substeps, GERDTS_START);
  Radau<GerdtsPendulum<1>, 5, 1, 2> index1(std::make_shared<GerdtsPendulum<1>>(), n, substeps, GERDTS_START);
  Radau<GerdtsPendulum<1>, 5, 1, 2> single(std::make_shared<GerdtsPendulum<1>>(), n, 1, GERDTS_START);
  Radau<GerdtsPendulum<0>, 4, 1, 2> index0(std::make_shared<GerdtsPendulum<0>>(), n, substeps, GERDTS_START.head<4>());
  ASSERT_TRUE(index3.succeeded && index2.succeeded && index1.succeeded && index0.succeeded);

  // Index 1 and 0 are the same discrete problem (RadauIIA eliminates y5 exactly at the stages)
  EXPECT_NEAR(index0.J, index1.J, 1.0e-8*index1.J);
  EXPECT_NEAR(index2.J, index1.J, 1.0e-5*index1.J);
  EXPECT_NEAR(index3.J, index1.J, 1.0e-5*index1.J);

  // The invariants hold exactly on the nodes
  EXPECT_LT(index2.invariants_violation(), 1.0e-12);
  EXPECT_LT(index1.invariants_violation(), 1.0e-12);
  EXPECT_LT(index0.invariants_violation(), 1.0e-12);

  // With one substep, multiple shooting and collocation (same mesh, controls, tableau, and node
  // corrections) are the same discrete problem
  ASSERT_TRUE(single.succeeded);
  Collocation<Real, 3, 5, 1, 2> collocation(std::make_shared<GerdtsPendulum<1>>(), RadauIIA5Tableau<Real>(), n);
  ASSERT_TRUE(collocation.simulate_guess(GERDTS_START, [] (Real) {return Vector1(0.0);}));
  Real J_collocation;
  ASSERT_EQ(ipopt_solve(collocation, options(), J_collocation), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_collocation, single.J, 1.0e-9*single.J);
}

/*\
 |  Weight of the invariants and inconsistent initial data
\*/

TEST(DAEOptimalControl, InvariantsWeight) {
  // From epsilon = 1 (continuity and invariants weighted equally) to epsilon = 0 (exact invariants)
  Real violation_previous{std::numeric_limits<Real>::infinity()}, J_exact{0.0};
  std::vector<Real> J;
  for (Real epsilon : {1.0, 1.0e-6, 0.0}) {
    Radau<GerdtsPendulum<1>, 5, 1, 2> solve(std::make_shared<GerdtsPendulum<1>>(), 30, 2, GERDTS_START, epsilon);
    ASSERT_TRUE(solve.succeeded) << "epsilon = " << epsilon;
    EXPECT_LE(solve.invariants_violation(), violation_previous) << "epsilon = " << epsilon;
    violation_previous = solve.invariants_violation();
    J.push_back(solve.J);
    J_exact = solve.J;
    if (epsilon == 1.0) {EXPECT_LT(solve.invariants_violation(), 1.0e-4);} // of the order of the local error
  }
  EXPECT_LT(violation_previous, 1.0e-12);
  EXPECT_NEAR(J[1], J_exact, 1.0e-8*J_exact);
}

TEST(DAEOptimalControl, InconsistentInitialData) {
  // The initial position is off the circle: the least-squares correction of the first node absorbs it
  auto problem = std::make_shared<GerdtsPendulum<1>>();
  problem->y1_ini = 1.01;
  Eigen::Vector<Real, 5> start{GERDTS_START};
  start(0) = 1.01;
  Radau<GerdtsPendulum<1>, 5, 1, 2> solve(problem, 30, 2, start);
  ASSERT_TRUE(solve.succeeded);
  EXPECT_LT(solve.invariants_violation(), 1.0e-12);
  EXPECT_NEAR(solve.nlp.x(solve.nlp.guess(), 0)(0), 1.0, 1.0e-12); // corrected onto the circle
  EXPECT_GT(solve.nlp.mu(solve.nlp.guess(), 0).norm(), 1.0e-3);     // by a visible correction
}

/*\
 |  Bounded swing-up (several swings), index 1 and index 0 with invariants, against the angle ODE
\*/

// Swing-up guess: the angle goes from 0 to pi with zero velocity at the ends, with consistent tension.
template <Integer N>
Eigen::Vector<Real, N> swing_up_guess(Real const t, Real const T)
{
  Real const q{M_PI*(1.0 - std::cos(M_PI*t/T))/2.0}, w{M_PI*M_PI/(2.0*T)*std::sin(M_PI*t/T)};
  Eigen::Vector<Real, N> x;
  x(X) = std::sin(q); x(Y) = -std::cos(q); x(VX) = w*std::cos(q); x(VY) = w*std::sin(q);
  if constexpr (N > LAMBDA) {x(LAMBDA) = (9.81*std::cos(q) + w*w)/2.0;}
  return x;
}

TEST(DAEOptimalControl, SwingUpBoundedControl) {
  PendulumParameters const p;
  Real const T{5.0};
  std::vector<Real> const continuation{UNBOUNDED, 6.0, 4.0};
  Radau<PendulumAngleOCP, 2, 1, 0> ode(std::make_shared<PendulumAngleOCP>(T), 80, 1, Eigen::Vector2d::Zero(), 0.0, continuation);
  Radau<PendulumSwingUpOCP<PendulumIndex1<2>>, 5, 1, 2> index1(std::make_shared<PendulumSwingUpOCP<PendulumIndex1<2>>>(T),
    80, 1, p.hanging_at_rest<5>(), 0.0, continuation, [T] (Real t) {return swing_up_guess<5>(t, T);});
  Radau<PendulumSwingUpOCP<PendulumIndex0<2>>, 4, 1, 2> index0(std::make_shared<PendulumSwingUpOCP<PendulumIndex0<2>>>(T),
    80, 1, p.hanging_at_rest<4>(), 0.0, continuation, [T] (Real t) {return swing_up_guess<4>(t, T);});
  ASSERT_TRUE(ode.succeeded && index1.succeeded && index0.succeeded);
  EXPECT_NEAR(index1.J, ode.J, 1.0e-3*ode.J);
  EXPECT_NEAR(index0.J, index1.J, 1.0e-6*index1.J);
  EXPECT_NEAR(index1.nlp.u_mesh(index1.nlp.guess()).cwiseAbs().maxCoeff(), 4.0, 1.0e-6); // the bound is active
  EXPECT_LT(index1.invariants_violation(), 1.0e-12);
}

/*\
 |  Index-1 DAE against the equivalent ODE
\*/

TEST(DAEOptimalControl, Index1DAEMatchesODE) {
  Radau<SemiLinearDAEOCP, 2, 1, 0> dae(std::make_shared<SemiLinearDAEOCP>(), 20, 1, Eigen::Vector2d(0.5, 0.25));
  Radau<SemiLinearODEOCP, 1, 1, 0> ode(std::make_shared<SemiLinearODEOCP>(), 20, 1, Vector1(0.5));
  ASSERT_TRUE(dae.succeeded && ode.succeeded);
  EXPECT_NEAR(dae.J, ode.J, 1.0e-10);
  EXPECT_LT((dae.nlp.x_mesh(dae.nlp.guess()).row(0) - ode.nlp.x_mesh(ode.nlp.guess()).row(0)).norm(), 1.0e-9);
}

/*\
 |  Attitude slew with the unit quaternion invariant
\*/

TEST(DAEOptimalControl, AttitudeSlew) {
  auto problem = std::make_shared<AttitudeSlew>(1.0, 2.0);
  Eigen::Vector<Real, 7> start(Eigen::Vector<Real, 7>::Zero());
  start(0) = 1.0;
  Radau<AttitudeSlew, 7, 3, 1> coarse(problem, 10, 1, start);
  Radau<AttitudeSlew, 7, 3, 1> fine(problem, 20, 1, start);
  ASSERT_TRUE(coarse.succeeded && fine.succeeded);
  EXPECT_LT(fine.invariants_violation(), 1.0e-12);

  // Piecewise-constant torques converge to the analytical optimum as O(h^2)
  Real const error_coarse{std::abs(coarse.J - problem->optimal_cost())};
  Real const error_fine{std::abs(fine.J - problem->optimal_cost())};
  EXPECT_LT(error_fine, 1.0e-2*problem->optimal_cost());
  EXPECT_GT(error_coarse/error_fine, 3.5);

  // Rotation about z only
  for (Integer k{0}; k < 20; ++k) {
    EXPECT_LT(fine.nlp.u(fine.nlp.guess(), k).head<2>().norm(), 1.0e-8);
  }
}

/*\
 |  Kepler orbit over 2 periods: invariants exact on the nodes, drift without them
\*/

TEST(DAEOptimalControl, KeplerOrbitInvariants) {
  Real const e{0.5};
  auto with_invariants = std::make_shared<KeplerOrbit<2>>(e, 2);
  Eigen::Vector4d const start{with_invariants->periapsis()};
  ShootingSolve<RK4Tableau<Real>, 4, 1, 2, KeplerOrbit<2>> with(with_invariants, 60, 4, start);
  ShootingSolve<RK4Tableau<Real>, 4, 1, 0, KeplerOrbit<0>> without(std::make_shared<KeplerOrbit<0>>(e, 2), 60, 4, start);
  ASSERT_TRUE(with.succeeded && without.succeeded);

  // Without the invariants, energy and angular momentum drift; with them, they are exact
  Real drift{0.0};
  for (Integer k{0}; k <= 60; ++k) {
    drift = std::max(drift, with_invariants->system()->h(without.nlp.x(without.nlp.guess(), k), 0.0).lpNorm<Eigen::Infinity>());
  }
  EXPECT_GT(drift, 1.0e-6);
  EXPECT_LT(with.invariants_violation(), 1.0e-12);

  // After two periods the orbit is back at the periapsis (up to the phase error of the integrator)
  EXPECT_LT((with.nlp.x(with.nlp.guess(), 60) - start).norm(), 5.0e-2);

  // Gauss-Legendre preserves the (quadratic) angular momentum by itself, not the energy
  ShootingSolve<GaussLegendre4Tableau<Real>, 4, 1, 0, KeplerOrbit<0>> gauss(std::make_shared<KeplerOrbit<0>>(e, 2), 60, 4, start);
  ASSERT_TRUE(gauss.succeeded);
  Real momentum_drift{0.0}, energy_drift{0.0};
  for (Integer k{0}; k <= 60; ++k) {
    Eigen::Vector2d const h{with_invariants->system()->h(gauss.nlp.x(gauss.nlp.guess(), k), 0.0)};
    energy_drift = std::max(energy_drift, std::abs(h(0)));
    momentum_drift = std::max(momentum_drift, std::abs(h(1)));
  }
  EXPECT_LT(momentum_drift, 1.0e-11);
  EXPECT_GT(energy_drift, 1.0e-8);
}

TEST(DAEOptimalControl, DependentInvariantsNeedPositiveWeight) {
  // On the circular orbit energy and angular momentum are dependent: a small positive weight keeps the
  // node least squares regular, with invariants violations of the order of the weight. Near the
  // circular orbit the invariants are degenerate (with L = 1, the energy is 1/2 v_r^2 + 1/2 (1/r - 1)^2),
  // so a violation delta only bounds the radius error by sqrt(delta)
  auto circular = std::make_shared<KeplerOrbit<2>>(0.0, 2);
  ShootingSolve<RK4Tableau<Real>, 4, 1, 2, KeplerOrbit<2>> solve(circular, 40, 3, circular->periapsis(), 1.0e-10);
  ASSERT_TRUE(solve.succeeded);
  Real const violation{solve.invariants_violation()};
  EXPECT_LT(violation, 1.0e-9);
  EXPECT_NEAR(solve.nlp.x(solve.nlp.guess(), 40).head<2>().norm(), 1.0, 10.0*std::sqrt(violation));
}

/*\
 |  Double pendulum: Cartesian index 1 with 4 invariants, against minimal coordinates
\*/

TEST(DAEOptimalControl, DoublePendulumMechanism) {
  DoublePendulumParameters const p;
  Eigen::Vector<Real, 10> cartesian(Eigen::Vector<Real, 10>::Zero());
  cartesian(1) = -p.L1;
  cartesian(3) = -p.L1 - p.L2;
  cartesian(8) = (p.m1 + p.m2)*p.g/(2.0*p.L1);
  cartesian(9) = p.m2*p.g/(2.0*p.L2);
  Radau<DoublePendulumCartesian, 10, 1, 4> dae(std::make_shared<DoublePendulumCartesian>(), 20, 1, cartesian);
  Radau<DoublePendulumAngles, 4, 1, 0> ode(std::make_shared<DoublePendulumAngles>(), 20, 1, Eigen::Vector4d::Zero());
  ASSERT_TRUE(dae.succeeded && ode.succeeded);
  EXPECT_NEAR(dae.J, ode.J, 1.0e-3*ode.J);
  EXPECT_LT(dae.invariants_violation(), 1.0e-12);
}

/*\
 |  Substeps where the dynamics is stiff (fixed mesh)
\*/

TEST(DAEOptimalControl, SubstepsForStiffness) {
  Vector1 const start(0.0);
  Radau<StiffTracking, 1, 1, 0> reference(std::make_shared<StiffTracking>(), 10, 10, start);
  ASSERT_TRUE(reference.succeeded);

  // RK4 with one step per interval is unstable (h lambda = 5), with 10 substeps it is accurate
  ShootingSolve<RK4Tableau<Real>, 1, 1, 0, StiffTracking> unstable(std::make_shared<StiffTracking>(), 10, 1, start);
  ShootingSolve<RK4Tableau<Real>, 1, 1, 0, StiffTracking> stable(std::make_shared<StiffTracking>(), 10, 10, start);
  ASSERT_TRUE(stable.succeeded);
  EXPECT_NEAR(stable.J, reference.J, 1.0e-3*reference.J);
  EXPECT_TRUE(!unstable.succeeded || std::abs(unstable.J - reference.J) > 1.0e-2*reference.J);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
