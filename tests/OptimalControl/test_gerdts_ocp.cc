/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// The DAE optimal control examples of M. Gerdts, Optimal Control of ODEs and DAEs (2nd edition, De
// Gruyter, 2023), solved by direct multiple shooting with IPOPT (tolerance 1e-12) in the formulation
// of the book and index-reduced with invariants, against the analytical optima and the numerical
// results reported in the book. The minimum energy (1.0.1, 3.1.11) and pendulum (1.1.9) examples are
// in test_dae_ocp.cc.

#include <gtest/gtest.h>

#include "Sandals/Ipopt.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"

#include "GerdtsProblemsAD.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;

// Quiet IPOPT with the default tolerance (1e-12); SANDALS_IPOPT_DEBUG=1 prints the iterations.
IpoptOptions options()
{
  IpoptOptions opts;
  opts.integers["print_level"] = std::getenv("SANDALS_IPOPT_DEBUG") ? 5 : 0;
  opts.strings["sb"] = "yes";
  opts.integers["max_iter"] = 500;
  return opts;
}

// A problem transcribed by multiple shooting with the Radau IIA (order 5) method, and solved by IPOPT.
template <Integer N, Integer U, Integer M>
struct Shooting
{
  using FunctionU = std::function<Eigen::Vector<Real, U>(Real)>;

  MultipleShooting<Real, 3, N, U, M> nlp;
  Real J{0.0};
  VectorX multipliers;

  template <typename Problem>
  Shooting(std::shared_ptr<Problem> problem, Integer const intervals, Integer const substeps = 1)
    : nlp(problem, RadauIIA5Tableau<Real>(), intervals, substeps) {}

  // Start from a forward simulation (zero controls by default), and solve
  bool solve(Eigen::Vector<Real, N> const & x_ini, FunctionU const & u_guess = nullptr, IpoptOptions const & opts = options())
  {
    auto const zero{[] (Real) {return Eigen::Vector<Real, U>::Zero().eval();}};
    return nlp.simulate_guess(x_ini, u_guess ? u_guess : zero) && resolve(opts);
  }

  // Solve again from the current solution (e.g., after changing the bounds)
  bool resolve(IpoptOptions const & opts = options())
  {return ipopt_solve(nlp, opts, J, nullptr, &multipliers) == Ipopt::Solve_Succeeded;}

  Eigen::Vector<Real, N> x(Integer const k) const {return nlp.x(nlp.guess(), k);}
  Eigen::Vector<Real, U> u(Integer const k) const {return nlp.u(nlp.guess(), k);}
  Real t(Integer const k) const {return nlp.t_mesh()(k);}
  Integer intervals() const {return nlp.intervals();}

  // The piecewise constant controls of the solution, as a function of time (to warm-start another
  // transcription on the same mesh)
  FunctionU controls() const
  {
    return [this] (Real const t) {
      Integer k{0};
      while (k + 1 < this->intervals() && this->t(k + 1) <= t) {++k;}
      return this->u(k);
    };
  }

  // Multiplier of the boundary condition i
  Real boundary_multiplier(Integer const i) const {return multipliers(nlp.boundary_row() + i);}
};

// Richardson extrapolation of second-order approximations on meshes of n and 2n intervals.
Real richardson(Real const coarse, Real const fine) {return fine + (fine - coarse)/3.0;}

/*\
 |  Example 3.1.13: cost on the algebraic state. The control is constant, so the piecewise constant
 |  controls contain the exact solution, J* = a^2/(2 + tf).
\*/

TEST(GerdtsExamples, AlgebraicCost) {
  auto index2 = std::make_shared<AlgebraicCost<2>>(1.0);
  auto index1 = std::make_shared<AlgebraicCost<1>>(1.0);
  Eigen::Vector3d const x_ini(1.0, 0.0, 0.0);

  Shooting<3, 1, 0> s2(index2, 10);
  Shooting<3, 1, 1> s1(index1, 10);
  ASSERT_TRUE(s2.solve(x_ini));
  ASSERT_TRUE(s1.solve(x_ini));
  EXPECT_NEAR(s2.J, index2->optimal_cost(), 1.0e-12);
  EXPECT_NEAR(s1.J, index2->optimal_cost(), 1.0e-12);
  for (Integer k{0}; k < 10; ++k) {
    EXPECT_NEAR(s2.u(k)(0), index2->optimal_control(), 1.0e-10);
    EXPECT_NEAR(s1.u(k)(0), index2->optimal_control(), 1.0e-10);
    EXPECT_NEAR(s2.x(k+1)(2), index2->optimal_control(), 1.0e-10); // y = u
    EXPECT_NEAR(s1.x(k+1)(1), 0.0, 1.0e-12);                       // invariant x2 = 0
  }
}

/*\
 |  Example 7.1.13: bang-singular-bang control, J* = 1/3. With a mesh multiple of three intervals the
 |  exact solution is in the discrete space; the bounds are not relaxed (IPOPT's bound_relax_factor),
 |  otherwise the bang arcs would exceed them by 1e-8.
\*/

TEST(GerdtsExamples, BangSingularBang) {
  IpoptOptions opts{options()};
  opts.numerics["bound_relax_factor"] = 0.0;
  Eigen::Vector3d const x_ini(1.0, 1.0, 0.0);

  Shooting<3, 1, 0> index2(std::make_shared<BangSingularBang<2>>(), 30);
  Shooting<3, 1, 1> index1(std::make_shared<BangSingularBang<1>>(), 30);
  ASSERT_TRUE(index2.solve(x_ini, nullptr, opts));
  ASSERT_TRUE(index1.solve(x_ini, nullptr, opts));
  EXPECT_NEAR(index2.J, 1.0/3.0, 1.0e-10);
  EXPECT_NEAR(index1.J, 1.0/3.0, 1.0e-10);
  for (Integer k{0}; k < 30; ++k) {
    Real const u{BangSingularBang<2>::optimal_control((index2.t(k) + index2.t(k+1))/2.0)};
    EXPECT_NEAR(index2.u(k)(0), u, 1.0e-7); // the singular control is determined at second order only
    EXPECT_NEAR(index1.u(k)(0), u, 1.0e-7);
    EXPECT_NEAR(index1.x(k)(0), index1.x(k)(1), 1.0e-12);
  }
}

/*\
 |  Example 3.2.13: the minimum energy problem (3.1.11) with the state constraint y1 <= l (a boundary
 |  arc on [3 l, 1 - 3 l] for l <= 1/6), J* = 4/(9 l). The state constraint is imposed at the nodes, and
 |  the cost converges with order two.
\*/

TEST(GerdtsExamples, MinimumEnergyStateConstraint) {
  Real const l{1.0/9.0};
  Eigen::Vector4d const y_ini(0.0, 1.0, 0.0, -3.0);
  Eigen::Vector4d upper(Eigen::Vector4d::Constant(INF));
  upper(0) = l;

  std::vector<Real> J2, J1;
  for (Integer const n : {60, 120}) {
    auto index2 = std::make_shared<MinimumEnergyIndex2>();
    auto index1 = std::make_shared<MinimumEnergyIndex1>();
    index2->x_bounds(Eigen::Vector4d::Constant(-INF), upper);
    index1->x_bounds(Eigen::Vector4d::Constant(-INF), upper);
    Shooting<4, 1, 0> s2(index2, n);
    Shooting<4, 1, 1> s1(index1, n);
    ASSERT_TRUE(s2.solve(y_ini));
    ASSERT_TRUE(s1.solve(y_ini));
    J2.push_back(s2.J);
    J1.push_back(s1.J);
    EXPECT_NEAR(s2.J, s1.J, 1.0e-10);

    // The boundary arc y1 = l on [3 l, 1 - 3 l]: the discrete solution touches the constraint at some
    // nodes and stays slightly below it between them (by 3e-5 with 60 intervals); IPOPT relaxes the
    // bounds by 1e-8
    for (Integer k{0}; k <= n; ++k) {
      EXPECT_LE(s2.x(k)(0), l + 2.0e-8);
      if (s2.t(k) > 3.0*l && s2.t(k) < 1.0 - 3.0*l) {EXPECT_NEAR(s2.x(k)(0), l, 1.0e-4);}
    }
  }
  EXPECT_NEAR(richardson(J2[0], J2[1]), 4.0/(9.0*l), 1.0e-5);
  EXPECT_NEAR(richardson(J1[0], J1[1]), 4.0/(9.0*l), 1.0e-5);
}

/*\
 |  Example 4.6.5 (and 1.1.19): the pendulum in the GGL formulation, on [0, 3]. The book solves the
 |  boundary value problem of the minimum principle by single shooting, and reports the multipliers of
 |  the initial conditions; those of y2(0) and y4(0) are unique (sigma_2 = -13.10630040572912, sigma_4 =
 |  -0.27699589145618; the others share the redundant information with the multipliers of the algebraic
 |  constraints at t = 0). The multipliers of the discrete boundary conditions converge to them with
 |  order two, in the GGL formulation and in the index-reduced ones with invariants (where the node
 |  least-squares multipliers take the redundant information).
\*/

TEST(GerdtsExamples, PendulumGGLMultipliers) {
  Real const sigma2{-13.10630040572912}, sigma4{-0.27699589145618};
  Eigen::Vector<Real, 6> const ggl_start{(Eigen::Vector<Real, 6>() << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0).finished()};
  Eigen::Vector<Real, 5> const start{ggl_start.head<5>()};

  std::vector<Real> ggl_sigma2, ggl_sigma4, index1_sigma2, index1_sigma4, ggl_y2;
  for (Integer const n : {30, 60}) {
    Shooting<6, 1, 0> ggl(std::make_shared<GerdtsPendulumGGL>(3.0), n, 2);
    Shooting<5, 1, 2> index1(std::make_shared<GerdtsPendulum<1>>(3.0), n, 2);
    ASSERT_TRUE(ggl.solve(ggl_start));
    ASSERT_TRUE(index1.solve(start));
    EXPECT_NEAR(ggl.J, index1.J, 1.0e-5*ggl.J);
    ggl_sigma2.push_back(ggl.boundary_multiplier(1));
    ggl_sigma4.push_back(ggl.boundary_multiplier(3));
    index1_sigma2.push_back(index1.boundary_multiplier(1));
    index1_sigma4.push_back(index1.boundary_multiplier(3));

    // The additional GGL multiplier vanishes on the exact solution (on the discrete one, it is of the
    // order of the discretization error)
    Real y2{0.0};
    for (Integer k{1}; k <= n; ++k) {y2 = std::max(y2, std::abs(ggl.x(k)(5)));}
    ggl_y2.push_back(y2);

    if (n == 30) {
      Shooting<5, 1, 1> index2(std::make_shared<GerdtsPendulum<2>>(3.0), n, 2);
      ASSERT_TRUE(index2.solve(start));
      EXPECT_NEAR(ggl.J, index2.J, 1.0e-5*ggl.J);
    }
  }
  EXPECT_LT(ggl_y2[0], 1.0e-3);
  EXPECT_LT(ggl_y2[1], ggl_y2[0]/4.0);
  EXPECT_NEAR(richardson(ggl_sigma2[0], ggl_sigma2[1]), sigma2, 1.0e-4);
  EXPECT_NEAR(richardson(ggl_sigma4[0], ggl_sigma4[1]), sigma4, 3.0e-5);
  EXPECT_NEAR(richardson(index1_sigma2[0], index1_sigma2[1]), sigma2, 1.0e-4);
  EXPECT_NEAR(richardson(index1_sigma4[0], index1_sigma4[1]), sigma4, 3.0e-5);
}

/*\
 |  Example 6.1.9: the GGL pendulum from 45 degrees below the horizontal, on [0, 4.5], with the control
 |  bounds |u| <= 1. Some control arcs saturate; the index-1 formulation with the invariants agrees.
\*/

TEST(GerdtsExamples, PendulumControlBounds) {
  Real const c{1.0/std::sqrt(2.0)};
  Eigen::Vector<Real, 1> const u_max(1.0);

  auto ggl = std::make_shared<GerdtsPendulumGGL>(4.5);
  ggl->x1_ini = c;
  ggl->x2_ini = -c;
  ggl->u_bounds(-u_max, u_max);
  auto index1 = std::make_shared<GerdtsPendulum<1>>(4.5);
  index1->y1_ini = c;
  index1->y2_ini = -c;
  index1->u_bounds(-u_max, u_max);

  Integer const n{60};
  Shooting<6, 1, 0> s_ggl(ggl, n, 2);
  Shooting<5, 1, 2> s_index1(index1, n, 2);
  ASSERT_TRUE(s_ggl.solve((Eigen::Vector<Real, 6>() << c, -c, 0.0, 0.0, 0.0, 0.0).finished()));
  ASSERT_TRUE(s_index1.solve((Eigen::Vector<Real, 5>() << c, -c, 0.0, 0.0, 0.0).finished()));
  EXPECT_NEAR(s_ggl.J, s_index1.J, 1.0e-6*s_ggl.J);

  Integer saturated{0};
  for (Integer k{0}; k < n; ++k) {
    if (std::abs(s_ggl.u(k)(0)) > 1.0 - 1.0e-6) {++saturated;}
    EXPECT_NEAR(s_ggl.u(k)(0), s_index1.u(k)(0), 1.0e-3);
  }
  EXPECT_GT(saturated, n/4);
  EXPECT_NEAR(s_ggl.x(n)(0), 0.0, 1.0e-10);
  EXPECT_LT(s_ggl.x(n)(1), -1.0 + 1.0e-10); // at rest at the bottom
}

/*\
 |  Example 5.7.3: the ski driver, minimum time, index 3. The steering angle may wind around (by 2 pi),
 |  which gives worse local minima: it is first bounded (|delta| <= 2.5), then released once the bound is
 |  found inactive. The book reports tf = 8.426448607129171 (with its own discretization on 401 points);
 |  here tf = 8.42724 (50 intervals), 8.42407 (100), 8.42318 (200), with the push force at its lower
 |  bound except for a short arc, as in the book. The index-3 formulation, started from the controls
 |  of the index-2 one, finds the same solution (the constraint is linear).
\*/

TEST(GerdtsExamples, SkiDriver) {
  Integer const n{50};
  auto index2 = std::make_shared<SkiDriver<2>>();
  auto index1 = std::make_shared<SkiDriver<1>>();
  auto index3 = std::make_shared<SkiDriver<3>>();
  auto const push{[] (Real) {return Eigen::Vector2d(0.0, -200.0);}};
  auto const bound_steering{[] (auto & problem, Real const delta_max) {
    Eigen::Vector<Real, 9> lower{problem->x_lower()}, upper{problem->x_upper()};
    lower(6) = -delta_max;
    upper(6) = delta_max;
    problem->x_bounds(lower, upper);
  }};

  Shooting<9, 2, 1> s2(index2, n);
  Shooting<9, 2, 2> s1(index1, n);
  bound_steering(index2, 2.5);
  bound_steering(index1, 2.5);
  ASSERT_TRUE(s2.solve(index2->initial_states(8.5), push));
  ASSERT_TRUE(s1.solve(index1->initial_states(8.5), push));
  Real steering{0.0};
  for (Integer k{0}; k <= n; ++k) {steering = std::max(steering, std::abs(s2.x(k)(6)));}
  EXPECT_LT(steering, 2.0); // the bound is inactive
  bound_steering(index2, INF);
  bound_steering(index1, INF);
  ASSERT_TRUE(s2.resolve());
  ASSERT_TRUE(s1.resolve());

  Real const T{s2.x(0)(8)};
  EXPECT_NEAR(s2.J, T, 1.0e-12);
  EXPECT_NEAR(T, SkiDriver<2>::REFERENCE_TIME, 1.0e-3);
  EXPECT_NEAR(s1.J, s2.J, 1.0e-9);

  // Index 3 from the controls of the index-2 solution
  Shooting<9, 2, 0> s3(index3, n);
  ASSERT_TRUE(s3.solve(index3->initial_states(T), s2.controls()));
  EXPECT_NEAR(s3.J, s2.J, 1.0e-9);

  // The final conditions, the plane, the constant normal force, and the push at its lower bound
  Integer pushing{0};
  for (Integer k{0}; k < n; ++k) {if (s2.u(k)(1) < -200.0 + 1.0e-6) {++pushing;}}
  EXPECT_GT(pushing, 3*n/4);
  EXPECT_NEAR(s3.x(n)(0), -100.0, 1.0e-10);
  EXPECT_NEAR(s3.x(n)(1), -20.0, 1.0e-10);
  EXPECT_NEAR(s3.x(n)(4), 0.0, 1.0e-10);
  for (Integer k{1}; k <= n; ++k) {
    EXPECT_NEAR(s2.x(k)(2), s2.x(k)(0)*index2->tan_psi, 1.0e-10);
    EXPECT_NEAR(s2.x(k)(7), index2->normal_force(), 1.0e-6*std::abs(index2->normal_force()));
  }
}

/*\
 |  Examples 1.1.17 and 5.7.4 (O. Kostyukova): the state constraint x3 = x1 >= 0 becomes active at about
 |  t = 2 (where the book observes oscillations of the controls). The book gives no cost; the index-2 and
 |  index-1 formulations agree, and the cost converges with order two to J* = 2.8 (by Richardson
 |  extrapolation: 2.80000000 from 120 and 240 intervals).
\*/

TEST(GerdtsExamples, Kostyukova) {
  Eigen::Vector<Real, 5> const x_ini{(Eigen::Vector<Real, 5>() << 1.0, 2.0/3.0, 1.0, 2.0/3.0, 0.0).finished()};
  std::vector<Real> J;
  for (Integer const n : {60, 120}) {
    Shooting<5, 2, 0> index2(std::make_shared<Kostyukova<2>>(), n);
    Shooting<5, 2, 1> index1(std::make_shared<Kostyukova<1>>(), n);
    ASSERT_TRUE(index2.solve(x_ini));
    ASSERT_TRUE(index1.solve(x_ini));
    EXPECT_NEAR(index2.J, index1.J, 1.0e-10);
    J.push_back(index2.J);
    for (Integer k{0}; k <= n; ++k) {
      EXPECT_NEAR(index2.x(k)(2), index2.x(k)(0), 1.0e-10); // x3 = x1
      if (index2.t(k) > 2.2) {EXPECT_NEAR(index2.x(k)(2), 0.0, 2.0e-8);} // on the constraint (relaxed by 1e-8)
    }
  }
  EXPECT_NEAR(richardson(J[0], J[1]), 2.8, 1.0e-6);
}

/*\
 |  Example 8.2.8: semi-batch chemical reactor, an index-1 DAE whose states range from 0 to 1.5e5. The
 |  book shows the solution in figures only; the cost converges with order two (to about -477.22), and
 |  the temperature satisfies the energy balance at the nodes.
\*/

TEST(GerdtsExamples, ChemicalReactor) {
  std::vector<Real> J;
  for (Integer const n : {40, 80, 160}) {
    Shooting<6, 2, 0> reactor(std::make_shared<ChemicalReactor>(), n, 2);
    ASSERT_TRUE(reactor.solve(ChemicalReactor::initial_states()));
    J.push_back(reactor.J);
    for (Integer k{0}; k <= n; ++k) {
      Eigen::Vector<Real, 6> const x{reactor.x(k)};
      EXPECT_NEAR(x(5), ChemicalReactor::temperature(x.head<4>(), x(4)), 1.0e-9*x(5));
    }
  }
  EXPECT_LT(J[2], J[1]);
  EXPECT_LT(J[1], J[0]);
  EXPECT_NEAR(richardson(J[0], J[1]), richardson(J[1], J[2]), 1.0e-2);
}

/*\
 |  Examples 1.1.10 and 6.1.10: emergency landing of a hypersonic flight system with maximum range, with
 |  the dynamic pressure constraint. The book reports tf = 726.57 and the cost -0.765 (200 intervals),
 |  with the pressure constraint active at the normalized time 0.185-0.19; here, with 50 intervals,
 |  tf = 726.21 and -0.76518 (with 200 intervals: 726.574 and -0.76488, active at 0.185). The states
 |  range from 1 (angles) to 3e4 (altitude): the NLP is scaled by their nominal magnitudes.
\*/

TEST(GerdtsExamples, EmergencyLanding) {
  Integer const n{50};
  auto problem = std::make_shared<EmergencyLanding>();
  Shooting<7, 2, 0> landing(problem, n);
  ASSERT_TRUE(landing.solve(EmergencyLanding::initial_states(700.0), [] (Real) {return Eigen::Vector2d(0.1, 0.0);}));
  EXPECT_NEAR(landing.J, EmergencyLanding::REFERENCE_COST, 1.0e-3);
  EXPECT_NEAR(landing.x(0)(6), EmergencyLanding::REFERENCE_TIME, 1.0);
  EXPECT_NEAR(landing.x(n)(3), EmergencyLanding::H_END, 1.0e-8);

  // The dynamic pressure constraint is active once, at the normalized time ~0.19
  Real q_max{0.0}, t_max{0.0};
  for (Integer k{0}; k <= n; ++k) {
    Real const q{EmergencyLanding::pressure(landing.x(k)(0), landing.x(k)(3))};
    if (q > q_max) {q_max = q; t_max = landing.t(k);}
  }
  EXPECT_NEAR(q_max, EmergencyLanding::Q_MAX, 1.0e-6*EmergencyLanding::Q_MAX);
  EXPECT_NEAR(t_max, 0.185, 0.01);
}

/*\
 |  Example 1.1.10 on a prescribed path: the bank angle is an algebraic state of index 3. The index-2 and
 |  index-1 formulations with invariants converge (tolerance 1e-12) to the same solution, also on the
 |  finer mesh. The index-3 formulation is not tested: its stage Jacobians have condition ~h^-3, and the
 |  round-off of the integrated states (~1e-10) keeps IPOPT's optimality error above 1e-8 (it stops at
 |  its acceptable level); this is what the index reduction with invariants removes. The initial
 |  azimuth is the consistent one of the book, chi(0) = -0.9955.
\*/

TEST(GerdtsExamples, EmergencyLandingPrescribedPath) {
  Integer const n{50};
  using Path1 = EmergencyLandingPath<1>;
  auto const lift{[] (Real) {return Eigen::Vector<Real, 1>(0.1);}};
  Shooting<8, 1, 2> index1(std::make_shared<Path1>(), n);
  Shooting<8, 1, 1> index2(std::make_shared<EmergencyLandingPath<2>>(), n);
  ASSERT_TRUE(index1.solve(Path1::initial_states(700.0), lift));
  ASSERT_TRUE(index2.solve(Path1::initial_states(700.0), lift));
  EXPECT_NEAR(index1.J, index2.J, 1.0e-10);
  EXPECT_NEAR(index1.x(0)(7), index2.x(0)(7), 1.0e-6);
  EXPECT_NEAR(index1.x(0)(2), Path1::REFERENCE_CHI0, 1.0e-4);

  // Twice as many intervals, from the simulation as well
  Shooting<8, 1, 2> fine1(std::make_shared<Path1>(), 2*n);
  Shooting<8, 1, 1> fine2(std::make_shared<EmergencyLandingPath<2>>(), 2*n);
  ASSERT_TRUE(fine1.solve(Path1::initial_states(700.0), lift));
  ASSERT_TRUE(fine2.solve(Path1::initial_states(700.0), lift));
  EXPECT_NEAR(fine1.J, fine2.J, 1.0e-10);

  // On the path, at the final altitude, and within the dynamic pressure limit
  for (Integer k{0}; k <= n; ++k) {
    Eigen::Vector<Real, 8> const x{index1.x(k)};
    EXPECT_NEAR(x(4) + x(5), Path1::LAMBDA0 + Path1::THETA0, 1.0e-12);
    EXPECT_NEAR(std::sin(x(2))*std::cos(x(4)) + std::cos(x(2)), 0.0, 1.0e-12);
    EXPECT_LT(Path1::pressure(x(0), x(3)), Path1::Q_MAX*(1.0 + 1.0e-8));
  }
  EXPECT_NEAR(index1.x(n)(3), Path1::H_END, 1.0e-8);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
