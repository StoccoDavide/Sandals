/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Application examples of M. Gerdts, Optimal Control of ODEs and DAEs (2nd edition, De Gruyter, 2023):
// the chemical reactor, the emergency landing (free and on a prescribed path) and the Stokes flow,
// solved by multiple shooting with IPOPT (tolerance 1e-12).

#include <gtest/gtest.h>

#include "GerdtsProblemsAD.hh"
#include "Gerdts/Shooting.hh"

using namespace TestOCP;

/*\
 |  Example 8.2.8: semi-batch chemical reactor, an index-1 DAE whose states range from 0 to 1.5e5. The
 |  book shows the solution in figures only; the cost converges with order two (to about -477.22), and
 |  the temperature satisfies the energy balance at the nodes.
\*/

TEST(GerdtsApplications, ChemicalReactor) {
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

TEST(GerdtsApplications, EmergencyLanding) {
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

TEST(GerdtsApplications, EmergencyLandingPrescribedPath) {
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

/*\
 |  Examples 1.0.4, 1.1.11 and 3.1.12: distributed control of the Stokes flow on the unit square,
 |  discretized on a 4 x 4 grid (18 velocities, 9 pressures, 18 controls). Linear-quadratic: one IPOPT
 |  iteration, the index-2 and index-1 formulations agree, and the flow stays divergence free. The
 |  desired flow, sampled on the grid, is not discretely divergence free: the part of it along the
 |  discrete gradients cannot be tracked, which bounds the cost from below (2.4825); the optimum is
 |  above it by the cost of the controls (2.8755, 2.8546, 2.8520 with 20, 40, 80 intervals).
\*/

TEST(GerdtsApplications, Stokes) {
  using Grid = StokesGrid<4>;
  Integer const n{20};
  auto index2 = std::make_shared<Stokes<4, 2>>();
  Shooting<Grid::N, Grid::NV, 0> s2(index2, n);
  Shooting<Grid::N, Grid::NV, Grid::NP> s1(std::make_shared<Stokes<4, 1>>(), n);
  Eigen::Vector<Real, Grid::N> const rest(Eigen::Vector<Real, Grid::N>::Zero());
  ASSERT_TRUE(s2.nlp.simulate_guess(rest, [] (Real) {return Eigen::Vector<Real, Grid::NV>::Zero().eval();}));
  Real const J_uncontrolled{s2.nlp.objective(s2.nlp.guess())};
  ASSERT_TRUE(s2.solve(rest));
  ASSERT_TRUE(s1.solve(rest));
  EXPECT_NEAR(s2.J, s1.J, 1.0e-10*s2.J);
  EXPECT_LT(s2.J, 0.3*J_uncontrolled);

  // Lower bound: half the time integral of the untrackable part (Simpson's rule)
  Integer const m{2000};
  Real bound{0.0};
  for (Integer k{0}; k <= m; ++k) {bound += (k == 0 || k == m ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0))*index2->unreachable(2.0*k/m);}
  bound *= 0.5*(2.0/m)/3.0;
  EXPECT_NEAR(bound, 2.4825, 1.0e-4);
  EXPECT_GT(s2.J, bound);
  EXPECT_LT(s2.J, 1.2*bound);
  for (Integer k{1}; k <= n; ++k) {
    Eigen::Vector<Real, Grid::NV> const v{s2.x(k).head<Grid::NV>()};
    EXPECT_LT(index2->divergence<Real>(v).lpNorm<Eigen::Infinity>(), 1.0e-10);
    EXPECT_LT((v - s1.x(k).head<Grid::NV>()).lpNorm<Eigen::Infinity>(), 1.0e-8);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
