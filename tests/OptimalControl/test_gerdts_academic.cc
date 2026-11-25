/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Small DAE optimal control examples of M. Gerdts, Optimal Control of ODEs and DAEs (2nd edition, De
// Gruyter, 2023), solved by multiple shooting with IPOPT (tolerance 1e-12), against exact optima.

#include <gtest/gtest.h>

#include "GerdtsProblemsAD.hh"
#include "Gerdts/Shooting.hh"

using namespace TestOCP;

/*\
 |  Example 3.1.13: cost on the algebraic state. The control is constant, so the piecewise constant
 |  controls contain the exact solution, J* = a^2/(2 + tf).
\*/

TEST(GerdtsAcademic, AlgebraicCost) {
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

TEST(GerdtsAcademic, BangSingularBang) {
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

TEST(GerdtsAcademic, MinimumEnergyStateConstraint) {
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
 |  Examples 1.1.17 and 5.7.4 (O. Kostyukova): the state constraint x3 = x1 >= 0 becomes active at about
 |  t = 2 (where the book observes oscillations of the controls). The book gives no cost; the index-2 and
 |  index-1 formulations agree, and the cost converges with order two to J* = 2.8 (by Richardson
 |  extrapolation: 2.80000000 from 120 and 240 intervals).
\*/

TEST(GerdtsAcademic, Kostyukova) {
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
 |  Example 7.1.6: every feasible control is optimal (J = 0), for alpha = 1 and for the indefinite
 |  case alpha = -1 of the book. IPOPT converges (tolerance 1e-12) although the reduced Hessian
 |  vanishes in the control directions, and keeps the starting controls.
\*/

TEST(GerdtsAcademic, Degenerate) {
  auto const wave{[] (Real const t) {return Eigen::Vector<Real, 1>(0.5*std::sin(6.0*t));}};
  for (Real const alpha : {1.0, -1.0}) {
    auto index2 = std::make_shared<Degenerate<2>>();
    auto index1 = std::make_shared<Degenerate<1>>();
    index2->alpha = index1->alpha = alpha;
    Shooting<2, 1, 0> s2(index2, 20);
    Shooting<2, 1, 1> s1(index1, 20);
    ASSERT_TRUE(s2.solve(Eigen::Vector2d::Zero(), wave));
    ASSERT_TRUE(s1.solve(Eigen::Vector2d::Zero(), wave));
    EXPECT_NEAR(s2.J, 0.0, 1.0e-12);
    EXPECT_NEAR(s1.J, 0.0, 1.0e-12);
    for (Integer k{0}; k < 20; ++k) {
      EXPECT_NEAR(s2.u(k)(0), s1.u(k)(0), 1.0e-6);
      EXPECT_NEAR(s2.x(k+1)(0), 0.0, 1.0e-12);              // x = 0
      EXPECT_NEAR(s2.x(k+1)(1), s2.u(k)(0), 1.0e-10);       // y = u
    }
  }
}

/*\
 |  Example 1.1.17 as first stated, with the control u2 in the algebraic constraint: the transcription
 |  with piecewise constant controls converges, at first order, to about 8.963 instead of the optimum
 |  2.8 of the equivalent smoothed problem (Kostyukova test above; see the problem for the reason).
\*/

TEST(GerdtsAcademic, KostyukovaUnsmoothed) {
  Eigen::Vector4d const x_ini(1.0, 2.0/3.0, 1.0, 0.0);
  auto const guess{[] (Real) {return Eigen::Vector2d(0.0, 2.0/3.0);}};
  std::vector<Real> J;
  for (Integer const n : {60, 120, 240}) {
    Shooting<4, 2, 0> s(std::make_shared<KostyukovaUnsmoothed>(), n);
    ASSERT_TRUE(s.solve(x_ini, guess));
    J.push_back(s.J);
    for (Integer k{0}; k < n; ++k) {EXPECT_NEAR(s.x(k+1)(1), s.u(k)(1), 1.0e-12);} // x2 = u2 at the nodes
  }
  EXPECT_NEAR((J[2] - J[1])/(J[1] - J[0]), 0.5, 0.02);   // first order
  EXPECT_NEAR(2.0*J[2] - J[1], 8.963, 2.0e-3);           // first-order extrapolation
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
