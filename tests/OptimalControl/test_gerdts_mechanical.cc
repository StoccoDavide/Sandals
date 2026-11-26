/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Mechanical DAE examples of M. Gerdts, Optimal Control of ODEs and DAEs (2nd edition, De Gruyter,
// 2023): the pendulum (index 2 and 3), the ski driver and the half-car, solved by multiple shooting
// with IPOPT (tolerance 1e-12) in the formulation of the book and index-reduced with invariants.

#include <gtest/gtest.h>

#include "GerdtsProblemsAD.hh"
#include "Gerdts/Shooting.hh"
#include "Sandals/ShootingInterval.hh"
#include "Sandals/RungeKutta/RadauIIA3.hh"
#include "Sandals/RungeKutta/RK4.hh"
#include "Sandals/RungeKutta/GaussLegendre4.hh"

using namespace TestOCP;

/*\
 |  Example 4.6.5 (and 1.1.19): the pendulum in the GGL formulation, on [0, 3]. The book solves the
 |  boundary value problem of the minimum principle by single shooting, and reports the multipliers of
 |  the initial conditions; those of y2(0) and y4(0) are unique (sigma_2 = -13.10630040572912, sigma_4 =
 |  -0.27699589145618; the others share the redundant information with the multipliers of the algebraic
 |  constraints at t = 0). The multipliers of the discrete boundary conditions converge to them with
 |  order two, in the GGL formulation and in the index-reduced ones with invariants (where the node
 |  least-squares multipliers take the redundant information).
\*/

TEST(GerdtsMechanical, PendulumGGLMultipliers) {
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

TEST(GerdtsMechanical, PendulumControlBounds) {
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

TEST(GerdtsMechanical, SkiDriver) {
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
 |  Example 5.7.6: semi-active chassis control of a half-car over a pothole (index 3). All the
 |  formulations converge; index 2 and 3 agree to 1e-4 relative. Index 1 integrates the road
 |  acceleration, which jumps at the edges of the pothole, inside the intervals: its error is larger
 |  (1% with 100 intervals) and decreases with substeps (0.3% with four). The damping is used to
 |  reduce the chassis accelerations and the wheel load variations, against the constant damping 1/2.
\*/

TEST(GerdtsMechanical, HalfCar) {
  Integer const n{100};
  auto const half{[] (Real) {return Eigen::Vector2d(0.5, 0.5);}};
  Shooting<10, 2, 0> index3(std::make_shared<HalfCar<3>>(), n);
  Shooting<10, 2, 2> index2(std::make_shared<HalfCar<2>>(), n);
  Shooting<10, 2, 4> index1(std::make_shared<HalfCar<1>>(), n, 4);
  ASSERT_TRUE(index3.nlp.simulate_guess(HalfCar<3>::initial_states(), half));
  Real const J_passive{index3.nlp.objective(index3.nlp.guess())};
  ASSERT_TRUE(index3.solve(HalfCar<3>::initial_states(), half));
  ASSERT_TRUE(index2.solve(HalfCar<2>::initial_states(), half));
  ASSERT_TRUE(index1.solve(HalfCar<1>::initial_states(), half));
  EXPECT_NEAR(index3.J, index2.J, 1.0e-4*index2.J);
  EXPECT_NEAR(index1.J, index2.J, 5.0e-3*index2.J);
  EXPECT_LT(index3.J, 0.95*J_passive);

  // The wheels follow the road; the contact forces are those of the static load before the pothole
  for (Integer k{0}; k <= n; ++k) {
    Eigen::Vector<Real, 10> const x{index2.x(k)};
    EXPECT_NEAR(x(1), HalfCar<2>::road(index2.t(k))(0), 1.0e-12);
    if (index2.t(k) > 0.0 && index2.t(k) < 0.29) {EXPECT_NEAR(x(8), HalfCar<2>::initial_states()(8), 1.0e-6);}
  }
  for (Integer k{0}; k < n; ++k) {
    EXPECT_GE(index3.u(k).minCoeff(), -1.0e-8);
    EXPECT_LE(index3.u(k).maxCoeff(), 1.0 + 1.0e-8);
  }
}

/*\
 |  Example 4.3.4: the 2-stage Radau IIA method (order 3) on the pendulum from (1, 0, 0, 0, 0), with
 |  the steps h = 1/N on [0, 1], against a reference with the 3-stage method. The orders are those of
 |  Hairer, Lubich and Roche, as in Tables 4.1 and 4.2 of the book: index 2, order 3 for the states and
 |  2 for the multiplier; index 3, order 3 for the positions, 2 for the velocities and 1 for the
 |  multiplier. The index-3 errors of the velocities and of the multiplier also match the book (its
 |  method is linearized, one Newton step per step: its position errors are larger).
\*/

// The largest errors (positions, velocities, multiplier) on the nodes, with N steps
template <Integer Index>
Eigen::Vector3d pendulum_errors(Integer const N)
{
  constexpr Integer M{GerdtsPendulumSizes<Index>::M};
  GerdtsPendulum<Index> const pendulum(1.0);
  auto const integrate = [&pendulum] (auto const & tableau, Integer const steps, Integer const substeps) {
    constexpr Integer S{std::decay_t<decltype(tableau)>::Vector::RowsAtCompileTime};
    std::vector<Eigen::Vector<Real, 5>> x{(Eigen::Vector<Real, 5>() << 1.0, 0.0, 0.0, 0.0, 0.0).finished()};
    for (Integer k{0}; k < steps; ++k) {
      ShootingInterval<Real, S, 5, 1, M> step(pendulum, tableau, Real(k)/steps, Real(k + 1)/steps, substeps);
      EXPECT_TRUE(step.integrate(x.back(), Eigen::Vector<Real, 1>::Zero()));
      x.push_back(step.end().template head<5>());
    }
    return x;
  };
  auto const x{integrate(RadauIIA3Tableau<Real>(), N, 1)};
  auto const reference{integrate(RadauIIA5Tableau<Real>(), N, 5120/N)};
  Eigen::Vector3d error(Eigen::Vector3d::Zero());
  for (Integer k{1}; k <= N; ++k) {
    Eigen::Vector<Real, 5> const e{(x[k] - reference[k]).cwiseAbs()};
    error = error.cwiseMax(Eigen::Vector3d(e.head<2>().maxCoeff(), e.segment<2>(2).maxCoeff(), e(4)));
  }
  return error;
}

TEST(GerdtsMechanical, PendulumIntegrationOrder) {
  Eigen::Vector3d const index2_orders(3.0, 3.0, 2.0), index3_orders(3.0, 2.0, 1.0);
  Eigen::Vector3d const e2_coarse{pendulum_errors<2>(160)}, e2_fine{pendulum_errors<2>(320)};
  Eigen::Vector3d const e3_coarse{pendulum_errors<3>(160)}, e3_fine{pendulum_errors<3>(320)};
  for (Integer i{0}; i < 3; ++i) {
    EXPECT_NEAR(std::log2(e2_coarse(i)/e2_fine(i)), index2_orders(i), 0.05) << "index 2, component " << i;
    EXPECT_NEAR(std::log2(e3_coarse(i)/e3_fine(i)), index3_orders(i), 0.05) << "index 3, component " << i;
  }
  // Table 4.2, N = 320: velocities 0.767e-4, multiplier 0.872e-1
  EXPECT_NEAR(e3_fine(1), 0.767e-4, 0.03*0.767e-4);
  EXPECT_NEAR(e3_fine(2), 0.872e-1, 0.05*0.872e-1);
}

/*\
 |  Only stiffly accurate methods (Radau IIA) integrate DAEs: the transcription refuses the others on a
 |  system with algebraic states, and accepts them on an ODE (index 0, here with invariants).
\*/

TEST(GerdtsMechanical, DAEsNeedStifflyAccurateMethods) {
  Eigen::Vector<Real, 1> const zero(Eigen::Vector<Real, 1>::Zero());
  Eigen::Vector<Real, 5> const start{(Eigen::Vector<Real, 5>() << 1.0, 0.0, 0.0, 0.0, 0.0).finished()};
  GerdtsPendulum<2> const dae(1.0);
  GerdtsPendulum<0> const ode(1.0);
  ShootingInterval<Real, 4, 5, 1, 1> rk4(dae, RK4Tableau<Real>(), 0.0, 0.1, 1);
  ShootingInterval<Real, 2, 5, 1, 1> gauss(dae, GaussLegendre4Tableau<Real>(), 0.0, 0.1, 1);
  ShootingInterval<Real, 3, 5, 1, 1> radau(dae, RadauIIA5Tableau<Real>(), 0.0, 0.1, 1);
  ShootingInterval<Real, 4, 4, 1, 2> explicit_ode(ode, RK4Tableau<Real>(), 0.0, 0.1, 1);
  EXPECT_THROW(rk4.integrate(start, zero), std::runtime_error);
  EXPECT_THROW(gauss.integrate(start, zero), std::runtime_error);
  EXPECT_TRUE(radau.integrate(start, zero));
  EXPECT_TRUE(explicit_ode.integrate(start.head<4>(), zero));
}

/*\
 |  Example 1.1.20: docking of two bodies, a DAE whose index changes from one (free flight, lambda = 0)
 |  to three (docked) at tc = 2 pi. The bodies meet with the same position and velocity; the consistent
 |  contact forces are those of the book, lambda(tc) = -m1 m2/((m1 + m2) r^2) |v1|^2 p1(tc); the docked
 |  bodies keep their angular momentum (central force); index 3, 2 and 1 (with invariants) agree.
\*/

// The states on the nodes of a simulation from the given states (no optimization)
template <Integer Phase, Integer M>
std::vector<Eigen::Vector<Real, 10>> simulate_docking(Eigen::Vector<Real, 10> const & start, Integer const n)
{
  auto problem = std::make_shared<Docking<Phase>>();
  problem->start = start;
  Shooting<10, 1, M> s(problem, n, 2);
  EXPECT_TRUE(s.nlp.simulate_guess(start, [] (Real) {return Eigen::Vector<Real, 1>::Zero().eval();}));
  std::vector<Eigen::Vector<Real, 10>> x;
  for (Integer k{0}; k <= n; ++k) {x.push_back(s.x(k));}
  return x;
}

TEST(GerdtsMechanical, Docking) {
  using Model = DockingModel;
  Integer const n{40};

  // Free flight: the bodies meet at (r, 0) with the same velocity
  Eigen::Vector<Real, 10> const docking{simulate_docking<0, 0>(Model::initial_states(), n).back()};
  EXPECT_NEAR((docking.segment<2>(0) - docking.segment<2>(2)).norm(), 0.0, 1.0e-7);
  EXPECT_NEAR((docking.segment<2>(4) - docking.segment<2>(6)).norm(), 0.0, 1.0e-7);
  EXPECT_NEAR((docking.segment<2>(0) - Eigen::Vector2d(Model::R, 0.0)).norm(), 0.0, 1.0e-7);

  // The consistent contact forces of the book
  Eigen::Vector<Real, 10> start{docking};
  start.tail<2>() = Model::coupling(docking);
  EXPECT_NEAR(start(8), -Model::M1*Model::M2/((Model::M1 + Model::M2)*Model::R*Model::R)*Model::V*Model::V*Model::R, 1.0e-7);
  EXPECT_NEAR((start - Model::docked_states()).norm(), 0.0, 1.0e-7);

  // Docked motion in the three formulations
  auto const index3{simulate_docking<3, 0>(start, n)};
  auto const index2{simulate_docking<2, 2>(start, n)};
  auto const index1{simulate_docking<1, 4>(start, n)};
  Real const L{Model::angular_momentum(start)};
  for (Integer k{1}; k <= n; ++k) {
    for (auto const * x : {&index3[k], &index2[k], &index1[k]}) {
      EXPECT_NEAR(Model::angular_momentum(*x), L, 1.0e-8);
    }
    EXPECT_NEAR((index3[k].segment<2>(0) - index3[k].segment<2>(2)).norm(), 0.0, 1.0e-7);
    EXPECT_NEAR((index2[k].segment<2>(0) - index2[k].segment<2>(2)).norm(), 0.0, 1.0e-12);
    EXPECT_NEAR((index1[k].segment<2>(0) - index1[k].segment<2>(2)).norm(), 0.0, 1.0e-12);
    EXPECT_NEAR((index1[k].segment<2>(4) - index1[k].segment<2>(6)).norm(), 0.0, 1.0e-12);
    EXPECT_NEAR((index2[k].tail<2>() - Model::coupling(index2[k])).norm(), 0.0, 1.0e-10);
    EXPECT_NEAR((index1[k].tail<2>() - Model::coupling(index1[k])).norm(), 0.0, 1.0e-10);
    EXPECT_NEAR((index3[k].tail<2>() - Model::coupling(index3[k])).norm(), 0.0, 1.0e-4);
    EXPECT_NEAR((index3[k].head<4>() - index1[k].head<4>()).norm(), 0.0, 1.0e-6);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
