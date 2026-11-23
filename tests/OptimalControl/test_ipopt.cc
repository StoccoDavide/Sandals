/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <fstream>
#include <memory>

#include <gtest/gtest.h>

#include "Sandals/Ipopt.hh"
#include "Sandals/RungeKutta/RadauIIA3.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/GaussLegendre4.hh"

#include "OptimalControlProblems.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using Vector1 = Eigen::Vector<Real, 1>;
using Vector2 = Eigen::Vector<Real, 2>;

// Quiet and tight IPOPT options (exact Hessian unless overridden)
IpoptOptions quiet_options()
{
  IpoptOptions opts;
  opts.integers["print_level"] = 0;
  opts.strings["sb"] = "yes";
  opts.numerics["tol"] = 1.0e-12;
  return opts;
}

// Rest-to-rest double integrator: optimal u = 6 - 12 t, J* = 6.
template <typename TableauType, Integer S = TableauType::Vector::RowsAtCompileTime>
std::unique_ptr<Collocation<Real, S, 2, 1>> rest_to_rest(Integer const n, Real & J)
{
  auto ocp = std::make_shared<DoubleIntegratorOCP>(Vector2(0.0, 0.0), Vector2(1.0, 0.0));
  auto nlp = std::make_unique<Collocation<Real, S, 2, 1>>(ocp, TableauType(), n);
  nlp->guess([] (Real t) {return Vector2(t, 1.0);}, [] (Real) {return Vector1(0.0);});
  EXPECT_EQ(ipopt_solve(*nlp, quiet_options(), J), Ipopt::Solve_Succeeded);
  return nlp;
}

TEST(Ipopt, DoubleIntegratorRestToRest) {
  Real J20, J10, J40, J40_gl;
  auto nlp{rest_to_rest<RadauIIA5Tableau<Real>>(20, J20)};
  EXPECT_NEAR(J20, 6.0, 2.0e-2);

  // Boundary conditions and the analytical control (at the intervals midpoints)
  VectorX const & z{nlp->guess()};
  EXPECT_LT((nlp->x(z, 0) - Vector2(0.0, 0.0)).norm(), 1.0e-8);
  EXPECT_LT((nlp->x(z, 20) - Vector2(1.0, 0.0)).norm(), 1.0e-8);
  for (Integer k{0}; k < 20; ++k) {
    Real const t_mid{(nlp->t_mesh()(k) + nlp->t_mesh()(k+1))/2.0};
    EXPECT_NEAR(nlp->u(z, k)(0), 6.0 - 12.0*t_mid, 5.0e-2) << "interval " << k;
  }

  // The objective converges to the analytical optimum with mesh refinement
  rest_to_rest<RadauIIA5Tableau<Real>>(10, J10);
  rest_to_rest<RadauIIA5Tableau<Real>>(40, J40);
  rest_to_rest<GaussLegendre4Tableau<Real>>(40, J40_gl);
  EXPECT_LT(std::abs(J40 - 6.0), std::abs(J10 - 6.0)/8.0);
  EXPECT_NEAR(J40_gl, 6.0, 1.0e-2);
}

TEST(Ipopt, BrysonDenham) {
  // State constrained x0 <= 1/9, J* = 4/(9 * (1/9)) = 4
  auto ocp = std::make_shared<DoubleIntegratorOCP>(Vector2(0.0, 1.0), Vector2(0.0, -1.0));
  Real const inf{std::numeric_limits<Real>::infinity()};
  ocp->x_bounds(Vector2(-inf, -inf), Vector2(1.0/9.0, inf));
  Collocation<Real, 3, 2, 1> nlp(ocp, RadauIIA5Tableau<Real>(), 60);
  nlp.guess([] (Real t) {return Vector2(0.0, 1.0 - 2.0*t);}, [] (Real) {return Vector1(-2.0);});
  Real J;
  ASSERT_EQ(ipopt_solve(nlp, quiet_options(), J), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J, 4.0, 2.0e-2);
  auto const X{nlp.x_mesh(nlp.guess())};
  EXPECT_LE(X.row(0).maxCoeff(), 1.0/9.0 + 1.0e-8);
  EXPECT_NEAR(X.row(0).maxCoeff(), 1.0/9.0, 1.0e-6); // the constraint is active
}

TEST(Ipopt, DAEMatchesEquivalentODE) {
  // The index-1 DAE (x, z = x^2) and the ODE x' = -x + x^2 + u give the same optimum with RadauIIA
  Integer const n{20};
  auto dae = std::make_shared<SemiLinearDAEOCP>(0.5, 2.0);
  auto ode = std::make_shared<SemiLinearODEOCP>(0.5, 2.0);
  Collocation<Real, 3, 2, 1> nlp_dae(dae, RadauIIA5Tableau<Real>(), n);
  Collocation<Real, 3, 1, 1> nlp_ode(ode, RadauIIA5Tableau<Real>(), n);
  ASSERT_TRUE(nlp_dae.simulate_guess(Vector2(0.5, 0.25), [] (Real) {return Vector1(0.0);}));
  ASSERT_TRUE(nlp_ode.simulate_guess(Vector1(0.5), [] (Real) {return Vector1(0.0);}));

  Real J_dae, J_ode;
  IpoptOptions dbg{quiet_options()};
  if (std::getenv("SANDALS_IPOPT_DEBUG")) {dbg.integers["print_level"] = 5;}
  ASSERT_EQ(ipopt_solve(nlp_dae, dbg, J_dae), Ipopt::Solve_Succeeded);
  ASSERT_EQ(ipopt_solve(nlp_ode, quiet_options(), J_ode), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_dae, J_ode, 1.0e-8);
  EXPECT_LT(J_dae, 0.5*0.5*2.0); // better than doing nothing with x = 0.5 constant (rough bound)

  auto const X{nlp_dae.x_mesh(nlp_dae.guess())};
  auto const Y{nlp_ode.x_mesh(nlp_ode.guess())};
  EXPECT_LT((X.row(0) - Y.row(0)).norm(), 1.0e-6);
  EXPECT_LT((X.row(1) - X.row(0).cwiseProduct(X.row(0))).norm(), 1.0e-8); // algebraic constraint at nodes
  EXPECT_LT((nlp_dae.u_mesh(nlp_dae.guess()) - nlp_ode.u_mesh(nlp_ode.guess())).norm(), 1.0e-6);
}

TEST(Ipopt, VanDerPolConstraints) {
  auto ocp = std::make_shared<VanDerPolOCP>(5.0);
  Collocation<Real, 3, 2, 1> nlp(ocp, RadauIIA5Tableau<Real>(), 50);
  ASSERT_TRUE(nlp.simulate_guess(Vector2(1.0, 0.0), [] (Real) {return Vector1(0.0);}));
  Real J;
  ASSERT_EQ(ipopt_solve(nlp, quiet_options(), J), Ipopt::Solve_Succeeded);

  VectorX const & z{nlp.guess()};
  VectorX z_l, z_u, c_l, c_u;
  nlp.bounds(z_l, z_u, c_l, c_u);
  VectorX const c{nlp.constraints(z)};
  EXPECT_LT((c - c.cwiseMax(c_l).cwiseMin(c_u)).lpNorm<Eigen::Infinity>(), 1.0e-8); // feasibility
  EXPECT_LT((z - z.cwiseMax(z_l).cwiseMin(z_u)).lpNorm<Eigen::Infinity>(), 1.0e-8);
  EXPECT_LT((nlp.x(z, 0) - Vector2(1.0, 0.0)).norm(), 1.0e-8);
  EXPECT_NEAR(J, nlp.objective(z), 1.0e-12);

  // Values-only problem (finite differences derivatives) and limited-memory Hessian reach the same optimum
  auto ocp_fd = std::make_shared<VanDerPolOCPValuesOnly>(5.0);
  ocp_fd->u_bounds(Vector1(-1.5), Vector1(1.5));
  Collocation<Real, 3, 2, 1> nlp_fd(ocp_fd, RadauIIA5Tableau<Real>(), 50);
  ASSERT_TRUE(nlp_fd.simulate_guess(Vector2(1.0, 0.0), [] (Real) {return Vector1(0.0);}));
  Real J_fd;
  IpoptOptions opts_fd{quiet_options()};
  opts_fd.numerics["tol"] = 1.0e-10; // finite differences derivatives
  EXPECT_EQ(ipopt_solve(nlp_fd, opts_fd, J_fd), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_fd, J, 1.0e-5);

  Collocation<Real, 3, 2, 1> nlp_lbfgs(ocp, RadauIIA5Tableau<Real>(), 50);
  ASSERT_TRUE(nlp_lbfgs.simulate_guess(Vector2(1.0, 0.0), [] (Real) {return Vector1(0.0);}));
  IpoptOptions opts_lbfgs{quiet_options()};
  opts_lbfgs.strings["hessian_approximation"] = "limited-memory";
  opts_lbfgs.numerics["tol"] = 1.0e-10; // approximate Hessian
  opts_lbfgs.integers["max_iter"] = 1000;
  Real J_lbfgs;
  EXPECT_EQ(ipopt_solve(nlp_lbfgs, opts_lbfgs, J_lbfgs), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_lbfgs, J, 1.0e-5);
}

/*\
 |  Index-3 DAE: swing-up of the pendulum in Cartesian coordinates, validated against the angle ODE
\*/

using Vector5 = Eigen::Vector<Real, 5>;

// Swing-up guess: theta from 0 to pi with zero velocity at the ends, consistent tension.
struct SwingUpGuess
{
  Real T;
  Real theta(Real t) const {return M_PI*(1.0 - std::cos(M_PI*t/T))/2.0;}
  Real omega(Real t) const {return M_PI*M_PI/(2.0*T)*std::sin(M_PI*t/T);}
  Vector2 angle(Real t) const {return Vector2(this->theta(t), this->omega(t));}
  Vector5 cartesian(Real t) const
  {
    Real const q{this->theta(t)}, w{this->omega(t)};
    Vector5 x;
    x << std::sin(q), -std::cos(q), w*std::cos(q), w*std::sin(q), (9.81*std::cos(q) + w*w)/2.0;
    return x;
  }
};

// Replay piecewise-constant controls on the angle ODE with fine RK4 steps: returns (theta(T), omega(T)).
Vector2 replay_on_angle_ode(Eigen::Matrix<Real, 1, Eigen::Dynamic> const & U, VectorX const & mesh)
{
  PendulumAngle p;
  Vector2 x(0.0, 0.0);
  for (Integer k{0}; k < U.cols(); ++k) {
    Vector1 const u(U(0, k));
    Real const h{(mesh(k+1) - mesh(k))/200.0};
    Real t{mesh(k)};
    for (Integer s{0}; s < 200; ++s, t += h) {
      Vector2 const k1{p.f(x, u, t)}, k2{p.f(x + h/2.0*k1, u, t)}, k3{p.f(x + h/2.0*k2, u, t)}, k4{p.f(x + h*k3, u, t)};
      x += h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
    }
  }
  return x;
}

// Angle error modulo full turns (the final condition x(T) = 0 does not count revolutions).
Real angle_error(Real const theta, Real const target) {return std::abs(std::remainder(theta - target, 2.0*M_PI));}

TEST(Ipopt, Index3PendulumUnbounded) {
  // Unbounded control: the index-3 DAE solved directly converges to the ODE optimum
  Real const T{2.0};
  SwingUpGuess const g{T};
  auto dae = std::make_shared<PendulumIndex3OCP>(T);
  auto ode = std::make_shared<PendulumAngleOCP>(T);
  Real drift_prev{0.0};
  for (Integer n : {20, 40, 80}) {
    Collocation<Real, 3, 5, 1> nlp_dae(dae, RadauIIA5Tableau<Real>(), n);
    Collocation<Real, 3, 2, 1> nlp_ode(ode, RadauIIA5Tableau<Real>(), n);
    nlp_dae.guess([&g] (Real t) {return g.cartesian(t);}, [] (Real) {return Vector1(0.0);});
    nlp_ode.guess([&g] (Real t) {return g.angle(t);}, [] (Real) {return Vector1(0.0);});
    Real J_dae, J_ode;
    ASSERT_EQ(ipopt_solve(nlp_dae, quiet_options(), J_dae), Ipopt::Solve_Succeeded) << "n = " << n;
    ASSERT_EQ(ipopt_solve(nlp_ode, quiet_options(), J_ode), Ipopt::Solve_Succeeded) << "n = " << n;
    if (n >= 40) {EXPECT_NEAR(J_dae, J_ode, 1.0e-4*J_ode) << "n = " << n;}

    // Position constraint exact at the nodes (k >= 1), hidden velocity constraint converges as O(h^3)
    auto const X{nlp_dae.x_mesh(nlp_dae.guess())};
    Real drift{0.0};
    for (Integer k{1}; k <= n; ++k) {
      EXPECT_LT(std::abs(X(0, k)*X(0, k) + X(1, k)*X(1, k) - 1.0), 1.0e-9);
      drift = std::max(drift, std::abs(X(0, k)*X(2, k) + X(1, k)*X(3, k)));
    }
    if (drift_prev > 0.0) {EXPECT_GT(drift_prev/drift, 6.0) << "n = " << n;}
    drift_prev = drift;
    EXPECT_NEAR(X(1, n), 1.0, 1.0e-9); // upright
  }
}

TEST(Ipopt, Index3PendulumBoundedContinuation) {
  // Bounded control |u| <= 4 needs several swings; continuation on the bound from the unbounded
  // solution, warm-starting each solve. The control is validated by replaying it on the angle ODE.
  Real const T{5.0};
  Integer const n{80};
  SwingUpGuess const g{T};
  auto dae = std::make_shared<PendulumIndex3OCP>(T);
  auto ode = std::make_shared<PendulumAngleOCP>(T);
  Collocation<Real, 3, 5, 1> nlp_dae(dae, RadauIIA5Tableau<Real>(), n);
  Collocation<Real, 3, 2, 1> nlp_ode(ode, RadauIIA5Tableau<Real>(), n);
  nlp_dae.guess([&g] (Real t) {return g.cartesian(t);}, [] (Real) {return Vector1(0.0);});
  nlp_ode.guess([&g] (Real t) {return g.angle(t);}, [] (Real) {return Vector1(0.0);});

  Real J_dae{0.0}, J_ode{0.0};
  for (Real u_max : {std::numeric_limits<Real>::infinity(), 6.0, 4.0}) {
    dae->u_bounds(Vector1(-u_max), Vector1(u_max));
    ode->u_bounds(Vector1(-u_max), Vector1(u_max));
    ASSERT_EQ(ipopt_solve(nlp_dae, quiet_options(), J_dae), Ipopt::Solve_Succeeded) << "u_max = " << u_max;
    ASSERT_EQ(ipopt_solve(nlp_ode, quiet_options(), J_ode), Ipopt::Solve_Succeeded) << "u_max = " << u_max;
    EXPECT_NEAR(J_dae, J_ode, 1.0e-3*J_ode) << "u_max = " << u_max;
  }

  // The bound is active, and the control swings the real pendulum up
  auto const U{nlp_dae.u_mesh(nlp_dae.guess())};
  EXPECT_NEAR(U.cwiseAbs().maxCoeff(), 4.0, 1.0e-6);
  Vector2 const x_end{replay_on_angle_ode(U, nlp_dae.t_mesh())};
  EXPECT_LT(angle_error(x_end(0), M_PI), 5.0e-2);
  EXPECT_LT(std::abs(x_end(1)), 1.0e-1);
}

TEST(Ipopt, DerivativeChecker) {
  // IPOPT's own second-order derivative checker on the transcription (exact derivatives)
  auto ocp = std::make_shared<VanDerPolOCP>(1.0);
  Collocation<Real, 2, 2, 1> nlp(ocp, RadauIIA3Tableau<Real>(), 5);
  ASSERT_TRUE(nlp.simulate_guess(Vector2(1.0, 0.0), [] (Real t) {return Vector1(0.2*t);}));
  // Move away from zero components, where forward differences of quadratic terms are inaccurate
  VectorX z{nlp.guess()};
  for (Integer i{0}; i < z.size(); ++i) {z(i) += 0.1 + 0.05*std::sin(Real(i));}
  nlp.guess(z);
  IpoptOptions opts{quiet_options()};
  opts.strings["derivative_test"] = "second-order";
  // Forward differences: evaluate at the guess (no random perturbation) with a matching tolerance
  opts.numerics["point_perturbation_radius"] = 0.0;
  opts.numerics["derivative_test_perturbation"] = 1.0e-7;
  opts.numerics["derivative_test_tol"] = 1.0e-4;
  opts.integers["max_iter"] = 0;
  opts.strings["output_file"] = "test_ipopt_derivative_check.txt"; // the checker reports to the journal
  opts.integers["file_print_level"] = 5;
  Real J;
  ipopt_solve(nlp, opts, J);
  std::ifstream log("test_ipopt_derivative_check.txt");
  ASSERT_TRUE(log.good());
  std::string text((std::istreambuf_iterator<char>(log)), std::istreambuf_iterator<char>());
  EXPECT_NE(text.find("No errors detected by derivative checker"), std::string::npos)
    << text.substr(0, std::min<std::size_t>(text.size(), 4000));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
