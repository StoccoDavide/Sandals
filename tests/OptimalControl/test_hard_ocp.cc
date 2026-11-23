/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Harder optimal control problems solved by multiple shooting with IPOPT: free final time, singular
// arcs, orbital dynamics, a nonconvex mechanism with active state constraints, boundary layers on a
// non-uniform mesh, and an invariant with a singular start.

#include <gtest/gtest.h>

#include "Sandals/Ipopt.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"

#include "HardProblemsAD.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using Vector1 = Eigen::Vector<Real, 1>;

IpoptOptions options()
{
  IpoptOptions opts;
  opts.integers["print_level"] = std::getenv("SANDALS_IPOPT_DEBUG") ? 5 : 0;
  opts.strings["sb"] = "yes";
  opts.integers["max_iter"] = 1000;
  return opts;
}

// Solve by multiple shooting (RadauIIA5) on a given mesh, from a forward simulation with the given
// controls (or from a states guess).
template <Integer N, Integer M, typename Problem>
struct Solve
{
  using NLP = MultipleShooting<Real, 3, N, 1, M>;
  using FunctionX = std::function<Eigen::Vector<Real, N>(Real)>;
  using FunctionU = std::function<Vector1(Real)>;

  std::shared_ptr<Problem> problem;
  NLP nlp;
  bool succeeded{false};
  Real J{0.0};
  Integer iterations{0};

  Solve(std::shared_ptr<Problem> p, VectorX const & mesh, Integer const substeps, Eigen::Vector<Real, N> const & x_ini,
    FunctionU const & u_guess, FunctionX const & x_guess = nullptr)
    : problem(p), nlp(p, RadauIIA5Tableau<Real>(), mesh, std::vector<Integer>(mesh.size() - 1, substeps))
  {
    if (x_guess) {nlp.guess(x_guess, u_guess);}
    else if (!nlp.simulate_guess(x_ini, u_guess)) {return;}
    std::vector<IpoptIteration> history;
    succeeded = ipopt_solve(nlp, options(), J, &history) == Ipopt::Solve_Succeeded;
    iterations = history.empty() ? 0 : history.back().iteration;
    std::printf("  %-22s iterations %4d  J = %.10f\n", p->name().c_str(), int(iterations), J);
  }

  Eigen::Vector<Real, N> x(Integer const k) const {return nlp.x(nlp.guess(), k);}
  Eigen::Vector<Real, N> x_end() const {return this->x(nlp.intervals());}
};

VectorX uniform(Real const t_end, Integer const n) {return VectorX::LinSpaced(n + 1, 0.0, t_end);}

/*\
 |  Goddard rocket: free final time and a singular arc
\*/

TEST(HardOptimalControl, GoddardRocket) {
  // Starting point: the constant thrust that burns the fuel in the guessed duration 0.2
  auto rocket = std::make_shared<GoddardRocket>();
  Eigen::Vector4d const start(1.0, 0.0, 1.0, 0.2);
  Solve<4, 0, GoddardRocket> solve(rocket, uniform(1.0, 60), 1, start, [] (Real) {return Vector1(1.0);});
  ASSERT_TRUE(solve.succeeded);
  EXPECT_NEAR(solve.x_end()(0), GoddardRocket::OPTIMAL_ALTITUDE, 2.0e-5);
  EXPECT_NEAR(solve.x_end()(2), GoddardRocket::FINAL_MASS, 1.0e-12);

  // Bang-singular-bang: full thrust at the start, an intermediate (singular) arc, then no thrust
  auto const U{solve.nlp.u_mesh(solve.nlp.guess())};
  EXPECT_NEAR(U(0, 0), GoddardRocket::U_MAX, 1.0e-6);
  EXPECT_NEAR(U(0, U.cols() - 1), 0.0, 1.0e-6);
  Integer singular{0};
  for (Integer k{0}; k < U.cols(); ++k) {if (U(0, k) > 0.05 && U(0, k) < GoddardRocket::U_MAX - 0.05) {++singular;}}
  EXPECT_GE(singular, 5);
}

/*\
 |  Low-thrust orbit raising: polar and Cartesian formulations
\*/

TEST(HardOptimalControl, OrbitRaising) {
  Real const T{OrbitRaisingData::DURATION};
  Eigen::Vector4d const polar_start(1.0, 0.0, 0.0, 1.0), cartesian_start(1.0, 0.0, 0.0, 1.0);
  // Starting point: tangential thrust (in the polar frame, phi = 0; in the inertial frame, psi = theta + pi/2)
  // The two formulations parametrize the steering differently (local and inertial angle), so on a
  // finite mesh they are different discrete problems: they converge to the same optimum as O(h^2)
  auto zero = [] (Real) {return Vector1(0.0);};
  auto tangential = [] (Real t) {return Vector1(t + M_PI/2.0);};
  Solve<4, 0, OrbitRaisingPolar> polar_coarse(std::make_shared<OrbitRaisingPolar>(), uniform(T, 50), 2, polar_start, zero);
  Solve<4, 0, OrbitRaisingPolar> polar(std::make_shared<OrbitRaisingPolar>(), uniform(T, 100), 2, polar_start, zero);
  Solve<4, 0, OrbitRaisingCartesian> cartesian_coarse(std::make_shared<OrbitRaisingCartesian>(), uniform(T, 50), 2, cartesian_start, tangential);
  Solve<4, 0, OrbitRaisingCartesian> cartesian(std::make_shared<OrbitRaisingCartesian>(), uniform(T, 100), 2, cartesian_start, tangential);
  ASSERT_TRUE(polar_coarse.succeeded && polar.succeeded && cartesian_coarse.succeeded && cartesian.succeeded);

  Real const r_polar{polar.x_end()(0)}, r_cartesian{cartesian.x_end().head<2>().norm()};
  Real const r_polar_coarse{polar_coarse.x_end()(0)}, r_cartesian_coarse{cartesian_coarse.x_end().head<2>().norm()};
  EXPECT_LT(std::abs(r_cartesian - r_polar), std::abs(r_cartesian_coarse - r_polar_coarse));

  // Richardson extrapolation (error O(h^2)): both formulations converge to the same radius, which
  // matches the reference to its three digits
  Real const r_polar_limit{r_polar + (r_polar - r_polar_coarse)/3.0};
  Real const r_cartesian_limit{r_cartesian + (r_cartesian - r_cartesian_coarse)/3.0};
  EXPECT_NEAR(r_polar_limit, r_cartesian_limit, 1.0e-5);
  EXPECT_NEAR(r_polar_limit, OrbitRaisingData::OPTIMAL_RADIUS, 5.0e-4);

  // The final orbit is circular
  Eigen::Vector4d const xf{cartesian.x_end()};
  EXPECT_LT(xf.head<2>().norm(), 2.0); // sanity
  EXPECT_NEAR(xf.head<2>().dot(xf.tail<2>()), 0.0, 1.0e-10);
  EXPECT_NEAR(xf.tail<2>().squaredNorm(), 1.0/r_cartesian, 1.0e-10);
}

/*\
 |  Cart-pole swing-up with an active track limit: Cartesian DAE against minimal coordinates
\*/

// Conversion between the minimal coordinates (p, theta, p', theta') and the Cartesian states
// (p, x, y, p', vx, vy, lambda) of the cart-pole.
Eigen::Vector<Real, 7> cart_pole_cartesian(Eigen::Vector4d const & q)
{
  using P = CartPoleParameters;
  Real const p{q(0)}, th{q(1)}, v{q(2)}, w{q(3)};
  Eigen::Vector<Real, 7> s;
  s << p, p + P::LENGTH*std::sin(th), -P::LENGTH*std::cos(th), v, v + P::LENGTH*w*std::cos(th), P::LENGTH*w*std::sin(th),
       P::POLE*(P::GRAVITY*std::cos(th) + P::LENGTH*w*w)/(2.0*P::LENGTH);
  return s;
}

Eigen::Vector4d cart_pole_angles(Eigen::Vector<Real, 7> const & s, Real const theta_previous)
{
  using P = CartPoleParameters;
  using S = CartPoleCartesian;
  Real const dx{s(S::X) - s(S::P)}, dvx{s(S::VX) - s(S::VP)};
  Real theta{std::atan2(dx, -s(S::Y))};
  theta += 2.0*M_PI*std::round((theta_previous - theta)/(2.0*M_PI)); // continuous angle
  return Eigen::Vector4d(s(S::P), theta, s(S::VP), (-s(S::Y)*dvx + dx*s(S::VY))/(P::LENGTH*P::LENGTH));
}

TEST(HardOptimalControl, CartPoleSwingUp) {
  using P = CartPoleParameters;
  using S = CartPoleCartesian;
  Real const T{P::DURATION};
  Integer const n{40};
  auto theta = [T] (Real t) {return M_PI*(1.0 - std::cos(M_PI*t/T))/2.0;};
  auto omega = [T] (Real t) {return M_PI*M_PI/(2.0*T)*std::sin(M_PI*t/T);};
  auto zero = [] (Real) {return Vector1(0.0);};
  auto swing_up = [&] (Real t) {return Eigen::Vector4d(0.0, theta(t), 0.0, omega(t));};

  // Both formulations from the same swing-up guess: the problem is nonconvex, and they end in
  // different local optima
  Solve<4, 0, CartPoleAngles> angles(std::make_shared<CartPoleAngles>(), uniform(T, n), 2, Eigen::Vector4d::Zero(), zero, swing_up);
  Solve<7, 2, CartPoleCartesian> cartesian(std::make_shared<CartPoleCartesian>(), uniform(T, n), 2, Eigen::Vector<Real, 7>::Zero(), zero,
    [&] (Real t) {return cart_pole_cartesian(swing_up(t));});
  ASSERT_TRUE(angles.succeeded && cartesian.succeeded);

  // Warm-started from each other's solution, each formulation reproduces the other's optimum
  MultipleShooting<Real, 3, 7, 1, 2> cartesian_warm(std::make_shared<CartPoleCartesian>(), RadauIIA5Tableau<Real>(), n, 2);
  VectorX z(VectorX::Zero(cartesian_warm.variables_number()));
  for (Integer k{0}; k <= n; ++k) {z.segment<7>(cartesian_warm.x_index(k)) = cart_pole_cartesian(angles.x(k));}
  for (Integer k{0}; k < n; ++k) {z(cartesian_warm.u_index(k)) = angles.nlp.u(angles.nlp.guess(), k)(0);}
  z.segment<7>(cartesian_warm.X0_index()) = z.segment<7>(cartesian_warm.x_index(0));
  cartesian_warm.guess(z);
  Real J_cartesian_warm;
  ASSERT_EQ(ipopt_solve(cartesian_warm, options(), J_cartesian_warm), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_cartesian_warm, angles.J, 1.0e-5*angles.J);

  MultipleShooting<Real, 3, 4, 1, 0> angles_warm(std::make_shared<CartPoleAngles>(), RadauIIA5Tableau<Real>(), n, 2);
  VectorX y(VectorX::Zero(angles_warm.variables_number()));
  Real theta_previous{0.0};
  for (Integer k{0}; k <= n; ++k) {
    Eigen::Vector4d const q{cart_pole_angles(cartesian.x(k), theta_previous)};
    theta_previous = q(1);
    y.segment<4>(angles_warm.x_index(k)) = q;
  }
  for (Integer k{0}; k < n; ++k) {y(angles_warm.u_index(k)) = cartesian.nlp.u(cartesian.nlp.guess(), k)(0);}
  angles_warm.guess(y);
  Real J_angles_warm;
  ASSERT_EQ(ipopt_solve(angles_warm, options(), J_angles_warm), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J_angles_warm, cartesian.J, 1.0e-5*cartesian.J);

  // The track limit is active, the pole ends upright, and the rod constraints hold on every node
  Real p_max{0.0};
  for (Integer k{0}; k <= n; ++k) {
    Eigen::Vector<Real, 7> const s{cartesian.x(k)};
    p_max = std::max(p_max, std::abs(s(S::P)));
    EXPECT_LT(cartesian.problem->system()->h(s, 0.0).lpNorm<Eigen::Infinity>(), 1.0e-12);
  }
  EXPECT_NEAR(p_max, P::TRACK, 1.0e-6);
  EXPECT_NEAR(cartesian.x_end()(S::Y), P::LENGTH, 1.0e-12);
}

/*\
 |  Hypersensitive problem: boundary layers on a mesh clustered at the ends
\*/

TEST(HardOptimalControl, HypersensitiveBoundaryLayers) {
  Real const T{40.0};
  auto clustered = [T] (Integer const n) {
    VectorX s{VectorX::LinSpaced(n + 1, 0.0, 1.0)};
    return VectorX((T/2.0)*(1.0 - (M_PI*s.array()).cos()));
  };
  auto zero = [] (Real) {return Vector1(0.0);};
  auto line = [T] (Real t) {return Vector1(1.5 - 0.5*t/T);};
  Solve<1, 0, Hypersensitive> coarse(std::make_shared<Hypersensitive>(T), clustered(40), 2, Vector1(1.5), zero, line);
  Solve<1, 0, Hypersensitive> fine(std::make_shared<Hypersensitive>(T), clustered(80), 2, Vector1(1.5), zero, line);
  ASSERT_TRUE(coarse.succeeded && fine.succeeded);

  // Piecewise-constant controls converge to the asymptotic optimum as O(h^2)
  Real const J_star{Hypersensitive::asymptotic_cost()};
  Real const error_coarse{std::abs(coarse.J - J_star)}, error_fine{std::abs(fine.J - J_star)};
  EXPECT_LT(error_fine, 5.0e-3*J_star);
  EXPECT_GT(error_coarse/error_fine, 3.5);

  // The state rests at zero between the layers
  Integer const middle{40};
  EXPECT_LT(std::abs(fine.x(middle)(0)), 1.0e-6);
}

/*\
 |  Brachistochrone: free final time, energy invariant, singular start at rest
\*/

TEST(HardOptimalControl, Brachistochrone) {
  auto problem = std::make_shared<Brachistochrone>();
  Eigen::Vector4d const start(0.0, 0.0, 0.0, 0.6);
  Solve<4, 1, Brachistochrone> solve(problem, uniform(1.0, 40), 2, start, [] (Real) {return Vector1(M_PI/4.0);});
  ASSERT_TRUE(solve.succeeded);
  EXPECT_NEAR(solve.J, Brachistochrone::optimal_time(), 1.0e-4);
  for (Integer k{0}; k <= 40; ++k) {
    EXPECT_LT(std::abs(problem->system()->h(solve.x(k), 0.0)(0)), 1.0e-12) << "node " << k;
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
