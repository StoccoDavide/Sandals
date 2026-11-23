/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <memory>

#include <gtest/gtest.h>

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/RungeKutta/RK4.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"

#include "BasicExplicit.hh"
#include "BasicImplicit.hh"
#include "BasicSemiExplicit.hh"
#include "SimplePendulumExplicit.hh"

using namespace Sandals;

using Real = double;
using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using Vector2 = Eigen::Vector<Real, 2>;

/*\
 |  Linear BVP x'' = 1, x(0) = 0, x'(1) = 0, with exact solution x = t^2/2 - t (RK4 is exact here)
\*/

template <typename Problem>
void check_linear_bvp(Problem & problem, bool const multiple)
{
  VectorX const t_mesh(VectorX::LinSpaced(11, 0.0, 1.0));
  Vector2 const guess(10.0, -10.0);
  ASSERT_TRUE(multiple ? problem.multiple_shooting(t_mesh, guess) : problem.single_shooting(t_mesh, guess));
  auto sol = problem.solution();
  ASSERT_EQ(sol->size(), t_mesh.size());
  for (Integer i{0}; i < sol->size(); ++i) {
    EXPECT_NEAR(sol->t(i), t_mesh(i), 1.0e-14);
    EXPECT_LT((sol->x.col(i) - problem.analytical_solution(sol->t(i))).norm(), 1.0e-10) << "at t = " << sol->t(i);
  }
}

TEST(Shooting, LinearSingle) {
  BasicExplicitProblem pe(std::make_shared<RK4<Real, 2, 0>>());
  BasicImplicitProblem pi(std::make_shared<RK4<Real, 2, 0>>());
  BasicSemiExplicitProblem ps(std::make_shared<RK4<Real, 2, 0>>());
  check_linear_bvp(pe, false);
  check_linear_bvp(pi, false);
  check_linear_bvp(ps, false);
}

TEST(Shooting, LinearMultiple) {
  BasicExplicitProblem pe(std::make_shared<RK4<Real, 2, 0>>());
  BasicImplicitProblem pi(std::make_shared<RK4<Real, 2, 0>>());
  BasicSemiExplicitProblem ps(std::make_shared<RK4<Real, 2, 0>>());
  check_linear_bvp(pe, true);
  check_linear_bvp(pi, true);
  check_linear_bvp(ps, true);
}

/*\
 |  Nonlinear BVP: pendulum theta'' = -g/l sin(theta), theta(0) = 1, theta(T) = -0.5
\*/

template <typename Integrator>
class PendulumBVP : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;

  PendulumBVP(std::shared_ptr<Integrator> rk)
    : Problem<Real, 2, 0, Integrator>("PendulumBVP", std::make_shared<SimplePendulumExplicit<Real>>(), rk) {}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {return VectorF(x_ini(0) - 1.0, x_end(0) + 0.5);}

  MatrixJF Jb_x_ini(VectorF const &, VectorF const &) const override
  {MatrixJF J(MatrixJF::Zero()); J(0, 0) = 1.0; return J;}

  MatrixJF Jb_x_end(VectorF const &, VectorF const &) const override
  {MatrixJF J(MatrixJF::Zero()); J(1, 0) = 1.0; return J;}
};

template <typename Integrator>
void check_pendulum_bvp(bool const multiple)
{
  PendulumBVP<Integrator> problem(std::make_shared<Integrator>());
  problem.tolerance(1.0e-12);
  problem.max_iterations(20);
  VectorX const t_mesh(VectorX::LinSpaced(101, 0.0, 0.8));
  Vector2 const guess(1.0, 0.0);

  ASSERT_TRUE(multiple ? problem.multiple_shooting(t_mesh, guess) : problem.single_shooting(t_mesh, guess));
  auto sol = problem.solution();
  Vector2 const x_ini(sol->x.col(0)), x_end(sol->x.col(sol->size() - 1));
  EXPECT_NEAR(x_ini(0), 1.0, 1.0e-9);
  EXPECT_NEAR(x_end(0), -0.5, 1.0e-9);

  // Re-integrating from the found initial state must reproduce the trajectory (continuity)
  Integrator rk(std::make_shared<SimplePendulumExplicit<Real>>());
  Solution<Real, 2, 0> check;
  ASSERT_TRUE(rk.solve(t_mesh, x_ini, check));
  EXPECT_LT((check.x - sol->x).norm(), 1.0e-8);
}

TEST(Shooting, PendulumSingle) {
  check_pendulum_bvp<RK4<Real, 2, 0>>(false);
  check_pendulum_bvp<RadauIIA5<Real, 2, 0>>(false);
}

TEST(Shooting, PendulumMultiple) {
  check_pendulum_bvp<RK4<Real, 2, 0>>(true);
  check_pendulum_bvp<RadauIIA5<Real, 2, 0>>(true);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
