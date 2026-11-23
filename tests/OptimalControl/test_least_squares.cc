/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Least-squares multiple shooting for DAE boundary value problems, solved with IPOPT (instead of a
// Gauss-Newton method), with the invariants enforced either in the least-squares sense or by
// projection.

#include <memory>

#include <gtest/gtest.h>

#include "Sandals/Ipopt.hh"
#include "Sandals/RungeKutta/RK4.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"

#include "MinimumEnergyProblems.hh"
#include "OscillatorExplicit.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

using Integrator8 = RK4<Real, 8, 2>;
using Shooting8 = LeastSquaresShooting<Real, 8, 2, Integrator8>;

// Quiet IPOPT with the default tight tolerance (1e-12)
IpoptOptions quiet_options()
{
  IpoptOptions options;
  options.integers["print_level"] = std::getenv("SANDALS_IPOPT_DEBUG") ? 5 : 0;
  options.strings["sb"] = "yes";
  return options;
}

// Dense constraints Jacobian from the sparse values (for testing).
template <typename NLP>
MatrixX dense_jacobian(NLP const & nlp, VectorX const & z)
{
  MatrixX J(MatrixX::Zero(nlp.constraints_number(), nlp.variables_number()));
  VectorX const values{nlp.jacobian_values(z)};
  for (Integer e{0}; e < nlp.jacobian_nnz(); ++e) {J(nlp.jacobian_rows()[e], nlp.jacobian_cols()[e]) += values(e);}
  return J;
}

// Central finite differences of the constraints.
template <typename NLP>
MatrixX finite_differences_jacobian(NLP const & nlp, VectorX const & z)
{
  Real const delta{1.0e-6};
  MatrixX J(nlp.constraints_number(), nlp.variables_number());
  for (Integer j{0}; j < z.size(); ++j) {
    VectorX z_p(z), z_m(z);
    z_p(j) += delta;
    z_m(j) -= delta;
    J.col(j) = (nlp.constraints(z_p) - nlp.constraints(z_m)) / (2.0*delta);
  }
  return J;
}

/*\
 |  Inconsistent overdetermined system: harmonic oscillator with its energy invariant
\*/

// x0' = x1, x1' = -x0 on [0, 2], with the energy invariant x0^2 + x1^2 = 1, exact solution
// (cos t, -sin t). RK4 does not preserve the energy, so with the invariants as rows the residual
// system is overdetermined and inconsistent (by the local error), and is solved in the least-squares
// sense. The boundary conditions are x1(2) = -sin(2) and either
//  - x0(0) = 1, which leaves x1(0) free: fine with the invariants as rows (h(x_0) = 0 is one of them),
//    but not with the projection, which never acts on the first node;
//  - x0(0)^2 + x1(0)^2 = 1, a consistent initial state, as required by the projection.
// With the second choice the problem has two solutions, (cos(t + p), -sin(t + p)) with p = 0 or
// p = pi - 4, so the initial guess (the exact solution perturbed by 10%) selects p = 0.
using Integrator2 = RK4<Real, 2, 1>;

class OscillatorBVP : public Problem<Real, 2, 1, Integrator2>
{
  bool m_consistent_start;

public:
  OscillatorBVP(bool const consistent_start = false)
    : Problem<Real, 2, 1, Integrator2>("OscillatorBVP", std::make_shared<OscillatorExplicit<Real>>(),
      std::make_shared<Integrator2>()), m_consistent_start(consistent_start) {}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    Real const start{this->m_consistent_start ? x_ini.squaredNorm() - 1.0 : x_ini(0) - 1.0};
    return VectorF(start, x_end(1) + std::sin(2.0));
  }

  MatrixJF Jb_x_ini(VectorF const & x_ini, VectorF const &) const override
  {
    MatrixJF J(MatrixJF::Zero());
    if (this->m_consistent_start) {J.row(0) = 2.0*x_ini.transpose();} else {J(0, 0) = 1.0;}
    return J;
  }

  MatrixJF Jb_x_end(VectorF const &, VectorF const &) const override
  {
    MatrixJF J(MatrixJF::Zero());
    J(1, 1) = 1.0;
    return J;
  }
};

using Shooting2 = LeastSquaresShooting<Real, 2, 1, Integrator2>;

/*\
 |  Nonlinear program structure and derivatives
\*/

TEST(LeastSquaresShooting, Sizes) {
  VectorX const mesh{VectorX::LinSpaced(6, 0.0, 1.0)}; // m = 5 intervals
  auto least_squares = std::make_shared<MinimumEnergyIndex0BVP<Integrator8>>(std::make_shared<Integrator8>());
  Shooting8 nlp_ls(least_squares, mesh, 2, Shooting8::Invariants::LEAST_SQUARES);
  EXPECT_EQ(nlp_ls.residuals_number(), 5*8 + 6*2 + 8); // continuity, invariants on every node, conditions
  EXPECT_EQ(nlp_ls.variables_number(), 6*8 + nlp_ls.residuals_number());

  auto projection = std::make_shared<MinimumEnergyIndex0BVP<Integrator8>>(std::make_shared<Integrator8>());
  Shooting8 nlp_pr(projection, mesh, 2, Shooting8::Invariants::PROJECTION);
  EXPECT_EQ(nlp_pr.residuals_number(), 5*8 + 8); // square: no invariants rows
  EXPECT_TRUE(projection->integrator()->projection());
  EXPECT_FALSE(least_squares->integrator()->projection());
}

TEST(LeastSquaresShooting, JacobianMatchesFiniteDifferences) {
  // The Jacobian chains the integration (and projection) sensitivities over the substeps
  VectorX const mesh{VectorX::LinSpaced(4, 0.0, 1.0)};
  for (auto invariants : {Shooting8::Invariants::LEAST_SQUARES, Shooting8::Invariants::PROJECTION}) {
    auto problem = std::make_shared<MinimumEnergyIndex0BVP<Integrator8>>(std::make_shared<Integrator8>());
    Shooting8 nlp(problem, mesh, 3, invariants);
    std::srand(11);
    VectorX const z{VectorX::Random(nlp.variables_number())};
    MatrixX const J{dense_jacobian(nlp, z)};
    MatrixX const J_fd{finite_differences_jacobian(nlp, z)};
    EXPECT_LT((J - J_fd).norm() / J_fd.norm(), 1.0e-7);
  }

  // A nonlinear invariant (oscillator energy)
  for (auto invariants : {Shooting2::Invariants::LEAST_SQUARES, Shooting2::Invariants::PROJECTION}) {
    Shooting2 nlp(std::make_shared<OscillatorBVP>(), VectorX::LinSpaced(4, 0.0, 1.5), 4, invariants);
    std::srand(13);
    VectorX z{VectorX::Random(nlp.variables_number())};
    z.head(8) += VectorX::Constant(8, 1.0); // keep the states away from the origin
    EXPECT_LT((dense_jacobian(nlp, z) - finite_differences_jacobian(nlp, z)).norm(), 1.0e-7);
  }
}

/*\
 |  Minimum energy problem (Gerdts, Example 1.0.1), solved from the zero initial guess
\*/

// Solve the index-0 boundary value problem with invariants by least-squares shooting with IPOPT.
// Returns the shooting nonlinear program (with the solution as guess).
std::unique_ptr<Shooting8> solve_minimum_energy(Shooting8::Invariants const invariants, Integer const intervals)
{
  auto problem = std::make_shared<MinimumEnergyIndex0BVP<Integrator8>>(std::make_shared<Integrator8>());
  auto nlp = std::make_unique<Shooting8>(problem, VectorX::LinSpaced(intervals + 1, 0.0, 1.0), 1, invariants);
  nlp->guess([] (Real) {return Eigen::Vector<Real, 8>::Zero().eval();});
  Real objective;
  EXPECT_EQ(ipopt_solve(*nlp, quiet_options(), objective), Ipopt::Solve_Succeeded);
  return nlp;
}

void check_minimum_energy_solution(Shooting8 const & nlp)
{
  VectorX const & z{nlp.guess()};
  for (Integer k{0}; k <= nlp.intervals(); ++k) {
    Real const t{nlp.t_mesh()(k)};
    EXPECT_NEAR(MinimumEnergyIndex0::control(nlp.x(z, k)), -2.0, 1.0e-9) << "t = " << t;
    EXPECT_LT((nlp.x(z, k) - MinimumEnergyIndex0::solution(t)).norm(), 1.0e-9) << "t = " << t;
  }
  // The solution is quadratic in time, so RK4 is exact and the overdetermined system is consistent
  EXPECT_LT(nlp.residuals(z).norm(), 1.0e-10);
}

TEST(LeastSquaresShooting, MinimumEnergyLeastSquaresInvariants) {
  auto nlp{solve_minimum_energy(Shooting8::Invariants::LEAST_SQUARES, 20)};
  check_minimum_energy_solution(*nlp);
}

TEST(LeastSquaresShooting, MinimumEnergyProjectedInvariants) {
  auto nlp{solve_minimum_energy(Shooting8::Invariants::PROJECTION, 20)};
  check_minimum_energy_solution(*nlp);
}

TEST(LeastSquaresShooting, MinimumEnergyDirectApproachAgrees) {
  // The direct approach (index-2 DAE, RadauIIA5 collocation) gives the same optimum: J* = 2, u = -2
  auto ocp = std::make_shared<MinimumEnergyIndex2OCP>();
  Collocation<Real, 3, 4, 1> nlp(ocp, RadauIIA5Tableau<Real>(), 20);
  Real J;
  ASSERT_EQ(ipopt_solve(nlp, quiet_options(), J), Ipopt::Solve_Succeeded);
  EXPECT_NEAR(J, 2.0, 1.0e-10);
  EXPECT_LT((nlp.u_mesh(nlp.guess()).array() + 2.0).abs().maxCoeff(), 1.0e-9);
  for (Integer k{1}; k <= 20; ++k) {
    Real const t{nlp.t_mesh()(k)};
    EXPECT_LT((nlp.x(nlp.guess(), k) - MinimumEnergyIndex0::solution(t).head<4>()).norm(), 1.0e-9) << "t = " << t;
  }
}

/*\
 |  Inconsistent overdetermined system (oscillator defined above)
\*/

// Largest error on the nodes against the exact solution, after solving with IPOPT.
Real oscillator_error(Shooting2::Invariants const invariants, Integer const intervals, Real & residual_norm)
{
  bool const consistent_start{invariants == Shooting2::Invariants::PROJECTION};
  Shooting2 nlp(std::make_shared<OscillatorBVP>(consistent_start), VectorX::LinSpaced(intervals + 1, 0.0, 2.0), 1, invariants);
  nlp.guess([] (Real t) {return Eigen::Vector<Real, 2>(1.1*std::cos(t), -1.1*std::sin(t));});
  Real objective;
  EXPECT_EQ(ipopt_solve(nlp, quiet_options(), objective), Ipopt::Solve_Succeeded);
  residual_norm = nlp.residuals(nlp.guess()).norm();
  Real error{0.0};
  for (Integer k{0}; k <= intervals; ++k) {
    Real const t{nlp.t_mesh()(k)};
    error = std::max(error, (nlp.x(nlp.guess(), k) - Eigen::Vector<Real, 2>(std::cos(t), -std::sin(t))).norm());
  }
  return error;
}

TEST(LeastSquaresShooting, InconsistentOverdeterminedSystem) {
  // Least squares: nonzero (but local-error small) residual, fourth-order convergence
  Real residual_coarse, residual_fine;
  Real const error_coarse{oscillator_error(Shooting2::Invariants::LEAST_SQUARES, 10, residual_coarse)};
  Real const error_fine{oscillator_error(Shooting2::Invariants::LEAST_SQUARES, 20, residual_fine)};
  EXPECT_GT(residual_coarse, 1.0e-8);
  EXPECT_LT(residual_coarse, 1.0e-3);
  EXPECT_GT(residual_coarse/residual_fine, 12.0);
  EXPECT_GT(error_coarse/error_fine, 12.0);

  // Projection: square and consistent system, zero residual, same order
  Real projected_coarse, projected_fine;
  Real const projected_error_coarse{oscillator_error(Shooting2::Invariants::PROJECTION, 10, projected_coarse)};
  Real const projected_error_fine{oscillator_error(Shooting2::Invariants::PROJECTION, 20, projected_fine)};
  EXPECT_LT(projected_coarse, 1.0e-10);
  EXPECT_LT(projected_fine, 1.0e-10);
  EXPECT_GT(projected_error_coarse/projected_error_fine, 12.0);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
