/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <set>

#include <gtest/gtest.h>

#include "Sandals/Collocation.hh"
#include "Sandals/RungeKutta/RadauIIA3.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/GaussLegendre4.hh"
#include "Sandals/RungeKutta/RK4.hh"

#include "OptimalControlProblems.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

// Central finite differences of a vector function (dynamic sizes).
template <typename Function>
MatrixX fd(Function && g, VectorX const & z, Real const delta = 1.0e-6)
{
  VectorX const g0{g(z)};
  MatrixX J(g0.size(), z.size());
  for (Integer j{0}; j < z.size(); ++j) {
    VectorX z_p(z), z_m(z);
    z_p(j) += delta;
    z_m(j) -= delta;
    J.col(j) = (g(z_p) - g(z_m)) / (2.0*delta);
  }
  return J;
}

Real rel_error(MatrixX const & a, MatrixX const & b) {return (a - b).norm() / std::max(Real(1.0), b.norm());}

// Check sizes, sparsity patterns, and first/second derivatives of the transcription at a random point.
template <Integer S, Integer N, Integer U>
void check_transcription(std::shared_ptr<OptimalControlProblem<Real, N, U>> ocp, Tableau<Real, S> const & tableau,
  Integer const intervals, Real const tol)
{
  SCOPED_TRACE(ocp->name() + " with " + tableau.name);
  Collocation<Real, S, N, U> nlp(ocp, tableau, intervals);

  Integer const nz{nlp.variables_number()}, nc{nlp.constraints_number()};
  EXPECT_EQ(nz, intervals*(N + N*S + U) + N);
  EXPECT_EQ(nc, intervals*(N*S + N) + ocp->boundary_conditions_number() + (intervals + 1)*ocp->path_constraints_number());

  // Jacobian pattern: in range and without duplicates
  std::set<std::pair<Integer, Integer>> jac_set;
  for (Integer e{0}; e < nlp.jacobian_nnz(); ++e) {
    Integer const r{nlp.jacobian_rows()[e]}, c{nlp.jacobian_cols()[e]};
    ASSERT_TRUE(r >= 0 && r < nc && c >= 0 && c < nz);
    EXPECT_TRUE(jac_set.emplace(r, c).second) << "duplicate Jacobian entry (" << r << ", " << c << ")";
  }

  // Hessian pattern: lower triangle, in range and without duplicates
  std::set<std::pair<Integer, Integer>> hess_set;
  for (Integer e{0}; e < nlp.hessian_nnz(); ++e) {
    Integer const r{nlp.hessian_rows()[e]}, c{nlp.hessian_cols()[e]};
    ASSERT_TRUE(r >= c && r < nz && c >= 0);
    EXPECT_TRUE(hess_set.emplace(r, c).second) << "duplicate Hessian entry (" << r << ", " << c << ")";
  }

  // Random evaluation point and multipliers
  std::srand(7);
  VectorX const z{VectorX::Random(nz)*0.5};
  VectorX const lambda{VectorX::Random(nc)};
  Real const sigma{0.7};

  // Gradient, Jacobian, and Hessian of the Lagrangian against finite differences
  VectorX const grad{nlp.objective_gradient(z)};
  EXPECT_LT(rel_error(grad.transpose(), fd([&] (VectorX const & w)
    {return VectorX::Constant(1, nlp.objective(w));}, z)), tol);

  MatrixX const J{nlp.jacobian_dense(z)};
  EXPECT_LT(rel_error(J, fd([&] (VectorX const & w) {return nlp.constraints(w);}, z)), tol);

  MatrixX const H{nlp.hessian_dense(z, sigma, lambda)};
  MatrixX const H_fd{fd([&] (VectorX const & w)
    {return VectorX(sigma*nlp.objective_gradient(w) + nlp.jacobian_dense(w).transpose()*lambda);}, z)};
  EXPECT_LT(rel_error(H, H_fd), tol);

  // Every nonzero of the dense derivatives must be in the sparsity patterns
  for (Integer r{0}; r < nc; ++r) {
    for (Integer c{0}; c < nz; ++c) {
      if (J(r, c) != 0.0) {EXPECT_TRUE(jac_set.count({r, c})) << "Jacobian entry out of pattern (" << r << ", " << c << ")";}
    }
  }
  for (Integer r{0}; r < nz; ++r) {
    for (Integer c{0}; c <= r; ++c) {
      if (std::abs(H_fd(r, c)) > 1.0e-6) {EXPECT_TRUE(hess_set.count({r, c})) << "Hessian entry out of pattern (" << r << ", " << c << ")";}
    }
  }
}

TEST(Collocation, VanDerPolDerivatives) {
  auto ocp = std::make_shared<VanDerPolOCP>(1.0);
  check_transcription<2, 2, 1>(ocp, RadauIIA3Tableau<Real>(), 4, 1.0e-6);
  check_transcription<3, 2, 1>(ocp, RadauIIA5Tableau<Real>(), 3, 1.0e-6);
  check_transcription<2, 2, 1>(ocp, GaussLegendre4Tableau<Real>(), 3, 1.0e-6);
  check_transcription<4, 2, 1>(ocp, RK4Tableau<Real>(), 3, 1.0e-6);
}

TEST(Collocation, VanDerPolFiniteDifferencesDefaults) {
  // All derivatives of the problem come from the finite differences defaults (looser tolerance)
  auto ocp = std::make_shared<VanDerPolOCPValuesOnly>(1.0);
  check_transcription<3, 2, 1>(ocp, RadauIIA5Tableau<Real>(), 3, 1.0e-4);
}

TEST(Collocation, DAEDerivatives) {
  auto ocp = std::make_shared<SemiLinearDAEOCP>();
  check_transcription<2, 2, 1>(ocp, RadauIIA3Tableau<Real>(), 4, 1.0e-6);
  check_transcription<3, 2, 1>(ocp, RadauIIA5Tableau<Real>(), 3, 1.0e-6);
}

TEST(Collocation, DefaultsMatchExactDerivatives) {
  // The finite differences defaults of the problem must agree with the exact derivatives
  VanDerPolOCP exact(1.0);
  VanDerPolOCPValuesOnly approx(1.0);
  Eigen::Vector<Real, 2> const x(0.3, -0.7), y(1.1, 0.2);
  Eigen::Vector<Real, 1> const u(0.4);
  VectorX const l{VectorX::Constant(2, 1.3)}, m{VectorX::Constant(1, -0.8)};
  EXPECT_LT(rel_error(approx.mayer_gradient(x, y), exact.mayer_gradient(x, y)), 1.0e-8);
  EXPECT_LT(rel_error(approx.mayer_hessian(x, y), exact.mayer_hessian(x, y)), 1.0e-4);
  EXPECT_LT(rel_error(approx.lagrange_gradient(x, u, 0.0), exact.lagrange_gradient(x, u, 0.0)), 1.0e-8);
  EXPECT_LT(rel_error(approx.lagrange_hessian(x, u, 0.0), exact.lagrange_hessian(x, u, 0.0)), 1.0e-4);
  EXPECT_LT(rel_error(approx.boundary_conditions_jacobian(x, y), exact.boundary_conditions_jacobian(x, y)), 1.0e-8);
  EXPECT_LT(rel_error(approx.boundary_conditions_hessian(x, y, l), exact.boundary_conditions_hessian(x, y, l)), 1.0e-4);
  EXPECT_LT(rel_error(approx.path_constraints_jacobian(x, u, 0.0), exact.path_constraints_jacobian(x, u, 0.0)), 1.0e-8);
  EXPECT_LT(rel_error(approx.path_constraints_hessian(x, u, 0.0, m), exact.path_constraints_hessian(x, u, 0.0, m)), 1.0e-4);

  // Controlled system Hessians defaults against the exact ones
  VanDerPol vdp;
  SemiLinearDAE dae;
  Eigen::Vector<Real, 2> const lam(0.6, -1.4), xd(0.2, 0.9);
  auto const & base = static_cast<Controlled<Real, 2, 1> const &>(vdp);
  EXPECT_LT(rel_error(vdp.Controlled<Real, 2, 1>::HF_x_x(x, xd, u, 0.0, lam), base.HF_x_x(x, xd, u, 0.0, lam)), 1.0e-8);
  EXPECT_LT(rel_error(dae.HF_x_x(x, xd, u, 0.0, lam), MatrixX((Eigen::Matrix<Real, 2, 2>() << -2.0*lam(1), 0.0, 0.0, 0.0).finished())), 1.0e-8);
  EXPECT_LT(dae.HF_x_x_dot(x, xd, u, 0.0, lam).norm(), 1.0e-8);
  EXPECT_LT(dae.HF_u_u(x, xd, u, 0.0, lam).norm(), 1.0e-8);
}

TEST(Collocation, SimulatedGuessSatisfiesDynamics) {
  // The forward simulation with the same Runge-Kutta method must satisfy stage and continuity rows
  auto check = [] (auto ocp, auto const & tableau, auto x0, Integer const n) {
    constexpr Integer S{std::decay_t<decltype(tableau)>::Vector::RowsAtCompileTime};
    using Ocp = typename decltype(ocp)::element_type;
    constexpr Integer N{Ocp::VectorF::RowsAtCompileTime}, U{Ocp::VectorU::RowsAtCompileTime};
    Collocation<Real, S, N, U> nlp(ocp, tableau, n);
    ASSERT_TRUE(nlp.simulate_guess(x0, [] (Real t) {return typename Ocp::VectorU(0.3*std::sin(t));}));
    VectorX const c{nlp.constraints(nlp.guess())};
    // Implicit stages are solved by the integrator Newton up to its tolerance
    EXPECT_LT(c.head(n*(N*S + N)).lpNorm<Eigen::Infinity>(), 1.0e-7) << tableau.name;
    EXPECT_LT((nlp.x(nlp.guess(), 0) - x0).norm(), 1.0e-15);
  };
  check(std::make_shared<VanDerPolOCP>(2.0), RadauIIA5Tableau<Real>(), Eigen::Vector<Real, 2>(1.0, 0.0), 10);
  check(std::make_shared<VanDerPolOCP>(2.0), RK4Tableau<Real>(), Eigen::Vector<Real, 2>(1.0, 0.0), 10);
  check(std::make_shared<VanDerPolOCP>(2.0), GaussLegendre4Tableau<Real>(), Eigen::Vector<Real, 2>(1.0, 0.0), 10);
  check(std::make_shared<SemiLinearDAEOCP>(), RadauIIA5Tableau<Real>(), Eigen::Vector<Real, 2>(0.5, 0.25), 10);
}

TEST(Collocation, Bounds) {
  auto ocp = std::make_shared<VanDerPolOCP>(1.0);
  Collocation<Real, 3, 2, 1> nlp(ocp, RadauIIA5Tableau<Real>(), 4);
  VectorX z_l, z_u, c_l, c_u;
  nlp.bounds(z_l, z_u, c_l, c_u);
  for (Integer k{0}; k < 4; ++k) {
    EXPECT_EQ(z_l(nlp.u_index(k)), -1.5);
    EXPECT_EQ(z_u(nlp.u_index(k)), 1.5);
    EXPECT_TRUE(std::isinf(z_l(nlp.K_index(k, 0))));
  }
  for (Integer k{0}; k <= 4; ++k) {
    EXPECT_TRUE(std::isinf(c_l(nlp.path_row(k))));
    EXPECT_EQ(c_u(nlp.path_row(k)), 4.0);
  }
  EXPECT_EQ(c_l.head(nlp.path_row(0)).norm(), 0.0);
  EXPECT_EQ(c_u.head(nlp.path_row(0)).norm(), 0.0);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
