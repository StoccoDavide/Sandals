/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Multiple shooting transcription: sizes, sparsity patterns, and first and second derivatives (through
// the Runge-Kutta substeps and the node least squares) against finite differences.

#include <set>

#include <gtest/gtest.h>

#include "Sandals/MultipleShooting.hh"
#include "Sandals/RungeKutta/RadauIIA3.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/SSPIRK33.hh"
#include "Sandals/RungeKutta/RK4.hh"

#include "OptimalControlProblems.hh"
#include "PendulumProblems.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

// Central finite differences of a vector function.
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

Real relative_error(MatrixX const & a, MatrixX const & b) {return (a - b).norm() / std::max(Real(1.0), b.norm());}

template <typename NLP>
MatrixX dense_jacobian(NLP const & nlp, VectorX const & z)
{
  MatrixX J(MatrixX::Zero(nlp.constraints_number(), nlp.variables_number()));
  VectorX const values{nlp.jacobian_values(z)};
  for (Integer e{0}; e < nlp.jacobian_nnz(); ++e) {J(nlp.jacobian_rows()[e], nlp.jacobian_cols()[e]) += values(e);}
  return J;
}

template <typename NLP>
MatrixX dense_hessian(NLP const & nlp, VectorX const & z, Real const sigma, VectorX const & lambda)
{
  MatrixX H(MatrixX::Zero(nlp.variables_number(), nlp.variables_number()));
  VectorX const values{nlp.hessian_values(z, sigma, lambda)};
  for (Integer e{0}; e < nlp.hessian_nnz(); ++e) {
    Integer const r{nlp.hessian_rows()[e]}, c{nlp.hessian_cols()[e]};
    H(r, c) += values(e);
    if (r != c) {H(c, r) += values(e);}
  }
  return H;
}

// Check the transcription of a problem at a point near its simulated starting point.
template <Integer S, Integer N, Integer U, Integer M, typename Problem>
void check_transcription(std::shared_ptr<Problem> problem, Tableau<Real, S> const & tableau, Integer const intervals,
  Integer const substeps, Real const epsilon, VectorX const & x_ini, Real const tolerance)
{
  std::shared_ptr<OptimalControlProblem<Real, N, U, M>> const ocp{problem};
  SCOPED_TRACE(ocp->name() + " with " + tableau.name + ", substeps " + std::to_string(substeps) +
    ", epsilon " + std::to_string(epsilon));
  MultipleShooting<Real, S, N, U, M> nlp(ocp, tableau, intervals, substeps, epsilon);
  Integer const nz{nlp.variables_number()}, nc{nlp.constraints_number()};

  // Sizes
  Integer const nodes{intervals*(N + U) + N};
  EXPECT_EQ(nz, M > 0 ? nodes + N + (intervals + 1)*M : nodes);
  EXPECT_EQ(nc, intervals*N + (M > 0 ? N + (intervals + 1)*M : 0) + ocp->boundary_conditions_number() +
    (intervals + 1)*ocp->path_constraints_number());

  // Patterns: in range, without duplicates, lower triangle for the Hessian
  std::set<std::pair<Integer, Integer>> jacobian_pattern, hessian_pattern;
  for (Integer e{0}; e < nlp.jacobian_nnz(); ++e) {
    Integer const r{nlp.jacobian_rows()[e]}, c{nlp.jacobian_cols()[e]};
    ASSERT_TRUE(r >= 0 && r < nc && c >= 0 && c < nz);
    EXPECT_TRUE(jacobian_pattern.emplace(r, c).second);
  }
  for (Integer e{0}; e < nlp.hessian_nnz(); ++e) {
    Integer const r{nlp.hessian_rows()[e]}, c{nlp.hessian_cols()[e]};
    ASSERT_TRUE(r >= c && r < nz && c >= 0);
    EXPECT_TRUE(hessian_pattern.emplace(r, c).second);
  }

  // Point: the simulated starting point, perturbed, with random multipliers
  ASSERT_TRUE(nlp.simulate_guess(x_ini, [] (Real t) {return Eigen::Vector<Real, U>::Constant(0.3*std::sin(2.0*t));}));
  std::srand(17);
  VectorX const z{nlp.guess() + 1.0e-2*VectorX::Random(nz)};
  VectorX const lambda{VectorX::Random(nc)};
  Real const sigma{0.8};

  // First and second derivatives against finite differences
  EXPECT_LT(relative_error(nlp.objective_gradient(z).transpose(),
    fd([&] (VectorX const & w) {return VectorX::Constant(1, nlp.objective(w));}, z)), tolerance);
  MatrixX const J{dense_jacobian(nlp, z)};
  EXPECT_LT(relative_error(J, fd([&] (VectorX const & w) {return nlp.constraints(w);}, z)), tolerance);
  MatrixX const H{dense_hessian(nlp, z, sigma, lambda)};
  MatrixX const H_fd{fd([&] (VectorX const & w) {
    return VectorX(sigma*nlp.objective_gradient(w) + dense_jacobian(nlp, w).transpose()*lambda);}, z)};
  EXPECT_LT(relative_error(H, H_fd), tolerance);

  // Every nonzero is in the patterns
  for (Integer r{0}; r < nc; ++r) {
    for (Integer c{0}; c < nz; ++c) {
      if (std::abs(J(r, c)) > 1.0e-12) {EXPECT_TRUE(jacobian_pattern.count({r, c})) << "Jacobian (" << r << ", " << c << ")";}
    }
  }
  for (Integer r{0}; r < nz; ++r) {
    for (Integer c{0}; c <= r; ++c) {
      if (std::abs(H_fd(r, c)) > 1.0e-6) {EXPECT_TRUE(hessian_pattern.count({r, c})) << "Hessian (" << r << ", " << c << ")";}
    }
  }

  // The equality constraints are independent (no algebraic gauges left free, no redundant rows);
  // inactive inequalities may be dependent without harming the constraint qualification
  VectorX z_l, z_u, c_l, c_u;
  nlp.bounds(z_l, z_u, c_l, c_u);
  std::vector<Integer> equalities;
  for (Integer r{0}; r < nc; ++r) {if (c_l(r) == c_u(r)) {equalities.push_back(r);}}
  Eigen::FullPivLU<MatrixX> lu(J(equalities, Eigen::all).transpose());
  lu.setThreshold(1.0e-10);
  EXPECT_EQ(lu.rank(), static_cast<Integer>(equalities.size()));

  // The simulated starting point satisfies the continuity and invariants rows
  VectorX const c0{nlp.constraints(nlp.guess())};
  EXPECT_LT(c0.head(nlp.boundary_row()).template lpNorm<Eigen::Infinity>(), 1.0e-11);
}

/*\
 |  Problems without invariants
\*/

TEST(MultipleShooting, VanDerPolDerivatives) {
  auto ocp = std::make_shared<VanDerPolOCP>(1.0);
  VectorX const x0{Eigen::Vector<Real, 2>(1.0, 0.0)};
  check_transcription<3, 2, 1, 0>(ocp, RadauIIA5Tableau<Real>(), 3, 1, 0.0, x0, 1.0e-6);
  check_transcription<3, 2, 1, 0>(ocp, RadauIIA5Tableau<Real>(), 3, 3, 0.0, x0, 1.0e-6);
  check_transcription<4, 2, 1, 0>(ocp, RK4Tableau<Real>(), 4, 2, 0.0, x0, 1.0e-6);
  check_transcription<3, 2, 1, 0>(ocp, SSPIRK33Tableau<Real>(), 3, 2, 0.0, x0, 1.0e-6);
}

TEST(MultipleShooting, Index1DAEDerivatives) {
  // Algebraic state z = x^2: its initial value is fixed by the boundary conditions
  auto ocp = std::make_shared<SemiLinearDAEOCP>();
  check_transcription<3, 2, 1, 0>(ocp, RadauIIA5Tableau<Real>(), 4, 2, 0.0, Eigen::Vector<Real, 2>(0.5, 0.25), 1.0e-6);
  check_transcription<2, 2, 1, 0>(ocp, RadauIIA3Tableau<Real>(), 4, 1, 0.0, Eigen::Vector<Real, 2>(0.5, 0.25), 1.0e-6);
}

/*\
 |  Index-reduced pendulums with invariants, for different weights of the invariants
\*/

TEST(MultipleShooting, InvariantsDerivatives) {
  PendulumParameters const p;
  for (Real epsilon : {0.0, 1.0e-3, 1.0}) {
    check_transcription<3, 5, 1, 1>(std::make_shared<PendulumSwingUpOCP<PendulumIndex2<1>>>(2.0),
      RadauIIA5Tableau<Real>(), 3, 2, epsilon, p.hanging_at_rest<5>(), 1.0e-6);
    check_transcription<3, 5, 1, 2>(std::make_shared<PendulumSwingUpOCP<PendulumIndex1<2>>>(2.0),
      RadauIIA5Tableau<Real>(), 3, 2, epsilon, p.hanging_at_rest<5>(), 1.0e-6);
  }
  // Index 0 (finite differences Hessians of the dynamics, looser tolerance)
  check_transcription<4, 4, 1, 2>(std::make_shared<PendulumSwingUpOCP<PendulumIndex0<2>>>(2.0),
    RK4Tableau<Real>(), 3, 3, 0.0, p.hanging_at_rest<4>(), 1.0e-5);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
