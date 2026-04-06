/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>

#include "Sandals.hh"
#include "Sandals/RungeKutta/GaussLegendre4.hh"
#include "TWPBVPC.hh"

using namespace Sandals;
using Real = double;

// Function to solve the problem
template <typename ProblemType, typename SystemType>
void test(const bool reverse) {
  // Set problem and options
  ProblemType problem;
  problem.verbose_mode(false);
  problem.integrator()->verbose_mode(false);
  problem.integrator()->reverse_mode(reverse);
  problem.tolerance(1.0e-8);
  problem.max_iterations(100);
  problem.subintervals(1);

  // Set time mesh
  constexpr Integer num_points{500};
  Eigen::Vector<Real, Eigen::Dynamic> time(
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(num_points,
                                                     problem.time_start(),
                                                     problem.time_end()));

  // Set initial guess
  auto guess(problem.guess(time));

  // Solve the problem with different lambda values
  std::vector<Real> lambda_vec{1.0e-0, 1.0e-1, 1.0e-2, 1.0e-3};
  for (Real lambda : lambda_vec) {
    static_cast<SystemType *>(problem.integrator()->system())->lambda(lambda);
    EXPECT_TRUE(problem.multiple_shooting(time, guess));
  }

  // Check the exact solutions
  if (problem.has_exact_solution()) {
    auto sol_n(problem.solution());
    auto sol_e(
        problem.exact_solution(reverse ? sol_n.t.reverse().eval() : sol_n.t));
    EXPECT_NEAR((sol_n.eigen_x(0) - sol_e).array().abs().maxCoeff(), 0.0, 0.5);
  }
}

#ifndef GENERATE_TEST
#define GENERATE_TEST(PROBLEM, INTEGRATOR)                                   \
  constexpr Integer PROBLEM##D{PROBLEM##Explicit<Real>::equations_number()}; \
  TEST(PROBLEM##Forward, INTEGRATOR##Explicit) {                             \
    test<PROBLEM##Problem<Real,                                              \
                          PROBLEM##Explicit<Real>,                           \
                          INTEGRATOR<Real, PROBLEM##D, 0>>,                  \
         PROBLEM##Explicit<Real>>(false);                                    \
  }                                                                          \
  TEST(PROBLEM##Forward, INTEGRATOR##Implicit) {                             \
    test<PROBLEM##Problem<Real,                                              \
                          PROBLEM##Implicit<Real>,                           \
                          INTEGRATOR<Real, PROBLEM##D, 0>>,                  \
         PROBLEM##Implicit<Real>>(false);                                    \
  }                                                                          \
  TEST(PROBLEM##Forward, INTEGRATOR##SemiExplicit) {                         \
    test<PROBLEM##Problem<Real,                                              \
                          PROBLEM##SemiExplicit<Real>,                       \
                          INTEGRATOR<Real, PROBLEM##D, 0>>,                  \
         PROBLEM##SemiExplicit<Real>>(false);                                \
  }                                                                          \
  TEST(PROBLEM##Reverse, INTEGRATOR##Explicit) {                             \
    test<PROBLEM##Problem<Real,                                              \
                          PROBLEM##Explicit<Real>,                           \
                          INTEGRATOR<Real, PROBLEM##D, 0>>,                  \
         PROBLEM##Explicit<Real>>(true);                                     \
  }                                                                          \
  TEST(PROBLEM##Reverse, INTEGRATOR##Implicit) {                             \
    test<PROBLEM##Problem<Real,                                              \
                          PROBLEM##Implicit<Real>,                           \
                          INTEGRATOR<Real, PROBLEM##D, 0>>,                  \
         PROBLEM##Implicit<Real>>(true);                                     \
  }                                                                          \
  TEST(PROBLEM##Reverse, INTEGRATOR##SemiExplicit) {                         \
    test<PROBLEM##Problem<Real,                                              \
                          PROBLEM##SemiExplicit<Real>,                       \
                          INTEGRATOR<Real, PROBLEM##D, 0>>,                  \
         PROBLEM##SemiExplicit<Real>>(true);                                 \
  }
#endif

// Generate the tests
GENERATE_TEST(BVPT1, GaussLegendre4)
GENERATE_TEST(BVPT2, GaussLegendre4)
GENERATE_TEST(BVPT3, GaussLegendre4)
GENERATE_TEST(BVPT4, GaussLegendre4)
GENERATE_TEST(BVPT5, GaussLegendre4)
GENERATE_TEST(BVPT6, GaussLegendre4)
GENERATE_TEST(BVPT7, GaussLegendre4)
GENERATE_TEST(BVPT8, GaussLegendre4)
GENERATE_TEST(BVPT9, GaussLegendre4)
GENERATE_TEST(BVPT10, GaussLegendre4)
GENERATE_TEST(BVPT11, GaussLegendre4)
GENERATE_TEST(BVPT12, GaussLegendre4)
GENERATE_TEST(BVPT13, GaussLegendre4)
GENERATE_TEST(BVPT14, GaussLegendre4)
GENERATE_TEST(BVPT15, GaussLegendre4)
GENERATE_TEST(BVPT16, GaussLegendre4)
GENERATE_TEST(BVPT17, GaussLegendre4)
GENERATE_TEST(BVPT18, GaussLegendre4)
GENERATE_TEST(BVPT19, GaussLegendre4)
GENERATE_TEST(BVPT20, GaussLegendre4)
GENERATE_TEST(BVPT21, GaussLegendre4)
GENERATE_TEST(BVPT22, GaussLegendre4)
GENERATE_TEST(BVPT23, GaussLegendre4)
GENERATE_TEST(BVPT24, GaussLegendre4)
GENERATE_TEST(BVPT25, GaussLegendre4)
GENERATE_TEST(BVPT26, GaussLegendre4)
GENERATE_TEST(BVPT27, GaussLegendre4)
GENERATE_TEST(BVPT28, GaussLegendre4)
GENERATE_TEST(BVPT29, GaussLegendre4)
GENERATE_TEST(BVPT30, GaussLegendre4)
GENERATE_TEST(BVPT31, GaussLegendre4)
GENERATE_TEST(BVPT32, GaussLegendre4)
GENERATE_TEST(BVPT33, GaussLegendre4)

// Run all the tests
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
