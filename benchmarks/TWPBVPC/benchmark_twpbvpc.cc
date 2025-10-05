/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>

#include "Sandals.hh"
#include "Sandals/RungeKutta/SSPRK104.hh"

#include "TWPBVPC.hh"

using namespace Sandals;
using Real = double;

// Function to solve the problem
template<typename ProblemType, typename SystemType>
void test(bool const reverse) {

  // Set problem and options
  ProblemType problem;
  problem.verbose_mode(false);
  problem.integrator()->verbose_mode(false);
  problem.integrator()->reverse_mode(reverse);
  problem.tolerance(1.0e-8);
  problem.max_iterations(100);
  problem.subintervals(4);

  // Set time mesh
  constexpr Integer num_points{100};
  Eigen::Vector<Real, Eigen::Dynamic> time(Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
    num_points, problem.time_start(), problem.time_end()
  ));

  // Set initial guess
  Eigen::Matrix<Real, 2, Eigen::Dynamic> guess(problem.guess(time));

  // Solve the problem with different lambda values
  std::vector<Real> lambda_vec{1.0e-0, 1.0e-1, 1.0e-2, 1.0e-3};
  for (Real lambda : lambda_vec) {
    static_cast<SystemType*>(problem.integrator()->system())->lambda(lambda);
    EXPECT_TRUE(problem.multiple_shooting(time, guess));
  }

    // Get the numerical and analytical solutions
    auto sol_n{problem.solution()};
    auto sol_a{problem.exact_solution(reverse ? sol_n.t.reverse().eval() : sol_n.t)};

    // Check the solutions
    EXPECT_NEAR((sol_n.eigen_x(0) - sol_a).array().abs().maxCoeff(), 0.0, 0.5);
}

#ifndef GENERATE_TEST
#define GENERATE_TEST(PROBLEM, INTEGRATOR) \
TEST(PROBLEM##Forward, INTEGRATOR##Explicit) \
{test<PROBLEM##Problem<Real, PROBLEM##Explicit<Real>, INTEGRATOR<Real, 2>>, PROBLEM##Explicit<Real>>(false);} \
TEST(PROBLEM##Forward, INTEGRATOR##Implicit) \
{test<PROBLEM##Problem<Real, PROBLEM##Implicit<Real>, INTEGRATOR<Real, 2>>, PROBLEM##Implicit<Real>>(false);} \
TEST(PROBLEM##Forward, INTEGRATOR##SemiExplicit) \
{test<PROBLEM##Problem<Real, PROBLEM##SemiExplicit<Real>, INTEGRATOR<Real, 2>>, PROBLEM##SemiExplicit<Real>>(false);} \
TEST(PROBLEM##Reverse, INTEGRATOR##Explicit) \
{test<PROBLEM##Problem<Real, PROBLEM##Explicit<Real>, INTEGRATOR<Real, 2>>, PROBLEM##Explicit<Real>>(true);} \
TEST(PROBLEM##Reverse, INTEGRATOR##Implicit) \
{test<PROBLEM##Problem<Real, PROBLEM##Implicit<Real>, INTEGRATOR<Real, 2>>, PROBLEM##Implicit<Real>>(true);} \
TEST(PROBLEM##Reverse, INTEGRATOR##SemiExplicit) \
{test<PROBLEM##Problem<Real, PROBLEM##SemiExplicit<Real>, INTEGRATOR<Real, 2>>, PROBLEM##SemiExplicit<Real>>(true);}
#endif

// Generate the tests
GENERATE_TEST(BVPT1, SSPRK104)
GENERATE_TEST(BVPT2, SSPRK104)
GENERATE_TEST(BVPT3, SSPRK104)
GENERATE_TEST(BVPT4, SSPRK104)
// GENERATE_TEST(BVPT5, SSPRK104)
// GENERATE_TEST(BVPT6, SSPRK104)
// GENERATE_TEST(BVPT7, SSPRK104)
// GENERATE_TEST(BVPT8, SSPRK104)
// GENERATE_TEST(BVPT9, SSPRK104)
// GENERATE_TEST(BVPT10, SSPRK104)
// GENERATE_TEST(BVPT11, SSPRK104)
// GENERATE_TEST(BVPT12, SSPRK104)
// GENERATE_TEST(BVPT13, SSPRK104)
// GENERATE_TEST(BVPT14, SSPRK104)
// GENERATE_TEST(BVPT15, SSPRK104)
// GENERATE_TEST(BVPT16, SSPRK104)
// GENERATE_TEST(BVPT17, SSPRK104)
// GENERATE_TEST(BVPT18, SSPRK104)
// GENERATE_TEST(BVPT19, SSPRK104)
// GENERATE_TEST(BVPT20, SSPRK104)
// GENERATE_TEST(BVPT21, SSPRK104)
// GENERATE_TEST(BVPT22, SSPRK104)
// GENERATE_TEST(BVPT23, SSPRK104)
// GENERATE_TEST(BVPT24, SSPRK104)
// GENERATE_TEST(BVPT25, SSPRK104)
// GENERATE_TEST(BVPT26, SSPRK104)
// GENERATE_TEST(BVPT27, SSPRK104)
// GENERATE_TEST(BVPT28, SSPRK104)
// GENERATE_TEST(BVPT29, SSPRK104)
// GENERATE_TEST(BVPT30, SSPRK104)
// GENERATE_TEST(BVPT31, SSPRK104)
// GENERATE_TEST(BVPT32, SSPRK104)
// GENERATE_TEST(BVPT33, SSPRK104)

// Run all the tests
int main(int argc, char ** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
