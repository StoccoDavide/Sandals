/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Astro project is distributed under the GNU GPLv3.                     *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"
#include <gtest/gtest.h>

using namespace Sandals;

// Funzione di Booth.
TEST(Newton, Booth) {
  // Non-linear solver.
  Sandals::Newton<2> nlsolver;
  // Starting entries.
  vec2 x_ini(0.0, 0.0);
  vec2 x_out(QUIET_NAN, QUIET_NAN);
  auto fun = [](vec2 const & X, vec2 & F) {
    F << X(0) + 2.0*X(1) - 7.0, 2.0*X(0) + X(1) - 5.0;
  };
  auto jac = [](vec2 const & /*X*/, mat2 & JF) {
    JF << 1.0, 2.0, 2.0, 1.0;
  };
  // Solve without damping.
  nlsolver.solve(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - vec2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
  // Solve with damping.
  nlsolver.solve_dumped(fun, jac, x_ini, x_out);
  EXPECT_LE((x_out - vec2(1.0, 3.0)).maxCoeff(), EPSILON_LOW);
  EXPECT_TRUE(nlsolver.converged());
}

// 2D Rosenbrock function.
TEST(Newton, Rosenbrock2D) {
  // Non-linear solver.
  Sandals::Newton<2> nlsolver;
  for (real a = 1.0; a <= 100.0; a += 1.0) {
    for (real b = 1.0; b <= 100.0; b += 1.0) {
      // Starting entries.
      vec2 x_ini(0.0, 0.0);
      vec2 x_out(QUIET_NAN, QUIET_NAN);
      auto fun = [a, b](vec2 const & X, vec2 & F) {
        F << a*(1 - X(0)), b*(X(1) - X(0)*X(0));
      };
      auto jac = [a, b](vec2 const & /*X*/, mat2 & JF) {
        JF << -a, real(0.0), real(0.0), b;
      };
      // Solve without damping.
      nlsolver.solve(fun, jac, x_ini, x_out);
      EXPECT_LE((x_out - vec2(1.0, 1.0)).maxCoeff(), EPSILON_LOW);
      EXPECT_TRUE(nlsolver.converged());
      // Solve with damping.
      nlsolver.solve_dumped(fun, jac, x_ini, x_out);
      EXPECT_LE((x_out - vec2(1.0, 1.0)).maxCoeff(), EPSILON_LOW);
      EXPECT_TRUE(nlsolver.converged());
    }
  }
}

// 3D Rosenbrock function.
TEST(Newton, Rosenbrock3D) {
  // Non-linear solver.
  Sandals::Newton<3> nlsolver;
  for (real a = 1.0; a <= 100.0; a += 1.0) {
    for (real b = 1.0; b <= 100.0; b += 1.0) {
      for (real c = 1.0; c <= 100.0; c += 1.0) {
        // Starting entries.
        vec3 x_ini(0.0, 0.0, 0.0);
        vec3 x_out(QUIET_NAN, QUIET_NAN, QUIET_NAN);
        auto fun = [a, b, c](vec3 const & X, vec3 & F) {
          F << a*(1 - X(0)), b*(X(1) - X(0)*X(0)), c*(X(2) - X(1)*X(1));
        };
        auto jac = [a, b, c](vec3 const & X, mat3 & JF) {
          JF << -a, real(0.0), real(0.0), -2*a*X(0), b, real(0.0), real(0.0), -2*b*X(1), c;
        };
        // Solve without damping.
        nlsolver.solve(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - vec3(1.0, 1.0, 1.0)).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
        // Solve with damping.
        nlsolver.solve_dumped(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - vec3(1.0, 1.0, 1.0)).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
      }
    }
  }
}

template <typename F, std::size_t... Is>
void repeat_unrolled(F&& f, std::index_sequence<Is...>) {
  (f(std::integral_constant<std::size_t, Is>{}), ...);
}

template <std::size_t... Iterations, typename F>
void for_unrolled(F&& f) {
  repeat_unrolled(std::forward<F>(f), std::index_sequence<Iterations...>{});
}

// ND Rosenbrock function.
TEST(Newton, RosenbrockND) {
  for_unrolled<2, 5, 10, 20>([](auto d) {
  // Non-linear solver.
    Sandals::Newton<d> nlsolver;
    for (real a = 1.0; a <= 100.0; a += 1.0) {
      for (real b = 1.0; b <= 100.0; b += 1.0) {
        // Starting entries.
        Eigen::Vector<real, d> x_ini = Eigen::Vector<real, d>::Ones();
        Eigen::Vector<real, d> x_out = Eigen::Vector<real, d>::Zero();
        auto fun = [a, b, d](Eigen::Vector<real, d> const & X, Eigen::Vector<real, d> & F) {
          F(0) = a*(1 - X(0));
          for (unsigned i = 1; i < d; ++i) {
            F(i) = b*(X(i-1) - X(i-1)*X(i-1));
          }
        };
        auto jac = [a, b, d](Eigen::Vector<real, d> const & X, Eigen::Matrix<real, d, d> & JF) {
          JF(0, 0) = -a;
          for (unsigned i = 1; i < d; ++i) {
            JF(i, i-1) = -2*b*X(i-1);
            JF(i, i)   = b;
          }
        };
        // Solve without damping.
        nlsolver.solve(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - Eigen::Vector<real, d>::Ones()).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
        // Solve with damping.
        nlsolver.solve_dumped(fun, jac, x_ini, x_out);
        EXPECT_LE((x_out - Eigen::Vector<real, d>::Ones()).maxCoeff(), EPSILON_LOW);
        EXPECT_TRUE(nlsolver.converged());
      }
    }
  });
}

// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
