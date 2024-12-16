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

using namespace Sandals;

template class Sandals::Newton<2>;

int main(void) {

  std::cout << "TEST BROYDEN NONLINEAR SOLVER" << std::endl << std::endl;

  try {

    // Non-linear solver.
    Sandals::Newton<2> nlsolver;

    // Starting entries.
    real a{1.0}, b{10.0};
    vec2 x_ini(1.5, 2.0);
    vec2 x_out(QUIET_NAN, QUIET_NAN);
    auto fun = [a, b](vec2 const & X, vec2 & F) {F << a*(1 - X(0)), b*(X(1) - X(0)*X(0));};
    auto jac = [a, b](vec2 const & /*X*/, mat2 & JF) {JF << -a, real(0.0), real(0.0), b;};

    // Solve without damping.
    nlsolver.solve(fun, jac, x_ini, x_out);

    // Display the results on command line.
    std::cout
      << "SOLVER: " << nlsolver.name() << std::endl
      << "\tTolerance   = " << nlsolver.tolerance() << " " << std::endl
      << "\tConverged   = " << nlsolver.converged() << " " << std::endl
      << "\tSolution    = [ " << x_out.transpose() << " ]'" << std::endl
      << "\tFun. eval.s = " << nlsolver.function_evaluations() << std::endl
      << "\tJac. eval.s = " << nlsolver.jacobian_evaluations() << std::endl
      << "\tIterations  = " << nlsolver.iterations() << std::endl
      << "\tRelaxations = " << nlsolver.relaxations() << std::endl
      << std::endl;

    // Solve with damping.
    nlsolver.solve_dumped(fun, jac, x_ini, x_out);

    // Display the results on command line.
    std::cout
      << "SOLVER: " << nlsolver.name() << " (damped)" << std::endl
      << "\tTolerance   = " << nlsolver.tolerance() << " " << std::endl
      << "\tConverged   = " << nlsolver.converged() << " " << std::endl
      << "\tSolution    = [" << x_out.transpose() << "]'" << std::endl
      << "\tFun. eval.s = " << nlsolver.function_evaluations() << std::endl
      << "\tJac. eval.s = " << nlsolver.jacobian_evaluations() << std::endl
      << "\tIterations  = " << nlsolver.iterations() << std::endl
      << "\tRelaxations = " << nlsolver.relaxations() << std::endl
      << std::endl;

  } catch(const std::exception& e) {
    std::cout << "ERROR: " << e.what() << std::endl;
  }

  std::cout << "END OF TEST" << std::endl;

  return 0;
}
