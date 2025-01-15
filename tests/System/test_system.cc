/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"
#include "ODEs/ThreeBodyImplicit.hh"
#include "ODEs/ThreeBodyExplicit.hh"

#ifdef SANDALS_ENABLE_PLOTTING
  #include <matplot/matplot.h>
#endif

using namespace Sandals;
#ifdef SANDALS_ENABLE_PLOTTING
using namespace matplot;
#endif

int main() {
  RK4<12, 0>::Time time = Eigen::VectorXd::LinSpaced(1000, 0.0, 100.0);

  ThreeBodyImplicit system_implicit;
  RK4<12, 0> integrator_implicit(std::make_shared<ThreeBodyImplicit>(system_implicit));
  integrator_implicit.enable_projection();
  integrator_implicit.enable_reverse_mode();
  RK4<12, 0>::Solution solution_implicit;
  RK4<12, 0>::Solution solution_implicit_adaptive;
  integrator_implicit.solve(time, system_implicit.ics(), solution_implicit);
  integrator_implicit.adaptive_solve(time, system_implicit.ics(), solution_implicit_adaptive);

  #ifdef SANDALS_ENABLE_PLOTTING
  auto fig_1a = figure();
  figure(fig_1a);
  plot(
    solution_implicit.std_t(), solution_implicit.std_x(0),
    solution_implicit.std_t(), solution_implicit.std_x(1)
  );
  plot(
    solution_implicit_adaptive.std_t(), solution_implicit_adaptive.std_x(0),
    solution_implicit_adaptive.std_t(), solution_implicit_adaptive.std_x(1)
  );

  auto fig_1b = figure();
  figure(fig_1b);
  plot(
    solution_implicit.std_t(), solution_implicit.std_h(0)
  );
  plot(
    solution_implicit_adaptive.std_t(), solution_implicit_adaptive.std_h(0)
  );
  #endif

  ThreeBodyExplicit system_explicit;
  RK4<12, 0> integrator_explicit(std::make_shared<ThreeBodyExplicit>(system_explicit));
  integrator_explicit.enable_projection();
  integrator_explicit.enable_reverse_mode();
  RK4<12, 0>::Solution solution_explicit;
  RK4<12, 0>::Solution solution_explicit_adaptive;
  integrator_explicit.solve(time, system_explicit.ics(), solution_explicit);
  integrator_explicit.adaptive_solve(time, system_explicit.ics(), solution_explicit_adaptive);

  #ifdef SANDALS_ENABLE_PLOTTING
  auto fig_2a = figure();
  figure(fig_2a);
  plot(
    solution_explicit.std_t(), solution_explicit.std_x(0),
    solution_explicit.std_t(), solution_explicit.std_x(1)
  );
  plot(
    solution_explicit_adaptive.std_t(), solution_explicit_adaptive.std_x(0),
    solution_explicit_adaptive.std_t(), solution_explicit_adaptive.std_x(1)
  );

  auto fig_2b = figure();
  figure(fig_2b);
  plot(
    solution_explicit.std_t(), solution_explicit.std_h(0)
  );
  plot(
    solution_explicit_adaptive.std_t(), solution_explicit_adaptive.std_h(0)
  );

  show();
  #endif

  return 0;
}

