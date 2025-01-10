/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"
#include "OscillatorImplicit.hh"
#include "OscillatorExplicit.hh"
#include <matplot/matplot.h>

using namespace Sandals;
using namespace matplot;

int main() {
  RK4<2,1>::Time time = Eigen::VectorXd::LinSpaced(100, 0.0, 100.0);

  OscillatorImplicit system_implicit;
  RK4<2,1> integrator_implicit(std::make_shared<OscillatorImplicit>(system_implicit));
  integrator_implicit.enable_projection();
  RK4<2,1>::Solution solution_implicit;
  RK4<2,1>::Solution solution_implicit_adaptive;
  integrator_implicit.solve(time, system_implicit.ics(), solution_implicit);
  integrator_implicit.adaptive_solve(time, system_implicit.ics(), solution_implicit_adaptive);

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

  OscillatorExplicit system_explicit;
  RK4<2,1> integrator_explicit(std::make_shared<OscillatorExplicit>(system_explicit));
  integrator_explicit.enable_projection();
  RK4<2,1>::Solution solution_explicit;
  RK4<2,1>::Solution solution_explicit_adaptive;
  integrator_explicit.solve(time, system_explicit.ics(), solution_explicit);
  integrator_explicit.adaptive_solve(time, system_explicit.ics(), solution_explicit_adaptive);

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

  return 0;
}

