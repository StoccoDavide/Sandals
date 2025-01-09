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
#include "SimplePendulumExplicit.hh"
#include <matplot/matplot.h>

using namespace Sandals;
using namespace matplot;

int main(void) {
  SimplePendulumExplicit system;
  RK4<2> integrator(std::make_shared<SimplePendulumExplicit>(system));
  RK4<2>::Time time = Eigen::VectorXd::LinSpaced(200, 0.0, 100.0);
  RK4<2>::Solution sol;
  RK4<2>::Solution sol_adaptive;
  integrator.solve(time, system.ics(), sol);
  integrator.adaptive_solve(time, system.ics(), sol_adaptive);

  plot(
    sol.std_time(), sol.std_state(0),
    sol.std_time(), sol.std_state(1)
  ); show();
  plot(
    sol_adaptive.std_time(), sol_adaptive.std_state(0),
    sol_adaptive.std_time(), sol_adaptive.std_state(1)
  ); show();

  return 0;
}