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
#include "SimplePendulumImplicit.hh"
#include <matplot/matplot.h>

using namespace Sandals;
using namespace matplot;

int main(void) {
  SimplePendulumImplicit system;
  RK4<2> integrator(std::make_shared<SimplePendulumImplicit>(system));
  RK4<2>::time time = Eigen::VectorXd::LinSpaced(200, 0.0, 100.0);
  RK4<2>::solution sol;

  // integrator.solve(time, system.ics(), sol);
  // plot(sol.std_time(), sol.std_state(0), sol.std_time(), sol.std_state(1));
  // show();

  integrator.adaptive_solve(time, system.ics(), sol);
  plot(sol.std_time(), sol.std_state(0), sol.std_time(), sol.std_state(1));
  show();

  std::cout << "Hi, I'm finding the sandal size for you!" << std::endl;
  return 0;
}
