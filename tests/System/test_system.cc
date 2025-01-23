/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"

#include "OscillatorImplicit.hh"
#include "OscillatorExplicit.hh"
#include "OscillatorSemiExplicit.hh"

#include "SimplePendulumImplicit.hh"
#include "SimplePendulumExplicit.hh"
#include "SimplePendulumSemiExplicit.hh"

#include "SinCosImplicit.hh"
#include "SinCosExplicit.hh"
#include "SinCosSemiExplicit.hh"

#include "ThreeBodyImplicit.hh"
#include "ThreeBodyExplicit.hh"
#include "ThreeBodySemiExplicit.hh"

#ifdef SANDALS_ENABLE_PLOTTING
  #include <matplot/matplot.h>
#endif

using namespace Sandals;
#ifdef SANDALS_ENABLE_PLOTTING
using namespace matplot;
#endif

int main() {
  RK4<2, 1>::Time time = Eigen::VectorXd::LinSpaced(1000, 0.0, 100.0);

  OscillatorImplicit system_implicit;
  RK4<2, 1> integrator_implicit(std::make_shared<OscillatorImplicit>(system_implicit));
  integrator_implicit.enable_projection();
  integrator_implicit.enable_reverse_mode();
  Solution<2, 1> solution_implicit;
  Solution<2, 1> solution_implicit_adaptive;
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


  OscillatorExplicit system_explicit;
  RK4<2,1> integrator_explicit(std::make_shared<OscillatorExplicit>(system_explicit));
  integrator_explicit.enable_projection();
  integrator_explicit.enable_reverse_mode();
  Solution<2,1> solution_explicit;
  Solution<2,1> solution_explicit_adaptive;
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
  #endif

  OscillatorSemiExplicit system_semiexplicit;
  RK4<2, 1> integrator_semiexplicit(std::make_shared<OscillatorSemiExplicit>(system_semiexplicit));
  integrator_semiexplicit.enable_projection();
  integrator_semiexplicit.enable_reverse_mode();
  Solution<2, 1> solution_semiexplicit;
  Solution<2, 1> solution_semiexplicit_adaptive;
  integrator_semiexplicit.solve(time, system_semiexplicit.ics(), solution_semiexplicit);
  integrator_semiexplicit.adaptive_solve(time, system_semiexplicit.ics(), solution_semiexplicit_adaptive);

  #ifdef SANDALS_ENABLE_PLOTTING
  auto fig_3a = figure();
  figure(fig_3a);
  plot(
    solution_semiexplicit.std_t(), solution_semiexplicit.std_x(0),
    solution_semiexplicit.std_t(), solution_semiexplicit.std_x(1)
  );
  plot(
    solution_semiexplicit_adaptive.std_t(), solution_semiexplicit_adaptive.std_x(0),
    solution_semiexplicit_adaptive.std_t(), solution_semiexplicit_adaptive.std_x(1)
  );

  auto fig_3b = figure();
  figure(fig_3b);
  plot(
    solution_semiexplicit.std_t(), solution_semiexplicit.std_h(0)
  );
  plot(
    solution_semiexplicit_adaptive.std_t(), solution_semiexplicit_adaptive.std_h(0)
  );

  show();
  #endif

  return 0;
}

