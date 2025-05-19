/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_SANDALS_HH
#define INCLUDE_SANDALS_HH

// C++ standard libraries
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <chrono>

// Eigen library
#include <Eigen/Dense>

// Optimist library
#include <Optimist.hh>
#include <Optimist/RootFinder/Newton.hh>

// Print Sandals errors
#ifndef SANDALS_ERROR
#define SANDALS_ERROR(MSG)              \
  {                                     \
    std::ostringstream os;              \
    os << MSG;                          \
    throw std::runtime_error(os.str()); \
  }
#endif

// Assert for Sandals
#ifndef SANDALS_ASSERT
#define SANDALS_ASSERT(COND, MSG) \
  if (!(COND))                    \
  {                               \
    SANDALS_ERROR(MSG);           \
  }
#endif

// Warning for Sandals
#ifndef SANDALS_WARNING
#define SANDALS_WARNING(MSG)        \
  {                                 \
    std::cout << MSG << std::endl;  \
  }
#endif

// Warning assert for Sandals
#ifndef SANDALS_ASSERT_WARNING
#define SANDALS_ASSERT_WARNING(COND, MSG) \
  if (!(COND))                            \
  {                                       \
    SANDALS_WARNING(MSG);                 \
  }
#endif

// Define the basic constants for Sandals
#ifndef SANDALS_BASIC_CONSTANTS
#define SANDALS_BASIC_CONSTANTS(Real) \
  static constexpr Real EPSILON{std::numeric_limits<Real>::epsilon()}; /**< Machine epsilon epsilon static constant value. */ \
  static constexpr Real EPSILON_HIGH{1.0e-12}; /**< High precision epsilon static constant value. */ \
  static constexpr Real EPSILON_MEDIUM{1.0e-10}; /**< Medium precision epsilon static constant value. */ \
  static constexpr Real EPSILON_LOW{1.0e-08}; /**< Low precision epsilon static constant value. */ \
  static constexpr Real INFTY{std::numeric_limits<Real>::infinity()}; /**< Infinity static constant value. */ \
  static constexpr Real QUIET_NAN{std::numeric_limits<Real>::quiet_NaN()}; /**< Not-a-number static constant value. */
#endif

#ifndef SANDALS_DEFAULT_INTEGER_TYPE
#define SANDALS_DEFAULT_INTEGER_TYPE int
#endif

/**
* \brief The namespace for the Sandals library.
*
* The namespace contains all the classes and functions of the Sandals library.
*/
namespace Sandals
{

  /**
  * \brief The Integer type as used for the API.
  *
  * The Integer type, \c \#define the preprocessor symbol \c SANDALS_DEFAULT_INTEGER_TYPE. The default
  * value is \c int.
  */
  using Integer = SANDALS_DEFAULT_INTEGER_TYPE;

} // namespace Sandals

// Time measurement
#include "Sandals/TicToc.hxx"

// ODEs/DAEs system
#include "Sandals/System/Implicit.hxx"
#include "Sandals/System/Explicit.hxx"
#include "Sandals/System/SemiExplicit.hxx"
#include "Sandals/System/Linear.hxx"

// Runge-Kutta integrator
#include "Sandals/Tableau.hxx"
#include "Sandals/Solution.hxx"
#include "Sandals/RungeKutta.hxx"

// Butcher tableau for Runge-Kutta methods
#include "Sandals/Tableau/Chebyshev51.hxx"
#include "Sandals/Tableau/ExplicitEuler.hxx"
#include "Sandals/Tableau/Fehlberg45.hxx"
#include "Sandals/Tableau/GaussLegendre2.hxx"
#include "Sandals/Tableau/GaussLegendre4.hxx"
#include "Sandals/Tableau/GaussLegendre6.hxx"
#include "Sandals/Tableau/Heun2.hxx"
#include "Sandals/Tableau/Heun3.hxx"
#include "Sandals/Tableau/LobattoIIIA2.hxx"
#include "Sandals/Tableau/ImplicitEuler.hxx"
#include "Sandals/Tableau/MTE22.hxx"
#include "Sandals/Tableau/RadauIIA3.hxx"
#include "Sandals/Tableau/RadauIIA5.hxx"
#include "Sandals/Tableau/Ralston2.hxx"
#include "Sandals/Tableau/Ralston3.hxx"
#include "Sandals/Tableau/Ralston4.hxx"
#include "Sandals/Tableau/RK4.hxx"
#include "Sandals/Tableau/SSPIRK33.hxx"
#include "Sandals/Tableau/SSPRK22.hxx"
#include "Sandals/Tableau/SSPRK22star.hxx"
#include "Sandals/Tableau/SSPRK33.hxx"
#include "Sandals/Tableau/SSPRK42.hxx"
#include "Sandals/Tableau/SSPRK43.hxx"
#include "Sandals/Tableau/SSPRK93.hxx"
#include "Sandals/Tableau/SSPRK104.hxx"

#endif // INCLUDE_SANDALS_HH
