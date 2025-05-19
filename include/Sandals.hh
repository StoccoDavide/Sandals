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

#endif // INCLUDE_SANDALS_HH
