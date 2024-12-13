/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_SANDALS_HH
#define INCLUDE_SANDALS_HH

#include <iostream>
#include <tuple>
#include <string>
#include <cmath>
#include <Eigen/Dense>

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


namespace Sandals
{

  /*\
   |      _    _ _
   |     / \  | (_) __ _ ___  ___  ___
   |    / _ \ | | |/ _` / __|/ _ \/ __|
   |   / ___ \| | | (_| \__ \  __/\__ \
   |  /_/   \_\_|_|\__,_|___/\___||___/
   |
  \*/

  using real    = double;       //!< Real number type.
  using integer = int;          //!< Integer number type.
  using ostream = std::ostream; //!< Output stream type.

  using vec2 = Eigen::Matrix<real, 2, 1>; //!< 2x1 vector of real number type (column vector)
  using mat2 = Eigen::Matrix<real, 2, 2>; //!< 2x2 matrix of real number type
  using vec3 = Eigen::Matrix<real, 3, 1>; //!< 3x1 vector of real number type (column vector)
  using mat3 = Eigen::Matrix<real, 3, 3>; //!< 3x3 matrix of real number type
  using vec4 = Eigen::Matrix<real, 4, 1>; //!< 4x1 vector of real number type (column vector)
  using mat4 = Eigen::Matrix<real, 4, 4>; //!< 4x4 matrix of real number type
  using vec5 = Eigen::Matrix<real, 5, 1>; //!< 5x1 vector of real number type (column vector)
  using mat5 = Eigen::Matrix<real, 5, 5>; //!< 5x5 matrix of real number type
  using vec6 = Eigen::Matrix<real, 6, 1>; //!< 6x1 vector of real number type (column vector)
  using mat6 = Eigen::Matrix<real, 6, 6>; //!< 6x6 matrix of real number type
  using vec7 = Eigen::Matrix<real, 7, 1>; //!< 7x1 vector of real number type (column vector)
  using mat7 = Eigen::Matrix<real, 7, 7>; //!< 7x7 matrix of real number type
  using vec8 = Eigen::Matrix<real, 8, 1>; //!< 8x1 vector of real number type (column vector)
  using mat8 = Eigen::Matrix<real, 8, 8>; //!< 8x8 matrix of real number type
  using vec9 = Eigen::Matrix<real, 9, 1>; //!< 9x1 vector of real number type (column vector)
  using mat9 = Eigen::Matrix<real, 9, 9>; //!< 9x9 matrix of real number type

  using vecD = Eigen::Matrix<real, Eigen::Dynamic, 1>;              //!< Nx1 vector of real number type (column vector).
  using matD = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>; //!< NxN matrix of real number type.

  using scale     = Eigen::DiagonalMatrix<real, 3>;           //!< 3D scaling transformation type.
  using translate = Eigen::Translation<real, 3>;              //!< 3D translation transformation type.
  using angleaxis = Eigen::AngleAxis<real>;                   //!< 3D rotation transformation type.
  using affine    = Eigen::Transform<real, 3, Eigen::Affine>; //!< 3D affine transformation type.

  /*\
   |    ____                _              _
   |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
   |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
   |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
   |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
   |
  \*/

  static real const EPSILON        = std::numeric_limits<real>::epsilon();            //!< Machine epsilon epsilon static constant value.
  static real const SQRT_EPSILON   = std::sqrt(EPSILON);                              //!< Square root of machine epsilon epsilon static constant value.
  static real const CBRT_EPSILON   = std::cbrt(EPSILON);                              //!< Cubic root of machine epsilon epsilon static constant value.
  static real const EPSILON_HIGH   = real(1.0e-12);                                   //!< High precision epsilon static constant value.
  static real const EPSILON_MEDIUM = real(1.0e-10);                                   //!< Medium precision epsilon static constant value.
  static real const EPSILON_LOW    = real(1.0e-08);                                   //!< Low precision epsilon static constant value.
  static real const INFTY          = std::numeric_limits<real>::infinity();           //!< Infinity static constant value.
  static real const QUIET_NAN      = std::numeric_limits<real>::quiet_NaN();          //!< Not-a-Number static constant value.
  static real const PI             = real(3.141592653589793238462643383279502884197); //!< Pi static constant value.
  static real const PIMUL2         = real(6.283185307179586476925286766559005768394); //<! The value of \f$ 2\pi \f$.
  static real const PIDIV2         = real(1.570796326794896619231321691639751442098); //<! The value of \f$ \pi/2 \f$.
  static real const DEG2RAD        = real(0.017453292519943295769236907684886127134); //!< The value of \f$ \pi/180 \f$.
  static real const RAD2DEG        = real(57.29577951308232087679815481410517033240); //!< The value of \f$ 180/\pi \f$.

  static vec2 const UNITX_VEC2    = vec2::UnitX();             //!< X axis unit vec2 static constant object.
  static vec2 const UNITY_VEC2    = vec2::UnitY();             //!< Y axis unit vec2 static constant object.
  static vec2 const NAN_VEC2      = vec2::Constant(QUIET_NAN); //!< Not-a-Number vec2 static constant object.
  static mat2 const NAN_MAT2      = mat2::Constant(QUIET_NAN); //!< Not-a-Number mat2 static constant object.
  static vec2 const ZEROS_VEC2    = vec2::Constant(real(0.0)); //!< Zeros vec2 static constant object.
  static mat2 const ZEROS_MAT2    = mat2::Constant(real(0.0)); //!< Zeros mat2 static constant object.
  static vec2 const ONES_VEC2     = vec2::Constant(real(1.0)); //!< Ones vec2 static constant object.
  static mat2 const ONES_MAT2     = mat2::Constant(real(1.0)); //!< Ones mat2 static constant object.
  static mat2 const IDENTITY_MAT2 = mat2::Identity();          //!< Identity mat2 static constant object.

  static vec3 const UNITX_VEC3    = vec3::UnitX();             //!< X axis unit vec3 type.
  static vec3 const UNITY_VEC3    = vec3::UnitY();             //!< Y axis unit vec3 type.
  static vec3 const UNITZ_VEC3    = vec3::UnitZ();             //!< Z axis unit vec3 type.
  static vec3 const NAN_VEC3      = vec3::Constant(QUIET_NAN); //!< Not-a-Number vec3 type.
  static mat3 const NAN_MAT3      = mat3::Constant(QUIET_NAN); //!< Not-a-Number mat3 type.
  static vec3 const ZEROS_VEC3    = vec3::Constant(real(0.0)); //!< Zeros vec3 type.
  static mat3 const ZEROS_MAT3    = mat3::Constant(real(0.0)); //!< Zeros mat3 type.
  static vec3 const ONES_VEC3     = vec3::Constant(real(1.0)); //!< Ones vec3 type.
  static mat3 const ONES_MAT3     = mat3::Constant(real(1.0)); //!< Ones mat3 type.
  static mat3 const IDENTITY_MAT3 = mat3::Identity();          //!< Identity mat3 type.

  static vec4 const UNITX_VEC4    = vec4::UnitX();             //!< X axis unit vec4 type.
  static vec4 const UNITY_VEC4    = vec4::UnitY();             //!< Y axis unit vec4 type.
  static vec4 const UNITZ_VEC4    = vec4::UnitZ();             //!< Z axis unit vec4 type.
  static vec4 const UNITW_VEC4    = vec4::UnitW();             //!< W axis unit vec4 type.
  static vec4 const NAN_VEC4      = vec4::Constant(QUIET_NAN); //!< Not-a-Number vec4 type.
  static mat4 const NAN_MAT4      = mat4::Constant(QUIET_NAN); //!< Not-a-Number mat4 type.
  static vec4 const ZEROS_VEC4    = vec4::Constant(real(0.0)); //!< Zeros vec4 type.
  static mat4 const ZEROS_MAT4    = mat4::Constant(real(0.0)); //!< Zeros mat4 type.
  static vec4 const ONES_VEC4     = vec4::Constant(real(1.0)); //!< Ones vec4 type.
  static mat4 const ONES_MAT4     = mat4::Constant(real(1.0)); //!< Ones mat4 type.
  static mat4 const IDENTITY_MAT4 = mat4::Identity();          //!< Identity mat4 type.

  static vec2 DUMMY_VEC2(NAN_VEC2); //!< Dummy vec2 type static non-const object.
  static vec3 DUMMY_VEC3(NAN_VEC3); //!< Dummy vec3 type static non-const object.
  static vec4 DUMMY_VEC4(NAN_VEC4); //!< Dummy vec4 type static non-const object.
  static mat2 DUMMY_MAT2(NAN_MAT2); //!< Dummy mat2 type static non-const object.
  static mat3 DUMMY_MAT3(NAN_MAT3); //!< Dummy mat3 type static non-const object.
  static mat4 DUMMY_MAT4(NAN_MAT4); //!< Dummy mat4 type static non-const object.

} // namespace Sandals

// Finite difference derivative
//#include "Sandals/FiniteDifference.hxx"

// Newton and Bronden non-linear systems of equations solvers
#include "Sandals/NonlinearSolver/NonlinearSolver.hxx"
#include "Sandals/NonlinearSolver/Newton.hxx"
#include "Sandals/NonlinearSolver/Broyden.hxx"

// Runge-Kutta integrators
#include "Sandals/RungeKutta/System/Implicit.hxx"
#include "Sandals/RungeKutta/System/Explicit.hxx"
//#include "Sandals/RungeKutta/System/SemiExplicit.hxx"
//#include "Sandals/RungeKutta/System/Linear.hxx"
#include "Sandals/RungeKutta/RungeKutta.hxx"
#include "Sandals/RungeKutta/Tableau/RK4.hxx"

#endif // INCLUDE_SANDALS_HH
