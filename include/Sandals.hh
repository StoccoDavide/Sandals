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
#include <map>
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

  using Real = double;   //!< Real number type.
  using Size = unsigned; //!< Size number type.

  using Vector2 = Eigen::Vector<Real, 2>;    //!< 2x1 vector of Real number type (column vector)
  using Matrix2 = Eigen::Matrix<Real, 2, 2>; //!< 2x2 matrix of Real number type
  using Vector3 = Eigen::Vector<Real, 3>;    //!< 3x1 vector of Real number type (column vector)
  using Matrix3 = Eigen::Matrix<Real, 3, 3>; //!< 3x3 matrix of Real number type
  using Vector4 = Eigen::Vector<Real, 4>;    //!< 4x1 vector of Real number type (column vector)
  using Matrix4 = Eigen::Matrix<Real, 4, 4>; //!< 4x4 matrix of Real number type
  using Vector5 = Eigen::Vector<Real, 5>;    //!< 5x1 vector of Real number type (column vector)
  using Matrix5 = Eigen::Matrix<Real, 5, 5>; //!< 5x5 matrix of Real number type
  using Vector6 = Eigen::Vector<Real, 6>;    //!< 6x1 vector of Real number type (column vector)
  using Matrix6 = Eigen::Matrix<Real, 6, 6>; //!< 6x6 matrix of Real number type
  using Vector7 = Eigen::Vector<Real, 7>;    //!< 7x1 vector of Real number type (column vector)
  using Matrix7 = Eigen::Matrix<Real, 7, 7>; //!< 7x7 matrix of Real number type
  using Vector8 = Eigen::Vector<Real, 8>;    //!< 8x1 vector of Real number type (column vector)
  using Matrix8 = Eigen::Matrix<Real, 8, 8>; //!< 8x8 matrix of Real number type
  using Vector9 = Eigen::Vector<Real, 9>;    //!< 9x1 vector of Real number type (column vector)
  using Matrix9 = Eigen::Matrix<Real, 9, 9>; //!< 9x9 matrix of Real number type

  using VectorD = Eigen::Vector<Real, Eigen::Dynamic>;                 //!< Nx1 vector of Real number type (column vector).
  using MatrixD = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; //!< NxN matrix of Real number type.

  /*\
   |    ____                _              _
   |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
   |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
   |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
   |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
   |
  \*/

  static Real const EPSILON        = std::numeric_limits<Real>::epsilon();            //!< Machine epsilon epsilon static constant value.
  static Real const SQRT_EPSILON   = std::sqrt(EPSILON);                              //!< Square root of machine epsilon epsilon static constant value.
  static Real const CBRT_EPSILON   = std::cbrt(EPSILON);                              //!< Cubic root of machine epsilon epsilon static constant value.
  static Real const EPSILON_HIGH   = Real(1.0e-12);                                   //!< High precision epsilon static constant value.
  static Real const EPSILON_MEDIUM = Real(1.0e-10);                                   //!< Medium precision epsilon static constant value.
  static Real const EPSILON_LOW    = Real(1.0e-08);                                   //!< Low precision epsilon static constant value.
  static Real const INFTY          = std::numeric_limits<Real>::infinity();           //!< Infinity static constant value.
  static Real const QUIET_NAN      = std::numeric_limits<Real>::quiet_NaN();          //!< Not-a-Number static constant value.
  static Real const PI             = Real(3.141592653589793238462643383279502884197); //!< Pi static constant value.
  static Real const PIMUL2         = Real(6.283185307179586476925286766559005768394); //<! The value of \f$ 2\pi \f$.
  static Real const PIDIV2         = Real(1.570796326794896619231321691639751442098); //<! The value of \f$ \pi/2 \f$.
  static Real const DEG2RAD        = Real(0.017453292519943295769236907684886127134); //!< The value of \f$ \pi/180 \f$.
  static Real const RAD2DEG        = Real(57.29577951308232087679815481410517033240); //!< The value of \f$ 180/\pi \f$.

  static Vector2 const UNITX_VEC2    = Vector2::UnitX();             //!< X axis unit vec2 static constant object.
  static Vector2 const UNITY_VEC2    = Vector2::UnitY();             //!< Y axis unit vec2 static constant object.
  static Vector2 const NAN_VEC2      = Vector2::Constant(QUIET_NAN); //!< Not-a-Number vec2 static constant object.
  static Matrix2 const NAN_MAT2      = Matrix2::Constant(QUIET_NAN); //!< Not-a-Number mat2 static constant object.
  static Vector2 const ZEROS_VEC2    = Vector2::Constant(Real(0.0)); //!< Zeros vec2 static constant object.
  static Matrix2 const ZEROS_MAT2    = Matrix2::Constant(Real(0.0)); //!< Zeros mat2 static constant object.
  static Vector2 const ONES_VEC2     = Vector2::Constant(Real(1.0)); //!< Ones vec2 static constant object.
  static Matrix2 const ONES_MAT2     = Matrix2::Constant(Real(1.0)); //!< Ones mat2 static constant object.
  static Matrix2 const IDENTITY_MAT2 = Matrix2::Identity();          //!< Identity mat2 static constant object.

  static Vector3 const UNITX_VEC3    = Vector3::UnitX();             //!< X axis unit vec3 type.
  static Vector3 const UNITY_VEC3    = Vector3::UnitY();             //!< Y axis unit vec3 type.
  static Vector3 const UNITZ_VEC3    = Vector3::UnitZ();             //!< Z axis unit vec3 type.
  static Vector3 const NAN_VEC3      = Vector3::Constant(QUIET_NAN); //!< Not-a-Number vec3 type.
  static Matrix3 const NAN_MAT3      = Matrix3::Constant(QUIET_NAN); //!< Not-a-Number mat3 type.
  static Vector3 const ZEROS_VEC3    = Vector3::Constant(Real(0.0)); //!< Zeros vec3 type.
  static Matrix3 const ZEROS_MAT3    = Matrix3::Constant(Real(0.0)); //!< Zeros mat3 type.
  static Vector3 const ONES_VEC3     = Vector3::Constant(Real(1.0)); //!< Ones vec3 type.
  static Matrix3 const ONES_MAT3     = Matrix3::Constant(Real(1.0)); //!< Ones mat3 type.
  static Matrix3 const IDENTITY_MAT3 = Matrix3::Identity();          //!< Identity mat3 type.

  static Vector4 const UNITX_VEC4    = Vector4::UnitX();             //!< X axis unit vec4 type.
  static Vector4 const UNITY_VEC4    = Vector4::UnitY();             //!< Y axis unit vec4 type.
  static Vector4 const UNITZ_VEC4    = Vector4::UnitZ();             //!< Z axis unit vec4 type.
  static Vector4 const UNITW_VEC4    = Vector4::UnitW();             //!< W axis unit vec4 type.
  static Vector4 const NAN_VEC4      = Vector4::Constant(QUIET_NAN); //!< Not-a-Number vec4 type.
  static Matrix4 const NAN_MAT4      = Matrix4::Constant(QUIET_NAN); //!< Not-a-Number mat4 type.
  static Vector4 const ZEROS_VEC4    = Vector4::Constant(Real(0.0)); //!< Zeros vec4 type.
  static Matrix4 const ZEROS_MAT4    = Matrix4::Constant(Real(0.0)); //!< Zeros mat4 type.
  static Vector4 const ONES_VEC4     = Vector4::Constant(Real(1.0)); //!< Ones vec4 type.
  static Matrix4 const ONES_MAT4     = Matrix4::Constant(Real(1.0)); //!< Ones mat4 type.
  static Matrix4 const IDENTITY_MAT4 = Matrix4::Identity();          //!< Identity mat4 type.

} // namespace Sandals

// Finite difference derivative
//#include "Sandals/FiniteDifference.hxx"

// Newton and Bronden nonlinear systems of equations solvers
#include "Sandals/NonlinearSolver/NonlinearSolver.hxx"
#include "Sandals/NonlinearSolver/Newton.hxx"
#include "Sandals/NonlinearSolver/Broyden.hxx"

// Runge-Kutta integrators
#include "Sandals/RungeKutta/System/Implicit.hxx"
#include "Sandals/RungeKutta/System/Explicit.hxx"
#include "Sandals/RungeKutta/System/SemiExplicit.hxx"
//#include "Sandals/RungeKutta/System/Linear.hxx"
#include "Sandals/RungeKutta/RungeKutta.hxx"
#include "Sandals/RungeKutta/Tableau/RK4.hxx"

#endif // INCLUDE_SANDALS_HH
