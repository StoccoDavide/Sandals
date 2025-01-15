/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef INCLUDE_SANDALS_HH
#define INCLUDE_SANDALS_HH

// C++ standard libraries
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <numeric>
#include <chrono>

// Eigen library
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

  using Real = double; //!< Real number type.
  using Size = int;    //!< Size number type.

  using Vector0 = Eigen::Vector<Real, 0>;    //!< \f$ 0 \times 1 \f$ vector of Real number type (column vector).
  using Matrix0 = Eigen::Matrix<Real, 0, 0>; //!< \f$ 0 \times 0 \f$ matrix of Real number type.
  using Vector1 = Eigen::Vector<Real, 1>;    //!< \f$ 1 \times 1 \f$ vector of Real number type (column vector).
  using Matrix1 = Eigen::Matrix<Real, 1, 1>; //!< \f$ 1 \times 1 \f$ matrix of Real number type.
  using Vector2 = Eigen::Vector<Real, 2>;    //!< \f$ 2 \times 1 \f$ vector of Real number type (column vector).
  using Matrix2 = Eigen::Matrix<Real, 2, 2>; //!< \f$ 2 \times 2 \f$ matrix of Real number type.
  using Vector3 = Eigen::Vector<Real, 3>;    //!< \f$ 3 \times 1 \f$ vector of Real number type (column vector).
  using Matrix3 = Eigen::Matrix<Real, 3, 3>; //!< \f$ 3 \times 3 \f$ matrix of Real number type.
  using Vector4 = Eigen::Vector<Real, 4>;    //!< \f$ 4 \times 1 \f$ vector of Real number type (column vector).
  using Matrix4 = Eigen::Matrix<Real, 4, 4>; //!< \f$ 4 \times 4 \f$ matrix of Real number type.
  using Vector5 = Eigen::Vector<Real, 5>;    //!< \f$ 5 \times 1 \f$ vector of Real number type (column vector).
  using Matrix5 = Eigen::Matrix<Real, 5, 5>; //!< \f$ 5 \times 5 \f$ matrix of Real number type.
  using Vector6 = Eigen::Vector<Real, 6>;    //!< \f$ 6 \times 1 \f$ vector of Real number type (column vector).
  using Matrix6 = Eigen::Matrix<Real, 6, 6>; //!< \f$ 6 \times 6 \f$ matrix of Real number type.
  using Vector7 = Eigen::Vector<Real, 7>;    //!< \f$ 7 \times 1 \f$ vector of Real number type (column vector).
  using Matrix7 = Eigen::Matrix<Real, 7, 7>; //!< \f$ 7 \times 7 \f$ matrix of Real number type.
  using Vector8 = Eigen::Vector<Real, 8>;    //!< \f$ 8 \times 1 \f$ vector of Real number type (column vector).
  using Matrix8 = Eigen::Matrix<Real, 8, 8>; //!< \f$ 8 \times 8 \f$ matrix of Real number type.
  using Vector9 = Eigen::Vector<Real, 9>;    //!< \f$ 9 \times 1 \f$ vector of Real number type (column vector).
  using Matrix9 = Eigen::Matrix<Real, 9, 9>; //!< \f$ 9 \times 9 \f$ matrix of Real number type.

  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;                 //!< \f$ N \times 1 \f$ vector of Real number type (column vector).
  using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; //!< \f$ N \times N \f$ matrix of Real number type.

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
  static Real const QUIET_NAN      = std::numeric_limits<Real>::quiet_NaN();          //!< Not-a-number static constant value.
  static Real const PI             = Real(3.141592653589793238462643383279502884197); //!< Pi static constant value.
  static Real const PIMUL2         = Real(6.283185307179586476925286766559005768394); //!< The value of \f$ 2\pi \f$.
  static Real const PIDIV2         = Real(1.570796326794896619231321691639751442098); //!< The value of \f$ \pi/2 \f$.
  static Real const DEG2RAD        = Real(0.017453292519943295769236907684886127134); //!< The value of \f$ \pi/180 \f$.
  static Real const RAD2DEG        = Real(57.29577951308232087679815481410517033240); //!< The value of \f$ 180/\pi \f$.

  static Vector1 const NAN_VEC1      = Vector1::Constant(QUIET_NAN); //!< Not-a-number \f$ 1 \times 1 \f$ vector static constant object.
  static Matrix1 const NAN_MAT1      = Matrix1::Constant(QUIET_NAN); //!< Not-a-number \f$ 1 \times 1 \f$ matrix static constant object.
  static Vector1 const ZEROS_VEC1    = Vector1::Zero();              //!< Zeros \f$ 1 \times 1 \f$ vector static constant object.
  static Matrix1 const ZEROS_MAT1    = Matrix1::Zero();              //!< Zeros \f$ 1 \times 1 \f$ matrix static constant object.
  static Vector1 const ONES_VEC1     = Vector1::Ones();              //!< Ones \f$ 1 \times 1 \f$ vector static constant object.
  static Matrix1 const ONES_MAT1     = Matrix1::Ones();              //!< Ones \f$ 1 \times 1 \f$ matrix static constant object.
  static Matrix1 const IDENTITY_MAT1 = Matrix1::Identity();          //!< Identity \f$ 1 \times 1 \f$ matrix static constant object.

  static Vector2 const NAN_VEC2      = Vector2::Constant(QUIET_NAN); //!< Not-a-number \f$ 2 \times 1 \f$ vector static constant object.
  static Matrix2 const NAN_MAT2      = Matrix2::Constant(QUIET_NAN); //!< Not-a-number \f$ 2 \times 2 \f$ matrix static constant object.
  static Vector2 const ZEROS_VEC2    = Vector2::Zero();              //!< Zeros \f$ 2 \times 1 \f$ vector static constant object.
  static Matrix2 const ZEROS_MAT2    = Matrix2::Zero();              //!< Zeros \f$ 2 \times 2 \f$ matrix static constant object.
  static Vector2 const ONES_VEC2     = Vector2::Ones();              //!< Ones \f$ 2 \times 1 \f$ vector static constant object.
  static Matrix2 const ONES_MAT2     = Matrix2::Ones();              //!< Ones \f$ 2 \times 2 \f$ matrix static constant object.
  static Matrix2 const IDENTITY_MAT2 = Matrix2::Identity();          //!< Identity \f$ 2 \times 2 \f$ matrix static constant object.

  static Vector3 const NAN_VEC3      = Vector3::Constant(QUIET_NAN); //!< Not-a-number \f$ 3 \times 1 \f$ vector static constant object.
  static Matrix3 const NAN_MAT3      = Matrix3::Constant(QUIET_NAN); //!< Not-a-number \f$ 3 \times 3 \f$ matrix static constant object.
  static Vector3 const ZEROS_VEC3    = Vector3::Zero();              //!< Zeros \f$ 3 \times 1 \f$ vector static constant object.
  static Matrix3 const ZEROS_MAT3    = Matrix3::Zero();              //!< Zeros \f$ 3 \times 3 \f$ matrix static constant object.
  static Vector3 const ONES_VEC3     = Vector3::Ones();              //!< Ones \f$ 3 \times 1 \f$ vector static constant object.
  static Matrix3 const ONES_MAT3     = Matrix3::Ones();              //!< Ones \f$ 3 \times 3 \f$ matrix static constant object.
  static Matrix3 const IDENTITY_MAT3 = Matrix3::Identity();          //!< Identity \f$ 3 \times 3 \f$ matrix static constant object.

  static Vector4 const NAN_VEC4      = Vector4::Constant(QUIET_NAN); //!< Not-a-number \f$ 4 \times 1 \f$ vector static constant object.
  static Matrix4 const NAN_MAT4      = Matrix4::Constant(QUIET_NAN); //!< Not-a-number \f$ 4 \times 4 \f$ matrix static constant object.
  static Vector4 const ZEROS_VEC4    = Vector4::Zero();              //!< Zeros \f$ 4 \times 1 \f$ vector static constant object.
  static Matrix4 const ZEROS_MAT4    = Matrix4::Zero();              //!< Zeros \f$ 4 \times 4 \f$ matrix static constant object.
  static Vector4 const ONES_VEC4     = Vector4::Ones();              //!< Ones \f$ 4 \times 1 \f$ vector static constant object.
  static Matrix4 const ONES_MAT4     = Matrix4::Ones();              //!< Ones \f$ 4 \times 4 \f$ matrix static constant object.
  static Matrix4 const IDENTITY_MAT4 = Matrix4::Identity();          //!< Identity \f$ 4 \times 4 \f$ matrix static constant object.

  static Vector5 const NAN_VEC5      = Vector5::Constant(QUIET_NAN); //!< Not-a-number \f$ 5 \times 1 \f$ vector static constant object.
  static Matrix5 const NAN_MAT5      = Matrix5::Constant(QUIET_NAN); //!< Not-a-number \f$ 5 \times 5 \f$ matrix static constant object.
  static Vector5 const ZEROS_VEC5    = Vector5::Zero();              //!< Zeros \f$ 5 \times 1 \f$ vector static constant object.
  static Matrix5 const ZEROS_MAT5    = Matrix5::Zero();              //!< Zeros \f$ 5 \times 5 \f$ matrix static constant object.
  static Vector5 const ONES_VEC5     = Vector5::Ones();              //!< Ones \f$ 5 \times 1 \f$ vector static constant object.
  static Matrix5 const ONES_MAT5     = Matrix5::Ones();              //!< Ones \f$ 5 \times 5 \f$ matrix static constant object.
  static Matrix5 const IDENTITY_MAT5 = Matrix5::Identity();          //!< Identity \f$ 5 \times 5 \f$ matrix static constant object.

  static Vector6 const NAN_VEC6      = Vector6::Constant(QUIET_NAN); //!< Not-a-number \f$ 6 \times 1 \f$ vector static constant object.
  static Matrix6 const NAN_MAT6      = Matrix6::Constant(QUIET_NAN); //!< Not-a-number \f$ 6 \times 6 \f$ matrix static constant object.
  static Vector6 const ZEROS_VEC6    = Vector6::Zero();              //!< Zeros \f$ 6 \times 1 \f$ vector static constant object.
  static Matrix6 const ZEROS_MAT6    = Matrix6::Zero();              //!< Zeros \f$ 6 \times 6 \f$ matrix static constant object.
  static Vector6 const ONES_VEC6     = Vector6::Ones();              //!< Ones \f$ 6 \times 1 \f$ vector static constant object.
  static Matrix6 const ONES_MAT6     = Matrix6::Ones();              //!< Ones \f$ 6 \times 6 \f$ matrix static constant object.
  static Matrix6 const IDENTITY_MAT6 = Matrix6::Identity();          //!< Identity \f$ 6 \times 6 \f$ matrix static constant object.

  static Vector7 const NAN_VEC7      = Vector7::Constant(QUIET_NAN); //!< Not-a-number \f$ 7 \times 1 \f$ vector static constant object.
  static Matrix7 const NAN_MAT7      = Matrix7::Constant(QUIET_NAN); //!< Not-a-number \f$ 7 \times 7 \f$ matrix static constant object.
  static Vector7 const ZEROS_VEC7    = Vector7::Zero();              //!< Zeros \f$ 7 \times 1 \f$ vector static constant object.
  static Matrix7 const ZEROS_MAT7    = Matrix7::Zero();              //!< Zeros \f$ 7 \times 7 \f$ matrix static constant object.
  static Vector7 const ONES_VEC7     = Vector7::Ones();              //!< Ones \f$ 7 \times 1 \f$ vector static constant object.
  static Matrix7 const ONES_MAT7     = Matrix7::Ones();              //!< Ones \f$ 7 \times 7 \f$ matrix static constant object.
  static Matrix7 const IDENTITY_MAT7 = Matrix7::Identity();          //!< Identity \f$ 7 \times 7 \f$ matrix static constant object.

  static Vector8 const NAN_VEC8      = Vector8::Constant(QUIET_NAN); //!< Not-a-number \f$ 8 \times 1 \f$ vector static constant object.
  static Matrix8 const NAN_MAT8      = Matrix8::Constant(QUIET_NAN); //!< Not-a-number \f$ 8 \times 8 \f$ matrix static constant object.
  static Vector8 const ZEROS_VEC8    = Vector8::Zero();              //!< Zeros \f$ 8 \times 1 \f$ vector static constant object.
  static Matrix8 const ZEROS_MAT8    = Matrix8::Zero();              //!< Zeros \f$ 8 \times 8 \f$ matrix static constant object.
  static Vector8 const ONES_VEC8     = Vector8::Ones();              //!< Ones \f$ 8 \times 1 \f$ vector static constant object.
  static Matrix8 const ONES_MAT8     = Matrix8::Ones();              //!< Ones \f$ 8 \times 8 \f$ matrix static constant object.
  static Matrix8 const IDENTITY_MAT8 = Matrix8::Identity();          //!< Identity \f$ 8 \times 8 \f$ matrix static constant object.

  static Vector9 const NAN_VEC9      = Vector9::Constant(QUIET_NAN); //!< Not-a-number \f$ 9 \times 1 \f$ vector static constant object.
  static Matrix9 const NAN_MAT9      = Matrix9::Constant(QUIET_NAN); //!< Not-a-number \f$ 9 \times 9 \f$ matrix static constant object.
  static Vector9 const ZEROS_VEC9    = Vector9::Zero();              //!< Zeros \f$ 9 \times 1 \f$ vector static constant object.
  static Matrix9 const ZEROS_MAT9    = Matrix9::Zero();              //!< Zeros \f$ 9 \times 9 \f$ matrix static constant object.
  static Vector9 const ONES_VEC9     = Vector9::Ones();              //!< Ones \f$ 9 \times 1 \f$ vector static constant object.
  static Matrix9 const ONES_MAT9     = Matrix9::Ones();              //!< Ones \f$ 9 \times 9 \f$ matrix static constant object.
  static Matrix9 const IDENTITY_MAT9 = Matrix9::Identity();          //!< Identity \f$ 9 \times 9 \f$ matrix static constant object.

  //! Print Sandals library information.
  //! \param[in] os Output stream.
  void Info(std::ostream &os) {
    os
      << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl
      << "* Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *" << std::endl
      << "*                                                                           *" << std::endl
      << "* The Sandals project is distributed under the BSD 2-Clause License.        *" << std::endl
      << "*                                                                           *" << std::endl
      << "* Davide Stocco                                           Enrico Bertolazzi *" << std::endl
      << "* University of Trento                                 University of Trento *" << std::endl
      << "* e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *" << std::endl
      << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << std::endl;
  }

} // namespace Sandals

// Finite difference derivative
//#include "Sandals/FiniteDifference.hxx"

// Time measurement
#include "Sandals/TicToc.hxx"

// Nonlinear systems of equations solvers
#include "Sandals/NonlinearSolver/NonlinearSolver.hxx"
#include "Sandals/NonlinearSolver/Newton.hxx"
#include "Sandals/NonlinearSolver/Broyden.hxx"

// ODE/DAE system of equations
#include "Sandals/System/Implicit.hxx"
#include "Sandals/System/Explicit.hxx"

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
