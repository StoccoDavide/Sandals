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

#ifndef TESTS_THREE_BODY_IMPLICIT_HH
#define TESTS_THREE_BODY_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;


/* Implementation of the Three Body Problem.

**Problem Statement**

The equations of motion are derived from the Lagrangian of the system

\f[
\mathcal{L} =
  \frac{1}{2} m_1 (u_1^2 + v_1^2) +
  \frac{1}{2} m_2 (u_2^2 + v_2^2) +
  \frac{1}{2} m_3 (u_3^2 + v_3^2) -
  \frac{G m_1 m_2}{\sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2}} -
  \frac{G m_1 m_3}{\sqrt{(x_1 - x_3)^2 + (y_1 - y_3)^2}} -
  \frac{G m_2 m_3}{\sqrt{(x_2 - x_3)^2 + (y_2 - y_3)^2}} \text{.}
\f]

The equations of motion are then derived by taking the partial derivatives of the Lagrangian
with respect to the generalized coordinates and setting them equal to zero. The equations of
motion written as a system of first order ODEs are then

\f[
\begin{cases}
  x_1^\prime - u_1 = 0 \\
  x_2^\prime - u_2 = 0 \\
  x_3^\prime - u_3 = 0 \\
  y_1^\prime - v_1 = 0 \\
  y_2^\prime - v_2 = 0 \\
  y_3^\prime - v_3 = 0 \\
  u_1^\prime + G m_2\frac{x_1 - x_2}{d_12^3} + G m_3\frac{x_1 - x_3}{d_13^3} = 0 \\
  u_2^\prime + G m_3\frac{x_2 - x_3}{d_23^3} + G m_1\frac{x_2 - x_1}{d_12^3} = 0 \\
  u_3^\prime + G m_1\frac{x_3 - x_1}{d_13^3} + G m_2\frac{x_3 - x_2}{d_23^3} = 0 \\
  v_1^\prime + G m_2\frac{y_1 - y_2}{d_12^3} + G m_3\frac{y_1 - y_3}{d_13^3} = 0 \\
  v_2^\prime + G m_3\frac{y_2 - y_3}{d_23^3} + G m_1\frac{y_2 - y_1}{d_12^3} = 0 \\
  v_3^\prime + G m_1\frac{y_3 - y_1}{d_13^3} + G m_2\frac{y_3 - y_2}{d_23^3} = 0  \text{,}
\end{cases}
\f]

where the distances are defined as

\f[
\begin{cases}
  d_12 = \sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2} \\
  d_13 = \sqrt{(x_1 - x_3)^2 + (y_1 - y_3)^2} \\
  d_23 = \sqrt{(x_2 - x_3)^2 + (y_2 - y_3)^2} \text{,}
\end{cases}
\f]

and \f$ G \f$ is the gravitational constant, \f$ m_1 \f$ is the mass of the first body, \f$ m_2
\f$ is the mass of the second body, and \f$ m_3 \f$ is the mass of the third mass. A feasible
set of initial conditions for creatinting a stable orbit (infinity shaped like) are

\f[
\begin{cases}
  x_1(0) =  0.97000436 \\
  x_2(0) = -0.97000436 \\
  x_3(0) =  0.0 \\
  y_1(0) = -0.24308753 \\
  y_2(0) =  0.24308753 \\
  y_3(0) =  0.0 \\
  u_1(0) =  \frac{0.93240737}{2} \\
  u_2(0) =  \frac{0.93240737}{2} \\
  u_3(0) = -0.93240737 \\
  v_1(0) =  \frac{0.86473146}{2} \\
  v_2(0) =  \frac{0.86473146}{2} \\
  v_3(0) = -0.86473146 \text{.}
\end{cases}
\f]
*/
template<typename Real = double>
class ThreeBodyImplicit : public Implicit<Real, 12, 0>
{
public:
  using VectorF  = typename Implicit<Real, 12, 0>::VectorF;
  using MatrixJF = typename Implicit<Real, 12, 0>::MatrixJF;
  using VectorH  = typename Implicit<Real, 12, 0>::VectorH;
  using MatrixJH = typename Implicit<Real, 12, 0>::MatrixJH;

private:
  Real m_g{1.0};   // Gravitational constant
  Real m_m_1{1.0}; // Body 1 mass
  Real m_m_2{1.0}; // Body 2 mass
  Real m_m_3{1.0}; // Body 3 mass

public:
  ThreeBodyImplicit() : Implicit<Real, 12, 0>("ThreeBodyImplicit") {}

  ~ThreeBodyImplicit() {}

  VectorF F(VectorF const &x, VectorF const &x_dot, Real /*t*/)  const override
  {
    #define CMD "Sandals::ThreeBodyImplicit::F(...): "

    // Compute the distances
    Real d_12{std::sqrt(std::pow(x(0) - x(1), 2.0) + std::pow(x(3) - x(4), 2.0))};
    Real d_13{std::sqrt(std::pow(x(0) - x(2), 2.0) + std::pow(x(3) - x(5), 2.0))};
    Real d_23{std::sqrt(std::pow(x(1) - x(2), 2.0) + std::pow(x(4) - x(5), 2.0))};

    // Compute the output
    VectorF F;
    F <<
      x_dot(0)  - x(6),
      x_dot(1)  - x(7),
      x_dot(2)  - x(8),
      x_dot(3)  - x(9),
      x_dot(4)  - x(10),
      x_dot(5)  - x(11),
      x_dot(6)  + this->m_g*this->m_m_2*(x(0) - x(1))/std::pow(d_12, 3.0) + this->m_g*this->m_m_3*(x(0) - x(2))/std::pow(d_13, 3.0),
      x_dot(7)  + this->m_g*this->m_m_3*(x(1) - x(2))/std::pow(d_23, 3.0) + this->m_g*this->m_m_1*(x(1) - x(0))/std::pow(d_12, 3.0),
      x_dot(8)  + this->m_g*this->m_m_1*(x(2) - x(0))/std::pow(d_13, 3.0) + this->m_g*this->m_m_2*(x(2) - x(1))/std::pow(d_23, 3.0),
      x_dot(9)  + this->m_g*this->m_m_2*(x(3) - x(4))/std::pow(d_12, 3.0) + this->m_g*this->m_m_3*(x(3) - x(5))/std::pow(d_13, 3.0),
      x_dot(10) + this->m_g*this->m_m_3*(x(4) - x(5))/std::pow(d_23, 3.0) + this->m_g*this->m_m_1*(x(4) - x(3))/std::pow(d_12, 3.0),
      x_dot(11) + this->m_g*this->m_m_1*(x(5) - x(3))/std::pow(d_13, 3.0) + this->m_g*this->m_m_2*(x(5) - x(4))/std::pow(d_23, 3.0);

    return F;

    #undef CMD
  }

  MatrixJF JF_x(VectorF const &x, VectorF const &/*x_dot*/, Real /*t*/) const override
  {
    #define CMD "Sandals::ThreeBodyImplicit::JF(...): "

    // Compute the distances
    Real d_12{std::sqrt(std::pow(x(0) - x(1), 2.0) + std::pow(x(3) - x(4), 2.0))};
    Real d_13{std::sqrt(std::pow(x(0) - x(2), 2.0) + std::pow(x(3) - x(5), 2.0))};
    Real d_23{std::sqrt(std::pow(x(1) - x(2), 2.0) + std::pow(x(4) - x(5), 2.0))};

    // Compute the Jacobian of the three body problem
    MatrixJF JF_x;
    JF_x <<
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
       this->m_g*this->m_m_2/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(0) - x(1))*(2*x(0) - 2*x(1))/std::pow(d_12, 5.0) + this->m_g*this->m_m_3/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(0) - x(2))*(2*x(0) - 2*x(2))/std::pow(d_13, 5.0),
      -this->m_g*this->m_m_2/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(0) - x(1))*(2*x(1) - 2*x(0))/std::pow(d_12, 5.0),
      -this->m_g*this->m_m_3/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(0) - x(2))*(2*x(2) - 2*x(0))/std::pow(d_13, 5.0),
      -3.0/2.0*this->m_g*this->m_m_2*(x(0) - x(1))*(2*x(3) - 2*x(4))/std::pow(d_12, 5.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(0) - x(2))*(2*x(3) - 2*x(5))/std::pow(d_13, 5.0),
      -3.0/2.0*this->m_g*this->m_m_2*(x(0) - x(1))*(2*x(4) - 2*x(3))/std::pow(d_12, 5.0),
      -3.0/2.0*this->m_g*this->m_m_3*(x(0) - x(2))*(2*x(5) - 2*x(3))/std::pow(d_13, 5.0),
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -this->m_g*this->m_m_1/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(1) - x(0))*(2*x(0) - 2*x(1))/std::pow(d_12, 5.0),
       this->m_g*this->m_m_3/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(1) - x(2))*(2*x(1) - 2*x(2))/std::pow(d_23, 5.0) + this->m_g*this->m_m_1/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(1) - x(0))*(2*x(1) - 2*x(0))/std::pow(d_12, 5.0),
      -this->m_g*this->m_m_3/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(1) - x(2))*(2*x(2) - 2*x(1))/std::pow(d_23, 5.0),
      -3.0/2.0*this->m_g*this->m_m_1*(x(1) - x(0))*(2*x(3) - 2*x(4))/std::pow(d_12, 5.0),
      -3.0/2.0*this->m_g*this->m_m_3*(x(1) - x(2))*(2*x(4) - 2*x(5))/std::pow(d_23, 5.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(1) - x(0))*(2*x(4) - 2*x(3))/std::pow(d_12, 5.0),
      -3.0/2.0*this->m_g*this->m_m_3*(x(1) - x(2))*(2*x(5) - 2*x(4))/std::pow(d_23, 5.0),
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -this->m_g*this->m_m_1/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(2) - x(0))*(2*x(0) - 2*x(2))/std::pow(d_13, 5.0),
      -this->m_g*this->m_m_2/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(2) - x(1))*(2*x(1) - 2*x(2))/std::pow(d_23, 5.0),
       this->m_g*this->m_m_1/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(2) - x(0))*(2*x(2) - 2*x(0))/std::pow(d_13, 5.0) + this->m_g*this->m_m_2/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(2) - x(1))*(2*x(2) - 2*x(1))/std::pow(d_23, 5.0),
      -3.0/2.0*this->m_g*this->m_m_1*(x(2) - x(0))*(2*x(3) - 2*x(5))/std::pow(d_13, 5.0),
      -3.0/2.0*this->m_g*this->m_m_2*(x(2) - x(1))*(2*x(4) - 2*x(5))/std::pow(d_23, 5.0),
      -3.0/2.0*this->m_g*this->m_m_1*(x(2) - x(0))*(2*x(5) - 2*x(3))/std::pow(d_13, 5.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(2) - x(1))*(2*x(5) - 2*x(4))/std::pow(d_23, 5.0),
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -3.0/2.0*this->m_g*this->m_m_2*(x(3) - x(4))*(2*x(0) - 2*x(1))/std::pow(d_12, 5.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(3) - x(5))*(2*x(0) - 2*x(2))/std::pow(d_13, 5.0),
      -3.0/2.0*this->m_g*this->m_m_2*(x(3) - x(4))*(2*x(1) - 2*x(0))/std::pow(d_12, 5.0),
      -3.0/2.0*this->m_g*this->m_m_3*(x(3) - x(5))*(2*x(2) - 2*x(0))/std::pow(d_13, 5.0),
       this->m_g*this->m_m_2/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(3) - x(4))*(2*x(3) - 2*x(4))/std::pow(d_12, 5.0) + this->m_g*this->m_m_3/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(3) - x(5))*(2*x(3) - 2*x(5))/std::pow(d_13, 5.0),
      -this->m_g*this->m_m_2/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(3) - x(4))*(2*x(4) - 2*x(3))/std::pow(d_12, 5.0),
      -this->m_g*this->m_m_3/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(3) - x(5))*(2*x(5) - 2*x(3))/std::pow(d_13, 5.0),
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -3.0/2.0*this->m_g*this->m_m_1*(x(4) - x(3))*(2*x(0) - 2*x(1))/std::pow(d_12, 5.0),
      -3.0/2.0*this->m_g*this->m_m_3*(x(4) - x(5))*(2*x(1) - 2*x(2))/std::pow(d_23, 5.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(4) - x(3))*(2*x(1) - 2*x(0))/std::pow(d_12, 5.0),
      -3.0/2.0*this->m_g*this->m_m_3*(x(4) - x(5))*(2*x(2) - 2*x(1))/std::pow(d_23, 5.0),
      -this->m_g*this->m_m_1/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(4) - x(3))*(2*x(3) - 2*x(4))/std::pow(d_12, 5.0),
       this->m_g*this->m_m_3/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(4) - x(5))*(2*x(4) - 2*x(5))/std::pow(d_23, 5.0) + this->m_g*this->m_m_1/std::pow(d_12, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(4) - x(3))*(2*x(4) - 2*x(3))/std::pow(d_12, 5.0),
      -this->m_g*this->m_m_3/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_3*(x(4) - x(5))*(2*x(5) - 2*x(4))/std::pow(d_23, 5.0),
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -3.0/2.0*this->m_g*this->m_m_1*(x(5) - x(3))*(2*x(0) - 2*x(2))/std::pow(d_13, 5.0),
      -3.0/2.0*this->m_g*this->m_m_2*(x(5) - x(4))*(2*x(1) - 2*x(2))/std::pow(d_23, 5.0),
      -3.0/2.0*this->m_g*this->m_m_1*(x(5) - x(3))*(2*x(2) - 2*x(0))/std::pow(d_13, 5.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(5) - x(4))*(2*x(2) - 2*x(1))/std::pow(d_23, 5.0),
      -this->m_g*this->m_m_1/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(5) - x(3))*(2*x(3) - 2*x(5))/std::pow(d_13, 5.0),
      -this->m_g*this->m_m_2/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(5) - x(4))*(2*x(4) - 2*x(5))/std::pow(d_23, 5.0),
       this->m_g*this->m_m_1/std::pow(d_13, 3.0) - 3.0/2.0*this->m_g*this->m_m_1*(x(5) - x(3))*(2*x(5) - 2*x(3))/std::pow(d_13, 5.0) + this->m_g*this->m_m_2/std::pow(d_23, 3.0) - 3.0/2.0*this->m_g*this->m_m_2*(x(5) - x(4))*(2*x(5) - 2*x(4))/std::pow(d_23, 5.0),
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    return JF_x;

    #undef CMD
  }

  MatrixJF JF_x_dot(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const override
  {return MatrixJF::Identity();}

  VectorH h(VectorF const &/*x*/, Real /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const &/*x*/, Real /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}

  VectorF ics() const
  {
    VectorF ics;
    ics <<
      0.97000436,
     -0.97000436,
      0.0,
     -0.24308753,
      0.24308753,
      0.0,
      0.93240737/2.0,
      0.93240737/2.0,
     -0.93240737,
      0.86473146/2.0,
      0.86473146/2.0,
     -0.86473146;
    return ics;
  }

};

#endif // TESTS_THREE_BODY_IMPLICIT_HH
