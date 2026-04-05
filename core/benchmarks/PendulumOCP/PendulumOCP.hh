/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PENDULUM_OCP_HH
#define TESTS_PENDULUM_OCP_HH

#include "Sandals/System/BoundaryValueProblem.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;
using std::pow;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class PendulumOCPindex3 : public Implicit<Real, 10, 0> {
 public:
  constexpr static Integer N = 10;
  constexpr static Integer M = 0;
  using typename Implicit<Real, N, M>::VectorF;
  using typename Implicit<Real, N, M>::MatrixJF;
  using typename Implicit<Real, N, M>::VectorH;
  using typename Implicit<Real, N, M>::MatrixJH;

 private:
  Real m_g{9.81};
  Real m_l{1.0};

 public:
  PendulumOCPindex3() : Implicit<Real, N, M>("PendulumOCPindex3") {}

  PendulumOCPindex3(Real g, Real l)
      : Implicit<Real, N, M>("PendulumOCPindex3"), m_g(g), m_l(l) {}

  ~PendulumOCPindex3() {}

  void g(const Real g) {
    this->m_g = g;
  }

  Real g() const {
    return this->m_g;
  }

  void l(const Real l) {
    this->m_l = l;
  }

  Real l() const {
    return this->m_l;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    const Real g = this->m_g;
    const Real l = this->m_l;
    VectorF F;
    F[0] = x_dot[0] - x[2];
    F[1] = x_dot[1] - x[3];
    F[2] = x_dot[2] + (-l * l * x[8] * x[1] + 4 * x[4]) * x[0] / 2 +
           l * l * x[1] * x[1] * x[7] / 2;
    F[3] = (x[0] * x[0] * x[8] - x[0] * x[1] * x[7]) * l * l / 2 +
           2 * x[1] * x[4] + g + x_dot[3];
    F[4] = -l * l + x[0] * x[0] + x[1] * x[1];
    F[5] = (-l * l * x[8] * x[8] - 4 * x[9]) * x[0] / 2 +
           (l * l * x[8] * x[1] - 4 * x[4]) * x[7] / 2 + x_dot[5];
    F[6] = (-l * l * x[7] * x[7] - 4 * x[9]) * x[1] / 2 +
           l * l * x[7] * x[0] * x[8] / 2 - 2 * x[4] * x[8] + x_dot[6];
    F[7] = x_dot[7] + x[5];
    F[8] = x_dot[8] + x[6];
    F[9] = -2 * x[0] * x[7] - 2 * x[1] * x[8] + x_dot[9];
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    const Real l = this->m_l;
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 2) = -1;
    JF_x(1, 3) = -1;
    JF_x(2, 0) = -l * l * x[8] * x[1] / 2 + 2 * x[4];
    JF_x(2, 1) = -l * l * x[8] * x[0] / 2 + l * l * x[1] * x[7];
    JF_x(2, 4) = 2 * x[0];
    JF_x(2, 7) = l * l * x[1] * x[1] / 2;
    JF_x(2, 8) = -l * l * x[1] * x[0] / 2;
    JF_x(3, 0) = l * l * x[8] * x[0] - l * l * x[1] * x[7] / 2;
    JF_x(3, 1) = -x[0] * x[7] * l * l / 2 + 2 * x[4];
    JF_x(3, 4) = 2 * x[1];
    JF_x(3, 7) = -l * l * x[1] * x[0] / 2;
    JF_x(3, 8) = x[0] * x[0] * l * l / 2;
    JF_x(4, 0) = 2 * x[0];
    JF_x(4, 1) = 2 * x[1];
    JF_x(5, 0) = -l * l * x[8] * x[8] / 2 - 2 * x[9];
    JF_x(5, 1) = l * l * x[8] * x[7] / 2;
    JF_x(5, 4) = -2 * x[7];
    JF_x(5, 7) = l * l * x[8] * x[1] / 2 - 2 * x[4];
    JF_x(5, 8) = -l * l * (x[0] * x[8] - x[1] * x[7] / 2);
    JF_x(5, 9) = -2 * x[0];
    JF_x(6, 0) = l * l * x[8] * x[7] / 2;
    JF_x(6, 1) = -l * l * x[7] * x[7] / 2 - 2 * x[9];
    JF_x(6, 4) = -2 * x[8];
    JF_x(6, 7) = -l * l * x[1] * x[7] + l * l * x[8] * x[0] / 2;
    JF_x(6, 8) = x[0] * x[7] * l * l / 2 - 2 * x[4];
    JF_x(6, 9) = -2 * x[1];
    JF_x(7, 5) = 1;
    JF_x(8, 6) = 1;
    JF_x(9, 0) = -2 * x[7];
    JF_x(9, 1) = -2 * x[8];
    JF_x(9, 7) = -2 * x[0];
    JF_x(9, 8) = -2 * x[1];
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    MatrixJF JF_x_dot;
    JF_x_dot.setZero();
    JF_x_dot(0, 0) = 1;
    JF_x_dot(1, 1) = 1;
    JF_x_dot(2, 2) = 1;
    JF_x_dot(3, 3) = 1;
    JF_x_dot(5, 5) = 1;
    JF_x_dot(6, 6) = 1;
    JF_x_dot(7, 7) = 1;
    JF_x_dot(8, 8) = 1;
    JF_x_dot(9, 9) = 1;
    return JF_x_dot;
  }

  VectorH h(const VectorF & /*x*/, const Real /*t*/) const override {
    return VectorH::Zero();
  }

  MatrixJH Jh_x(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixJH::Zero();
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }

};  // PendulumOCPindex3

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class PendulumOCPindex0 : public Implicit<Real, 10, 2> {
 public:
  constexpr static Integer N = 10;
  constexpr static Integer M = 2;
  using typename Implicit<Real, N, M>::VectorF;
  using typename Implicit<Real, N, M>::MatrixJF;
  using typename Implicit<Real, N, M>::VectorH;
  using typename Implicit<Real, N, M>::MatrixJH;

 private:
  Real m_g{9.81};
  Real m_l{1.0};

 public:
  PendulumOCPindex0() : Implicit<Real, N, M>("PendulumOCPindex0") {}

  PendulumOCPindex0(Real g, Real l)
      : Implicit<Real, N, M>("PendulumOCPindex0"), m_g(g), m_l(l) {}

  ~PendulumOCPindex0() {}

  void g(const Real g) {
    this->m_g = g;
  }

  Real g() const {
    return this->m_g;
  }

  void l(const Real l) {
    this->m_l = l;
  }

  Real l() const {
    return this->m_l;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    const Real g = this->m_g;
    const Real l = this->m_l;
    VectorF F;
    F[0] = x_dot[0] - x[2];
    F[1] = x_dot[1] - x[3];
    F[2] = (pow(x[1], 4) * x[7] * l * l - l * l * x[0] * pow(x[1], 3) * x[8] +
            (l * l * x[0] * x[0] * x[7] + l * l * x[2] * x[9] +
             4 * x[0] * x[4] + 2 * x_dot[2]) *
                x[1] * x[1] -
            l * l * x[0] * (x[8] * x[0] * x[0] + x[3] * x[9]) * x[1] +
            4 * pow(x[0], 3) * x[4] + 2 * x[0] * x[0] * x_dot[2]) /
           (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    F[3] = (pow(x[0], 4) * x[8] * l * l - x[1] * x[7] * pow(x[0], 3) * l * l +
            (x[8] * l * l * x[1] * x[1] + x[9] * l * l * x[3] +
             4 * x[1] * x[4] + 2 * g + 2 * x_dot[3]) *
                x[0] * x[0] -
            l * l * x[1] * (x[7] * x[1] * x[1] + x[2] * x[9]) * x[0] +
            2 * x[1] * x[1] * (2 * x[1] * x[4] + g + x_dot[3])) /
           (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    F[4] =
        ((l * l * x[3] * x[8] + 2 * x_dot[4]) * pow(x[0], 4) +
         (-l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] + 8 * x[2] * x[4]) *
             pow(x[0], 3) +
         (((x[2] * x[7] + x[3] * x[8]) * l * l + 4 * x_dot[4]) * x[1] * x[1] +
          8 * x[4] * x[3] * x[1] + x[9] * l * l * x[3] * x[3] + 3 * g * x[3]) *
             x[0] * x[0] -
         2 *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] / 2 -
              4 * x[2] * x[4] * x[1] + x[9] * l * l * x[2] * x[3]) *
             x[1] * x[0] +
         (l * l * x[2] * x[7] + 2 * x_dot[4]) * pow(x[1], 4) +
         8 * pow(x[1], 3) * x[3] * x[4] +
         (l * l * x[2] * x[2] * x[9] + 3 * g * x[3]) * x[1] * x[1]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    F[5] =
        (pow(x[1], 7) * x[7] * x[8] * l * l +
         (-x[8] * x[8] * l * l * x[0] - 4 * x[4] * x[7] + 2 * x_dot[5]) *
             pow(x[1], 6) +
         (3 * x[7] * x[8] * x[0] * x[0] + x[9] * (x[8] * x[2] + x[7] * x[3])) *
             l * l * pow(x[1], 5) +
         (-3 * x[8] * x[8] * l * l * pow(x[0], 3) +
          (-12 * x[4] * x[7] + 6 * x_dot[5]) * x[0] * x[0] +
          2 * x[9] * l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] -
          8 * x[9] * x[2] * x[4]) *
             pow(x[1], 4) +
         (3 * x[7] * x[8] * l * l * pow(x[0], 4) +
          x[9] * x[9] * l * l * x[2] * x[3] + 16 * x[9] * x[3] * x[0] * x[4]) *
             pow(x[1], 3) +
         6 * x[0] *
             (-x[8] * x[8] * l * l * pow(x[0], 4) / 2 +
              (-2 * x[4] * x[7] + x_dot[5]) * pow(x[0], 3) +
              x[9] * l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] * x[0] / 3 +
              ((x[2] * x[2] - x[3] * x[3] / 2) * l * l * x[9] / 3 + g * x[3]) *
                  x[9]) *
             x[1] * x[1] -
         3 * x[0] * x[0] *
             (-x[7] * x[8] * l * l * pow(x[0], 4) / 3 +
              x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) * x[0] * x[0] / 3 -
              0.16e2 / 0.3e1 * x[9] * x[3] * x[0] * x[4] +
              x[9] * x[9] * l * l * x[2] * x[3]) *
             x[1] -
         pow(x[0], 7) * x[8] * x[8] * l * l +
         (-4 * x[4] * x[7] + 2 * x_dot[5]) * pow(x[0], 6) +
         8 * pow(x[0], 4) * x[2] * x[9] * x[4] +
         (x[9] * x[9] * l * l * x[3] * x[3] + 6 * g * x[9] * x[3]) *
             pow(x[0], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    F[6] =
        (x[7] * pow(x[0], 7) * x[8] * l * l +
         (-x[7] * x[7] * l * l * x[1] - 4 * x[4] * x[8] + 2 * x_dot[6]) *
             pow(x[0], 6) +
         (3 * x[7] * x[8] * x[1] * x[1] + x[9] * (x[8] * x[2] + x[7] * x[3])) *
             l * l * pow(x[0], 5) +
         (-3 * x[7] * x[7] * l * l * pow(x[1], 3) +
          (-12 * x[4] * x[8] + 6 * x_dot[6]) * x[1] * x[1] -
          2 * x[9] * l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] -
          8 * x[9] * x[3] * x[4]) *
             pow(x[0], 4) +
         (3 * x[7] * x[8] * l * l * pow(x[1], 4) +
          x[9] * x[9] * l * l * x[2] * x[3] + 16 * x[9] * x[2] * x[1] * x[4]) *
             pow(x[0], 3) +
         6 *
             (-x[7] * x[7] * l * l * pow(x[1], 4) / 2 +
              (-2 * x[4] * x[8] + x_dot[6]) * pow(x[1], 3) -
              x[9] * l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] / 3 +
              x[9] * (-l * l * (x[2] * x[2] - 2 * x[3] * x[3]) * x[9] / 6 +
                      g * x[3])) *
             x[1] * x[0] * x[0] -
         3 *
             (-x[7] * x[8] * l * l * pow(x[1], 4) / 3 +
              x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] / 3 -
              0.16e2 / 0.3e1 * x[9] * x[2] * x[1] * x[4] +
              x[9] * x[9] * l * l * x[2] * x[3]) *
             x[1] * x[1] * x[0] -
         pow(x[1], 7) * x[7] * x[7] * l * l +
         (-4 * x[4] * x[8] + 2 * x_dot[6]) * pow(x[1], 6) +
         8 * pow(x[1], 4) * x[3] * x[9] * x[4] +
         (l * l * x[2] * x[2] * x[9] * x[9] + 6 * g * x[9] * x[3]) *
             pow(x[1], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    F[7] = ((-l * l * x[7] * x[9] + 2 * x[5] + 2 * x_dot[7]) * pow(x[1], 4) +
            l * l * x[0] * pow(x[1], 3) * x[8] * x[9] +
            ((-l * l * x[7] * x[9] + 4 * x[5] + 4 * x_dot[7]) * x[0] * x[0] -
             8 * x[9] * x[4] * x[0] - x[2] * l * l * x[9] * x[9]) *
                x[1] * x[1] +
            l * l * x[9] * x[0] * (x[8] * x[0] * x[0] + x[3] * x[9]) * x[1] +
            (2 * x_dot[7] + 2 * x[5]) * pow(x[0], 4) -
            8 * pow(x[0], 3) * x[4] * x[9]) *
           pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    F[8] =
        ((-l * l * x[8] * x[9] + 2 * x[6] + 2 * x_dot[8]) * pow(x[0], 4) +
         l * l * pow(x[0], 3) * x[1] * x[7] * x[9] +
         ((-l * l * x[8] * x[9] + 4 * x[6] + 4 * x_dot[8]) * x[1] * x[1] -
          8 * x[9] * x[1] * x[4] - x[9] * x[9] * l * l * x[3] - 3 * g * x[9]) *
             x[0] * x[0] +
         x[9] * l * l * x[1] * (x[7] * x[1] * x[1] + x[2] * x[9]) * x[0] +
         (2 * x_dot[8] + 2 * x[6]) * pow(x[1], 4) -
         8 * pow(x[1], 3) * x[4] * x[9] - 3 * g * x[1] * x[1] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    F[9] = (-2 * pow(x[1], 3) * x[8] +
            (-2 * x[0] * x[7] + x_dot[9]) * x[1] * x[1] +
            (-2 * x[8] * x[0] * x[0] - 4 * x[3] * x[9]) * x[1] -
            2 * pow(x[0], 3) * x[7] + x[0] * x[0] * x_dot[9] -
            4 * x[0] * x[2] * x[9]) /
           (x[0] * x[0] + x[1] * x[1]);
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    const Real g = this->m_g;
    const Real l = this->m_l;
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 2) = -1;
    JF_x(1, 3) = -1;
    JF_x(2, 0) =
        (-l * l * pow(x[0], 4) * x[1] * x[8] -
         2 * l * l * x[0] * x[0] * pow(x[1], 3) * x[8] -
         l * l * pow(x[1], 5) * x[8] +
         l * l * x[0] * x[0] * x[1] * x[3] * x[9] -
         2 * l * l * x[0] * x[1] * x[1] * x[2] * x[9] -
         l * l * pow(x[1], 3) * x[3] * x[9] + 4 * pow(x[0], 4) * x[4] +
         8 * x[0] * x[0] * x[1] * x[1] * x[4] + 4 * pow(x[1], 4) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(2, 1) =
        -l * l *
        (pow(x[0], 5) * x[8] - 2 * pow(x[0], 4) * x[1] * x[7] +
         2 * pow(x[0], 3) * x[1] * x[1] * x[8] -
         4 * x[0] * x[0] * pow(x[1], 3) * x[7] + x[0] * pow(x[1], 4) * x[8] -
         2 * pow(x[1], 5) * x[7] + pow(x[0], 3) * x[3] * x[9] -
         2 * x[0] * x[0] * x[1] * x[2] * x[9] -
         x[0] * x[1] * x[1] * x[3] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(2, 2) =
        x[1] * x[1] * l * l * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(2, 3) =
        -x[1] * l * l * x[0] * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(2, 4) = 2 * x[0];
    JF_x(2, 7) = l * l * x[1] * x[1] / 2;
    JF_x(2, 8) = -l * l * x[1] * x[0] / 2;
    JF_x(2, 9) = -x[1] * l * l * (x[0] * x[3] - x[1] * x[2]) /
                 (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(3, 0) =
        l * l *
        (2 * pow(x[0], 5) * x[8] - pow(x[0], 4) * x[1] * x[7] +
         4 * pow(x[0], 3) * x[1] * x[1] * x[8] -
         2 * x[0] * x[0] * pow(x[1], 3) * x[7] +
         2 * x[0] * pow(x[1], 4) * x[8] - pow(x[1], 5) * x[7] +
         x[0] * x[0] * x[1] * x[2] * x[9] +
         2 * x[0] * x[1] * x[1] * x[3] * x[9] - pow(x[1], 3) * x[2] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(3, 1) =
        (-l * l * pow(x[0], 5) * x[7] -
         2 * l * l * pow(x[0], 3) * x[1] * x[1] * x[7] -
         l * l * x[0] * pow(x[1], 4) * x[7] -
         l * l * pow(x[0], 3) * x[2] * x[9] -
         2 * l * l * x[0] * x[0] * x[1] * x[3] * x[9] +
         l * l * x[0] * x[1] * x[1] * x[2] * x[9] + 4 * pow(x[0], 4) * x[4] +
         8 * x[0] * x[0] * x[1] * x[1] * x[4] + 4 * pow(x[1], 4) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(3, 2) =
        -x[1] * l * l * x[0] * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(3, 3) =
        l * l * x[0] * x[0] * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(3, 4) = 2 * x[1];
    JF_x(3, 7) = -l * l * x[1] * x[0] / 2;
    JF_x(3, 8) = l * l * x[0] * x[0] / 2;
    JF_x(3, 9) = l * l * x[0] * (x[0] * x[3] - x[1] * x[2]) /
                 (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(4, 0) =
        (-l * l * (x[8] * x[2] + x[7] * x[3]) * pow(x[1], 5) +
         (-2 * l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] + 8 * x[2] * x[4]) *
             pow(x[1], 4) -
         2 * x[3] * (l * l * x[2] * x[9] + 8 * x[0] * x[4]) * pow(x[1], 3) -
         6 * x[0] *
             (l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] * x[0] / 3 +
              0.2e1 / 0.3e1 * (x[2] * x[2] - x[3] * x[3] / 2) * l * l * x[9] +
              g * x[3]) *
             x[1] * x[1] +
         6 * x[0] * x[0] *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[0] * x[0] / 6 -
              0.8e1 / 0.3e1 * x[3] * x[0] * x[4] + x[9] * l * l * x[2] * x[3]) *
             x[1] -
         8 * pow(x[0], 4) * x[2] * x[4] +
         (-2 * x[9] * l * l * x[3] * x[3] - 6 * g * x[3]) * pow(x[0], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(4, 1) =
        (-l * l * (x[8] * x[2] + x[7] * x[3]) * pow(x[0], 5) +
         (2 * l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] + 8 * x[4] * x[3]) *
             pow(x[0], 4) -
         2 * x[2] * (x[9] * l * l * x[3] + 8 * x[1] * x[4]) * pow(x[0], 3) -
         6 *
             (-l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] / 3 -
              l * l * (x[2] * x[2] - 2 * x[3] * x[3]) * x[9] / 3 + g * x[3]) *
             x[1] * x[0] * x[0] +
         6 *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] / 6 -
              0.8e1 / 0.3e1 * x[2] * x[4] * x[1] + x[9] * l * l * x[2] * x[3]) *
             x[1] * x[1] * x[0] -
         8 * pow(x[1], 4) * x[3] * x[4] +
         (-2 * l * l * x[2] * x[2] * x[9] - 6 * g * x[3]) * pow(x[1], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(4, 2) =
        (pow(x[1], 4) * x[7] * l * l - l * l * x[0] * pow(x[1], 3) * x[8] +
         (l * l * x[0] * x[0] * x[7] + 2 * l * l * x[2] * x[9] +
          8 * x[0] * x[4]) *
             x[1] * x[1] -
         l * l * x[0] * (x[8] * x[0] * x[0] + 2 * x[3] * x[9]) * x[1] +
         8 * pow(x[0], 3) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(4, 3) =
        (pow(x[0], 4) * x[8] * l * l - x[1] * x[7] * pow(x[0], 3) * l * l +
         (x[8] * l * l * x[1] * x[1] + 2 * x[9] * l * l * x[3] +
          8 * x[1] * x[4] + 3 * g) *
             x[0] * x[0] +
         (-l * l * pow(x[1], 3) * x[7] - 2 * l * l * x[1] * x[2] * x[9]) *
             x[0] +
         8 * pow(x[1], 3) * x[4] + 3 * g * x[1] * x[1]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(4, 4) =
        (4 * x[0] * x[2] + 4 * x[1] * x[3]) / (x[0] * x[0] + x[1] * x[1]);
    JF_x(4, 7) = -x[1] * l * l * (x[0] * x[3] - x[1] * x[2]) /
                 (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(4, 8) = l * l * x[0] * (x[0] * x[3] - x[1] * x[2]) /
                 (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(4, 9) = l * l * pow(x[0] * x[3] - x[1] * x[2], 2) *
                 pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(5, 0) =
        (-l * l * pow(x[1], 8) * x[8] * x[8] +
         2 * l * l *
             ((x[2] * x[7] - x[3] * x[8]) * x[9] -
              2 * x[0] * x[0] * x[8] * x[8]) *
             pow(x[1], 6) -
         6 *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[0] -
              0.8e1 / 0.3e1 * x[4] * x[3]) *
             x[9] * pow(x[1], 5) +
         (l * l * (2 * x[2] * x[2] - x[3] * x[3]) * x[9] * x[9] +
          (-4 * l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] * x[0] +
           48 * x[4] * x[2] * x[0] + 6 * g * x[3]) *
              x[9] -
          6 * x[8] * x[8] * l * l * pow(x[0], 4)) *
             pow(x[1], 4) -
         12 * x[0] *
             (x[9] * l * l * x[2] * x[3] +
              (l * l * (x[8] * x[2] + x[7] * x[3]) * x[0] + 8 * x[4] * x[3]) *
                  x[0] / 3) *
             x[9] * pow(x[1], 3) -
         12 *
             (0.5e1 / 0.6e1 * (x[2] * x[2] - 0.4e1 / 0.5e1 * x[3] * x[3]) * l *
                  l * x[9] * x[9] +
              (l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] * x[0] / 2 -
               0.8e1 / 0.3e1 * x[4] * x[2] * x[0] + g * x[3]) *
                  x[9] +
              x[8] * x[8] * l * l * pow(x[0], 4) / 3) *
             x[0] * x[0] * x[1] * x[1] +
         12 *
             (x[9] * l * l * x[2] * x[3] +
              x[0] *
                  (l * l * (x[8] * x[2] + x[7] * x[3]) * x[0] -
                   24 * x[4] * x[3]) /
                  6) *
             pow(x[0], 3) * x[9] * x[1] -
         3 * l * l * pow(x[0], 4) * x[3] * x[3] * x[9] * x[9] +
         (-16 * pow(x[0], 5) * x[2] * x[4] - 18 * g * pow(x[0], 4) * x[3]) *
             x[9] -
         l * l * pow(x[0], 8) * x[8] * x[8]) *
        pow(x[0] * x[0] + x[1] * x[1], -4) / 2;
    JF_x(5, 1) =
        (l * l * pow(x[0], 8) * x[7] * x[8] -
         l * l *
             (-4 * x[7] * x[8] * x[1] * x[1] +
              x[9] * (x[8] * x[2] + x[7] * x[3])) *
             pow(x[0], 6) +
         4 * (l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] + 4 * x[4] * x[3]) *
             x[9] * pow(x[0], 5) +
         (6 * x[7] * x[8] * l * l * pow(x[1], 4) +
          5 * x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] -
          48 * x[9] * x[2] * x[1] * x[4] -
          3 * x[9] * x[9] * l * l * x[2] * x[3]) *
             pow(x[0], 4) -
         24 *
             (0.4e1 / 0.3e1 * x[4] * x[3] * x[1] -
              l * l * (x[2] * x[2] - 2 * x[3] * x[3]) * x[9] / 6 + g * x[3]) *
             x[1] * x[9] * pow(x[0], 3) +
         18 *
             (0.2e1 / 0.9e1 * x[7] * x[8] * l * l * pow(x[1], 4) +
              0.5e1 / 0.18e2 * x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) *
                  x[1] * x[1] -
              0.16e2 / 0.9e1 * x[9] * x[2] * x[1] * x[4] +
              x[9] * x[9] * l * l * x[2] * x[3]) *
             x[1] * x[1] * x[0] * x[0] -
         24 *
             (l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] / 6 +
              2 * x[4] * x[3] * x[1] +
              (x[2] * x[2] - x[3] * x[3] / 2) * l * l * x[9] / 3 + g * x[3]) *
             pow(x[1], 3) * x[9] * x[0] -
         3 *
             (-x[7] * x[8] * l * l * pow(x[1], 4) / 3 +
              x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] / 3 -
              0.16e2 / 0.3e1 * x[9] * x[2] * x[1] * x[4] +
              x[9] * x[9] * l * l * x[2] * x[3]) *
             pow(x[1], 4)) *
        pow(x[0] * x[0] + x[1] * x[1], -4) / 2;
    JF_x(5, 2) =
        -x[9] *
        (-l * l * pow(x[1], 5) * x[8] +
         (-2 * l * l * x[0] * x[7] + 8 * x[4]) * pow(x[1], 4) -
         l * l * pow(x[1], 3) * x[3] * x[9] -
         2 * l * l * x[0] * (x[0] * x[0] * x[7] + 2 * x[2] * x[9]) * x[1] *
             x[1] +
         l * l * x[0] * x[0] * (x[8] * x[0] * x[0] + 3 * x[3] * x[9]) * x[1] -
         8 * pow(x[0], 4) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(5, 3) =
        3 *
        (-l * l * pow(x[0], 4) * x[1] * x[7] / 6 +
         (-x[8] * l * l * x[1] * x[1] / 3 + x[9] * l * l * x[3] / 3 +
          0.8e1 / 0.3e1 * x[1] * x[4] + g) *
             pow(x[0], 3) -
         l * l * x[0] * x[0] * x[1] * x[2] * x[9] / 2 +
         x[1] * x[1] *
             (-x[8] * l * l * x[1] * x[1] / 3 - x[9] * l * l * x[3] / 3 +
              0.8e1 / 0.3e1 * x[1] * x[4] + g) *
             x[0] +
         l * l * pow(x[1], 3) * (x[7] * x[1] * x[1] + x[2] * x[9]) / 6) *
        x[9] * pow(x[0] * x[0] + x[1] * x[1], -3);
    JF_x(5, 4) = (-2 * pow(x[1], 4) * x[7] +
                  (-4 * x[0] * x[0] * x[7] - 4 * x[2] * x[9]) * x[1] * x[1] +
                  8 * x[3] * x[0] * x[1] * x[9] - 2 * pow(x[0], 4) * x[7] +
                  4 * x[0] * x[0] * x[2] * x[9]) *
                 pow(x[0] * x[0] + x[1] * x[1], -2);
    JF_x(5, 7) =
        (l * l * pow(x[0], 4) * x[1] * x[8] +
         2 * l * l * x[0] * x[0] * pow(x[1], 3) * x[8] +
         l * l * pow(x[1], 5) * x[8] -
         l * l * x[0] * x[0] * x[1] * x[3] * x[9] +
         2 * l * l * x[0] * x[1] * x[1] * x[2] * x[9] +
         l * l * pow(x[1], 3) * x[3] * x[9] - 4 * pow(x[0], 4) * x[4] -
         8 * x[0] * x[0] * x[1] * x[1] * x[4] - 4 * pow(x[1], 4) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(5, 8) =
        -l * l *
        (2 * pow(x[0], 5) * x[8] - pow(x[0], 4) * x[1] * x[7] +
         4 * pow(x[0], 3) * x[1] * x[1] * x[8] -
         2 * x[0] * x[0] * pow(x[1], 3) * x[7] +
         2 * x[0] * pow(x[1], 4) * x[8] - pow(x[1], 5) * x[7] +
         x[0] * x[0] * x[1] * x[2] * x[9] +
         2 * x[0] * x[1] * x[1] * x[3] * x[9] - pow(x[1], 3) * x[2] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(5, 9) =
        (l * l * (x[8] * x[2] + x[7] * x[3]) * pow(x[1], 5) +
         (2 * l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] - 8 * x[2] * x[4]) *
             pow(x[1], 4) +
         2 * x[3] * (l * l * x[2] * x[9] + 8 * x[0] * x[4]) * pow(x[1], 3) +
         6 * x[0] *
             (l * l * (x[2] * x[7] - x[3] * x[8]) * x[0] * x[0] / 3 +
              0.2e1 / 0.3e1 * (x[2] * x[2] - x[3] * x[3] / 2) * l * l * x[9] +
              g * x[3]) *
             x[1] * x[1] -
         6 * x[0] * x[0] *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[0] * x[0] / 6 -
              0.8e1 / 0.3e1 * x[3] * x[0] * x[4] + x[9] * l * l * x[2] * x[3]) *
             x[1] +
         8 * pow(x[0], 4) * x[2] * x[4] +
         (2 * x[9] * l * l * x[3] * x[3] + 6 * g * x[3]) * pow(x[0], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(6, 0) =
        (l * l * pow(x[0], 8) * x[7] * x[8] -
         l * l *
             (-4 * x[7] * x[8] * x[1] * x[1] +
              x[9] * (x[8] * x[2] + x[7] * x[3])) *
             pow(x[0], 6) +
         4 * (l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] + 4 * x[4] * x[3]) *
             x[9] * pow(x[0], 5) +
         (6 * x[7] * x[8] * l * l * pow(x[1], 4) +
          5 * x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] -
          48 * x[9] * x[2] * x[1] * x[4] -
          3 * x[9] * x[9] * l * l * x[2] * x[3]) *
             pow(x[0], 4) -
         24 *
             (0.4e1 / 0.3e1 * x[4] * x[3] * x[1] -
              l * l * (x[2] * x[2] - 2 * x[3] * x[3]) * x[9] / 6 + g * x[3]) *
             x[1] * x[9] * pow(x[0], 3) +
         18 *
             (0.2e1 / 0.9e1 * x[7] * x[8] * l * l * pow(x[1], 4) +
              0.5e1 / 0.18e2 * x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) *
                  x[1] * x[1] -
              0.16e2 / 0.9e1 * x[9] * x[2] * x[1] * x[4] +
              x[9] * x[9] * l * l * x[2] * x[3]) *
             x[1] * x[1] * x[0] * x[0] -
         24 *
             (l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] / 6 +
              2 * x[4] * x[3] * x[1] +
              (x[2] * x[2] - x[3] * x[3] / 2) * l * l * x[9] / 3 + g * x[3]) *
             pow(x[1], 3) * x[9] * x[0] -
         3 *
             (-x[7] * x[8] * l * l * pow(x[1], 4) / 3 +
              x[9] * l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] / 3 -
              0.16e2 / 0.3e1 * x[9] * x[2] * x[1] * x[4] +
              x[9] * x[9] * l * l * x[2] * x[3]) *
             pow(x[1], 4)) *
        pow(x[0] * x[0] + x[1] * x[1], -4) / 2;
    JF_x(6, 1) =
        (-l * l * pow(x[0], 8) * x[7] * x[7] -
         2 *
             ((x[2] * x[7] - x[3] * x[8]) * x[9] +
              2 * x[1] * x[1] * x[7] * x[7]) *
             l * l * pow(x[0], 6) -
         6 *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] -
              0.8e1 / 0.3e1 * x[2] * x[4]) *
             x[9] * pow(x[0], 5) +
         (-l * l * (x[2] * x[2] - 2 * x[3] * x[3]) * x[9] * x[9] +
          (4 * l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] +
           48 * x[4] * x[3] * x[1] + 6 * g * x[3]) *
              x[9] -
          6 * x[7] * x[7] * l * l * pow(x[1], 4)) *
             pow(x[0], 4) -
         12 *
             (x[9] * l * l * x[2] * x[3] +
              (l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] + 8 * x[2] * x[4]) *
                  x[1] / 3) *
             x[1] * x[9] * pow(x[0], 3) -
         12 *
             (-0.2e1 / 0.3e1 * (x[2] * x[2] - 0.5e1 / 0.4e1 * x[3] * x[3]) * l *
                  l * x[9] * x[9] +
              (-l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] / 2 -
               0.8e1 / 0.3e1 * x[4] * x[3] * x[1] + g * x[3]) *
                  x[9] +
              x[7] * x[7] * l * l * pow(x[1], 4) / 3) *
             x[1] * x[1] * x[0] * x[0] +
         12 *
             (x[9] * l * l * x[2] * x[3] +
              (l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] - 24 * x[2] * x[4]) *
                  x[1] / 6) *
             pow(x[1], 3) * x[9] * x[0] -
         18 *
             (l * l * x[2] * x[2] * x[9] * x[9] / 6 +
              x[3] * (0.8e1 / 0.9e1 * x[1] * x[4] + g) * x[9] +
              x[7] * x[7] * l * l * pow(x[1], 4) / 18) *
             pow(x[1], 4)) *
        pow(x[0] * x[0] + x[1] * x[1], -4) / 2;
    JF_x(6, 2) =
        x[9] *
        (l * l * pow(x[0], 5) * x[8] - 2 * l * l * pow(x[0], 4) * x[1] * x[7] +
         (x[9] * l * l * x[3] + 16 * x[1] * x[4]) * pow(x[0], 3) -
         2 * l * l * x[1] * (x[7] * x[1] * x[1] + x[2] * x[9]) * x[0] * x[0] +
         (-l * l * pow(x[1], 4) * x[8] - 3 * x[9] * l * l * x[3] * x[1] * x[1] +
          16 * pow(x[1], 3) * x[4]) *
             x[0] +
         2 * l * l * pow(x[1], 3) * x[2] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(6, 3) =
        x[9] *
        (l * l * pow(x[0], 5) * x[7] + 2 * l * l * pow(x[0], 4) * x[1] * x[8] +
         2 * l * l * x[0] * x[0] * pow(x[1], 3) * x[8] -
         l * l * x[0] * pow(x[1], 4) * x[7] +
         l * l * pow(x[0], 3) * x[2] * x[9] +
         4 * l * l * x[0] * x[0] * x[1] * x[3] * x[9] -
         3 * l * l * x[0] * x[1] * x[1] * x[2] * x[9] -
         8 * pow(x[0], 4) * x[4] + 8 * pow(x[1], 4) * x[4] +
         6 * g * x[0] * x[0] * x[1] + 6 * g * pow(x[1], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(6, 4) = (-2 * pow(x[1], 4) * x[8] +
                  (-4 * x[8] * x[0] * x[0] + 4 * x[3] * x[9]) * x[1] * x[1] +
                  8 * x[0] * x[1] * x[2] * x[9] - 2 * pow(x[0], 4) * x[8] -
                  4 * x[0] * x[0] * x[3] * x[9]) *
                 pow(x[0] * x[0] + x[1] * x[1], -2);
    JF_x(6, 7) =
        l * l *
        (pow(x[0], 5) * x[8] - 2 * pow(x[0], 4) * x[1] * x[7] +
         2 * pow(x[0], 3) * x[1] * x[1] * x[8] -
         4 * x[0] * x[0] * pow(x[1], 3) * x[7] + x[0] * pow(x[1], 4) * x[8] -
         2 * pow(x[1], 5) * x[7] + pow(x[0], 3) * x[3] * x[9] -
         2 * x[0] * x[0] * x[1] * x[2] * x[9] -
         x[0] * x[1] * x[1] * x[3] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(6, 8) =
        (l * l * pow(x[0], 5) * x[7] +
         2 * l * l * pow(x[0], 3) * x[1] * x[1] * x[7] +
         l * l * x[0] * pow(x[1], 4) * x[7] +
         l * l * pow(x[0], 3) * x[2] * x[9] +
         2 * l * l * x[0] * x[0] * x[1] * x[3] * x[9] -
         l * l * x[0] * x[1] * x[1] * x[2] * x[9] - 4 * pow(x[0], 4) * x[4] -
         8 * x[0] * x[0] * x[1] * x[1] * x[4] - 4 * pow(x[1], 4) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(6, 9) =
        (l * l * (x[8] * x[2] + x[7] * x[3]) * pow(x[0], 5) +
         (-2 * l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] - 8 * x[4] * x[3]) *
             pow(x[0], 4) +
         2 * x[2] * (x[9] * l * l * x[3] + 8 * x[1] * x[4]) * pow(x[0], 3) +
         6 *
             (-l * l * (x[2] * x[7] - x[3] * x[8]) * x[1] * x[1] / 3 -
              l * l * (x[2] * x[2] - 2 * x[3] * x[3]) * x[9] / 3 + g * x[3]) *
             x[1] * x[0] * x[0] -
         6 *
             (l * l * (x[8] * x[2] + x[7] * x[3]) * x[1] * x[1] / 6 -
              0.8e1 / 0.3e1 * x[2] * x[4] * x[1] + x[9] * l * l * x[2] * x[3]) *
             x[1] * x[1] * x[0] +
         8 * pow(x[1], 4) * x[3] * x[4] +
         (2 * l * l * x[2] * x[2] * x[9] + 6 * g * x[3]) * pow(x[1], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(7, 0) =
        -x[9] *
        (-l * l * pow(x[1], 5) * x[8] +
         (-2 * l * l * x[0] * x[7] + 8 * x[4]) * pow(x[1], 4) -
         l * l * pow(x[1], 3) * x[3] * x[9] -
         2 * l * l * x[0] * (x[0] * x[0] * x[7] + 2 * x[2] * x[9]) * x[1] *
             x[1] +
         l * l * x[0] * x[0] * (x[8] * x[0] * x[0] + 3 * x[3] * x[9]) * x[1] -
         8 * pow(x[0], 4) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(7, 1) =
        x[9] *
        (l * l * pow(x[0], 5) * x[8] - 2 * l * l * pow(x[0], 4) * x[1] * x[7] +
         (x[9] * l * l * x[3] + 16 * x[1] * x[4]) * pow(x[0], 3) -
         2 * l * l * x[1] * (x[7] * x[1] * x[1] + x[2] * x[9]) * x[0] * x[0] +
         (-l * l * pow(x[1], 4) * x[8] - 3 * x[9] * l * l * x[3] * x[1] * x[1] +
          16 * pow(x[1], 3) * x[4]) *
             x[0] +
         2 * l * l * pow(x[1], 3) * x[2] * x[9]) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(7, 2) = -l * l * x[1] * x[1] * x[9] * x[9] *
                 pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(7, 3) = l * l * x[1] * x[9] * x[9] * x[0] *
                 pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(7, 4) = -4 * x[0] * x[9] / (x[0] * x[0] + x[1] * x[1]);
    JF_x(7, 5) = 1;
    JF_x(7, 7) =
        -x[1] * x[1] * l * l * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(7, 8) =
        x[1] * l * l * x[0] * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(7, 9) =
        (-pow(x[1], 4) * x[7] * l * l + l * l * x[0] * pow(x[1], 3) * x[8] +
         (-l * l * x[0] * x[0] * x[7] - 2 * l * l * x[2] * x[9] -
          8 * x[0] * x[4]) *
             x[1] * x[1] +
         l * l * x[0] * (x[8] * x[0] * x[0] + 2 * x[3] * x[9]) * x[1] -
         8 * pow(x[0], 3) * x[4]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(8, 0) =
        3 *
        (-l * l * pow(x[0], 4) * x[1] * x[7] / 6 +
         (-x[8] * l * l * x[1] * x[1] / 3 + x[9] * l * l * x[3] / 3 +
          0.8e1 / 0.3e1 * x[1] * x[4] + g) *
             pow(x[0], 3) -
         l * l * x[0] * x[0] * x[1] * x[2] * x[9] / 2 +
         x[1] * x[1] *
             (-x[8] * l * l * x[1] * x[1] / 3 - x[9] * l * l * x[3] / 3 +
              0.8e1 / 0.3e1 * x[1] * x[4] + g) *
             x[0] +
         l * l * pow(x[1], 3) * (x[7] * x[1] * x[1] + x[2] * x[9]) / 6) *
        x[9] * pow(x[0] * x[0] + x[1] * x[1], -3);
    JF_x(8, 1) =
        x[9] *
        (l * l * pow(x[0], 5) * x[7] + 2 * l * l * pow(x[0], 4) * x[1] * x[8] +
         2 * l * l * x[0] * x[0] * pow(x[1], 3) * x[8] -
         l * l * x[0] * pow(x[1], 4) * x[7] +
         l * l * pow(x[0], 3) * x[2] * x[9] +
         4 * l * l * x[0] * x[0] * x[1] * x[3] * x[9] -
         3 * l * l * x[0] * x[1] * x[1] * x[2] * x[9] -
         8 * pow(x[0], 4) * x[4] + 8 * pow(x[1], 4) * x[4] +
         6 * g * x[0] * x[0] * x[1] + 6 * g * pow(x[1], 3)) *
        pow(x[0] * x[0] + x[1] * x[1], -3) / 2;
    JF_x(8, 2) = l * l * x[1] * x[9] * x[9] * x[0] *
                 pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(8, 3) = -x[9] * x[9] * l * l * x[0] * x[0] *
                 pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(8, 4) = -4 * x[1] * x[9] / (x[0] * x[0] + x[1] * x[1]);
    JF_x(8, 6) = 1;
    JF_x(8, 7) =
        x[1] * l * l * x[0] * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(8, 8) =
        -l * l * x[0] * x[0] * x[9] / (2 * x[0] * x[0] + 2 * x[1] * x[1]);
    JF_x(8, 9) =
        (-pow(x[0], 4) * x[8] * l * l + x[1] * x[7] * pow(x[0], 3) * l * l +
         (-x[8] * l * l * x[1] * x[1] - 2 * x[9] * l * l * x[3] -
          8 * x[1] * x[4] - 3 * g) *
             x[0] * x[0] +
         (l * l * pow(x[1], 3) * x[7] + 2 * l * l * x[1] * x[2] * x[9]) * x[0] -
         8 * pow(x[1], 3) * x[4] - 3 * g * x[1] * x[1]) *
        pow(x[0] * x[0] + x[1] * x[1], -2) / 2;
    JF_x(9, 0) = (-2 * pow(x[1], 4) * x[7] +
                  (-4 * x[0] * x[0] * x[7] - 4 * x[2] * x[9]) * x[1] * x[1] +
                  8 * x[3] * x[0] * x[1] * x[9] - 2 * pow(x[0], 4) * x[7] +
                  4 * x[0] * x[0] * x[2] * x[9]) *
                 pow(x[0] * x[0] + x[1] * x[1], -2);
    JF_x(9, 1) = (-2 * pow(x[1], 4) * x[8] +
                  (-4 * x[8] * x[0] * x[0] + 4 * x[3] * x[9]) * x[1] * x[1] +
                  8 * x[0] * x[1] * x[2] * x[9] - 2 * pow(x[0], 4) * x[8] -
                  4 * x[0] * x[0] * x[3] * x[9]) *
                 pow(x[0] * x[0] + x[1] * x[1], -2);
    JF_x(9, 2) = -4 * x[0] * x[9] / (x[0] * x[0] + x[1] * x[1]);
    JF_x(9, 3) = -4 * x[1] * x[9] / (x[0] * x[0] + x[1] * x[1]);
    JF_x(9, 7) = -2 * x[0];
    JF_x(9, 8) = -2 * x[1];
    JF_x(9, 9) =
        (-4 * x[0] * x[2] - 4 * x[1] * x[3]) / (x[0] * x[0] + x[1] * x[1]);
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    MatrixJF JF_x_dot;
    JF_x_dot.setIdentity();
    return JF_x_dot;
  }

  VectorH h(const VectorF &x, const Real /*t*/) const override {
    const Real g = this->m_g;
    VectorH h;
    h[0] = -2 * x[0] * x[2] - 2 * x[1] * x[3];
    h[1] = (4 * x[0] * x[0] + 4 * x[1] * x[1]) * x[4] + 2 * g * x[1] -
           2 * x[2] * x[2] - 2 * x[3] * x[3];
    return h;
  }

  MatrixJH Jh_x(const VectorF &x, const Real /*t*/) const override {
    const Real g = this->m_g;
    MatrixJH Jh_x;
    Jh_x.setZero();
    Jh_x(0, 0) = -2 * x[2];
    Jh_x(0, 1) = -2 * x[3];
    Jh_x(0, 2) = -2 * x[0];
    Jh_x(0, 3) = -2 * x[1];
    Jh_x(1, 0) = 8 * x[0] * x[4];
    Jh_x(1, 1) = 8 * x[1] * x[4] + 2 * g;
    Jh_x(1, 2) = -4 * x[2];
    Jh_x(1, 3) = -4 * x[3];
    Jh_x(1, 4) = 4 * x[0] * x[0] + 4 * x[1] * x[1];
    return Jh_x;
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }

};  // PendulumOCPindex0

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class PendulumOCPProblem
    : public BoundaryValueProblem<Real, System::N, System::M, Integrator> {
 public:
  using typename BoundaryValueProblem<Real, System::N, System::M, Integrator>::
      SystemPtr;
  using typename BoundaryValueProblem<Real, System::N, System::M, Integrator>::
      IntegratorPtr;
  using typename BoundaryValueProblem<Real, System::N, System::M, Integrator>::
      SolutionPtr;
  using typename BoundaryValueProblem<Real, System::N, System::M, Integrator>::
      VectorF;
  using typename BoundaryValueProblem<Real, System::N, System::M, Integrator>::
      MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, System::N, Eigen::Dynamic>;

  PendulumOCPProblem()
      : BoundaryValueProblem<Real, System::N, System::M, Integrator>(
            "PendulumOCPProblem",
            std::make_unique<System>(),
            std::make_unique<Integrator>()) {}

  ~PendulumOCPProblem() {}

  Real g() const {
    return static_cast<const System *>(this->integrator()->system())->g();
  }

  Real l() const {
    return static_cast<const System *>(this->integrator()->system())->l();
  }

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 1.0;
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    const Real l =
        static_cast<const System *>(this->integrator()->system())->l();
    VectorF b;
    b[0] = x_ini[0] - l;
    b[1] = x_ini[1];
    b[2] = x_end[0] + l;
    b[3] = x_end[1];
    b[4] = x_ini[7];
    b[5] = x_ini[8];
    b[6] = x_ini[9];
    b[7] = -x_end[7];
    b[8] = -x_end[8];
    b[9] = -x_end[9];
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini;
    Jb_x_ini.setZero();
    Jb_x_ini(0, 0) = 1;
    Jb_x_ini(1, 1) = 1;
    Jb_x_ini(4, 7) = 1;
    Jb_x_ini(5, 8) = 1;
    Jb_x_ini(6, 9) = 1;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end;
    Jb_x_end.setZero();
    Jb_x_end(2, 0) = 1;
    Jb_x_end(3, 1) = 1;
    Jb_x_end(7, 7) = -1;
    Jb_x_end(8, 8) = -1;
    Jb_x_end(9, 9) = -1;
    return Jb_x_end;
  }

  static VectorF guess(const Real t) {
    Real l = 1.0;
    Real g = 9.81;

    // Angle theta is defined such that theta=0 corresponds to the initial
    // position (l, 0) and theta=-pi corresponds to the final position (-l, 0).
    const Real theta{-M_PI * t / (time_end() - time_start())};
    const Real theta_dot{-M_PI / (time_end() - time_start())};

    // Cartesian coordinates
    const Real x{l * std::cos(theta)};
    const Real y{l * std::sin(theta)};
    const Real u{-l * std::sin(theta) * theta_dot};
    const Real v{l * std::cos(theta) * theta_dot};
    const Real lambda{-(-u * u + y * g - v * v) / (2 * (x * x + y * y))};

    // Adjoints lagrangian multipliers

    // Assemble the guess vector
    VectorF guess;
    guess << x, y, u, v, lambda, 0.0, 0.0, 0.0, 0.0, 0.0;
    return guess;
  }

  static MatrixX guess(const VectorX &t) {
    MatrixX out(System::N, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      out.col(i) = guess(t[i]);
    }
    return out;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PENDULUM_OCP_HH
