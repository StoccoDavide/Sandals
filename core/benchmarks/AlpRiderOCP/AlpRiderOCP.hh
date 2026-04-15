/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_BENCHMARKS_ALP_RIDER_OCP_HH
#define SANDALS_BENCHMARKS_ALP_RIDER_OCP_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class AlpRiderOCPImplicit : public Implicit<Real, 8, 0> {
 public:
  using typename Implicit<Real, 8, 0>::VectorF;
  using typename Implicit<Real, 8, 0>::MatrixJF;
  using typename Implicit<Real, 8, 0>::VectorH;
  using typename Implicit<Real, 8, 0>::MatrixJH;

 private:
  Real m_mu{0.1};

 public:
  AlpRiderOCPImplicit() : Implicit<Real, 8, 0>("AlpRiderOCPImplicit") {}

  AlpRiderOCPImplicit(const Real mu)
      : Implicit<Real, 8, 0>("AlpRiderOCPImplicit"), m_mu(mu) {}

  ~AlpRiderOCPImplicit() {}

  void mu(const Real mu) {
    this->m_mu = mu;
  }

  Real mu() const {
    return this->m_mu;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real t) const override {
    const Real &mu{this->m_mu};
    VectorF F;
    F[0] = x_dot[0] + 10 * x[0] + 100 * x[4] + 150 * x[5] + 200 * x[7];
    F[1] =
        x_dot[1] + 2 * x[1] + 150 * x[4] + 250 * x[5] - 50 * x[6] + 350 * x[7];
    F[2] = x_dot[2] + 3 * x[2] - 5 * x[3] - 50 * x[5] + 100 * x[6] - 100 * x[7];
    F[3] = x_dot[3] - 5 * x[2] + 3 * x[3] + 200 * x[4] + 350 * x[5] -
           100 * x[6] + 500 * x[7];
    F[4] =
        ((160000 * x[0] - 8000 * x[4] + 800 * x_dot[4]) *
             exp(-4 * pow(t - 15, 2)) +
         (60000 * x[0] - 3000 * x[4] + 300 * x_dot[4]) *
             exp(-6 * pow(t - 10, 2)) +
         (60000 * x[0] - 3000 * x[4] + 300 * x_dot[4]) *
             exp(-10 * pow(t - 6, 2)) +
         (60000 * x[0] - 3000 * x[4] + 300 * x_dot[4]) *
             exp(-12 * pow(t - 3, 2)) -
         20000 * pow(x[0], 3) + (1000 * x[4] - 100 * x_dot[4]) * x[0] * x[0] +
         (-20000 * x[1] * x[1] - 20000 * x[2] * x[2] - 20000 * x[3] * x[3] +
          200 * mu + 200) *
             x[0] +
         1000 * (x[4] - x_dot[4] / 10) *
             (x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 0.1e1 / 0.100e3)) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    F[5] = ((160000 * x[1] - 1600 * x[5] + 800 * x_dot[5]) *
                exp(-4 * pow(t - 15, 2)) +
            (60000 * x[1] - 600 * x[5] + 300 * x_dot[5]) *
                exp(-6 * pow(t - 10, 2)) +
            (60000 * x[1] - 600 * x[5] + 300 * x_dot[5]) *
                exp(-10 * pow(t - 6, 2)) +
            (60000 * x[1] - 600 * x[5] + 300 * x_dot[5]) *
                exp(-12 * pow(t - 3, 2)) -
            20000 * pow(x[1], 3) + (200 * x[5] - 100 * x_dot[5]) * x[1] * x[1] +
            (-20000 * x[0] * x[0] - 20000 * x[2] * x[2] - 20000 * x[3] * x[3] +
             200 * mu + 200) *
                x[1] +
            200 * (x[5] - x_dot[5] / 2) *
                (x[0] * x[0] + x[2] * x[2] + x[3] * x[3] - 0.1e1 / 0.100e3)) /
           (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
            100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
            300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
            800 * exp(-4 * pow(t - 15, 2)) + 1);
    F[6] = ((160000 * x[2] - 2400 * x[6] + 4000 * x[7] + 800 * x_dot[6]) *
                exp(-4 * pow(t - 15, 2)) +
            (60000 * x[2] - 900 * x[6] + 1500 * x[7] + 300 * x_dot[6]) *
                exp(-6 * pow(t - 10, 2)) +
            (60000 * x[2] - 900 * x[6] + 1500 * x[7] + 300 * x_dot[6]) *
                exp(-10 * pow(t - 6, 2)) +
            (60000 * x[2] - 900 * x[6] + 1500 * x[7] + 300 * x_dot[6]) *
                exp(-12 * pow(t - 3, 2)) -
            20000 * pow(x[2], 3) +
            (300 * x[6] - 500 * x[7] - 100 * x_dot[6]) * x[2] * x[2] +
            (-20000 * x[0] * x[0] - 20000 * x[1] * x[1] - 20000 * x[3] * x[3] +
             200 * mu + 200) *
                x[2] +
            300 * (x[6] - 0.5e1 / 0.3e1 * x[7] - x_dot[6] / 3) *
                (x[0] * x[0] + x[1] * x[1] + x[3] * x[3] - 0.1e1 / 0.100e3)) /
           (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
            100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
            300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
            800 * exp(-4 * pow(t - 15, 2)) + 1);
    F[7] = ((160000 * x[3] + 4000 * x[6] - 2400 * x[7] + 800 * x_dot[7]) *
                exp(-4 * pow(t - 15, 2)) +
            (60000 * x[3] + 1500 * x[6] - 900 * x[7] + 300 * x_dot[7]) *
                exp(-6 * pow(t - 10, 2)) +
            (60000 * x[3] + 1500 * x[6] - 900 * x[7] + 300 * x_dot[7]) *
                exp(-10 * pow(t - 6, 2)) +
            (60000 * x[3] + 1500 * x[6] - 900 * x[7] + 300 * x_dot[7]) *
                exp(-12 * pow(t - 3, 2)) -
            20000 * pow(x[3], 3) +
            (-500 * x[6] + 300 * x[7] - 100 * x_dot[7]) * x[3] * x[3] +
            (-20000 * x[0] * x[0] - 20000 * x[1] * x[1] - 20000 * x[2] * x[2] +
             200 * mu + 200) *
                x[3] -
            500 * (x[6] - 0.3e1 / 0.5e1 * x[7] + x_dot[7] / 5) *
                (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 0.1e1 / 0.100e3)) /
           (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
            100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
            300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
            800 * exp(-4 * pow(t - 15, 2)) + 1);
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF &x_dot,
                const Real t) const override {
    const Real &mu{this->m_mu};
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 0) = 10;
    JF_x(0, 4) = 100;
    JF_x(0, 5) = 150;
    JF_x(0, 7) = 200;
    JF_x(1, 1) = 2;
    JF_x(1, 4) = 150;
    JF_x(1, 5) = 250;
    JF_x(1, 6) = -50;
    JF_x(1, 7) = 350;
    JF_x(2, 2) = 3;
    JF_x(2, 3) = -5;
    JF_x(2, 5) = -50;
    JF_x(2, 6) = 100;
    JF_x(2, 7) = -100;
    JF_x(3, 2) = -5;
    JF_x(3, 3) = 3;
    JF_x(3, 4) = 200;
    JF_x(3, 5) = 350;
    JF_x(3, 6) = -100;
    JF_x(3, 7) = 500;
    JF_x(4, 0) =
        (-160000 * exp(-4 * pow(t - 15, 2)) - 60000 * exp(-6 * pow(t - 10, 2)) -
         60000 * exp(-10 * pow(t - 6, 2)) - 60000 * exp(-12 * pow(t - 3, 2)) +
         60000 * x[0] * x[0] + (-2000 * x[4] + 200 * x_dot[4]) * x[0] +
         20000 * x[1] * x[1] + 20000 * x[2] * x[2] + 20000 * x[3] * x[3] -
         200 * mu - 200) /
            (100 * x[0] * x[0] + 100 * x[1] * x[1] + 100 * x[2] * x[2] +
             100 * x[3] * x[3] - 800 * exp(-4 * pow(t - 15, 2)) -
             300 * exp(-6 * pow(t - 10, 2)) - 300 * exp(-10 * pow(t - 6, 2)) -
             300 * exp(-12 * pow(t - 3, 2)) - 1) +
        4 * x[0] *
            ((800 * x[0] + 4 * x_dot[4] - 40 * x[4]) *
                 exp(-4 * pow(t - 15, 2)) +
             (300 * x[0] + 0.3e1 / 0.2e1 * x_dot[4] - 15 * x[4]) *
                 exp(-6 * pow(t - 10, 2)) +
             (300 * x[0] + 0.3e1 / 0.2e1 * x_dot[4] - 15 * x[4]) *
                 exp(-10 * pow(t - 6, 2)) +
             (300 * x[0] + 0.3e1 / 0.2e1 * x_dot[4] - 15 * x[4]) *
                 exp(-12 * pow(t - 3, 2)) -
             100 * pow(x[0], 3) + (5 * x[4] - x_dot[4] / 2) * x[0] * x[0] +
             (-100 * x[1] * x[1] - 100 * x[2] * x[2] - 100 * x[3] * x[3] + mu +
              1) *
                 x[0] +
             5 * (x[4] - x_dot[4] / 10) *
                 (x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 0.1e1 / 0.100e3)) *
            pow(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] -
                    3 * exp(-12 * pow(t - 3, 2)) -
                    3 * exp(-10 * pow(t - 6, 2)) -
                    3 * exp(-6 * pow(t - 10, 2)) -
                    8 * exp(-4 * pow(t - 15, 2)) - 0.1e1 / 0.100e3,
                -2);
    JF_x(4, 1) =
        40000 * x[1] * mu * x[0] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(4, 2) =
        40000 * x[2] * mu * x[0] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(4, 3) =
        40000 * x[3] * mu * x[0] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(4, 4) =
        (-8000 * exp(-4 * pow(t - 15, 2)) - 3000 * exp(-6 * pow(t - 10, 2)) -
         3000 * exp(-10 * pow(t - 6, 2)) - 3000 * exp(-12 * pow(t - 3, 2)) +
         1000 * x[0] * x[0] + 1000 * x[1] * x[1] + 1000 * x[2] * x[2] +
         1000 * x[3] * x[3] - 10) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    JF_x(5, 0) =
        40000 * x[1] * mu * x[0] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(5, 1) =
        (-160000 * exp(-4 * pow(t - 15, 2)) - 60000 * exp(-6 * pow(t - 10, 2)) -
         60000 * exp(-10 * pow(t - 6, 2)) - 60000 * exp(-12 * pow(t - 3, 2)) +
         60000 * x[1] * x[1] + (-400 * x[5] + 200 * x_dot[5]) * x[1] +
         20000 * x[0] * x[0] + 20000 * x[2] * x[2] + 20000 * x[3] * x[3] -
         200 * mu - 200) /
            (100 * x[0] * x[0] + 100 * x[1] * x[1] + 100 * x[2] * x[2] +
             100 * x[3] * x[3] - 800 * exp(-4 * pow(t - 15, 2)) -
             300 * exp(-6 * pow(t - 10, 2)) - 300 * exp(-10 * pow(t - 6, 2)) -
             300 * exp(-12 * pow(t - 3, 2)) - 1) +
        4 * x[1] *
            ((800 * x[1] + 4 * x_dot[5] - 8 * x[5]) * exp(-4 * pow(t - 15, 2)) +
             (300 * x[1] + 0.3e1 / 0.2e1 * x_dot[5] - 3 * x[5]) *
                 exp(-6 * pow(t - 10, 2)) +
             (300 * x[1] + 0.3e1 / 0.2e1 * x_dot[5] - 3 * x[5]) *
                 exp(-10 * pow(t - 6, 2)) +
             (300 * x[1] + 0.3e1 / 0.2e1 * x_dot[5] - 3 * x[5]) *
                 exp(-12 * pow(t - 3, 2)) -
             100 * pow(x[1], 3) + (x[5] - x_dot[5] / 2) * x[1] * x[1] +
             (-100 * x[0] * x[0] - 100 * x[2] * x[2] - 100 * x[3] * x[3] + mu +
              1) *
                 x[1] +
             (x[5] - x_dot[5] / 2) *
                 (x[0] * x[0] + x[2] * x[2] + x[3] * x[3] - 0.1e1 / 0.100e3)) *
            pow(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] -
                    3 * exp(-12 * pow(t - 3, 2)) -
                    3 * exp(-10 * pow(t - 6, 2)) -
                    3 * exp(-6 * pow(t - 10, 2)) -
                    8 * exp(-4 * pow(t - 15, 2)) - 0.1e1 / 0.100e3,
                -2);
    JF_x(5, 2) =
        40000 * x[2] * mu * x[1] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(5, 3) =
        40000 * x[3] * mu * x[1] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(5, 5) =
        (-1600 * exp(-4 * pow(t - 15, 2)) - 600 * exp(-6 * pow(t - 10, 2)) -
         600 * exp(-10 * pow(t - 6, 2)) - 600 * exp(-12 * pow(t - 3, 2)) +
         200 * x[1] * x[1] + 200 * x[0] * x[0] + 200 * x[2] * x[2] +
         200 * x[3] * x[3] - 2) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    JF_x(6, 0) =
        40000 * x[2] * mu * x[0] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(6, 1) =
        40000 * x[2] * mu * x[1] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(6, 2) =
        (-160000 * exp(-4 * pow(t - 15, 2)) - 60000 * exp(-6 * pow(t - 10, 2)) -
         60000 * exp(-10 * pow(t - 6, 2)) - 60000 * exp(-12 * pow(t - 3, 2)) +
         60000 * x[2] * x[2] +
         (-600 * x[6] + 1000 * x[7] + 200 * x_dot[6]) * x[2] +
         20000 * x[0] * x[0] + 20000 * x[1] * x[1] + 20000 * x[3] * x[3] -
         200 * mu - 200) /
            (100 * x[0] * x[0] + 100 * x[1] * x[1] + 100 * x[2] * x[2] +
             100 * x[3] * x[3] - 800 * exp(-4 * pow(t - 15, 2)) -
             300 * exp(-6 * pow(t - 10, 2)) - 300 * exp(-10 * pow(t - 6, 2)) -
             300 * exp(-12 * pow(t - 3, 2)) - 1) +
        4 *
            ((800 * x[2] - 12 * x[6] + 20 * x[7] + 4 * x_dot[6]) *
                 exp(-4 * pow(t - 15, 2)) +
             (300 * x[2] - 0.9e1 / 0.2e1 * x[6] + 0.15e2 / 0.2e1 * x[7] +
              0.3e1 / 0.2e1 * x_dot[6]) *
                 exp(-6 * pow(t - 10, 2)) +
             (300 * x[2] - 0.9e1 / 0.2e1 * x[6] + 0.15e2 / 0.2e1 * x[7] +
              0.3e1 / 0.2e1 * x_dot[6]) *
                 exp(-10 * pow(t - 6, 2)) +
             (300 * x[2] - 0.9e1 / 0.2e1 * x[6] + 0.15e2 / 0.2e1 * x[7] +
              0.3e1 / 0.2e1 * x_dot[6]) *
                 exp(-12 * pow(t - 3, 2)) -
             100 * pow(x[2], 3) +
             (0.3e1 / 0.2e1 * x[6] - 0.5e1 / 0.2e1 * x[7] - x_dot[6] / 2) *
                 x[2] * x[2] +
             (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[3] * x[3] + mu +
              1) *
                 x[2] +
             0.3e1 / 0.2e1 * (x[6] - 0.5e1 / 0.3e1 * x[7] - x_dot[6] / 3) *
                 (x[0] * x[0] + x[1] * x[1] + x[3] * x[3] - 0.1e1 / 0.100e3)) *
            x[2] *
            pow(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] -
                    3 * exp(-12 * pow(t - 3, 2)) -
                    3 * exp(-10 * pow(t - 6, 2)) -
                    3 * exp(-6 * pow(t - 10, 2)) -
                    8 * exp(-4 * pow(t - 15, 2)) - 0.1e1 / 0.100e3,
                -2);
    JF_x(6, 3) =
        40000 * x[3] * mu * x[2] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(6, 6) =
        (-2400 * exp(-4 * pow(t - 15, 2)) - 900 * exp(-6 * pow(t - 10, 2)) -
         900 * exp(-10 * pow(t - 6, 2)) - 900 * exp(-12 * pow(t - 3, 2)) +
         300 * x[2] * x[2] + 300 * x[0] * x[0] + 300 * x[1] * x[1] +
         300 * x[3] * x[3] - 3) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    JF_x(6, 7) =
        (4000 * exp(-4 * pow(t - 15, 2)) + 1500 * exp(-6 * pow(t - 10, 2)) +
         1500 * exp(-10 * pow(t - 6, 2)) + 1500 * exp(-12 * pow(t - 3, 2)) -
         500 * x[2] * x[2] - 500 * x[0] * x[0] - 500 * x[1] * x[1] -
         500 * x[3] * x[3] + 5) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    JF_x(7, 0) =
        40000 * x[3] * mu * x[0] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(7, 1) =
        40000 * x[3] * mu * x[1] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(7, 2) =
        40000 * x[3] * mu * x[2] *
        pow(-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
                100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
                300 * exp(-10 * pow(t - 6, 2)) +
                300 * exp(-6 * pow(t - 10, 2)) +
                800 * exp(-4 * pow(t - 15, 2)) + 1,
            -2);
    JF_x(7, 3) =
        (-160000 * exp(-4 * pow(t - 15, 2)) - 60000 * exp(-6 * pow(t - 10, 2)) -
         60000 * exp(-10 * pow(t - 6, 2)) - 60000 * exp(-12 * pow(t - 3, 2)) +
         60000 * x[3] * x[3] +
         (1000 * x[6] - 600 * x[7] + 200 * x_dot[7]) * x[3] +
         20000 * x[0] * x[0] + 20000 * x[1] * x[1] + 20000 * x[2] * x[2] -
         200 * mu - 200) /
            (100 * x[0] * x[0] + 100 * x[1] * x[1] + 100 * x[2] * x[2] +
             100 * x[3] * x[3] - 800 * exp(-4 * pow(t - 15, 2)) -
             300 * exp(-6 * pow(t - 10, 2)) - 300 * exp(-10 * pow(t - 6, 2)) -
             300 * exp(-12 * pow(t - 3, 2)) - 1) +
        4 *
            ((800 * x[3] + 20 * x[6] - 12 * x[7] + 4 * x_dot[7]) *
                 exp(-4 * pow(t - 15, 2)) +
             (300 * x[3] + 0.15e2 / 0.2e1 * x[6] - 0.9e1 / 0.2e1 * x[7] +
              0.3e1 / 0.2e1 * x_dot[7]) *
                 exp(-6 * pow(t - 10, 2)) +
             (300 * x[3] + 0.15e2 / 0.2e1 * x[6] - 0.9e1 / 0.2e1 * x[7] +
              0.3e1 / 0.2e1 * x_dot[7]) *
                 exp(-10 * pow(t - 6, 2)) +
             (300 * x[3] + 0.15e2 / 0.2e1 * x[6] - 0.9e1 / 0.2e1 * x[7] +
              0.3e1 / 0.2e1 * x_dot[7]) *
                 exp(-12 * pow(t - 3, 2)) -
             100 * pow(x[3], 3) +
             (-0.5e1 / 0.2e1 * x[6] + 0.3e1 / 0.2e1 * x[7] - x_dot[7] / 2) *
                 x[3] * x[3] +
             (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] + mu +
              1) *
                 x[3] -
             0.5e1 / 0.2e1 * (x[6] - 0.3e1 / 0.5e1 * x[7] + x_dot[7] / 5) *
                 (x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 0.1e1 / 0.100e3)) *
            x[3] *
            pow(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] -
                    3 * exp(-12 * pow(t - 3, 2)) -
                    3 * exp(-10 * pow(t - 6, 2)) -
                    3 * exp(-6 * pow(t - 10, 2)) -
                    8 * exp(-4 * pow(t - 15, 2)) - 0.1e1 / 0.100e3,
                -2);
    JF_x(7, 6) =
        (4000 * exp(-4 * pow(t - 15, 2)) + 1500 * exp(-6 * pow(t - 10, 2)) +
         1500 * exp(-10 * pow(t - 6, 2)) + 1500 * exp(-12 * pow(t - 3, 2)) -
         500 * x[2] * x[2] - 500 * x[0] * x[0] - 500 * x[1] * x[1] -
         500 * x[3] * x[3] + 5) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    JF_x(7, 7) =
        (-2400 * exp(-4 * pow(t - 15, 2)) - 900 * exp(-6 * pow(t - 10, 2)) -
         900 * exp(-10 * pow(t - 6, 2)) - 900 * exp(-12 * pow(t - 3, 2)) +
         300 * x[2] * x[2] + 300 * x[0] * x[0] + 300 * x[1] * x[1] +
         300 * x[3] * x[3] - 3) /
        (-100 * x[0] * x[0] - 100 * x[1] * x[1] - 100 * x[2] * x[2] -
         100 * x[3] * x[3] + 300 * exp(-12 * pow(t - 3, 2)) +
         300 * exp(-10 * pow(t - 6, 2)) + 300 * exp(-6 * pow(t - 10, 2)) +
         800 * exp(-4 * pow(t - 15, 2)) + 1);
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
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
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class AlpRiderOCPProblem : public BVP<Real, 8, 0, Integrator> {
 public:
  using typename BVP<Real, 8, 0, Integrator>::SystemPtr;
  using typename BVP<Real, 8, 0, Integrator>::IntegratorPtr;
  using typename BVP<Real, 8, 0, Integrator>::SolutionPtr;
  using typename BVP<Real, 8, 0, Integrator>::VectorF;
  using typename BVP<Real, 8, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 8, Eigen::Dynamic>;

  AlpRiderOCPProblem()
      : BVP<Real, 8, 0, Integrator>("AlpRiderOCPProblem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}

  ~AlpRiderOCPProblem() {}

  Real mu() const {
    return static_cast<const System *>(this->integrator()->system())->mu();
  }

  void mu(const Real mu) {
    static_cast<System *>(this->integrator()->system())->mu(mu);
  }

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 20.0;
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b[0] = x_ini[0] - 2;
    b[1] = x_ini[1] - 1;
    b[2] = x_ini[2] - 2;
    b[3] = x_ini[3] - 1;
    b[4] = x_end[0] - 2;
    b[5] = x_end[1] - 3;
    b[6] = x_end[2] - 1;
    b[7] = x_end[3] + 2;
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini;
    Jb_x_ini.setZero();
    Jb_x_ini(0, 0) = 1;
    Jb_x_ini(1, 1) = 1;
    Jb_x_ini(2, 2) = 1;
    Jb_x_ini(3, 3) = 1;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end;
    Jb_x_end.setZero();
    Jb_x_end(4, 0) = 1;
    Jb_x_end(5, 1) = 1;
    Jb_x_end(6, 2) = 1;
    Jb_x_end(7, 3) = 1;
    return Jb_x_end;
  }

  VectorF guess(const Real /*t*/) {
    return VectorF::Zero();
  }

  MatrixX guess(const VectorX &t) {
    MatrixX guess_vec(8, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      guess_vec.col(i) = this->guess(t[i]);
    }
    return guess_vec;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // SANDALS_BENCHMARKS_ALP_RIDER_OCP_HH
