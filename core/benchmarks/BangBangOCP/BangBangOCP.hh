/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef BENCHMARK_BANG_BANG_OCP_HH
#define BENCHMARK_BANG_BANG_OCP_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class BangBangOCPImplicit : public Implicit<Real, 6, 0> {
 public:
  using typename Implicit<Real, 6, 0>::VectorF;
  using typename Implicit<Real, 6, 0>::MatrixJF;
  using typename Implicit<Real, 6, 0>::VectorH;
  using typename Implicit<Real, 6, 0>::MatrixJH;

 private:
  Real m_mu{0.1};
  Real m_h{0.1};

 public:
  BangBangOCPImplicit() : Implicit<Real, 6, 0>("BangBangOCPImplicit") {}

  ~BangBangOCPImplicit() {}

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    const Real &mu{this->m_mu};
    const Real &h{this->m_h};
    VectorF F;
    F[0] = -x[2] * x[1] + x_dot[0];
    F[1] =
        (x_dot[1] * sqrt((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2)) * h +
         x[2] * x[2] * x[4]) *
        pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), -0.1e1 / 0.2e1) /
        h;
    F[2] = x_dot[2];
    F[3] = x_dot[3];
    F[4] = x[2] * x[3] + x_dot[4];
    F[5] =
        -2 *
        pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), -0.1e1 / 0.2e1) *
        (((-x[1] * x[3] / 2 - x_dot[5] / 2) * x[2] * x[2] +
          (2 * x[1] * x[3] + mu + 2 * x_dot[5]) * x[2] - 2 * mu) *
             h * sqrt((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2)) +
         x[2] * x[2] * x[4] * x[4] * (-4 + x[2]) / 2) /
        h / x[2] / (-4 + x[2]);
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    const Real &mu{this->m_mu};
    const Real &h{this->m_h};
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 1) = -x[2];
    JF_x(0, 2) = -x[1];
    JF_x(1, 2) =
        pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), -0.3e1 / 0.2e1) *
        (x[2] * x[2] * x[4] * x[4] + 2 * h * h) * x[2] * x[4] * pow(h, -3);
    JF_x(1, 4) =
        x[2] * x[2] / h *
        pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), -0.3e1 / 0.2e1);
    JF_x(4, 2) = x[3];
    JF_x(4, 3) = x[2];
    JF_x(5, 1) = x[3];
    JF_x(5, 2) =
        2 *
        pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), -0.3e1 / 0.2e1) *
        (pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), 0.3e1 / 0.2e1) *
             (x[2] * x[2] - 4 * x[2] + 8) * mu * h -
         x[2] * x[2] * x[4] * x[4] * pow(-4 + x[2], 2) / 2) /
        h * pow(-4 + x[2], -2) * pow(x[2], -2);
    JF_x(5, 3) = x[1];
    JF_x(5, 4) =
        -pow((x[2] * x[2] * x[4] * x[4] + h * h) * pow(h, -2), -0.3e1 / 0.2e1) *
        (x[2] * x[2] * x[4] * x[4] + 2 * h * h) * x[2] * x[4] * pow(h, -3);
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF &x,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    MatrixJF JF_x_dot;
    JF_x_dot(0, 0) = 1;
    JF_x_dot(1, 1) = 1;
    JF_x_dot(2, 2) = 1;
    JF_x_dot(3, 3) = 1;
    JF_x_dot(4, 4) = 1;
    JF_x_dot(5, 5) = -2 * (-x[2] * x[2] / 2 + 2 * x[2]) / x[2] / (-4 + x[2]);
    return JF_x_dot;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class BangBangOCPProblem : public BVP<Real, 6, 0, Integrator> {
 public:
  using typename BVP<Real, 6, 0, Integrator>::VectorF;
  using typename BVP<Real, 6, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 6, Eigen::Dynamic>;

  BangBangOCPProblem()
      : BVP<Real, 6, 0, Integrator>("BangBangOCPProblem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}

  ~BangBangOCPProblem() {}

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 1.0;
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b[0] = x_ini[0];
    b[1] = x_ini[1];
    b[2] = x_end[0] - 1;
    b[3] = x_end[1];
    b[4] = x_ini[5];
    b[5] = -x_end[5];
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini;
    Jb_x_ini.setZero();
    Jb_x_ini(0, 0) = 1;
    Jb_x_ini(1, 1) = 1;
    Jb_x_ini(4, 5) = 1;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end;
    Jb_x_end.setZero();
    Jb_x_end(2, 0) = 1;
    Jb_x_end(3, 1) = 1;
    Jb_x_end(5, 5) = -1;
    return Jb_x_end;
  }

  VectorF guess(const Real t) {
    VectorF guess;
    guess << t, t, 2.0, 0.1, 0.1, 0.1;
    return guess;
  }

  MatrixX guess(const VectorX &t) {
    MatrixX guess_vec(6, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      guess_vec.col(i) = this->guess(t[i]);
    }
    return guess_vec;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // BENCHMARK_BANG_BANG_OCP_HH
