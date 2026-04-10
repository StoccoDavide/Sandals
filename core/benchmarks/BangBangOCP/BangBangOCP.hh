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

#include "Sandals/System/BoundaryValueProblem.hh"
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
  Real m_mu{10.0};
  Real m_h{0.001};

 public:
  BangBangOCPImplicit() : Implicit<Real, 6, 0>("BangBangOCPImplicit") {}

  ~BangBangOCPImplicit() {}

  Real u(const VectorF &x) const {
    return sin(atan(x[2] * x[5] / this->m_h));
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    const Real &mu{this->m_mu};
    VectorF F;
    F[0] = -x[2] * x[1] + x_dot[0];
    F[1] = x_dot[1] - x[2] * this->u(x);
    F[2] = x_dot[2];
    F[3] = x_dot[3];
    F[4] = x[2] * x[3] + x_dot[4];
    F[5] = x_dot[5] + x[1] * x[3] + this->u(x) * x[4] -
           2 * mu * (-2 + x[2]) / x[2] / (-4 + x[2]);
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    const Real &mu{this->m_mu};
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 1) = -x[2];
    JF_x(0, 2) = -x[1];
    JF_x(1, 2) = -u(x);
    JF_x(4, 2) = x[3];
    JF_x(4, 3) = x[2];
    JF_x(5, 1) = x[3];
    JF_x(5, 2) = 2 * mu * (x[2] * x[2] - 4 * x[2] + 8) * pow(x[2], -2) *
                 pow(-4 + x[2], -2);
    JF_x(5, 3) = x[1];
    JF_x(5, 4) = this->u(x);
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
class BangBangOCPProblem : public BoundaryValueProblem<Real, 6, 0, Integrator> {
 public:
  using typename BoundaryValueProblem<Real, 6, 0, Integrator>::SystemPtr;
  using typename BoundaryValueProblem<Real, 6, 0, Integrator>::IntegratorPtr;
  using typename BoundaryValueProblem<Real, 6, 0, Integrator>::SolutionPtr;
  using typename BoundaryValueProblem<Real, 6, 0, Integrator>::VectorF;
  using typename BoundaryValueProblem<Real, 6, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 6, Eigen::Dynamic>;

  BangBangOCPProblem()
      : BoundaryValueProblem<Real, 6, 0, Integrator>(
            "BangBangOCPProblem",
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
    guess << t, t, 2.0, 0.0, 0.0, 0.0;
    return guess;
  }

  MatrixX guess(const VectorX &t) {
    MatrixX guess_vec(6, t.size());
    for (int i = 0; i < t.size(); ++i) {
      guess_vec.col(i) = guess(t[i]);
    }
    return guess_vec;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // BENCHMARK_BANG_BANG_OCP_HH
