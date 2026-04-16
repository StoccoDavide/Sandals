/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_BENCHMARKS_RAYLEIGH_OCP_HH
#define SANDALS_BENCHMARKS_RAYLEIGH_OCP_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class RayleighOCPImplicit : public Implicit<Real, 4, 0> {
 public:
  using typename Implicit<Real, 4, 0>::VectorF;
  using typename Implicit<Real, 4, 0>::MatrixJF;
  using typename Implicit<Real, 4, 0>::VectorH;
  using typename Implicit<Real, 4, 0>::MatrixJH;

 private:
  Real m_h{0.01};
  Real m_p{0.14};

 public:
  RayleighOCPImplicit() : Implicit<Real, 4, 0>("RayleighOCPImplicit") {}

  ~RayleighOCPImplicit() {}

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    const Real &h{this->m_h};
    const Real &p{this->m_p};
    VectorF F;
    F[0] = x_dot[0] - x[1];
    F[1] = x_dot[1] + x[0] - 0.14e1 * x[1] + p * pow(x[1], 3) -
           2 * sqrt(4 * x[3] * x[3] + h - 4 * x[3] + 1) +
           2 * sqrt(4 * x[3] * x[3] + h + 4 * x[3] + 1);
    F[2] = x_dot[2] + 2 * x[0] - x[3];
    F[3] = x_dot[3] + x[2] + 0.14e1 * x[3] - 0.3e1 * x[3] * p * x[1] * x[1];
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    const Real &h{this->m_h};
    const Real &p{this->m_p};
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 1) = -1;
    JF_x(1, 0) = 1;
    JF_x(1, 1) = -0.14e1 + 3 * p * x[1] * x[1];
    JF_x(1, 3) = (-8 * x[3] + 4) *
                     pow(4 * x[3] * x[3] + h - 4 * x[3] + 1, -0.1e1 / 0.2e1) +
                 pow(4 * x[3] * x[3] + h + 4 * x[3] + 1, -0.1e1 / 0.2e1) *
                     (8 * x[3] + 4);
    JF_x(2, 0) = 2;
    JF_x(2, 3) = -1;
    JF_x(3, 1) = -0.6e1 * x[3] * p * x[1];
    JF_x(3, 2) = 1;
    JF_x(3, 3) = 0.14e1 - 0.3e1 * p * x[1] * x[1];
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class RayleighOCPProblem : public BVP<Real, 4, 0, Integrator> {
 public:
  using typename BVP<Real, 4, 0, Integrator>::VectorF;
  using typename BVP<Real, 4, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 4, Eigen::Dynamic>;

  RayleighOCPProblem()
      : BVP<Real, 4, 0, Integrator>("RayleighOCPProblem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}

  ~RayleighOCPProblem() {}

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 4.5;
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b[0] = x_ini[0] + 2;
    b[1] = x_ini[1] + 5;
    b[2] = x_end[0];
    b[3] = x_end[1];
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini;
    Jb_x_ini.setZero();
    Jb_x_ini(0, 0) = 1;
    Jb_x_ini(1, 1) = 1;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end;
    Jb_x_end.setZero();
    Jb_x_end(2, 0) = 1;
    Jb_x_end(3, 1) = 1;
    return Jb_x_end;
  }

  VectorF guess(const Real /*t*/) {
    return VectorF::Zero();
  }

  MatrixX guess(const VectorX &t) {
    MatrixX guess_vec(4, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      guess_vec.col(i) = this->guess(t[i]);
    }
    return guess_vec;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // SANDALS_BENCHMARKS_RAYLEIGH_OCP_HH
