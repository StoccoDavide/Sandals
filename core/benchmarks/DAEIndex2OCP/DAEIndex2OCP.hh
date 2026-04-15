/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_DAE_INDEX_2_OCP_HH
#define SANDALS_DAE_INDEX_2_OCP_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;
using std::pow;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class DAEIndex2OCPindex2 : public Implicit<Real, 8, 0> {
 public:
  constexpr static Integer N = 8;
  constexpr static Integer M = 0;
  using typename Implicit<Real, 8, 0>::VectorF;
  using typename Implicit<Real, 8, 0>::MatrixJF;
  using typename Implicit<Real, 8, 0>::VectorH;
  using typename Implicit<Real, 8, 0>::MatrixJH;

 private:
  Real m_g{9.81};

 public:
  DAEIndex2OCPindex2() : Implicit<Real, 8, 0>("DAEIndex2OCPindex2") {}

  ~DAEIndex2OCPindex2() {}

  void g(const Real g) {
    this->m_g = g;
  }

  Real g() const {
    return this->m_g;
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F[0] = x_dot[0] + x[4] + x[5] + x[3];
    F[1] = x_dot[1] + x[4] + x[5];
    F[2] = x_dot[2] + x[1];
    F[3] = -x[0] - x[2];
    F[4] = x_dot[4] + x[7];
    F[5] = x_dot[5] - x[6];
    F[6] = x_dot[6] + x[7];
    F[7] = x_dot[7] - x[4];
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 3) = 1;
    JF_x(0, 4) = 1;
    JF_x(0, 5) = 1;
    JF_x(1, 4) = 1;
    JF_x(1, 5) = 1;
    JF_x(2, 1) = 1;
    JF_x(3, 0) = -1;
    JF_x(3, 2) = -1;
    JF_x(4, 7) = 1;
    JF_x(5, 6) = -1;
    JF_x(6, 7) = 1;
    JF_x(7, 4) = -1;
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
    JF_x_dot(4, 4) = 1;
    JF_x_dot(5, 5) = 1;
    JF_x_dot(6, 6) = 1;
    JF_x_dot(7, 7) = 1;
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

};  // DAEIndex2OCPindex2

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class DAEIndex2OCPindex0 : public Implicit<Real, 8, 2> {
 public:
  constexpr static Integer N = 8;
  constexpr static Integer M = 2;
  using typename Implicit<Real, 8, 2>::VectorF;
  using typename Implicit<Real, 8, 2>::MatrixJF;
  using typename Implicit<Real, 8, 2>::VectorH;
  using typename Implicit<Real, 8, 2>::MatrixJH;

 public:
  DAEIndex2OCPindex0() : Implicit<Real, 8, 2>("DAEIndex2OCPindex0") {}

  ~DAEIndex2OCPindex0() {}

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F[0] = x_dot[0] + x[4] + x[5] - x[7] + x[3];
    F[1] = x_dot[1] + x[4] + x[5] - x[7];
    F[2] = x_dot[2] + x[1];
    F[3] = x_dot[3] + x[6] - 2 * x[4] - x[5] + x[7];
    F[4] = x_dot[4];
    F[5] = x_dot[5] - x[6];
    F[6] = x_dot[6];
    F[7] = x_dot[7] - x[4];
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 3) = 1;
    JF_x(0, 4) = 1;
    JF_x(0, 5) = 1;
    JF_x(0, 7) = -1;
    JF_x(1, 4) = 1;
    JF_x(1, 5) = 1;
    JF_x(1, 7) = -1;
    JF_x(2, 1) = 1;
    JF_x(3, 4) = -2;
    JF_x(3, 5) = -1;
    JF_x(3, 6) = 1;
    JF_x(3, 7) = 1;
    JF_x(5, 6) = -1;
    JF_x(7, 4) = -1;
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }

  VectorH h(const VectorF &x, const Real /*t*/) const override {
    VectorH h;
    h[0] = x[0] + x[2];
    h[1] = -x[4] - x[5] + x[7] - x[3] - x[1];
    return h;
  }

  MatrixJH Jh_x(const VectorF & /*x*/, const Real /*t*/) const override {
    MatrixJH Jh_x;
    Jh_x.setZero();
    Jh_x(0, 0) = 1;
    Jh_x(0, 2) = 1;
    Jh_x(1, 1) = -1;
    Jh_x(1, 3) = -1;
    Jh_x(1, 4) = -1;
    Jh_x(1, 5) = -1;
    Jh_x(1, 7) = 1;
    return Jh_x;
  }

  bool in_domain(const VectorF & /*x*/, const Real /*t*/) const override {
    return true;
  }

};  // DAEIndex2OCPindex0

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class DAEIndex2OCPProblem : public BVP<Real, System::N, System::M, Integrator> {
 public:
  using typename BVP<Real, System::N, System::M, Integrator>::VectorF;
  using typename BVP<Real, System::N, System::M, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, System::N, Eigen::Dynamic>;

  DAEIndex2OCPProblem()
      : BVP<Real, System::N, System::M, Integrator>(
            "DAEIndex2OCPProblem",
            std::make_unique<System>(),
            std::make_unique<Integrator>()) {}

  ~DAEIndex2OCPProblem() {}

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 1.0;
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b[0] = x_ini[0];
    b[1] = x_end[0];
    b[2] = x_ini[1] - 1;
    b[3] = x_end[1] + 1;
    b[4] = x_ini[2];
    b[5] = x_ini[7];
    b[6] = -x_end[6];
    b[7] = -x_end[7];
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini;
    Jb_x_ini.setZero();
    Jb_x_ini(0, 0) = 1;
    Jb_x_ini(2, 1) = 1;
    Jb_x_ini(4, 2) = 1;
    Jb_x_ini(5, 7) = 1;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end;
    Jb_x_end.setZero();
    Jb_x_end(1, 0) = 1;
    Jb_x_end(3, 1) = 1;
    Jb_x_end(6, 6) = -1;
    Jb_x_end(7, 7) = -1;
    return Jb_x_end;
  }

  VectorF guess(const Real /*t*/) {
    return VectorF::Constant(0.0);
  }

  MatrixX guess(const VectorX &t) {
    MatrixX guess_vec(System::N, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      guess_vec.col(i) = this->guess(t[i]);
    }
    return guess_vec;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // SANDALS_DAE_INDEX_2_OCP_HH
