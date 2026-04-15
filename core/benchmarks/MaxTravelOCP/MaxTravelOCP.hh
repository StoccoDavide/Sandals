/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_BENCHMARKS_MAX_TRAVEL_OCP_HH
#define SANDALS_BENCHMARKS_MAX_TRAVEL_OCP_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class MaxTravelOCPImplicit : public Implicit<Real, 4, 0> {
 public:
  using typename Implicit<Real, 4, 0>::VectorF;
  using typename Implicit<Real, 4, 0>::MatrixJF;
  using typename Implicit<Real, 4, 0>::VectorH;
  using typename Implicit<Real, 4, 0>::MatrixJH;

 private:
  Real m_m{1.0};
  Real m_h{0.001};

 public:
  MaxTravelOCPImplicit() : Implicit<Real, 4, 0>("MaxTravelOCPImplicit") {}

  ~MaxTravelOCPImplicit() {}

  Real f(const VectorF &x) const {
    // return (x[3] < 0.0) ? 1.0 : -1.0;
    return -sin(atan(x[3] / this->m_h));
  }

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    const Real &m{this->m_m};
    VectorF F;
    F[0] = x_dot[0] - x[1];
    F[1] = x_dot[1] - this->f(x) / m;
    F[2] = x_dot[2];
    F[3] = x_dot[3] + x[2];
    return F;
  }

  MatrixJF JF_x(const VectorF & /*x*/,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x;
    JF_x.setZero();
    JF_x(0, 1) = -1;
    JF_x(3, 2) = 1;
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
class MaxTravelOCPProblem : public BVP<Real, 4, 0, Integrator> {
 public:
  using typename BVP<Real, 4, 0, Integrator>::VectorF;
  using typename BVP<Real, 4, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 4, Eigen::Dynamic>;

  MaxTravelOCPProblem()
      : BVP<Real, 4, 0, Integrator>("MaxTravelOCPProblem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}

  ~MaxTravelOCPProblem() {}

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
    b[2] = x_end[1];
    b[3] = -1 - x_end[2];
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
    Jb_x_end(2, 1) = 1;
    Jb_x_end(3, 2) = -1;
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

#endif  // SANDALS_BENCHMARKS_MAX_TRAVEL_OCP_HH
