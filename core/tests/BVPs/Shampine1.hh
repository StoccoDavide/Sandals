/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef TESTS_PROBLEMS_SHAMPINE1_HH
#define TESTS_PROBLEMS_SHAMPINE1_HH

#include "Sandals/System/BVP.hh"
#include "Sandals/System/Explicit.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class Shampine1Explicit : public Explicit<Real, 5, 0> {
 public:
  using typename Explicit<Real, 5, 0>::VectorF;
  using typename Explicit<Real, 5, 0>::MatrixJF;
  using typename Explicit<Real, 5, 0>::VectorH;
  using typename Explicit<Real, 5, 0>::MatrixJH;

  Shampine1Explicit() : Explicit<Real, 5, 0>("Shampine1Explicit") {}

  ~Shampine1Explicit() {}

  VectorF f(const VectorF &x, const Real /*t*/) const override {
    VectorF f;
    f << 0.5 * x(0) * (x(2) - x(0)) / x(1), -0.5 * (x(2) - x(0)),
        (0.9 - 1000 * (x(2) - x(4)) - 0.5 * x(2) * (x(2) - x(0))) / x(3),
        0.5 * (x(2) - x(0)), 100 * (x(2) - x(4));
    return f;
  }

  MatrixJF Jf_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 0) = 0.5 * (x(2) - 2.0 * x(0)) / x(1);
    Jf_x(0, 1) = -0.5 * x(0) * (x(2) - x(0)) / std::pow(x(1), 2.0);
    Jf_x(0, 2) = 0.5 * x(0) / x(1);
    Jf_x(1, 0) = 0.5;
    Jf_x(1, 2) = -0.5;
    Real N{0.9 - 1000.0 * (x(2) - x(4)) - 0.5 * std::pow(x(2), 2.0) +
           0.5 * x(0) * x(2)};
    Jf_x(2, 0) = 0.5 * x(2) / x(3);
    Jf_x(2, 2) = (-1000.0 - x(2) + 0.5 * x(0)) / x(3);
    Jf_x(2, 3) = -N / std::pow(x(3), 2.0);
    Jf_x(2, 4) = 1000.0 / x(3);
    Jf_x(3, 0) = -0.5;
    Jf_x(3, 2) = 0.5;
    Jf_x(4, 2) = 100.0;
    Jf_x(4, 4) = -100.0;
    return Jf_x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class Shampine1Implicit : public Implicit<Real, 5, 0> {
 public:
  using typename Implicit<Real, 5, 0>::VectorF;
  using typename Implicit<Real, 5, 0>::MatrixJF;
  using typename Implicit<Real, 5, 0>::VectorH;
  using typename Implicit<Real, 5, 0>::MatrixJH;

  Shampine1Implicit() : Implicit<Real, 5, 0>("Shampine1Implicit") {}

  ~Shampine1Implicit() {}

  VectorF F(const VectorF &x,
            const VectorF &x_dot,
            const Real /*t*/) const override {
    VectorF F;
    F << x_dot(0) - 0.5 * x(0) * (x(2) - x(0)) / x(1),
        x_dot(1) + 0.5 * (x(2) - x(0)),
        x_dot(2) -
            (0.9 - 1000 * (x(2) - x(4)) - 0.5 * x(2) * (x(2) - x(0))) / x(3),
        x_dot(3) - 0.5 * (x(2) - x(0)), x_dot(4) - 100 * (x(2) - x(4));
    return F;
  }

  MatrixJF JF_x(const VectorF &x,
                const VectorF & /*x_dot*/,
                const Real /*t*/) const override {
    MatrixJF JF_x(MatrixJF::Zero());
    JF_x(0, 0) = -0.5 * (x(2) - 2.0 * x(0)) / x(1);
    JF_x(0, 1) = 0.5 * x(0) * (x(2) - x(0)) / std::pow(x(1), 2.0);
    JF_x(0, 2) = -0.5 * x(0) / x(1);
    JF_x(1, 0) = -0.5;
    JF_x(1, 2) = 0.5;
    Real N{0.9 - 1000.0 * (x(2) - x(4)) - 0.5 * std::pow(x(2), 2.0) +
           0.5 * x(0) * x(2)};
    JF_x(2, 0) = -0.5 * x(2) / x(3);
    JF_x(2, 2) = (1000.0 + x(2) - 0.5 * x(0)) / x(3);
    JF_x(2, 3) = N / std::pow(x(3), 2.0);
    JF_x(2, 4) = -1000.0 / x(3);
    JF_x(3, 0) = 0.5;
    JF_x(3, 2) = -0.5;
    JF_x(4, 2) = -100.0;
    JF_x(4, 4) = 100.0;
    return JF_x;
  }

  MatrixJF JF_x_dot(const VectorF & /*x*/,
                    const VectorF & /*x_dot*/,
                    const Real /*t*/) const override {
    return MatrixJF::Identity();
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real = double>
class Shampine1SemiExplicit : public SemiExplicit<Real, 5, 0> {
 public:
  using typename SemiExplicit<Real, 5, 0>::VectorF;
  using typename SemiExplicit<Real, 5, 0>::MatrixA;
  using typename SemiExplicit<Real, 5, 0>::TensorTA;
  using typename SemiExplicit<Real, 5, 0>::VectorB;
  using typename SemiExplicit<Real, 5, 0>::MatrixJB;
  using typename SemiExplicit<Real, 5, 0>::VectorH;
  using typename SemiExplicit<Real, 5, 0>::MatrixJH;

  Shampine1SemiExplicit() : SemiExplicit<Real, 5, 0>("Shampine1SemiExplicit") {}

  ~Shampine1SemiExplicit() {}

  MatrixA A(const VectorF & /*x*/, const Real /*t*/) const override {
    return MatrixA::Identity();
  }

  TensorTA TA_x(const VectorF & /*x*/, const Real /*t*/) const override {
    TensorTA TA_x;
    TA_x[0].setZero();
    TA_x[1].setZero();
    TA_x[2].setZero();
    TA_x[3].setZero();
    TA_x[4].setZero();
    return TA_x;
  }

  VectorB b(const VectorF &x, const Real /*t*/) const override {
    VectorB b;
    b << 0.5 * x(0) * (x(2) - x(0)) / x(1), -0.5 * (x(2) - x(0)),
        (0.9 - 1000 * (x(2) - x(4)) - 0.5 * x(2) * (x(2) - x(0))) / x(3),
        0.5 * (x(2) - x(0)), 100 * (x(2) - x(4));
    return b;
  }

  MatrixJB Jb_x(const VectorF &x, const Real /*t*/) const override {
    MatrixJB Jb_x(MatrixJB::Zero());
    Jb_x(0, 0) = 0.5 * (x(2) - 2.0 * x(0)) / x(1);
    Jb_x(0, 1) = -0.5 * x(0) * (x(2) - x(0)) / std::pow(x(1), 2.0);
    Jb_x(0, 2) = 0.5 * x(0) / x(1);
    Jb_x(1, 0) = 0.5;
    Jb_x(1, 2) = -0.5;
    Real N{0.9 - 1000.0 * (x(2) - x(4)) - 0.5 * std::pow(x(2), 2.0) +
           0.5 * x(0) * x(2)};
    Jb_x(2, 0) = 0.5 * x(2) / x(3);
    Jb_x(2, 2) = (-1000.0 - x(2) + 0.5 * x(0)) / x(3);
    Jb_x(2, 3) = -N / std::pow(x(3), 2.0);
    Jb_x(2, 4) = 1000.0 / x(3);
    Jb_x(3, 0) = -0.5;
    Jb_x(3, 2) = 0.5;
    Jb_x(4, 2) = 100.0;
    Jb_x(4, 4) = -100.0;
    return Jb_x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Real, typename System, typename Integrator>
class Shampine1Problem : public BVP<Real, 5, 0, Integrator> {
 public:
  using typename BVP<Real, 5, 0, Integrator>::VectorF;
  using typename BVP<Real, 5, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 5, Eigen::Dynamic>;

  Shampine1Problem()
      : BVP<Real, 5, 0, Integrator>("Shampine1Problem",
                                    std::make_unique<System>(),
                                    std::make_unique<Integrator>()) {}

  ~Shampine1Problem() {}

  static Real time_start() {
    return 0.0;
  }

  static Real time_end() {
    return 40.0;
  }

  VectorF b(const VectorF &x_ini, const VectorF &x_end) const override {
    VectorF b;
    b << x_ini(0) - 1.0, x_ini(1) - 1.0, x_ini(2) - 1.0, x_ini(3) + 10.0,
        x_end(2) - x_end(4);
    return b;
  }

  MatrixJF Jb_x_ini(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_ini(MatrixJF::Zero());
    Jb_x_ini(0, 0) = 1.0;
    Jb_x_ini(1, 1) = 1.0;
    Jb_x_ini(2, 2) = 1.0;
    Jb_x_ini(3, 3) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(const VectorF & /*x_ini*/,
                    const VectorF & /*x_end*/) const override {
    MatrixJF Jb_x_end(MatrixJF::Zero());
    Jb_x_end(4, 2) = 1.0;
    Jb_x_end(4, 4) = -1.0;
    return Jb_x_end;
  }

  static VectorF exact_solution(const Real /*t*/) {
    return VectorF::Zero();
  }

  static MatrixX exact_solution(const VectorX &t) {
    return MatrixX::Zero(5, t.size());
  }

  VectorF guess(const Real t) {
    VectorF x;
    x << 1.0, 1.0, -4.5 * t * t + 8.91 * t + 1.0, -10.0,
        -4.5 * t * t + 9.0 * t + 0.91;
    return x;
  }

  MatrixX guess(const VectorX &t) {
    MatrixX x(5, t.size());
    for (Integer i{0}; i < t.size(); ++i) {
      x.col(i) = Shampine1Problem<Real, System, Integrator>::guess(t(i));
    }
    return x;
  }
};

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif  // TESTS_PROBLEMS_SHAMPINE1_HH
