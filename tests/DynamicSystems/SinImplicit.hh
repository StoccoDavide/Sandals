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

#ifndef TESTS_SIN_IMPLICIT_HH
#define TESTS_SIN_IMPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Implicit.hh"

using namespace Sandals;
using namespace Eigen;

template<typename Real = double>
class SinImplicit : public Implicit<Real, 1, 0>
{
public:
  using VectorF  = typename Implicit<Real, 1, 0>::VectorF;
  using MatrixJF = typename Implicit<Real, 1, 0>::MatrixJF;
  using VectorH  = typename Implicit<Real, 1, 0>::VectorH;
  using MatrixJH = typename Implicit<Real, 1, 0>::MatrixJH;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX  = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

  SinImplicit() : Implicit<Real, 1, 0>("SinImplicit") {}

  ~SinImplicit() {}

  VectorF F(VectorF const &/*x*/, VectorF const &x_dot, Real t) const override
  {
    VectorF F;
    F << x_dot(0) - std::cos(t);
    return F;
  }

  MatrixJF JF_x(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const override
  {return MatrixJF::Zero();}

  MatrixJF JF_x_dot(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const override
  {return MatrixJF::Identity();}

  VectorH h(VectorF const &/*x*/, Real /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const &/*x*/, Real /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}

  VectorF ics() const {return VectorF::Zero();}

  VectorF analytical_solution(Real t) const {
    VectorF x;
    x << std::sin(t);
    return x;
  }

  MatrixX analytical_solution(VectorX const &t) const {
    MatrixX x(1, t.size());
    for (int i = 0; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }
};

#endif // TESTS_SIN_IMPLICIT_HH
