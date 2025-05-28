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

#ifndef TESTS_SIN_EXPLICIT_HH
#define TESTS_SIN_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;
using namespace Eigen;

template<typename Real = double>
class SinExplicit : public Explicit<Real, 1, 0>
{
public:
  using VectorF  = typename Explicit<Real, 1, 0>::VectorF;
  using MatrixJF = typename Explicit<Real, 1, 0>::MatrixJF;
  using VectorH  = typename Explicit<Real, 1, 0>::VectorH;
  using MatrixJH = typename Explicit<Real, 1, 0>::MatrixJH;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX  = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

  SinExplicit() : Explicit<Real, 1, 0>("SinExplicit") {}

  ~SinExplicit() {}

  VectorF f(VectorF const &/*x*/, Real t) const override
  {
    VectorF f;
    f << std::cos(t);
    return f;
  }

  MatrixJF Jf_x(VectorF const &/*x*/, Real /*t*/) const override
  {return MatrixJF::Zero();}

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

#endif // TESTS_SIN_EXPLICIT_HH
