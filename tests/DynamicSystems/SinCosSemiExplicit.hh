/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_SINCOS_SEMIEXPLICIT_HH
#define TESTS_SINCOS_SEMIEXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/System/SemiExplicit.hh"

using namespace Sandals;
using namespace Eigen;

template<typename Real = double>
class SinCosSemiExplicit : public SemiExplicit<Real, 2, 0>
{
public:
  using VectorF  = typename SemiExplicit<Real, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<Real, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<Real, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<Real, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<Real, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<Real, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<Real, 2, 0>::MatrixJH;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX  = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  SinCosSemiExplicit() : SemiExplicit<Real, 2, 0>("SinCosSemiExplicit") {}

  ~SinCosSemiExplicit() {}

  MatrixA A(VectorF const &/*x*/, Real /*t*/) const override
  {
    MatrixA A;
    A.setIdentity();
    return A;
  }

  TensorTA TA_x(VectorF const &/*x*/, Real /*t*/) const override
  {
    TensorTA TA_x(2);
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(VectorF const &/*x*/, Real t) const override
  {
    VectorB b;
    b << std::cos(t), std::sin(t);
    return b;
  }

  MatrixJB Jb_x(VectorF const &/*x*/, Real /*t*/) const override
  {return MatrixJB::Zero();}

  VectorH h(VectorF const &/*x*/, Real /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const &/*x*/, Real /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}

  VectorF ics() const {return VectorF::Zero();}

  VectorF analytical_solution(Real t) const {
    VectorF x;
    x << std::sin(t), Real(1.0)-std::cos(t);
    return x;
  }

  MatrixX analytical_solution(VectorX const &t) const {
    MatrixX x(2, t.size());
    for (int i = 0; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }
};

#endif // TESTS_SINCOS_SEMIEXPLICIT_HH
