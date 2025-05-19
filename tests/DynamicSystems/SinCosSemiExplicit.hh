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

class SinCosSemiExplicit : public SemiExplicit<double, 2, 0>
{
private:
  VectorF m_ics{0.0, 1.0}; // Initial conditions

public:
  using VectorF  = typename SemiExplicit<double, 2, 0>::VectorF;
  using MatrixA  = typename SemiExplicit<double, 2, 0>::MatrixA;
  using TensorTA = typename SemiExplicit<double, 2, 0>::TensorTA;
  using VectorB  = typename SemiExplicit<double, 2, 0>::VectorB;
  using MatrixJB = typename SemiExplicit<double, 2, 0>::MatrixJB;
  using VectorH  = typename SemiExplicit<double, 2, 0>::VectorH;
  using MatrixJH = typename SemiExplicit<double, 2, 0>::MatrixJH;


  SinCosSemiExplicit() : SemiExplicit<double, 2, 0>("SinCosSemiExplicit") {}

  ~SinCosSemiExplicit() {}

  MatrixA A(VectorF const &/*x*/, double /*t*/) const override
  {
    MatrixA A;
    A.setIdentity();
    return A;
  }

  TensorTA TA_x(VectorF const &/*x*/, double /*t*/) const override
  {
    TensorTA TA_x(2);
    TA_x[0].setZero();
    TA_x[1].setZero();
    return TA_x;
  }

  VectorB b(VectorF const &/*x*/, double t) const override
  {
    VectorB b;
    b << std::cos(t), std::sin(t);
    return b;
  }

  MatrixJB Jb_x(VectorF const &/*x*/, double /*t*/) const override
  {return MatrixJB::Zero();}

  VectorH h(VectorF const &/*x*/, double /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const &/*x*/, double /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const &/*x*/, double /*t*/) const override {return true;}

  VectorF const &ics() const {return this->m_ics;}

  VectorF analytical_solution(double t) const {
    VectorF x;
    x << std::sin(t), double(1.0)-std::cos(t);
    return x;
  }

  MatrixXd analytical_solution(VectorXd const &t) const {
    MatrixXd x(2, t.size());
    for (int i = 0; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }
};

#endif // TESTS_SINCOS_SEMIEXPLICIT_HH
