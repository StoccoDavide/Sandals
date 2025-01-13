/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_SINE_IMPLICIT_HH

using namespace Sandals;

class SineImplicit : public Implicit<2, 0>
{
public:
  using Matrix0 = Eigen::Matrix<Real, 0, 2>;

  SineImplicit() : Implicit<2, 0>("SineImplicit") {}

  ~SineImplicit() {}

  Vector2 F(Vector2 const &/*x*/, Vector2 const &x_dot, Real t) const override
  {
    Vector2 F;
    F << x_dot(0) - std::cos(t), x_dot(1) - std::sin(t);
    return F;
  }

  Matrix2 JF_x(Vector2 const &/*x*/, Vector2 const &/*x_dot*/, Real /*t*/) const override
  {return Matrix2::Zero();}

  Matrix2 JF_x_dot(Vector2 const &/*x*/, Vector2 const &/*x_dot*/, Real /*t*/) const override
  {return Matrix2::Identity();}

  Vector0 h(Vector2 const &/*x*/, Real /*t*/) const override {return Vector0::Zero();}

  Matrix0 Jh_x(Vector2 const &/*x*/, Real /*t*/) const override {return Matrix0::Zero();}

  bool in_domain(Vector2 const &/*x*/, Real /*t*/) const override {return true;}

  Vector2 ics() const {return Vector2::Zero();}

  Vector2 analytical_solution(Real t) const {
    Vector2 x;
    x << std::sin(t), std::cos(t);
    return x;
  }

  MatrixX analytical_solution(VectorX const &t) const {
    MatrixX x(2, t.size());
    for (int i = 0; i < t.size(); ++i) {
      x(0, i) = std::sin(t(i));
      x(1, i) = std::cos(t(i));
    }
    return x;
  }
};

#endif // TESTS_SINE_IMPLICIT_HH
