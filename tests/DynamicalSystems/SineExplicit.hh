/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_SINE_EXPLICIT_HH

using namespace Sandals;

class SineExplicit : public Explicit<2, 0>
{
public:
  using Matrix0 = Eigen::Matrix<Real, 0, 2>;

  SineExplicit() : Explicit<2, 0>("SineExplicit") {}

  ~SineExplicit() {}

  Vector2 f(Vector2 const &/*x*/, Real t) const override
  {
    Vector2 f;
    f << std::cos(t), std::sin(t);
    return f;
  }

  Matrix2 Jf_x(Vector2 const &/*x*/, Real /*t*/) const override
  {return Matrix2::Zero();}

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
      x.col(i) = this->analytical_solution(t(i));
    }
    return x;
  }
};

#endif // TESTS_SINE_EXPLICIT_HH
