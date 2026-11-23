/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Optimal control problems with automatic derivatives (CppAD): the recorded model reproduces the
// hand-coded derivatives, and the third derivatives of the invariants are exact.

#include <gtest/gtest.h>

#include "Sandals/OptimalControlAD.hh"

#include "PendulumProblems.hh"

using namespace Sandals;
using namespace TestOCP;

// Index-1 pendulum of PendulumProblems.hh, written once as templates.
class PendulumIndex1AD : public OptimalControlProblemAD<PendulumIndex1AD, 5, 1, 2, 2, 1>
{
public:
  Real m{1.0}, L{1.0}, g{9.81};

  PendulumIndex1AD() : OptimalControlProblemAD("PendulumIndex1AD", 0.0, 2.0) {}

  template <typename T>
  StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
  {
    StateVector<T> F;
    F(0) = x_dot(X) - x(VX);
    F(1) = x_dot(Y) - x(VY);
    F(2) = m*x_dot(VX) + 2.0*x(X)*x(LAMBDA) - u(0);
    F(3) = m*x_dot(VY) + 2.0*x(Y)*x(LAMBDA) + m*g;
    F(4) = x(VX)*x(VX) + x(VY)*x(VY) + x(X)*x_dot(VX) + x(Y)*x_dot(VY);
    return F;
  }

  template <typename T>
  InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
  {
    InvariantVector<T> h;
    h(0) = x(X)*x(X) + x(Y)*x(Y) - L*L;
    h(1) = x(X)*x(VX) + x(Y)*x(VY);
    return h;
  }

  template <typename T>
  T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const & t) const
  {return 0.5*u(0)*u(0) + sin(t)*x(X)*x(Y);}

  template <typename T>
  T terminal_cost(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
  {return x_end(VX)*x_end(VX) + x_ini(X)*x_end(Y)*x_end(Y);}

  template <typename T>
  BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const & x_end) const
  {return BoundaryVector<T>(x_ini(X)*x_ini(Y), exp(x_end(VY)));}

  template <typename T>
  PathVector<T> path(StateVector<T> const & x, ControlVector<T> const & u, T const &) const
  {return PathVector<T>(x(VX)*u(0)*u(0));}
};

// A system with a non-quadratic invariant, whose third derivatives do not vanish.
class CurvedInvariant : public OptimalControlProblemAD<CurvedInvariant, 2, 1, 1>
{
public:
  CurvedInvariant() : OptimalControlProblemAD("CurvedInvariant", 0.0, 1.0) {}

  template <typename T>
  StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const & t) const
  {return StateVector<T>(x_dot(0) - sin(x(1))*u(0), x_dot(1) - exp(x(0)*t));}

  template <typename T>
  InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
  {return InvariantVector<T>(sin(x(0))*x(1)*x(1) + exp(x(1)));}
};

using Vector5 = Eigen::Vector<Real, 5>;
using Vector1 = Eigen::Vector<Real, 1>;
using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;

Real error(Eigen::MatrixXd const & a, Eigen::MatrixXd const & b) {return (a - b).norm() / std::max(1.0, b.norm());}

TEST(AutoDiff, SystemMatchesHandCodedDerivatives) {
  PendulumIndex1AD ad;
  PendulumIndex1<2> hand;
  auto const & sys{*ad.system()};
  std::srand(1);
  Vector5 const x{Vector5::Random()}, x_dot{Vector5::Random()}, w{Vector5::Random()};
  Vector1 const u(0.3);
  Eigen::Vector2d const mu(0.7, -1.2);
  Real const t{0.1};

  EXPECT_LT(error(sys.F(x, x_dot, u, t), hand.F(x, x_dot, u, t)), 1.0e-15);
  EXPECT_LT(error(sys.JF_x(x, x_dot, u, t), hand.JF_x(x, x_dot, u, t)), 1.0e-15);
  EXPECT_LT(error(sys.JF_x_dot(x, x_dot, u, t), hand.JF_x_dot(x, x_dot, u, t)), 1.0e-15);
  EXPECT_LT(error(sys.JF_u(x, x_dot, u, t), hand.JF_u(x, x_dot, u, t)), 1.0e-15);
  EXPECT_LT(error(sys.HF_x_x(x, x_dot, u, t, w), hand.HF_x_x(x, x_dot, u, t, w)), 1.0e-15);
  EXPECT_LT(error(sys.HF_x_x_dot(x, x_dot, u, t, w), hand.HF_x_x_dot(x, x_dot, u, t, w)), 1.0e-15);
  EXPECT_LT(error(sys.HF_x_dot_x_dot(x, x_dot, u, t, w), hand.HF_x_dot_x_dot(x, x_dot, u, t, w)), 1.0e-15);
  EXPECT_LT(error(sys.HF_x_u(x, x_dot, u, t, w), hand.HF_x_u(x, x_dot, u, t, w)), 1.0e-15);
  EXPECT_LT(error(sys.HF_u_u(x, x_dot, u, t, w), hand.HF_u_u(x, x_dot, u, t, w)), 1.0e-15);
  EXPECT_LT(error(sys.h(x, t), hand.h(x, t)), 1.0e-15);
  EXPECT_LT(error(sys.Jh_x(x, t), hand.Jh_x(x, t)), 1.0e-15);
  EXPECT_LT(error(sys.Hh_x_x(x, t, mu), hand.Hh_x_x(x, t, mu)), 1.0e-15);
  EXPECT_LT(sys.Th_x_x(x, t, mu, w).norm(), 1.0e-15); // quadratic invariants
}

TEST(AutoDiff, ProblemTermsMatchFiniteDifferences) {
  // The finite differences defaults of OptimalControlProblem, applied to the recorded values
  PendulumIndex1AD ad;
  using Defaults = OptimalControlProblem<Real, 5, 1, 2>;
  std::srand(2);
  Vector5 const x{Vector5::Random()}, y{Vector5::Random()};
  Vector1 const u(-0.4);
  VectorX const l_b{VectorX::Random(2)}, l_g{VectorX::Random(1)};
  Real const t{0.7};

  EXPECT_LT(error(ad.mayer_gradient(x, y), ad.Defaults::mayer_gradient(x, y)), 1.0e-8);
  EXPECT_LT(error(ad.mayer_hessian(x, y), ad.Defaults::mayer_hessian(x, y)), 1.0e-6);
  EXPECT_LT(error(ad.lagrange_gradient(x, u, t), ad.Defaults::lagrange_gradient(x, u, t)), 1.0e-8);
  EXPECT_LT(error(ad.lagrange_hessian(x, u, t), ad.Defaults::lagrange_hessian(x, u, t)), 1.0e-6);
  EXPECT_LT(error(ad.boundary_conditions_jacobian(x, y), ad.Defaults::boundary_conditions_jacobian(x, y)), 1.0e-8);
  EXPECT_LT(error(ad.boundary_conditions_hessian(x, y, l_b), ad.Defaults::boundary_conditions_hessian(x, y, l_b)), 1.0e-6);
  EXPECT_LT(error(ad.path_constraints_jacobian(x, u, t), ad.Defaults::path_constraints_jacobian(x, u, t)), 1.0e-8);
  EXPECT_LT(error(ad.path_constraints_hessian(x, u, t, l_g), ad.Defaults::path_constraints_hessian(x, u, t, l_g)), 1.0e-6);
  EXPECT_EQ(ad.boundary_conditions_number(), 2);
  EXPECT_EQ(ad.path_constraints_number(), 1);
}

TEST(AutoDiff, ThirdDerivativesOfInvariants) {
  // Nested taping against finite differences of the exact Hessian of the invariants
  CurvedInvariant ad;
  auto const & sys{*ad.system()};
  Eigen::Vector2d const x(0.4, -0.7), w(0.3, 1.1);
  Vector1 const mu(1.7);
  Eigen::Matrix2d const T{sys.Th_x_x(x, 0.2, mu, w)};
  Eigen::Matrix2d const T_fd{sys.Controlled<Real, 2, 1, 1>::Th_x_x(x, 0.2, mu, w)};
  EXPECT_GT(T.norm(), 1.0);
  EXPECT_LT(error(T, T_fd), 1.0e-8);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
