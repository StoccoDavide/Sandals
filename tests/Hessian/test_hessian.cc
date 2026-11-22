/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <memory>

#include <gtest/gtest.h>

#include "Sandals.hh"
#include "Sandals/RungeKutta/RK4.hh"

#include "SimplePendulumExplicit.hh"
#include "SimplePendulumImplicit.hh"
#include "SimplePendulumSemiExplicit.hh"
#include "OscillatorExplicit.hh"
#include "OscillatorImplicit.hh"
#include "OscillatorSemiExplicit.hh"
#include "ThreeBodyExplicit.hh"
#include "ThreeBodyImplicit.hh"
#include "ThreeBodySemiExplicit.hh"
#include "SinCosSemiExplicit.hh"

using namespace Sandals;

using Real = double;

/*\
 |  Reference derivatives computed by central finite differences, independently of the library
\*/

// Jacobian of g(z) at z by central finite differences (column j = dg/dz_j).
template <Integer N, Integer R, typename Function>
Eigen::Matrix<Real, R, N> fd(Function && g, Eigen::Vector<Real, N> const & z, Real const delta = 1.0e-6)
{
  Eigen::Matrix<Real, R, N> jac;
  for (Integer j{0}; j < N; ++j) {
    Eigen::Vector<Real, N> z_p(z), z_m(z);
    z_p(j) += delta;
    z_m(j) -= delta;
    jac.col(j) = (g(z_p) - g(z_m)) / (2.0*delta);
  }
  return jac;
}

// Relative error between two matrices.
template <typename A, typename B>
Real rel_error(A const & a, B const & b) {return (a - b).norm() / std::max(Real(1.0), b.norm());}

// Check the first and second derivatives of the implicit form of a system at a random point.
template <Integer N, Integer M>
void check_system(std::shared_ptr<Implicit<Real, N, M>> sys, Real const tol = 1.0e-6)
{
  using VectorF = Eigen::Vector<Real, N>;
  using VectorH = Eigen::Vector<Real, M>;

  std::srand(42);
  VectorF const x{VectorF::Random()};
  VectorF const x_dot{VectorF::Random()};
  VectorF const lambda{VectorF::Random()};
  Real const t{0.3};

  SCOPED_TRACE(sys->name());

  // First derivatives against finite differences of F
  EXPECT_LT(rel_error(sys->JF_x(x, x_dot, t),
    fd<N, N>([&] (VectorF const & z) {return sys->F(z, x_dot, t);}, x)), tol);
  EXPECT_LT(rel_error(sys->JF_x_dot(x, x_dot, t),
    fd<N, N>([&] (VectorF const & z) {return sys->F(x, z, t);}, x_dot)), tol);

  // Second derivatives against finite differences of the weighted Jacobians
  auto grad_x     = [&] (VectorF const & xx, VectorF const & xd) {return VectorF(sys->JF_x(xx, xd, t).transpose()*lambda);};
  auto grad_x_dot = [&] (VectorF const & xx, VectorF const & xd) {return VectorF(sys->JF_x_dot(xx, xd, t).transpose()*lambda);};

  auto HF_x_x         = sys->HF_x_x(x, x_dot, t, lambda);
  auto HF_x_x_dot     = sys->HF_x_x_dot(x, x_dot, t, lambda);
  auto HF_x_dot_x_dot = sys->HF_x_dot_x_dot(x, x_dot, t, lambda);

  EXPECT_LT(rel_error(HF_x_x, fd<N, N>([&] (VectorF const & z) {return grad_x(z, x_dot);}, x)), tol);
  EXPECT_LT(rel_error(HF_x_x_dot, fd<N, N>([&] (VectorF const & z) {return grad_x(x, z);}, x_dot)), tol);
  EXPECT_LT(rel_error(HF_x_x_dot.transpose(), fd<N, N>([&] (VectorF const & z) {return grad_x_dot(z, x_dot);}, x)), tol);
  EXPECT_LT(rel_error(HF_x_dot_x_dot, fd<N, N>([&] (VectorF const & z) {return grad_x_dot(x, z);}, x_dot)), tol);

  // Symmetry
  EXPECT_LT(rel_error(HF_x_x, HF_x_x.transpose()), 1.0e-12);
  EXPECT_LT(rel_error(HF_x_dot_x_dot, HF_x_dot_x_dot.transpose()), 1.0e-12);

  // Time-reversed Hessians against finite differences of the time-reversed Jacobians
  auto grad_x_rev = [&] (VectorF const & xx, VectorF const & xd) {return VectorF(sys->JF_x_reverse(xx, xd, t).transpose()*lambda);};
  auto grad_x_dot_rev = [&] (VectorF const & xx, VectorF const & xd) {return VectorF(sys->JF_x_dot_reverse(xx, xd, t).transpose()*lambda);};
  EXPECT_LT(rel_error(sys->HF_x_x_reverse(x, x_dot, t, lambda),
    fd<N, N>([&] (VectorF const & z) {return grad_x_rev(z, x_dot);}, x)), tol);
  EXPECT_LT(rel_error(sys->HF_x_x_dot_reverse(x, x_dot, t, lambda),
    fd<N, N>([&] (VectorF const & z) {return grad_x_rev(x, z);}, x_dot)), tol);
  EXPECT_LT(rel_error(sys->HF_x_dot_x_dot_reverse(x, x_dot, t, lambda),
    fd<N, N>([&] (VectorF const & z) {return grad_x_dot_rev(x, z);}, x_dot)), tol);

  // Invariants Hessian
  if constexpr (M > 0) {
    VectorH const mu{VectorH::Random()};
    EXPECT_LT(rel_error(sys->Hh_x_x(x, t, mu), fd<N, N>([&] (VectorF const & z)
      {return VectorF(sys->Jh_x(z, t).transpose()*mu);}, x)), tol);
  } else {
    EXPECT_EQ(sys->Hh_x_x(x, t, VectorH()).norm(), 0.0);
  }
}

/*\
 |  Hessians of the test systems (explicit, implicit and semi-explicit forms)
\*/

TEST(Hessian, SimplePendulum) {
  check_system<2, 0>(std::make_shared<SimplePendulumExplicit<Real>>());
  check_system<2, 0>(std::make_shared<SimplePendulumImplicit<Real>>());
  check_system<2, 0>(std::make_shared<SimplePendulumSemiExplicit<Real>>());
}

TEST(Hessian, Oscillator) {
  check_system<2, 1>(std::make_shared<OscillatorExplicit<Real>>());
  check_system<2, 1>(std::make_shared<OscillatorImplicit<Real>>());
  check_system<2, 1>(std::make_shared<OscillatorSemiExplicit<Real>>());
}

TEST(Hessian, ThreeBody) {
  check_system<12, 0>(std::make_shared<ThreeBodyExplicit<Real>>(), 1.0e-5);
  check_system<12, 0>(std::make_shared<ThreeBodyImplicit<Real>>(), 1.0e-5);
  check_system<12, 0>(std::make_shared<ThreeBodySemiExplicit<Real>>(), 1.0e-5);
}

TEST(Hessian, SinCos) {
  check_system<2, 0>(std::make_shared<SinCosSemiExplicit<Real>>());
}

/*\
 |  Exact Hessians against the finite differences defaults
\*/

// Pendulum f = [x1, -g/l sin(x0)], so d2(lambda^T f)/dx2 = [[lambda1 g/l sin(x0), 0], [0, 0]].
TEST(Hessian, ExplicitExactVersusDefault) {
  using Wrapper = ExplicitWrapper<Real, 2, 0>;
  using VectorF = Wrapper::VectorF;
  using MatrixJF = Wrapper::MatrixJF;
  Real const g{9.81};
  auto f    = [g] (VectorF const & x, Real) {return VectorF(x(1), -g*std::sin(x(0)));};
  auto Jf_x = [g] (VectorF const & x, Real) {MatrixJF J; J << 0.0, 1.0, -g*std::cos(x(0)), 0.0; return J;};
  auto Hf_x_x = [g] (VectorF const & x, Real, VectorF const & l)
    {MatrixJF H(MatrixJF::Zero()); H(0, 0) = l(1)*g*std::sin(x(0)); return H;};

  Wrapper fd_sys("PendulumFD", f, Jf_x);
  Wrapper exact_sys("PendulumExact", f, Jf_x);
  exact_sys.Hf_x_x() = Hf_x_x;

  VectorF const x(0.7, -0.2), x_dot(0.1, 0.4), lambda(0.5, -1.5);
  MatrixJF const H_exact{Hf_x_x(x, 0.0, lambda)};
  EXPECT_LT(rel_error(exact_sys.Hf_x_x(x, 0.0, lambda), H_exact), 1.0e-15);
  EXPECT_LT(rel_error(fd_sys.Hf_x_x(x, 0.0, lambda), H_exact), 1.0e-8);
  EXPECT_LT(rel_error(exact_sys.HF_x_x(x, x_dot, 0.0, lambda), MatrixJF(-H_exact)), 1.0e-15);
  EXPECT_EQ(exact_sys.HF_x_x_dot(x, x_dot, 0.0, lambda).norm(), 0.0);
  EXPECT_EQ(exact_sys.HF_x_dot_x_dot(x, x_dot, 0.0, lambda).norm(), 0.0);
}

// Implicit F = [x0' x1' - x1^2, sin(x0) x1']: user-provided Hessians are dispatched by the wrapper.
TEST(Hessian, ImplicitWrapperDispatch) {
  using Wrapper = ImplicitWrapper<Real, 2, 0>;
  using VectorF = Wrapper::VectorF;
  using MatrixJF = Wrapper::MatrixJF;
  auto F = [] (VectorF const & x, VectorF const & xd, Real)
    {return VectorF(xd(0)*xd(1) - x(1)*x(1), std::sin(x(0))*xd(1));};
  auto JF_x = [] (VectorF const & x, VectorF const & xd, Real)
    {MatrixJF J; J << 0.0, -2.0*x(1), std::cos(x(0))*xd(1), 0.0; return J;};
  auto JF_x_dot = [] (VectorF const & x, VectorF const & xd, Real)
    {MatrixJF J; J << xd(1), xd(0), 0.0, std::sin(x(0)); return J;};
  auto HF_x_x = [] (VectorF const & x, VectorF const & xd, Real, VectorF const & l)
    {MatrixJF H; H << -l(1)*std::sin(x(0))*xd(1), 0.0, 0.0, -2.0*l(0); return H;};
  auto HF_x_x_dot = [] (VectorF const & x, VectorF const &, Real, VectorF const & l)
    {MatrixJF H; H << 0.0, l(1)*std::cos(x(0)), 0.0, 0.0; return H;};
  auto HF_x_dot_x_dot = [] (VectorF const &, VectorF const &, Real, VectorF const & l)
    {MatrixJF H; H << 0.0, l(0), l(0), 0.0; return H;};

  auto fd_sys = std::make_shared<Wrapper>("ImplicitFD", F, JF_x, JF_x_dot);
  auto exact_sys = std::make_shared<Wrapper>("ImplicitExact", F, JF_x, JF_x_dot);
  exact_sys->HF_x_x() = HF_x_x;
  exact_sys->HF_x_x_dot() = HF_x_x_dot;
  exact_sys->HF_x_dot_x_dot() = HF_x_dot_x_dot;

  check_system<2, 0>(fd_sys);
  check_system<2, 0>(exact_sys);

  VectorF const x(0.3, -0.8), xd(1.2, 0.5), l(-0.4, 0.9);
  EXPECT_LT(rel_error(fd_sys->HF_x_x(x, xd, 0.0, l), HF_x_x(x, xd, 0.0, l)), 1.0e-8);
  EXPECT_LT(rel_error(fd_sys->HF_x_x_dot(x, xd, 0.0, l), HF_x_x_dot(x, xd, 0.0, l)), 1.0e-8);
  EXPECT_LT(rel_error(fd_sys->HF_x_dot_x_dot(x, xd, 0.0, l), HF_x_dot_x_dot(x, xd, 0.0, l)), 1.0e-8);
}

// Semi-explicit with a state-dependent mass matrix A(x) = [[1 + x0^2, x1], [0, 2 + sin(x1)]] and
// b(x) = [x0 x1, cos(x0)]: checks JF_x, the mixed Hessian (from TA_x), and HA_x_x / Hb_x_x.
TEST(Hessian, SemiExplicitStateDependentMass) {
  using Wrapper = SemiExplicitWrapper<Real, 2, 0>;
  using VectorF = Wrapper::VectorF;
  using MatrixA = Wrapper::MatrixA;
  using TensorTA = Wrapper::TensorTA;
  auto A = [] (VectorF const & x, Real) {MatrixA A; A << 1.0 + x(0)*x(0), x(1), 0.0, 2.0 + std::sin(x(1)); return A;};
  auto TA_x = [] (VectorF const & x, Real) {
    TensorTA T(2);
    T[0] << 2.0*x(0), 0.0, 0.0, 0.0;
    T[1] << 0.0, 1.0, 0.0, std::cos(x(1));
    return T;
  };
  auto b = [] (VectorF const & x, Real) {return VectorF(x(0)*x(1), std::cos(x(0)));};
  auto Jb_x = [] (VectorF const & x, Real) {MatrixA J; J << x(1), x(0), -std::sin(x(0)), 0.0; return J;};

  auto sys = std::make_shared<Wrapper>("SemiExplicitStateDependentMass", A, TA_x, b, Jb_x);
  check_system<2, 0>(sys);

  // Exact second derivatives
  VectorF const x(0.4, -0.6), v(0.7, 1.1), l(1.3, -0.5);
  MatrixA HA; HA << 2.0*l(0)*v(0), 0.0, 0.0, -l(1)*std::sin(x(1))*v(1);
  MatrixA Hb; Hb << -l(1)*std::cos(x(0)), l(0), l(0), 0.0;
  EXPECT_LT(rel_error(sys->HA_x_x(x, v, 0.0, l), HA), 1.0e-8);
  EXPECT_LT(rel_error(sys->Hb_x_x(x, 0.0, l), Hb), 1.0e-8);
  EXPECT_LT(rel_error(sys->HF_x_x(x, v, 0.0, l), MatrixA(HA - Hb)), 1.0e-8);
}

// Linear systems have zero Hessians.
TEST(Hessian, LinearIsZero) {
  using Wrapper = LinearWrapper<Real, 2, 0>;
  using VectorF = Wrapper::VectorF;
  using MatrixE = Wrapper::MatrixE;
  auto E = [] (Real) {MatrixE E; E << 2.0, 0.0, 1.0, 1.0; return E;};
  auto A = [] (Real t) {MatrixE A; A << 0.0, 1.0, -1.0, t; return A;};
  auto b = [] (Real t) {return VectorF(std::sin(t), 1.0);};
  auto sys = std::make_shared<Wrapper>("Linear", E, A, b);
  check_system<2, 0>(sys);
  VectorF const x(1.0, 2.0), l(3.0, 4.0);
  EXPECT_EQ(sys->Hf_x_x(x, 0.1, l).norm(), 0.0);
  EXPECT_EQ(sys->HF_x_x(x, x, 0.1, l).norm(), 0.0);
}

// Energy invariant h = (m x1^2 + k x0^2)/2 - h0 has Hessian mu diag(k, m) (with k = m = 1).
TEST(Hessian, InvariantsExact) {
  OscillatorExplicit<Real> sys;
  Eigen::Vector<Real, 2> const x(0.3, -1.2);
  Eigen::Vector<Real, 1> const mu(2.5);
  EXPECT_LT(rel_error(sys.Hh_x_x(x, 0.0, mu), Eigen::Matrix<Real, 2, 2>(2.5*Eigen::Matrix<Real, 2, 2>::Identity())), 1.0e-8);
}

/*\
 |  Sensitivity through the projection on the invariants (uses Hh_x_x)
\*/

TEST(Hessian, ProjectionSensitivity) {
  using VectorF = Eigen::Vector<Real, 2>;
  using MatrixJF = Eigen::Matrix<Real, 2, 2>;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;

  auto sys = std::make_shared<OscillatorExplicit<Real>>();
  RK4<Real, 2, 1> rk(sys);
  rk.enable_projection();
  rk.projection_tolerance(1.0e-14);
  rk.max_projection_iterations(20);

  VectorX const t_mesh(VectorX::LinSpaced(21, 0.0, 2.0));
  VectorF const ics(1.0, 0.0);

  Solution<Real, 2, 1> sol;
  MatrixJF Jx;
  ASSERT_TRUE(rk.solve(t_mesh, ics, sol, Jx));

  // The projected flow keeps the energy fixed at h0: perturbations leave the trajectory on the
  // manifold, so Jx must map onto the tangent space at the final point.
  MatrixJF Jx_fd(fd<2, 2>([&] (VectorF const & z) {
    Solution<Real, 2, 1> s;
    EXPECT_TRUE(rk.solve(t_mesh, z, s));
    return VectorF(s.x.col(s.size() - 1));
  }, ics, 1.0e-7));
  EXPECT_LT(rel_error(Jx, Jx_fd), 1.0e-5) << "Jx =\n" << Jx << "\nJx_fd =\n" << Jx_fd;

  VectorF const x_end(sol.x.col(sol.size() - 1));
  EXPECT_LT((sys->Jh_x(x_end, 2.0) * Jx).norm(), 1.0e-6);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
