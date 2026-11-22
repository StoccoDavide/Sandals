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
#include "Sandals/RungeKutta/GaussLegendre4.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/SSPIRK33.hh"

#include "SimplePendulumExplicit.hh"
#include "SimplePendulumImplicit.hh"
#include "ThreeBodyExplicit.hh"
#include "ThreeBodyImplicit.hh"

using namespace Sandals;

using Real = double;

// Compare the propagated Jacobian dx(t_end)/dx(t_ini) with central finite differences of the
// final state with respect to the initial conditions.
template <typename Integrator, Integer N, typename System>
void check_sensitivity(System system, Eigen::Vector<Real, N> const & ics, Real const t_end,
  Integer const steps, Real const tolerance)
{
  using VectorN  = Eigen::Vector<Real, N>;
  using MatrixJX = Eigen::Matrix<Real, N, N>;
  using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;

  Integrator rk(system);
  rk.disable_projection();
  VectorX t_mesh(VectorX::LinSpaced(steps + 1, 0.0, t_end));

  Solution<Real, N, 0> sol;
  MatrixJX Jx;
  ASSERT_TRUE(rk.solve(t_mesh, ics, sol, Jx));

  MatrixJX Jx_fd;
  Real const delta{1.0e-6};
  for (Integer j{0}; j < N; ++j) {
    VectorN ics_p(ics), ics_m(ics);
    ics_p(j) += delta;
    ics_m(j) -= delta;
    Solution<Real, N, 0> sol_p, sol_m;
    ASSERT_TRUE(rk.solve(t_mesh, ics_p, sol_p));
    ASSERT_TRUE(rk.solve(t_mesh, ics_m, sol_m));
    Jx_fd.col(j) = (sol_p.x.col(sol_p.size() - 1) - sol_m.x.col(sol_m.size() - 1)) / (2.0*delta);
  }

  Real const error{(Jx - Jx_fd).norm() / std::max(Real(1.0), Jx_fd.norm())};
  EXPECT_LT(error, tolerance) << rk.name() << " on " << system->name() << std::endl
    << "Jx =" << std::endl << Jx << std::endl << "Jx_fd =" << std::endl << Jx_fd;
}

TEST(Sensitivity, PendulumExplicit) {
  auto sys = std::make_shared<SimplePendulumExplicit<Real>>();
  check_sensitivity<RK4<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
  check_sensitivity<GaussLegendre4<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
  check_sensitivity<RadauIIA5<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
  check_sensitivity<SSPIRK33<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
}

TEST(Sensitivity, PendulumImplicit) {
  auto sys = std::make_shared<SimplePendulumImplicit<Real>>();
  check_sensitivity<RK4<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
  check_sensitivity<GaussLegendre4<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
  check_sensitivity<RadauIIA5<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
  check_sensitivity<SSPIRK33<Real, 2, 0>, 2>(sys, sys->ics(), 2.0, 200, 1.0e-6);
}

TEST(Sensitivity, ThreeBodyExplicit) {
  auto sys = std::make_shared<ThreeBodyExplicit<Real>>();
  check_sensitivity<RK4<Real, 12, 0>, 12>(sys, sys->ics(), 0.5, 100, 1.0e-5);
  check_sensitivity<RadauIIA5<Real, 12, 0>, 12>(sys, sys->ics(), 0.5, 100, 1.0e-5);
}

TEST(Sensitivity, ThreeBodyImplicit) {
  auto sys = std::make_shared<ThreeBodyImplicit<Real>>();
  check_sensitivity<RK4<Real, 12, 0>, 12>(sys, sys->ics(), 0.5, 100, 1.0e-5);
  check_sensitivity<GaussLegendre4<Real, 12, 0>, 12>(sys, sys->ics(), 0.5, 100, 1.0e-5);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
