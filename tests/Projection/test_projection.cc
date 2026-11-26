/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Projection of the states on the invariants: of all the states (project), and of a selection of the
// states on a selection of the invariants (project_ics).

#include <gtest/gtest.h>

#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;
using Real = double;

// A rotation in the (x0, x1) plane and a decay of x2, with the invariants h0 = x0^2 + x1^2 - 1 (the
// unit circle) and h1 = x2 - x0 (to have a second, independent invariant)
using System = ExplicitWrapper<Real, 3, 2>;

std::shared_ptr<System> rotation()
{
  return std::make_shared<System>(
    [] (System::VectorF const & x, Real) {return System::VectorF(-x(1), x(0), -x(2));},
    [] (System::VectorF const &, Real) {System::MatrixJF J; J << 0, -1, 0, 1, 0, 0, 0, 0, -1; return J;},
    [] (System::VectorF const & x, Real) {return System::VectorH(x(0)*x(0) + x(1)*x(1) - 1.0, x(2) - x(0));},
    [] (System::VectorF const & x, Real) {System::MatrixJH J; J << 2*x(0), 2*x(1), 0, -1, 0, 1; return J;});
}

TEST(Projection, AllStatesOnAllInvariants) {
  RadauIIA5<Real, 3, 2> integrator(rotation());
  System::VectorF const x(1.1, 0.2, 0.5);
  System::VectorF x_projected;
  ASSERT_TRUE(integrator.project(x, 0.0, x_projected));
  EXPECT_NEAR(x_projected(0)*x_projected(0) + x_projected(1)*x_projected(1), 1.0, 1.0e-10);
  EXPECT_NEAR(x_projected(2), x_projected(0), 1.0e-10);
}

TEST(Projection, SelectedStatesOnSelectedInvariants) {
  RadauIIA5<Real, 3, 2> integrator(rotation());
  System::VectorF const x(1.1, 0.2, 0.5);
  System::VectorF x_projected;

  // Only (x0, x1) on the circle: x2 is not touched, and the closest point is along the radius
  ASSERT_TRUE(integrator.project_ics(x, 0.0, {0, 1}, {0}, x_projected));
  EXPECT_NEAR(x_projected(0)*x_projected(0) + x_projected(1)*x_projected(1), 1.0, 1.0e-10);
  EXPECT_EQ(x_projected(2), x(2));
  EXPECT_NEAR(x_projected(1)/x_projected(0), x(1)/x(0), 1.0e-10);

  // Only x2 on the second invariant: x2 = x0, the others are not touched
  ASSERT_TRUE(integrator.project_ics(x, 0.0, {2}, {1}, x_projected));
  EXPECT_EQ(x_projected.head<2>(), x.head<2>());
  EXPECT_NEAR(x_projected(2), x(0), 1.0e-12);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
