/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "Sandals.hh"
#include <gtest/gtest.h>

using namespace Sandals;

// Explicit methods
TEST(Tableau, ExplicitEuler)  {ExplicitEulerTableau tbl;  EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Heun2)          {Heun2Tableau tbl;          EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Heun3)          {Heun3Tableau tbl;          EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, RK4)            {RK4Tableau tbl;            EXPECT_TRUE(tbl.check(true));}

// Implicit methods
TEST(Tableau, ImplicitEuler)  {ImplicitEulerTableau tbl;  EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, RadauIIA3)      {RadauIIA3Tableau tbl;      EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, RadauIIA5)      {RadauIIA5Tableau tbl;      EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, GaussLegendre2) {GaussLegendre2Tableau tbl; EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, GaussLegendre4) {GaussLegendre4Tableau tbl; EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, GaussLegendre6) {GaussLegendre6Tableau tbl; EXPECT_TRUE(tbl.check(true));}

// Embedded explicit methods
TEST(Tableau, RKF45)          {RKF45Tableau tbl;          EXPECT_TRUE(tbl.check(true));}

// Embedded implicit methods


// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}