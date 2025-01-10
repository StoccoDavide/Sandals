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
TEST(Tableau, Chebyshev51)    {Chebyshev51Tableau tbl;    EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, ExplicitEuler)  {ExplicitEulerTableau tbl;  EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Fehlberg45)     {Fehlberg45Tableau tbl;     EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, GaussLegendre2) {GaussLegendre2Tableau tbl; EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, GaussLegendre4) {GaussLegendre4Tableau tbl; EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, GaussLegendre6) {GaussLegendre6Tableau tbl; EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Heun2)          {Heun2Tableau tbl;          EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Heun3)          {Heun3Tableau tbl;          EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, LobattoIIIA2)   {LobattoIIIA2Tableau tbl;   EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, ImplicitEuler)  {ImplicitEulerTableau tbl;  EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, MTE22)          {MTE22Tableau tbl;          EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, RadauIIA3)      {RadauIIA3Tableau tbl;      EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, RadauIIA5)      {RadauIIA5Tableau tbl;      EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Ralston2)       {Ralston2Tableau tbl;       EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Ralston3)       {Ralston3Tableau tbl;       EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, Ralston4)       {Ralston4Tableau tbl;       EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, RK4)            {RK4Tableau tbl;            EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPIRK33)       {SSPIRK33Tableau tbl;       EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK22)        {SSPRK22Tableau tbl;        EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK22star)    {SSPRK22starTableau tbl;    EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK33)        {SSPRK33Tableau tbl;        EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK42)        {SSPRK42Tableau tbl;        EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK43)        {SSPRK43Tableau tbl;        EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK93)        {SSPRK93Tableau tbl;        EXPECT_TRUE(tbl.check(true));}
TEST(Tableau, SSPRK104)       {SSPRK104Tableau tbl;       EXPECT_TRUE(tbl.check(true));}

// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}