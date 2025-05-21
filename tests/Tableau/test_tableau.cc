/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>

#include "Sandals.hh"

#include "Sandals/RungeKutta/Chebyshev51.hh"
#include "Sandals/RungeKutta/ExplicitEuler.hh"
#include "Sandals/RungeKutta/Fehlberg45.hh"
#include "Sandals/RungeKutta/GaussLegendre2.hh"
#include "Sandals/RungeKutta/GaussLegendre4.hh"
#include "Sandals/RungeKutta/GaussLegendre6.hh"
#include "Sandals/RungeKutta/Heun2.hh"
#include "Sandals/RungeKutta/Heun3.hh"
#include "Sandals/RungeKutta/LobattoIIIA2.hh"
#include "Sandals/RungeKutta/ImplicitEuler.hh"
#include "Sandals/RungeKutta/MTE22.hh"
#include "Sandals/RungeKutta/RadauIIA3.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/Ralston2.hh"
#include "Sandals/RungeKutta/Ralston3.hh"
#include "Sandals/RungeKutta/Ralston4.hh"
#include "Sandals/RungeKutta/RK4.hh"
#include "Sandals/RungeKutta/SSPIRK33.hh"
#include "Sandals/RungeKutta/SSPRK22.hh"
#include "Sandals/RungeKutta/SSPRK22star.hh"
#include "Sandals/RungeKutta/SSPRK33.hh"
#include "Sandals/RungeKutta/SSPRK42.hh"
#include "Sandals/RungeKutta/SSPRK43.hh"
#include "Sandals/RungeKutta/SSPRK93.hh"
#include "Sandals/RungeKutta/SSPRK104.hh"

using namespace Sandals;

using Real = double;

TEST(Tableau, Chebyshev51)    {Chebyshev51Tableau<Real> tbl;    EXPECT_TRUE(tbl.check());}
TEST(Tableau, ExplicitEuler)  {ExplicitEulerTableau<Real> tbl;  EXPECT_TRUE(tbl.check());}
TEST(Tableau, Fehlberg45)     {Fehlberg45Tableau<Real> tbl;     EXPECT_TRUE(tbl.check());}
TEST(Tableau, GaussLegendre2) {GaussLegendre2Tableau<Real> tbl; EXPECT_TRUE(tbl.check());}
TEST(Tableau, GaussLegendre4) {GaussLegendre4Tableau<Real> tbl; EXPECT_TRUE(tbl.check());}
TEST(Tableau, GaussLegendre6) {GaussLegendre6Tableau<Real> tbl; EXPECT_TRUE(tbl.check());}
TEST(Tableau, Heun2)          {Heun2Tableau<Real> tbl;          EXPECT_TRUE(tbl.check());}
TEST(Tableau, Heun3)          {Heun3Tableau<Real> tbl;          EXPECT_TRUE(tbl.check());}
TEST(Tableau, LobattoIIIA2)   {LobattoIIIA2Tableau<Real> tbl;   EXPECT_TRUE(tbl.check());}
TEST(Tableau, ImplicitEuler)  {ImplicitEulerTableau<Real> tbl;  EXPECT_TRUE(tbl.check());}
TEST(Tableau, RadauIIA3)      {RadauIIA3Tableau<Real> tbl;      EXPECT_TRUE(tbl.check());}
TEST(Tableau, RadauIIA5)      {RadauIIA5Tableau<Real> tbl;      EXPECT_TRUE(tbl.check());}
TEST(Tableau, Ralston2)       {Ralston2Tableau<Real> tbl;       EXPECT_TRUE(tbl.check());}
TEST(Tableau, Ralston3)       {Ralston3Tableau<Real> tbl;       EXPECT_TRUE(tbl.check());}
TEST(Tableau, Ralston4)       {Ralston4Tableau<Real> tbl;       EXPECT_TRUE(tbl.check());}
TEST(Tableau, RK4)            {RK4Tableau<Real> tbl;            EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPIRK33)       {SSPIRK33Tableau<Real> tbl;       EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK22)        {SSPRK22Tableau<Real> tbl;        EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK22star)    {SSPRK22starTableau<Real> tbl;    EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK33)        {SSPRK33Tableau<Real> tbl;        EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK42)        {SSPRK42Tableau<Real> tbl;        EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK43)        {SSPRK43Tableau<Real> tbl;        EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK93)        {SSPRK93Tableau<Real> tbl;        EXPECT_TRUE(tbl.check());}
TEST(Tableau, SSPRK104)       {SSPRK104Tableau<Real> tbl;       EXPECT_TRUE(tbl.check());}

// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
