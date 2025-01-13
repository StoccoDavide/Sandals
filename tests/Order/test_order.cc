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
#include "SineImplicit.hh"
#include "SineExplicit.hh"
#include <gtest/gtest.h>

using namespace Sandals;

static Real tol{0.05};
static std::shared_ptr<SineImplicit> im_sys;
static std::shared_ptr<SineExplicit> ex_sys;
static Vector2 ics;
static std::vector<Eigen::VectorXd> t;
static std::vector<Eigen::MatrixXd> x;

class GlobalTestEnvironment : public testing::Environment {
public:
  void SetUp() override {
    Real t_ini{0.0};
    Real t_end{100.0};
    t.resize(12);
    t[0]  = Eigen::VectorXd::LinSpaced(3000, t_ini, t_end);
    t[1]  = Eigen::VectorXd::LinSpaced(2750, t_ini, t_end);
    t[2]  = Eigen::VectorXd::LinSpaced(2500, t_ini, t_end);
    t[3]  = Eigen::VectorXd::LinSpaced(2250, t_ini, t_end);
    t[4]  = Eigen::VectorXd::LinSpaced(2000, t_ini, t_end);
    t[5]  = Eigen::VectorXd::LinSpaced(1750, t_ini, t_end);
    t[6]  = Eigen::VectorXd::LinSpaced(1500, t_ini, t_end);
    t[7]  = Eigen::VectorXd::LinSpaced(1250, t_ini, t_end);
    t[8]  = Eigen::VectorXd::LinSpaced(1000, t_ini, t_end);
    t[9]  = Eigen::VectorXd::LinSpaced(750,  t_ini, t_end);
    t[10] = Eigen::VectorXd::LinSpaced(500,  t_ini, t_end);
    t[11] = Eigen::VectorXd::LinSpaced(250,  t_ini, t_end);
    x.resize(12);
    for (Size i{0}; i < static_cast<Size>(t.size()); ++i) {x[i] = im_sys->analytical_solution(t[i]);}
    im_sys = std::make_shared<SineImplicit>();
    ex_sys = std::make_shared<SineExplicit>();
    SANDALS_ASSERT((ex_sys->ics() - im_sys->ics()).norm() < SQRT_EPSILON, "Initial conditions are not equal.");
    ics = im_sys->ics();
  }

  void TearDown() override {}
};

TEST(OrderImplicit, Chebyshev51)    {Chebyshev51<2> rk(im_sys);    EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, ExplicitEuler)  {ExplicitEuler<2> rk(im_sys);  EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, Fehlberg45)     {Fehlberg45<2> rk(im_sys);     EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, GaussLegendre2) {GaussLegendre2<2> rk(im_sys); EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, GaussLegendre4) {GaussLegendre4<2> rk(im_sys); EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, GaussLegendre6) {GaussLegendre6<2> rk(im_sys); EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, Heun2)          {Heun2<2> rk(im_sys);          EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, Heun3)          {Heun3<2> rk(im_sys);          EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, LobattoIIIA2)   {LobattoIIIA2<2> rk(im_sys);   EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, ImplicitEuler)  {ImplicitEuler<2> rk(im_sys);  EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, MTE22)          {MTE22<2> rk(im_sys);          EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, RadauIIA3)      {RadauIIA3<2> rk(im_sys);      EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, RadauIIA5)      {RadauIIA5<2> rk(im_sys);      EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, Ralston2)       {Ralston2<2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, Ralston3)       {Ralston3<2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, Ralston4)       {Ralston4<2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, RK4)            {RK4<2> rk(im_sys);            EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPIRK33)       {SSPIRK33<2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK22)        {SSPRK22<2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK22star)    {SSPRK22star<2> rk(im_sys);    EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK33)        {SSPRK33<2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK42)        {SSPRK42<2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK43)        {SSPRK43<2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK93)        {SSPRK93<2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderImplicit, SSPRK104)       {SSPRK104<2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}

TEST(OrderExplicit, Chebyshev51)    {Chebyshev51<2> rk(ex_sys);    EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, ExplicitEuler)  {ExplicitEuler<2> rk(ex_sys);  EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, Fehlberg45)     {Fehlberg45<2> rk(ex_sys);     EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, GaussLegendre2) {GaussLegendre2<2> rk(ex_sys); EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, GaussLegendre4) {GaussLegendre4<2> rk(ex_sys); EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, GaussLegendre6) {GaussLegendre6<2> rk(ex_sys); EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, Heun2)          {Heun2<2> rk(ex_sys);          EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, Heun3)          {Heun3<2> rk(ex_sys);          EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, LobattoIIIA2)   {LobattoIIIA2<2> rk(ex_sys);   EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, ImplicitEuler)  {ImplicitEuler<2> rk(ex_sys);  EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, MTE22)          {MTE22<2> rk(ex_sys);          EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, RadauIIA3)      {RadauIIA3<2> rk(ex_sys);      EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, RadauIIA5)      {RadauIIA5<2> rk(ex_sys);      EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, Ralston2)       {Ralston2<2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, Ralston3)       {Ralston3<2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, Ralston4)       {Ralston4<2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, RK4)            {RK4<2> rk(ex_sys);            EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPIRK33)       {SSPIRK33<2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK22)        {SSPRK22<2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK22star)    {SSPRK22star<2> rk(ex_sys);    EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK33)        {SSPRK33<2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK42)        {SSPRK42<2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK43)        {SSPRK43<2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK93)        {SSPRK93<2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}
TEST(OrderExplicit, SSPRK104)       {SSPRK104<2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, x), rk.order()-tol);}

// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::AddGlobalTestEnvironment(new GlobalTestEnvironment());
  return RUN_ALL_TESTS();
}
