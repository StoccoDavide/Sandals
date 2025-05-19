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

#include "SinCosImplicit.hh"
#include "SinCosExplicit.hh"

using namespace Sandals;

static const double tolerance{0.05};
static std::shared_ptr<SinCosImplicit> im_sys;
static std::shared_ptr<SinCosExplicit> ex_sys;
static Eigen::Vector2d ics;
static std::vector<Eigen::VectorXd> t;
static std::function<Eigen::MatrixXd(Eigen::VectorXd)> sol;
static const double SQRT_EPSILON{std::sqrt(std::numeric_limits<double>::epsilon())};

class GlobalTestEnvironment : public testing::Environment {
public:
  void SetUp() override {
    double t_ini{0.0}, t_end{10.0};
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
    im_sys = std::make_shared<SinCosImplicit>();
    ex_sys = std::make_shared<SinCosExplicit>();
    SANDALS_ASSERT((ex_sys->ics() - im_sys->ics()).norm() < SQRT_EPSILON, "Initial conditions are not equal.");
    ics = im_sys->ics();
    sol = [](Eigen::VectorXd t) -> Eigen::MatrixXd {return im_sys->analytical_solution(t);};
  }

  void TearDown() override {}
};

TEST(OrderImplicit, Chebyshev51)    {Chebyshev51<double, 2> rk(im_sys);    EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, ExplicitEuler)  {ExplicitEuler<double, 2> rk(im_sys);  EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, Fehlberg45)     {Fehlberg45<double, 2> rk(im_sys);     EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, GaussLegendre2) {GaussLegendre2<double, 2> rk(im_sys); EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, GaussLegendre4) {GaussLegendre4<double, 2> rk(im_sys); EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, GaussLegendre6) {GaussLegendre6<double, 2> rk(im_sys); EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, Heun2)          {Heun2<double, 2> rk(im_sys);          EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, Heun3)          {Heun3<double, 2> rk(im_sys);          EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, LobattoIIIA2)   {LobattoIIIA2<double, 2> rk(im_sys);   EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, ImplicitEuler)  {ImplicitEuler<double, 2> rk(im_sys);  EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, RadauIIA3)      {RadauIIA3<double, 2> rk(im_sys);      EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, RadauIIA5)      {RadauIIA5<double, 2> rk(im_sys);      EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, Ralston2)       {Ralston2<double, 2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, Ralston3)       {Ralston3<double, 2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, Ralston4)       {Ralston4<double, 2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, RK4)            {RK4<double, 2> rk(im_sys);            EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPIRK33)       {SSPIRK33<double, 2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK22)        {SSPRK22<double, 2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK22star)    {SSPRK22star<double, 2> rk(im_sys);    EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK33)        {SSPRK33<double, 2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK42)        {SSPRK42<double, 2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK43)        {SSPRK43<double, 2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK93)        {SSPRK93<double, 2> rk(im_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderImplicit, SSPRK104)       {SSPRK104<double, 2> rk(im_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}

TEST(OrderExplicit, Chebyshev51)    {Chebyshev51<double, 2> rk(ex_sys);    EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, ExplicitEuler)  {ExplicitEuler<double, 2> rk(ex_sys);  EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, Fehlberg45)     {Fehlberg45<double, 2> rk(ex_sys);     EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, GaussLegendre2) {GaussLegendre2<double, 2> rk(ex_sys); EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, GaussLegendre4) {GaussLegendre4<double, 2> rk(ex_sys); EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, GaussLegendre6) {GaussLegendre6<double, 2> rk(ex_sys); EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, Heun2)          {Heun2<double, 2> rk(ex_sys);          EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, Heun3)          {Heun3<double, 2> rk(ex_sys);          EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, LobattoIIIA2)   {LobattoIIIA2<double, 2> rk(ex_sys);   EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, ImplicitEuler)  {ImplicitEuler<double, 2> rk(ex_sys);  EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, RadauIIA3)      {RadauIIA3<double, 2> rk(ex_sys);      EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, RadauIIA5)      {RadauIIA5<double, 2> rk(ex_sys);      EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, Ralston2)       {Ralston2<double, 2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, Ralston3)       {Ralston3<double, 2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, Ralston4)       {Ralston4<double, 2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, RK4)            {RK4<double, 2> rk(ex_sys);            EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPIRK33)       {SSPIRK33<double, 2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK22)        {SSPRK22<double, 2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK22star)    {SSPRK22star<double, 2> rk(ex_sys);    EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK33)        {SSPRK33<double, 2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK42)        {SSPRK42<double, 2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK43)        {SSPRK43<double, 2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK93)        {SSPRK93<double, 2> rk(ex_sys);        EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}
TEST(OrderExplicit, SSPRK104)       {SSPRK104<double, 2> rk(ex_sys);       EXPECT_GE(rk.estimate_order(t, ics, sol), rk.order()-tolerance);}

// Run all the tests.
int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::AddGlobalTestEnvironment(new GlobalTestEnvironment());
  return RUN_ALL_TESTS();
}
