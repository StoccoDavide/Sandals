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
#include "SinCosSemiExplicit.hh"

using namespace Sandals;
using Real = double;

static const Real tolerance{0.05};
static Eigen::Vector2d ics;
static std::vector<Eigen::VectorXd> t;
static std::function<Eigen::MatrixXd(Eigen::VectorXd)> sol;
static const Real SQRT_EPSILON{std::sqrt(std::numeric_limits<Real>::epsilon())};

class GlobalTestEnvironment : public testing::Environment {
public:
  void SetUp() override {
    Real t_ini{0.0}, t_end{10.0};
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
    SANDALS_ASSERT(
      (SinCosImplicit<Real>::ics() - SinCosExplicit<Real>::ics()).norm() < SQRT_EPSILON &&
      (SinCosImplicit<Real>::ics() - SinCosSemiExplicit<Real>::ics()).norm() < SQRT_EPSILON,
      "Initial conditions are not equal.");
    sol = [] (Eigen::VectorXd t) -> Eigen::MatrixXd {return SinCosImplicit<Real>::analytical_solution(t);};
  }

  void TearDown() override {}
};

#ifndef ORDER_IMPLICIT_SYSTEM
#define ORDER_IMPLICIT_SYSTEM(RK) \
  RK<Real, 2, 0> rk(std::make_unique<SinCosImplicit<Real>>()); \
  EXPECT_GE(rk.estimate_order(t, SinCosImplicit<Real>::ics(), sol), rk.order() - tolerance);
#endif

#ifndef ORDER_EXPLICIT_SYSTEM
#define ORDER_EXPLICIT_SYSTEM(RK) \
  RK<Real, 2, 0> rk(std::make_unique<SinCosExplicit<Real>>()); \
  EXPECT_GE(rk.estimate_order(t, SinCosExplicit<Real>::ics(), sol), rk.order() - tolerance);
#endif

#ifndef ORDER_SEMIEXPLICIT_SYSTEM
#define ORDER_SEMIEXPLICIT_SYSTEM(RK) \
  RK<Real, 2, 0> rk(std::make_unique<SinCosSemiExplicit<Real>>()); \
  EXPECT_GE(rk.estimate_order(t, SinCosSemiExplicit<Real>::ics(), sol), rk.order() - tolerance);
#endif

TEST(OrderImplicit, Chebyshev51)    {ORDER_IMPLICIT_SYSTEM(Chebyshev51)}
TEST(OrderImplicit, ExplicitEuler)  {ORDER_IMPLICIT_SYSTEM(ExplicitEuler)}
TEST(OrderImplicit, Fehlberg45)     {ORDER_IMPLICIT_SYSTEM(Fehlberg45)}
TEST(OrderImplicit, GaussLegendre2) {ORDER_IMPLICIT_SYSTEM(GaussLegendre2)}
TEST(OrderImplicit, GaussLegendre4) {ORDER_IMPLICIT_SYSTEM(GaussLegendre4)}
TEST(OrderImplicit, GaussLegendre6) {ORDER_IMPLICIT_SYSTEM(GaussLegendre6)}
TEST(OrderImplicit, Heun2)          {ORDER_IMPLICIT_SYSTEM(Heun2)}
TEST(OrderImplicit, Heun3)          {ORDER_IMPLICIT_SYSTEM(Heun3)}
TEST(OrderImplicit, LobattoIIIA2)   {ORDER_IMPLICIT_SYSTEM(LobattoIIIA2)}
TEST(OrderImplicit, ImplicitEuler)  {ORDER_IMPLICIT_SYSTEM(ImplicitEuler)}
TEST(OrderImplicit, RadauIIA3)      {ORDER_IMPLICIT_SYSTEM(RadauIIA3)}
TEST(OrderImplicit, RadauIIA5)      {ORDER_IMPLICIT_SYSTEM(RadauIIA5)}
TEST(OrderImplicit, Ralston2)       {ORDER_IMPLICIT_SYSTEM(Ralston2)}
TEST(OrderImplicit, Ralston3)       {ORDER_IMPLICIT_SYSTEM(Ralston3)}
TEST(OrderImplicit, Ralston4)       {ORDER_IMPLICIT_SYSTEM(Ralston4)}
TEST(OrderImplicit, RK4)            {ORDER_IMPLICIT_SYSTEM(RK4)}
TEST(OrderImplicit, SSPIRK33)       {ORDER_IMPLICIT_SYSTEM(SSPIRK33)}
TEST(OrderImplicit, SSPRK22)        {ORDER_IMPLICIT_SYSTEM(SSPRK22)}
TEST(OrderImplicit, SSPRK22star)    {ORDER_IMPLICIT_SYSTEM(SSPRK22star)}
TEST(OrderImplicit, SSPRK33)        {ORDER_IMPLICIT_SYSTEM(SSPRK33)}
TEST(OrderImplicit, SSPRK42)        {ORDER_IMPLICIT_SYSTEM(SSPRK42)}
TEST(OrderImplicit, SSPRK43)        {ORDER_IMPLICIT_SYSTEM(SSPRK43)}
TEST(OrderImplicit, SSPRK93)        {ORDER_IMPLICIT_SYSTEM(SSPRK93)}
TEST(OrderImplicit, SSPRK104)       {ORDER_IMPLICIT_SYSTEM(SSPRK104)}

TEST(OrderExplicit, Chebyshev51)    {ORDER_EXPLICIT_SYSTEM(Chebyshev51)}
TEST(OrderExplicit, ExplicitEuler)  {ORDER_EXPLICIT_SYSTEM(ExplicitEuler)}
TEST(OrderExplicit, Fehlberg45)     {ORDER_EXPLICIT_SYSTEM(Fehlberg45)}
TEST(OrderExplicit, GaussLegendre2) {ORDER_EXPLICIT_SYSTEM(GaussLegendre2)}
TEST(OrderExplicit, GaussLegendre4) {ORDER_EXPLICIT_SYSTEM(GaussLegendre4)}
TEST(OrderExplicit, GaussLegendre6) {ORDER_EXPLICIT_SYSTEM(GaussLegendre6)}
TEST(OrderExplicit, Heun2)          {ORDER_EXPLICIT_SYSTEM(Heun2)}
TEST(OrderExplicit, Heun3)          {ORDER_EXPLICIT_SYSTEM(Heun3)}
TEST(OrderExplicit, LobattoIIIA2)   {ORDER_EXPLICIT_SYSTEM(LobattoIIIA2)}
TEST(OrderExplicit, ImplicitEuler)  {ORDER_EXPLICIT_SYSTEM(ImplicitEuler)}
TEST(OrderExplicit, RadauIIA3)      {ORDER_EXPLICIT_SYSTEM(RadauIIA3)}
TEST(OrderExplicit, RadauIIA5)      {ORDER_EXPLICIT_SYSTEM(RadauIIA5)}
TEST(OrderExplicit, Ralston2)       {ORDER_EXPLICIT_SYSTEM(Ralston2)}
TEST(OrderExplicit, Ralston3)       {ORDER_EXPLICIT_SYSTEM(Ralston3)}
TEST(OrderExplicit, Ralston4)       {ORDER_EXPLICIT_SYSTEM(Ralston4)}
TEST(OrderExplicit, RK4)            {ORDER_EXPLICIT_SYSTEM(RK4)}
TEST(OrderExplicit, SSPIRK33)       {ORDER_EXPLICIT_SYSTEM(SSPIRK33)}
TEST(OrderExplicit, SSPRK22)        {ORDER_EXPLICIT_SYSTEM(SSPRK22)}
TEST(OrderExplicit, SSPRK22star)    {ORDER_EXPLICIT_SYSTEM(SSPRK22star)}
TEST(OrderExplicit, SSPRK33)        {ORDER_EXPLICIT_SYSTEM(SSPRK33)}
TEST(OrderExplicit, SSPRK42)        {ORDER_EXPLICIT_SYSTEM(SSPRK42)}
TEST(OrderExplicit, SSPRK43)        {ORDER_EXPLICIT_SYSTEM(SSPRK43)}
TEST(OrderExplicit, SSPRK93)        {ORDER_EXPLICIT_SYSTEM(SSPRK93)}
TEST(OrderExplicit, SSPRK104)       {ORDER_EXPLICIT_SYSTEM(SSPRK104)}

TEST(OrderSemiExplicit, Chebyshev51)    {ORDER_SEMIEXPLICIT_SYSTEM(Chebyshev51)}
TEST(OrderSemiExplicit, ExplicitEuler)  {ORDER_SEMIEXPLICIT_SYSTEM(ExplicitEuler)}
TEST(OrderSemiExplicit, Fehlberg45)     {ORDER_SEMIEXPLICIT_SYSTEM(Fehlberg45)}
TEST(OrderSemiExplicit, GaussLegendre2) {ORDER_SEMIEXPLICIT_SYSTEM(GaussLegendre2)}
TEST(OrderSemiExplicit, GaussLegendre4) {ORDER_SEMIEXPLICIT_SYSTEM(GaussLegendre4)}
TEST(OrderSemiExplicit, GaussLegendre6) {ORDER_SEMIEXPLICIT_SYSTEM(GaussLegendre6)}
TEST(OrderSemiExplicit, Heun2)          {ORDER_SEMIEXPLICIT_SYSTEM(Heun2)}
TEST(OrderSemiExplicit, Heun3)          {ORDER_SEMIEXPLICIT_SYSTEM(Heun3)}
TEST(OrderSemiExplicit, LobattoIIIA2)   {ORDER_SEMIEXPLICIT_SYSTEM(LobattoIIIA2)}
TEST(OrderSemiExplicit, ImplicitEuler)  {ORDER_SEMIEXPLICIT_SYSTEM(ImplicitEuler)}
TEST(OrderSemiExplicit, RadauIIA3)      {ORDER_SEMIEXPLICIT_SYSTEM(RadauIIA3)}
TEST(OrderSemiExplicit, RadauIIA5)      {ORDER_SEMIEXPLICIT_SYSTEM(RadauIIA5)}
TEST(OrderSemiExplicit, Ralston2)       {ORDER_SEMIEXPLICIT_SYSTEM(Ralston2)}
TEST(OrderSemiExplicit, Ralston3)       {ORDER_SEMIEXPLICIT_SYSTEM(Ralston3)}
TEST(OrderSemiExplicit, Ralston4)       {ORDER_SEMIEXPLICIT_SYSTEM(Ralston4)}
TEST(OrderSemiExplicit, RK4)            {ORDER_SEMIEXPLICIT_SYSTEM(RK4)}
TEST(OrderSemiExplicit, SSPIRK33)       {ORDER_SEMIEXPLICIT_SYSTEM(SSPIRK33)}
TEST(OrderSemiExplicit, SSPRK22)        {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK22)}
TEST(OrderSemiExplicit, SSPRK22star)    {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK22star)}
TEST(OrderSemiExplicit, SSPRK33)        {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK33)}
TEST(OrderSemiExplicit, SSPRK42)        {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK42)}
TEST(OrderSemiExplicit, SSPRK43)        {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK43)}
TEST(OrderSemiExplicit, SSPRK93)        {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK93)}
TEST(OrderSemiExplicit, SSPRK104)       {ORDER_SEMIEXPLICIT_SYSTEM(SSPRK104)}

// Run all the tests
int main(int argc, char ** argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::AddGlobalTestEnvironment(new GlobalTestEnvironment());
  return RUN_ALL_TESTS();
}
