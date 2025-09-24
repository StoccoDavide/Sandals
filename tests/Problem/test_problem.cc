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

#include "Sandals.hh"

// #include "Sandals/RungeKutta/ExplicitEuler.hh" // ERK - OK
// #include "Sandals/RungeKutta/Heun2.hh" // ERK - OK
// #include "Sandals/RungeKutta/Heun3.hh" // ERK - OK
#include "Sandals/RungeKutta/RK4.hh" // ERK - OK
#include "Sandals/RungeKutta/SSPRK104.hh" // ERK - OK
// #include "Sandals/RungeKutta/LobattoIIIA2.hh" // DIRK - OK
// #include "Sandals/RungeKutta/SSPIRK33.hh" // DIRK - OK
#include "Sandals/RungeKutta/ImplicitEuler.hh" // IRK - OK
#include "Sandals/RungeKutta/GaussLegendre2.hh" // IRK
#include "Sandals/RungeKutta/GaussLegendre4.hh" // IRK
#include "Sandals/RungeKutta/GaussLegendre6.hh" // IRK
#include "Sandals/RungeKutta/RadauIIA3.hh" // IRK
#include "Sandals/RungeKutta/RadauIIA5.hh" // IRK

#include "Basic.hh"
#include "Poisson.hh"
#include "Shampine0.hh"
#include "Shampine1.hh"
#include "Shampine2.hh"
#include "Hypersensitive0.hh"
#include "Hypersensitive1.hh"

#ifdef SANDALS_ENABLE_PLOTTING
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TColor.h>
#include <TLegend.h>

template<typename VecX, typename VecY>
TGraph* to_TGraph(const VecX & x, const VecY & y) {
  if (x.size() != y.size()) {
    throw std::runtime_error("x and y vectors must have the same size");
  }
  return new TGraph(x.size(), x.data(), y.data());
}

std::vector<Color_t> matlab_lines_colormap() {
  std::vector<Color_t> colors;
  colors.push_back(TColor::GetColor(0, 114, 189));
  colors.push_back(TColor::GetColor(217, 83, 25));
  colors.push_back(TColor::GetColor(237, 177, 32));
  colors.push_back(TColor::GetColor(126, 47, 142));
  colors.push_back(TColor::GetColor(119, 172, 48));
  colors.push_back(TColor::GetColor(77, 190, 238));
  colors.push_back(TColor::GetColor(162, 20, 47));
  return colors;
}
#endif

using namespace Sandals;
using Real = double;

#ifndef PROBLEM_INIT
#define PROBLEM_INIT(PROBLEM, INTEGRATOR)                                                \
  constexpr Integer D{PROBLEM##Explicit<Real>::equations_number()};                      \
  PROBLEM##Problem<Real, PROBLEM##Explicit<Real>, INTEGRATOR<Real, D>> problem_explicit; \
  PROBLEM##Problem<Real, PROBLEM##Implicit<Real>, INTEGRATOR<Real, D>> problem_implicit; \
  PROBLEM##Problem<Real, PROBLEM##SemiExplicit<Real>, INTEGRATOR<Real, D>> problem_semiexplicit;
#endif

int main(int argc, char ** argv) {

#ifdef SANDALS_ENABLE_PLOTTING
  TApplication app("plot", &argc, argv);
#else
  // Prevent compiler warnings
  (void)argc;
  (void)argv;
#endif

  // Istantiate the problems
  PROBLEM_INIT(Poisson, GaussLegendre6)

  // Set verbose mode
  constexpr bool verbose{true};
  problem_explicit.verbose_mode(verbose);
  problem_implicit.verbose_mode(verbose);
  problem_semiexplicit.verbose_mode(verbose);
  problem_explicit.integrator()->verbose_mode(false);
  problem_implicit.integrator()->verbose_mode(false);
  problem_semiexplicit.integrator()->verbose_mode(false);

  // Set solver tolerance
  problem_explicit.tolerance(1.0e-8);
  problem_implicit.tolerance(1.0e-8);
  problem_semiexplicit.tolerance(1.0e-8);

  // Set solver maximum number of iterations
  problem_explicit.max_iterations(100);
  problem_implicit.max_iterations(100);
  problem_semiexplicit.max_iterations(100);

  // Set solution parameters
  constexpr Integer num_points{200};
  problem_explicit.subintervals(5);
  problem_implicit.subintervals(5);
  problem_semiexplicit.subintervals(5);

  // Set time mesh
  Eigen::Vector<Real, Eigen::Dynamic> time(Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
    num_points, problem_explicit.time_start(), problem_explicit.time_end()
  ));

  // Set initial guess
  Eigen::Matrix<Real, D, Eigen::Dynamic> guess(problem_explicit.guess(time));

  // Solve the problems with shooting
  std::cout << "\n=== Basic Explicit Problem ===\n" << std::endl;
  problem_explicit.multiple_shooting(time, guess);
  std::cout << "\n=== Basic Implicit Problem ===\n" << std::endl;
  problem_implicit.multiple_shooting(time, guess);
  std::cout << "\n=== Basic Semi-Explicit Problem ===\n" << std::endl;
  problem_semiexplicit.multiple_shooting(time, guess);

  #ifdef SANDALS_ENABLE_PLOTTING
  auto esol = problem_explicit.solution();
  auto isol = problem_implicit.solution();
  auto ssol = problem_semiexplicit.solution();
  auto asol = problem_explicit.analytical_solution(esol.t);

  auto colors = matlab_lines_colormap();

  TCanvas *canvas = new TCanvas("canvas", "Solution Comparison", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  TGraph *graph_ex = to_TGraph(esol.t, esol.eigen_x(0));
  TGraph *graph_ey = to_TGraph(esol.t, esol.eigen_x(1));
  TGraph *graph_ax = to_TGraph(esol.t, asol.row(0).transpose().eval());
  TGraph *graph_ay = to_TGraph(esol.t, asol.row(1).transpose().eval());
  graph_ex->SetTitle("Explicit vs Analytical");
  graph_ex->SetLineColor(colors[0]); graph_ex->Draw("AL");
  graph_ey->SetLineColor(colors[1]); graph_ey->Draw("L SAME");
  graph_ax->SetMarkerColor(colors[0]); graph_ax->SetMarkerStyle(8); graph_ax->Draw("P SAME");
  graph_ay->SetMarkerColor(colors[1]); graph_ay->SetMarkerStyle(8); graph_ay->Draw("P SAME");
  graph_ex->GetXaxis()->SetTitle("t (s)");
  graph_ex->GetYaxis()->SetTitle("x, y (-)");
  graph_ex->GetXaxis()->SetLimits(time(0), time(Eigen::last));
  graph_ex->GetYaxis()->SetRangeUser(
    std::min(asol.minCoeff(), esol.x.minCoeff())*1.1, std::max(asol.maxCoeff(), esol.x.maxCoeff())*1.1
  );
  TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg1->AddEntry(graph_ex, "x (explicit)", "l");
  leg1->AddEntry(graph_ey, "y (explicit)", "l");
  leg1->AddEntry(graph_ax, "x (analytical)", "p");
  leg1->AddEntry(graph_ay, "y (analytical)", "p");
  leg1->Draw();

  canvas->cd(2);
  TGraph *graph_ix = to_TGraph(isol.t, isol.eigen_x(0));
  TGraph *graph_iy = to_TGraph(isol.t, isol.eigen_x(1));
  graph_ix->SetTitle("Implicit vs Analytical");
  graph_ix->SetLineColor(colors[0]); graph_ix->SetLineStyle(9); graph_ix->Draw("AL");
  graph_iy->SetLineColor(colors[1]); graph_iy->SetLineStyle(9); graph_iy->Draw("L SAME");
  graph_ax->SetMarkerColor(colors[0]); graph_ax->Draw("P SAME");
  graph_ay->SetMarkerColor(colors[1]); graph_ay->Draw("P SAME");
  graph_ix->GetXaxis()->SetTitle("t (s)");
  graph_ix->GetYaxis()->SetTitle("x, y (-)");
  graph_ix->GetXaxis()->SetLimits(time(0), time(Eigen::last));
  graph_ix->GetYaxis()->SetRangeUser(
    std::min(asol.minCoeff(), esol.x.minCoeff())*1.1, std::max(asol.maxCoeff(), esol.x.maxCoeff())*1.1
  );
  TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg2->AddEntry(graph_ix, "x (implicit)", "l");
  leg2->AddEntry(graph_iy, "y (implicit)", "l");
  leg2->AddEntry(graph_ax, "x (analytical)", "p");
  leg2->AddEntry(graph_ay, "y (analytical)", "p");
  leg2->Draw();

  canvas->cd(3);
  TGraph *graph_sx = to_TGraph(ssol.t, ssol.eigen_x(0));
  TGraph *graph_sy = to_TGraph(ssol.t, ssol.eigen_x(1));
  graph_sx->SetTitle("Semi-Explicit vs Analytical");
  graph_sx->SetLineColor(colors[0]); graph_sx->SetLineStyle(7); graph_sx->Draw("AL");
  graph_sy->SetLineColor(colors[1]); graph_sy->SetLineStyle(7); graph_sy->Draw("L SAME");
  graph_ax->SetMarkerColor(colors[0]); graph_ax->Draw("P SAME");
  graph_ay->SetMarkerColor(colors[1]); graph_ay->Draw("P SAME");
  graph_sx->GetXaxis()->SetTitle("t (s)");
  graph_sx->GetYaxis()->SetTitle("x, y (-)");
  graph_sx->GetXaxis()->SetLimits(time(0), time(Eigen::last));
  graph_sx->GetYaxis()->SetRangeUser(
    std::min(asol.minCoeff(), esol.x.minCoeff())*1.1, std::max(asol.maxCoeff(), esol.x.maxCoeff())*1.1
  );
  TLegend *leg3 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg3->AddEntry(graph_sx, "x (semi-explicit)", "l");
  leg3->AddEntry(graph_sy, "y (semi-explicit)", "l");
  leg3->AddEntry(graph_ax, "x (analytical)", "p");
  leg3->AddEntry(graph_ay, "y (analytical)", "p");
  leg3->Draw();

  canvas->Update();
  app.Run();
  #endif

  return 0;
}
