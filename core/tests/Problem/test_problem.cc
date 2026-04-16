/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <memory>

#include "Sandals.hh"
#include "Sandals/RungeKutta/ExplicitEuler.hh"  // ERK - OK
// #include "Sandals/RungeKutta/Heun2.hh" // ERK - OK
// #include "Sandals/RungeKutta/Heun3.hh" // ERK - OK
#include "Basic.hh"
#include "FalknerSkan.hh"
#include "Hypersensitive0.hh"
#include "Hypersensitive1.hh"
#include "Poisson.hh"
#include "Sandals/RungeKutta/GaussLegendre2.hh"  // IRK
#include "Sandals/RungeKutta/GaussLegendre4.hh"  // IRK
#include "Sandals/RungeKutta/GaussLegendre6.hh"  // IRK
#include "Sandals/RungeKutta/ImplicitEuler.hh"   // IRK - OK
#include "Sandals/RungeKutta/LobattoIIIA2.hh"    // DIRK - OK
#include "Sandals/RungeKutta/RK4.hh"             // ERK - OK
#include "Sandals/RungeKutta/RadauIIA3.hh"       // IRK
#include "Sandals/RungeKutta/RadauIIA5.hh"       // IRK
#include "Sandals/RungeKutta/SSPIRK33.hh"        // DIRK - OK
#include "Sandals/RungeKutta/SSPRK104.hh"        // ERK - OK
#include "Shampine0.hh"
#include "Shampine1.hh"
#include "Shock.hh"
#include "Troesch.hh"

#ifdef SANDALS_ENABLE_PLOTTING
#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

template <typename VecX, typename VecY>
TGraph *to_TGraph(const VecX &x, const VecY &y) {
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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR)                                  \
  constexpr Integer D{PROBLEM##Explicit<Real>::equations()};               \
  PROBLEM##Problem<Real, PROBLEM##Explicit<Real>, INTEGRATOR<Real, D>>     \
      problem_explicit;                                                    \
  PROBLEM##Problem<Real, PROBLEM##Implicit<Real>, INTEGRATOR<Real, D>>     \
      problem_implicit;                                                    \
  PROBLEM##Problem<Real, PROBLEM##SemiExplicit<Real>, INTEGRATOR<Real, D>> \
      problem_semiexplicit;
#endif

int main(int argc, char **argv) {
#ifdef SANDALS_ENABLE_PLOTTING
  TApplication app("plot", &argc, argv);
#else
  // Prevent compiler warnings
  (void)argc;
  (void)argv;
#endif

  // Istantiate the problems
  PROBLEM_INIT(FalknerSkan, SSPRK104)

  // Set verbose mode
  constexpr bool verbose{true};
  problem_explicit.verbose_mode(verbose);
  problem_implicit.verbose_mode(verbose);
  problem_semiexplicit.verbose_mode(verbose);
  problem_explicit.integrator()->verbose_mode(false);
  problem_implicit.integrator()->verbose_mode(false);
  problem_semiexplicit.integrator()->verbose_mode(false);

  // Set reverse mode
  constexpr bool reverse{true};
  problem_explicit.integrator()->reverse_mode(reverse);
  problem_implicit.integrator()->reverse_mode(reverse);
  problem_semiexplicit.integrator()->reverse_mode(reverse);

  // Set solver tolerance
  problem_explicit.tolerance(1.0e-8);
  problem_implicit.tolerance(1.0e-8);
  problem_semiexplicit.tolerance(1.0e-8);

  // Set solver maximum number of iterations
  problem_explicit.max_iterations(50);
  problem_implicit.max_iterations(50);
  problem_semiexplicit.max_iterations(50);

  // Set solution parameters
  constexpr Integer num_subintervals{2};
  problem_explicit.subintervals(num_subintervals);
  problem_implicit.subintervals(num_subintervals);
  problem_semiexplicit.subintervals(num_subintervals);

  // Set time mesh
  constexpr Integer num_points{13};
  Eigen::Vector<Real, Eigen::Dynamic> time(
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
          num_points,
          problem_explicit.time_start(),
          problem_explicit.time_end()));

  // Set initial guess
  Eigen::Matrix<Real, D, Eigen::Dynamic> guess(problem_explicit.guess(time));

  // Solve the problems with shooting
  std::cout << std::endl << "=== Explicit Problem" << std::endl << std::endl;
  problem_explicit.multiple_shooting(time, guess);
  std::cout << std::endl << "=== Implicit Problem" << std::endl << std::endl;
  problem_implicit.multiple_shooting(time, guess);
  std::cout << std::endl
            << "=== Semi-Explicit Problem" << std::endl
            << std::endl;
  problem_semiexplicit.multiple_shooting(time, guess);

#ifdef SANDALS_ENABLE_PLOTTING
  auto sol_e = problem_explicit.solution();
  auto sol_i = problem_implicit.solution();
  auto sol_s = problem_semiexplicit.solution();
  auto sol_a = problem_explicit.exact_solution(
      reverse ? sol_e.t.reverse().eval() : sol_e.t);

  auto colors = matlab_lines_colormap();

  TCanvas *canvas = new TCanvas("canvas", "Solution Comparison", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  TGraph *graph_e0 = to_TGraph(sol_e.t, sol_e.eigen_x(0));
  TGraph *graph_e1 = to_TGraph(sol_e.t, sol_e.eigen_x(1));
  TGraph *graph_a0 = to_TGraph(sol_e.t, sol_a.row(0).transpose().eval());
  TGraph *graph_a1 = to_TGraph(sol_e.t, sol_a.row(1).transpose().eval());
  graph_e0->SetTitle("Explicit vs Exact");
  graph_e0->SetLineColor(colors[0]);
  graph_e0->Draw("AL");
  graph_a0->SetLineColor(colors[1]);
  graph_a0->Draw("L SAME");
  graph_e1->SetLineColor(colors[0]);
  graph_e1->SetLineStyle(2);
  graph_e1->Draw("L SAME");
  graph_a1->SetLineColor(colors[1]);
  graph_a1->SetLineStyle(2);
  graph_a1->Draw("L SAME");
  graph_e0->GetXaxis()->SetTitle("t");
  graph_e0->GetYaxis()->SetTitle("x");
  graph_e0->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  TLegend *legend_1 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend_1->AddEntry(graph_e0, "Numerical", "l");
  legend_1->AddEntry(graph_a0, "Exact", "l");
  legend_1->Draw();

  canvas->cd(2);
  TGraph *graph_i0 = to_TGraph(sol_i.t, sol_i.eigen_x(0));
  TGraph *graph_i1 = to_TGraph(sol_i.t, sol_i.eigen_x(1));
  graph_i0->SetTitle("Implicit vs Exact");
  graph_i0->SetLineColor(colors[0]);
  graph_i0->Draw("AL");
  graph_a0->SetLineColor(colors[1]);
  graph_a0->Draw("L SAME");
  graph_i1->SetLineColor(colors[0]);
  graph_i1->SetLineStyle(2);
  graph_i1->Draw("L SAME");
  graph_a1->SetLineColor(colors[1]);
  graph_a1->SetLineStyle(2);
  graph_a1->Draw("L SAME");
  graph_i0->GetXaxis()->SetTitle("t");
  graph_i0->GetYaxis()->SetTitle("x");
  graph_i0->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  TLegend *legend_2 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend_2->AddEntry(graph_i0, "Numerical", "l");
  legend_2->AddEntry(graph_a0, "Exact", "l");
  legend_2->Draw();

  canvas->cd(3);
  TGraph *graph_s0 = to_TGraph(sol_s.t, sol_s.eigen_x(0));
  TGraph *graph_s1 = to_TGraph(sol_s.t, sol_s.eigen_x(1));
  graph_s0->SetTitle("Semi-Explicit vs Exact");
  graph_s0->SetLineColor(colors[0]);
  graph_s0->Draw("AL");
  graph_a0->SetLineColor(colors[1]);
  graph_a0->Draw("L SAME");
  graph_s1->SetLineColor(colors[0]);
  graph_s1->SetLineStyle(2);
  graph_s1->Draw("L SAME");
  graph_a1->SetLineColor(colors[1]);
  graph_a1->SetLineStyle(2);
  graph_a1->Draw("L SAME");
  graph_s0->GetXaxis()->SetTitle("t");
  graph_s0->GetYaxis()->SetTitle("x");
  graph_s0->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  TLegend *legend_3 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend_3->AddEntry(graph_s0, "Numerical", "l");
  legend_3->AddEntry(graph_a0, "Exact", "l");
  legend_3->Draw();

  canvas->Update();
  app.Run();
#endif

  return 0;
}
