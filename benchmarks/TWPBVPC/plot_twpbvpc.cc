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
#include "Sandals/RungeKutta/GaussLegendre6.hh"

#include "TWPBVPC.hh"

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
    throw std::runtime_error("vectors must have the same size");
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
  PROBLEM_INIT(BVPT4, GaussLegendre6)

  // Set verbose mode
  constexpr bool verbose{false};
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
  problem_explicit.max_iterations(100);
  problem_implicit.max_iterations(100);
  problem_semiexplicit.max_iterations(100);

  // Set solution parameters
  constexpr Integer num_subintervals{2};
  problem_explicit.subintervals(num_subintervals);
  problem_implicit.subintervals(num_subintervals);
  problem_semiexplicit.subintervals(num_subintervals);

  // Set time mesh
  constexpr Integer num_points{100};
  Eigen::Vector<Real, Eigen::Dynamic> time(Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
    num_points, problem_explicit.time_start(), problem_explicit.time_end()
  ));

  // Set initial guess
  Eigen::Matrix<Real, D, Eigen::Dynamic> guess(problem_explicit.guess(time));

  // Solve the problems with shooting
  problem_explicit.multiple_shooting(time, guess);
  problem_implicit.multiple_shooting(time, guess);
  problem_semiexplicit.multiple_shooting(time, guess);

  #ifdef SANDALS_ENABLE_PLOTTING
  auto sol_e = problem_explicit.solution();
  auto sol_i = problem_implicit.solution();
  auto sol_s = problem_semiexplicit.solution();
  auto sol_a = problem_explicit.exact_solution(reverse ? sol_e.t.reverse().eval() : sol_e.t);

  auto colors = matlab_lines_colormap();

  TCanvas *canvas = new TCanvas("canvas", "Solution Comparison", 1200, 400);
  canvas->Divide(3, 1);

  canvas->cd(1);
  TGraph *graph_e = to_TGraph(sol_e.t, sol_e.eigen_x(0));
  TGraph *graph_a = to_TGraph(sol_e.t, sol_a);
  graph_e->SetTitle("Explicit vs Exact");
  graph_e->SetLineColor(colors[0]); graph_e->Draw("AL");
  graph_a->SetLineColor(colors[1]); graph_a->Draw("L SAME");
  graph_e->GetXaxis()->SetTitle("t");
  graph_e->GetYaxis()->SetTitle("x");
  graph_e->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  TLegend *legend_1 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend_1->AddEntry(graph_e, "Numerical", "l");
  legend_1->AddEntry(graph_a, "Exact", "l");
  legend_1->Draw();

  canvas->cd(2);
  TGraph *graph_i = to_TGraph(sol_i.t, sol_i.eigen_x(0));
  graph_i->SetTitle("Implicit vs Exact");
  graph_i->SetLineColor(colors[0]); graph_i->Draw("AL");
  graph_a->SetLineColor(colors[1]); graph_a->Draw("L SAME");
  graph_i->GetXaxis()->SetTitle("t");
  graph_i->GetYaxis()->SetTitle("x");
  graph_i->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  TLegend *legend_2 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend_2->AddEntry(graph_i, "Numerical", "l");
  legend_2->AddEntry(graph_a, "Exact", "l");
  legend_2->Draw();

  canvas->cd(3);
  TGraph *graph_s = to_TGraph(sol_s.t, sol_s.eigen_x(0));
  graph_s->SetTitle("Semi-Explicit vs Exact");
  graph_s->SetLineColor(colors[0]); graph_s->Draw("AL");
  graph_a->SetLineColor(colors[1]); graph_a->Draw("L SAME");
  graph_s->GetXaxis()->SetTitle("t");
  graph_s->GetYaxis()->SetTitle("x");
  graph_s->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  TLegend *legend_3 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend_3->AddEntry(graph_s, "Numerical", "l");
  legend_3->AddEntry(graph_a, "Exact", "l");
  legend_3->Draw();

  canvas->Update();
  app.Run();
  #endif

  return 0;
}
