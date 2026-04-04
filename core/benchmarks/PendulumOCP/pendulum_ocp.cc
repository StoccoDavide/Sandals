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

#include "PendulumOCP.hh"
#include "Sandals.hh"
#include "Sandals/RungeKutta/RK4.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"

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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR)                               \
  constexpr Integer N{PROBLEM##index0<Real>::equations_number()};       \
  constexpr Integer M{PROBLEM##index0<Real>::invariants_number()};      \
  PROBLEM##Problem<Real, PROBLEM##index3<Real>, INTEGRATOR<Real, N, 0>> \
      problem_index_3;                                                  \
  PROBLEM##Problem<Real, PROBLEM##index0<Real>, INTEGRATOR<Real, N, M>> \
      problem_index_0;
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
  PROBLEM_INIT(PendulumOCP, RadauIIA5)

  // Set verbose mode
  constexpr bool verbose{true};
  problem_index_3.verbose_mode(verbose);
  problem_index_0.verbose_mode(verbose);
  problem_index_3.integrator()->verbose_mode(false);
  problem_index_0.integrator()->verbose_mode(false);

  // Set projection mode
  problem_index_3.integrator()->projection_mode(true);
  problem_index_0.integrator()->projection_mode(false);

  // Set reverse mode
  constexpr bool reverse{false};
  problem_index_3.integrator()->reverse_mode(reverse);
  problem_index_0.integrator()->reverse_mode(reverse);

  // Set solver tolerance
  problem_index_3.tolerance(1.0e-8);
  problem_index_0.tolerance(1.0e-8);

  // Set solver maximum number of iterations
  problem_index_3.max_iterations(30);
  problem_index_0.max_iterations(30);

  // Set solution parameters
  constexpr Integer num_subintervals{1};
  problem_index_3.subintervals(num_subintervals);
  problem_index_0.subintervals(num_subintervals);

  // Set time mesh
  constexpr Integer num_points{300};
  Eigen::Vector<Real, Eigen::Dynamic> time(
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
          num_points,
          problem_index_3.time_start(),
          problem_index_3.time_end()));

  // Set initial guess
  Eigen::Matrix<Real, N, Eigen::Dynamic> guess(problem_index_3.guess(time));

  // Solve the problems with shooting
  problem_index_3.sigma(1.0);
  problem_index_0.sigma(1.0);
  problem_index_3.multiple_shooting(time, guess);
  problem_index_0.multiple_shooting(time, guess);

#ifdef SANDALS_ENABLE_PLOTTING
  auto sol_index_3 = problem_index_3.solution();
  auto sol_index_0 = problem_index_0.solution();

  auto colors = matlab_lines_colormap();

  TCanvas *canvas = new TCanvas("canvas", "Solution Comparison", 1200, 400);
  canvas->Divide(2, 1);

  canvas->cd(1);
  TGraph *graph_x_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(0));
  TGraph *graph_x_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(0));
  graph_x_3->SetLineColor(colors[0]);
  graph_x_3->Draw("AL");
  graph_x_3->GetXaxis()->SetTitle("t");
  graph_x_3->GetYaxis()->SetTitle("x");
  graph_x_0->SetLineColor(colors[1]);
  graph_x_0->Draw("L SAME");
  graph_x_0->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());

  canvas->cd(2);
  TGraph *graph_f_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(4));
  TGraph *graph_f_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(4));
  graph_f_3->SetLineColor(colors[1]);
  graph_f_3->Draw("AL");
  graph_f_0->SetLineColor(colors[0]);
  graph_f_0->Draw("L SAME");
  graph_f_3->GetXaxis()->SetTitle("t");
  graph_f_3->GetYaxis()->SetTitle("f");
  graph_f_3->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  graph_f_0->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());

  canvas->Update();
  app.Run();
#endif

  return 0;
}
