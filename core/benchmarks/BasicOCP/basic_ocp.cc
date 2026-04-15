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

#include "BasicOCP.hh"
#include "Sandals.hh"
#include "Sandals/RungeKutta/GaussLegendre6.hh"

#ifdef SANDALS_ENABLE_PLOTTING
#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMultiGraph.h>
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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR)                                 \
  PROBLEM##Problem<Real, PROBLEM##Implicit<Real>, INTEGRATOR<Real, 4, 0>> \
      problem;
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
  PROBLEM_INIT(BasicOCP, GaussLegendre6)

  // Set verbose mode
  problem.verbose_mode(true);
  problem.integrator()->verbose_mode(false);

  // Set reverse mode
  problem.integrator()->reverse_mode(false);

  // Set solver tolerance
  problem.tolerance(1.0e-14);

  // Set solver maximum number of iterations
  problem.max_iterations(100);

  // Set solution parameters
  problem.subintervals(1);

  // Set time mesh
  constexpr Integer num_points{100};
  Eigen::Vector<Real, Eigen::Dynamic> time(
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(num_points,
                                                     problem.time_start(),
                                                     problem.time_end()));

  // Set initial guess
  Eigen::Matrix<Real, 4, Eigen::Dynamic> guess(problem.guess(time));

  // Solve the problems with shooting
  Solution<Real, 4, 0> solution(time.size());
  try {
    std::cout << "Solving problem..." << std::endl;
    problem.sigma(1.0);
    problem.integrator()->projection_mode(false);
    problem.multiple_shooting(time, guess);
    solution = problem.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem: " << e.what() << std::endl;
  }

#ifdef SANDALS_ENABLE_PLOTTING

  auto colors = matlab_lines_colormap();

  TCanvas *canvas  = new TCanvas("canvas", "Solution Comparison", 1200, 600);
  TGraph *graph_x1 = to_TGraph(solution.t, solution.eigen_x(0));
  TGraph *graph_x2 = to_TGraph(solution.t, solution.eigen_x(1));
  TGraph *graph_l1 = to_TGraph(solution.t, solution.eigen_x(2));
  TGraph *graph_l2 = to_TGraph(solution.t, solution.eigen_x(3));

  graph_x1->SetLineColor(colors[0]);
  graph_x2->SetLineColor(colors[1]);
  graph_l1->SetLineColor(colors[2]);
  graph_l2->SetLineColor(colors[3]);

  canvas->Divide(2, 1);

  canvas->cd(1);
  gPad->SetGrid();
  auto mg1 = new TMultiGraph();
  mg1->SetTitle("States trajectories");
  mg1->Add(graph_x1);
  mg1->Add(graph_x2);
  mg1->Draw("AL");
  mg1->GetXaxis()->SetTitle("t");
  mg1->GetYaxis()->SetTitle("x1, x2");
  auto legend1 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend1->AddEntry(graph_x1, "x1", "l");
  legend1->AddEntry(graph_x2, "x2", "l");
  legend1->Draw();

  canvas->cd(2);
  gPad->SetGrid();
  auto mg2 = new TMultiGraph();
  mg2->SetTitle("Lagrange multipliers");
  mg2->Add(graph_l1);
  mg2->Add(graph_l2);
  mg2->Draw("AL");
  mg2->GetXaxis()->SetTitle("t");
  mg2->GetYaxis()->SetTitle("l1, l2");
  auto legend2 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend2->AddEntry(graph_l1, "l1", "l");
  legend2->AddEntry(graph_l2, "l2", "l");
  legend2->Draw();

  canvas->Update();
  app.Run();
#endif

  return 0;
}
