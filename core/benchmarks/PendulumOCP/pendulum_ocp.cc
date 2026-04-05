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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR)                               \
  constexpr Integer N{PROBLEM##index0<Real>::equations_number()};       \
  PROBLEM##Problem<Real, PROBLEM##index3<Real>, INTEGRATOR<Real, N, 0>> \
      problem_index_3;                                                  \
  PROBLEM##Problem<Real, PROBLEM##index0<Real>, INTEGRATOR<Real, N, 2>> \
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
  problem_index_3.max_iterations(100);
  problem_index_0.max_iterations(100);

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
  Solution<Real, N, 0> sol_index_3(time.size());
  try {
    std::cout << "Solving problem with index 3..." << std::endl;
    problem_index_3.sigma(1.0);
    problem_index_3.multiple_shooting(time, guess);
    sol_index_3 = problem_index_3.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 3: " << e.what()
              << std::endl;
  }

  Solution<Real, N, 2> sol_index_0(time.size());
  try {
    std::cout << "Solving problem with index 0..." << std::endl;
    problem_index_0.sigma(1.0);
    problem_index_0.multiple_shooting(time, guess);
    sol_index_0 = problem_index_0.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 0: " << e.what()
              << std::endl;
  }

#ifdef SANDALS_ENABLE_PLOTTING

  auto colors = matlab_lines_colormap();

  TCanvas *canvas    = new TCanvas("canvas", "Solution Comparison", 1200, 1200);
  TGraph *graph_x_3  = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(0));
  TGraph *graph_y_3  = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(1));
  TGraph *graph_u_3  = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(2));
  TGraph *graph_v_3  = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(3));
  TGraph *graph_l_3  = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(4));
  TGraph *graph_l1_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(5));
  TGraph *graph_l2_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(6));
  TGraph *graph_l3_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(7));
  TGraph *graph_l4_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(8));
  TGraph *graph_l5_3 = to_TGraph(sol_index_3.t, sol_index_3.eigen_x(9));

  TGraph *graph_x_0  = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(0));
  TGraph *graph_y_0  = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(1));
  TGraph *graph_u_0  = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(2));
  TGraph *graph_v_0  = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(3));
  TGraph *graph_l_0  = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(4));
  TGraph *graph_l1_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(5));
  TGraph *graph_l2_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(6));
  TGraph *graph_l3_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(7));
  TGraph *graph_l4_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(8));
  TGraph *graph_l5_0 = to_TGraph(sol_index_0.t, sol_index_0.eigen_x(9));

  graph_x_3->SetLineColor(colors[0]);
  graph_y_3->SetLineColor(colors[0]);
  graph_u_3->SetLineColor(colors[0]);
  graph_v_3->SetLineColor(colors[0]);
  graph_l_3->SetLineColor(colors[0]);
  graph_l1_3->SetLineColor(colors[1]);
  graph_l2_3->SetLineColor(colors[2]);
  graph_l3_3->SetLineColor(colors[3]);
  graph_l4_3->SetLineColor(colors[4]);
  graph_l5_3->SetLineColor(colors[5]);
  graph_x_3->SetLineWidth(1.0);
  graph_y_3->SetLineWidth(1.0);
  graph_u_3->SetLineWidth(1.0);
  graph_v_3->SetLineWidth(1.0);
  graph_l_3->SetLineWidth(1.0);
  graph_l1_3->SetLineWidth(1.0);
  graph_l2_3->SetLineWidth(1.0);
  graph_l3_3->SetLineWidth(1.0);
  graph_l4_3->SetLineWidth(1.0);
  graph_l5_3->SetLineWidth(1.0);

  graph_x_0->SetLineColor(colors[1]);
  graph_y_0->SetLineColor(colors[1]);
  graph_u_0->SetLineColor(colors[1]);
  graph_v_0->SetLineColor(colors[1]);
  graph_l_0->SetLineColor(colors[1]);
  graph_l1_0->SetLineColor(colors[1]);
  graph_l2_0->SetLineColor(colors[2]);
  graph_l3_0->SetLineColor(colors[3]);
  graph_l4_0->SetLineColor(colors[4]);
  graph_l5_0->SetLineColor(colors[5]);
  graph_x_0->SetLineWidth(1.0);
  graph_y_0->SetLineWidth(1.0);
  graph_u_0->SetLineWidth(1.0);
  graph_v_0->SetLineWidth(1.0);
  graph_l_0->SetLineWidth(1.0);
  graph_l1_0->SetLineWidth(1.0);
  graph_l2_0->SetLineWidth(1.0);
  graph_l3_0->SetLineWidth(1.0);
  graph_l4_0->SetLineWidth(1.0);
  graph_l5_0->SetLineWidth(1.0);
  graph_x_0->SetLineStyle(2);
  graph_y_0->SetLineStyle(2);
  graph_u_0->SetLineStyle(2);
  graph_v_0->SetLineStyle(2);
  graph_l_0->SetLineStyle(2);
  graph_l1_0->SetLineStyle(2);
  graph_l2_0->SetLineStyle(2);
  graph_l3_0->SetLineStyle(2);
  graph_l4_0->SetLineStyle(2);
  graph_l5_0->SetLineStyle(2);

  canvas->Divide(2, 2);

  canvas->cd(1);
  gPad->SetGrid();
  auto mg1 = new TMultiGraph();
  mg1->SetTitle("Cartesian coordinates");
  mg1->Add(graph_x_3);
  mg1->Add(graph_y_3);
  mg1->Add(graph_x_0);
  mg1->Add(graph_y_0);
  mg1->Draw("AL");
  mg1->GetXaxis()->SetTitle("t");
  mg1->GetYaxis()->SetTitle("x, y (m)");
  mg1->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->AddEntry(graph_x_3, "x (index 3)", "l");
  legend1->AddEntry(graph_y_3, "y (index 3)", "l");
  legend1->AddEntry(graph_x_0, "x (index 0)", "l");
  legend1->AddEntry(graph_y_0, "y (index 0)", "l");
  legend1->Draw();

  canvas->cd(2);
  gPad->SetGrid();
  auto mg2 = new TMultiGraph();
  mg2->SetTitle("Velocities");
  mg2->Add(graph_u_3);
  mg2->Add(graph_v_3);
  mg2->Add(graph_u_0);
  mg2->Add(graph_v_0);
  mg2->Draw("AL");
  mg2->GetXaxis()->SetTitle("t");
  mg2->GetYaxis()->SetTitle("u, v (m/s)");
  mg2->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2->AddEntry(graph_u_3, "u (index 3)", "l");
  legend2->AddEntry(graph_v_3, "v (index 3)", "l");
  legend2->AddEntry(graph_u_0, "u (index 0)", "l");
  legend2->AddEntry(graph_v_0, "v (index 0)", "l");
  legend2->Draw();

  canvas->cd(3);
  gPad->SetGrid();
  auto mg3 = new TMultiGraph();
  mg3->SetTitle("Lagrange multiplier (DAE)");
  mg3->Add(graph_l_3);
  mg3->Add(graph_l_0);
  mg3->Draw("AL");
  mg3->GetXaxis()->SetTitle("t");
  mg3->GetYaxis()->SetTitle("lambda (N)");
  mg3->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend3 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend3->AddEntry(graph_l_3, "l (index 3)", "l");
  legend3->AddEntry(graph_l_0, "l (index 0)", "l");
  legend3->Draw();

  canvas->cd(4);
  gPad->SetGrid();
  auto mg4 = new TMultiGraph();
  mg4->SetTitle("Lagrange multipliers (OCP)");
  mg4->Add(graph_l1_3);
  mg4->Add(graph_l2_3);
  mg4->Add(graph_l3_3);
  mg4->Add(graph_l4_3);
  mg4->Add(graph_l5_3);
  mg4->Add(graph_l1_0);
  mg4->Add(graph_l2_0);
  mg4->Add(graph_l3_0);
  mg4->Add(graph_l4_0);
  mg4->Add(graph_l5_0);
  mg4->Draw("AL");
  mg4->GetXaxis()->SetTitle("t");
  mg4->GetYaxis()->SetTitle("l1, l2, l3, l4, l5 (-)");
  mg4->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend4 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend4->AddEntry(graph_l1_3, "l1 (index 3)", "l");
  legend4->AddEntry(graph_l2_3, "l2 (index 3)", "l");
  legend4->AddEntry(graph_l3_3, "l3 (index 3)", "l");
  legend4->AddEntry(graph_l4_3, "l4 (index 3)", "l");
  legend4->AddEntry(graph_l5_3, "l5 (index 3)", "l");
  legend4->AddEntry(graph_l1_0, "l1 (index 0)", "l");
  legend4->AddEntry(graph_l2_0, "l2 (index 0)", "l");
  legend4->AddEntry(graph_l3_0, "l3 (index 0)", "l");
  legend4->AddEntry(graph_l4_0, "l4 (index 0)", "l");
  legend4->AddEntry(graph_l5_0, "l5 (index 0)", "l");
  legend4->Draw();

  canvas->Update();
  app.Run();
#endif

  return 0;
}
