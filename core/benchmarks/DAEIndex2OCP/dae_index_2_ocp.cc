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

#include "DAEIndex2OCP.hh"
#include "Sandals.hh"
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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR_2, INTEGRATOR_0)                 \
  PROBLEM##Problem<Real, PROBLEM##index2<Real>, INTEGRATOR_2<Real, 8, 0>> \
      problem_index_2;                                                    \
  PROBLEM##Problem<Real, PROBLEM##index0<Real>, INTEGRATOR_0<Real, 8, 2>> \
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
  PROBLEM_INIT(DAEIndex2OCP, RadauIIA5, RadauIIA5)

  // Set verbose mode
  constexpr bool verbose{true};
  problem_index_2.verbose_mode(verbose);
  problem_index_0.verbose_mode(verbose);
  problem_index_2.integrator()->verbose_mode(false);
  problem_index_0.integrator()->verbose_mode(false);

  // Set reverse mode
  constexpr bool reverse{false};
  problem_index_2.integrator()->reverse_mode(reverse);
  problem_index_0.integrator()->reverse_mode(reverse);

  // Set solver tolerance
  problem_index_2.tolerance(1.0e-12);
  problem_index_0.tolerance(1.0e-12);

  // Set solver maximum number of iterations
  problem_index_2.max_iterations(200);
  problem_index_0.max_iterations(200);

  // Set solution parameters
  constexpr Integer num_subintervals{1};
  problem_index_2.subintervals(num_subintervals);
  problem_index_0.subintervals(num_subintervals);

  // Set time mesh
  constexpr Integer num_points{100};
  Eigen::Vector<Real, Eigen::Dynamic> time(
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
          num_points,
          problem_index_2.time_start(),
          problem_index_2.time_end()));

  // Set initial guess
  Eigen::Matrix<Real, 8, Eigen::Dynamic> guess(problem_index_2.guess(time));

  // Solve the problems with shooting
  Solution<Real, 8, 0> sol_index_2(time.size());
  try {
    std::cout << "Solving problem with index 3..." << std::endl;
    problem_index_2.sigma(1.0);
    problem_index_0.lambda(0.0);
    problem_index_2.integrator()->projection_mode(true);
    problem_index_2.multiple_shooting(time, guess);
    sol_index_2 = problem_index_2.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 3: " << e.what()
              << std::endl;
  }

  Solution<Real, 8, 2> sol_index_0_lesq(time.size());
  try {
    std::cout << "Solving problem with index 0 (least squares)..." << std::endl;
    problem_index_0.sigma(1.0);
    problem_index_0.lambda(0.0);
    problem_index_0.integrator()->projection_mode(false);
    problem_index_0.multiple_shooting(time, guess);
    sol_index_0_lesq = problem_index_0.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 0 (least squares): "
              << e.what() << std::endl;
  }

  Solution<Real, 8, 2> sol_index_0_proj(time.size());
  try {
    std::cout << "Solving problem with index 0 (projection)..." << std::endl;
    problem_index_0.sigma(1.0);
    problem_index_0.lambda(0.0);
    problem_index_0.integrator()->projection_mode(true);
    problem_index_0.multiple_shooting(time, guess);
    sol_index_0_proj = problem_index_0.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 0 (projection): " << e.what()
              << std::endl;
  }

#ifdef SANDALS_ENABLE_PLOTTING

  auto colors = matlab_lines_colormap();

  // Compute the invariants violations for the index 2 solutions

  Solution<Real, 8, 2> sol_index_2_aug(time.size());
  sol_index_2_aug.t = sol_index_2.t;
  sol_index_2_aug.x = sol_index_2.x;
  for (Integer i{0}; i < sol_index_2.t.size(); ++i) {
    sol_index_2_aug.h.col(i) =
        problem_index_0.system()->h(sol_index_2.x.col(i), sol_index_2.t[i]);
  }

  TCanvas *canvas    = new TCanvas("canvas", "Solution Comparison", 1200, 1200);
  TGraph *graph_x1_2 = to_TGraph(sol_index_2.t, sol_index_2.eigen_x(0));
  TGraph *graph_x2_2 = to_TGraph(sol_index_2.t, sol_index_2.eigen_x(1));
  TGraph *graph_x3_2 = to_TGraph(sol_index_2.t, sol_index_2.eigen_x(2));
  TGraph *graph_x4_2 = to_TGraph(sol_index_2.t, sol_index_2.eigen_x(3));
  TGraph *graph_h1_2 = to_TGraph(sol_index_2_aug.t, sol_index_2_aug.eigen_h(0));
  TGraph *graph_h2_2 = to_TGraph(sol_index_2_aug.t, sol_index_2_aug.eigen_h(1));

  TGraph *graph_x1_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(0));
  TGraph *graph_x2_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(1));
  TGraph *graph_x3_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(2));
  TGraph *graph_x4_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(3));
  TGraph *graph_h1_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_h(0));
  TGraph *graph_h2_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_h(1));

  TGraph *graph_x1_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(0));
  TGraph *graph_x2_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(1));
  TGraph *graph_x3_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(2));
  TGraph *graph_x4_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(3));
  TGraph *graph_h1_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_h(0));
  TGraph *graph_h2_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_h(1));

  graph_x1_2->SetLineColor(colors[0]);
  graph_x2_2->SetLineColor(colors[0]);
  graph_x3_2->SetLineColor(colors[0]);
  graph_x4_2->SetLineColor(colors[0]);
  graph_x1_2->SetLineWidth(2.0);
  graph_x2_2->SetLineWidth(2.0);
  graph_x3_2->SetLineWidth(2.0);
  graph_x4_2->SetLineWidth(2.0);

  graph_x1_0_lesq->SetLineColor(colors[1]);
  graph_x2_0_lesq->SetLineColor(colors[1]);
  graph_x3_0_lesq->SetLineColor(colors[1]);
  graph_x4_0_lesq->SetLineColor(colors[1]);
  graph_h1_0_lesq->SetLineColor(colors[1]);
  graph_h2_0_lesq->SetLineColor(colors[1]);
  graph_x1_0_lesq->SetLineWidth(2.0);
  graph_x2_0_lesq->SetLineWidth(2.0);
  graph_x3_0_lesq->SetLineWidth(2.0);
  graph_x4_0_lesq->SetLineWidth(2.0);
  graph_h1_0_lesq->SetLineWidth(2.0);
  graph_h2_0_lesq->SetLineWidth(2.0);

  graph_x1_0_proj->SetLineColor(colors[2]);
  graph_x2_0_proj->SetLineColor(colors[2]);
  graph_x3_0_proj->SetLineColor(colors[2]);
  graph_x4_0_proj->SetLineColor(colors[2]);
  graph_h1_0_proj->SetLineColor(colors[2]);
  graph_h2_0_proj->SetLineColor(colors[2]);
  graph_x1_0_proj->SetLineWidth(2.0);
  graph_x2_0_proj->SetLineWidth(2.0);
  graph_x3_0_proj->SetLineWidth(2.0);
  graph_x4_0_proj->SetLineWidth(2.0);
  graph_h1_0_proj->SetLineWidth(2.0);
  graph_h2_0_proj->SetLineWidth(2.0);
  graph_x1_0_lesq->SetLineWidth(2.0);
  graph_x2_0_lesq->SetLineWidth(2.0);
  graph_x3_0_lesq->SetLineWidth(2.0);
  graph_x4_0_lesq->SetLineWidth(2.0);
  graph_h1_0_lesq->SetLineWidth(2.0);
  graph_h2_0_lesq->SetLineWidth(2.0);

  canvas->Divide(2, 3);

  canvas->cd(1);
  gPad->SetGrid();
  auto mg1 = new TMultiGraph();
  mg1->SetTitle("State 1");
  mg1->Add(graph_x1_2);
  mg1->Add(graph_x1_0_lesq);
  mg1->Add(graph_x1_0_proj);
  mg1->Draw("AL");
  mg1->GetXaxis()->SetTitle("t");
  mg1->GetYaxis()->SetTitle("x1");
  mg1->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->AddEntry(graph_x1_2, "x1 (index 2)", "l");
  legend1->AddEntry(graph_x1_0_lesq, "x1 (index 0, least squares)", "l");
  legend1->AddEntry(graph_x1_0_proj, "x1 (index 0, projection)", "l");
  legend1->Draw();

  canvas->cd(2);
  gPad->SetGrid();
  auto mg2 = new TMultiGraph();
  mg2->SetTitle("State 2");
  mg2->Add(graph_x2_2);
  mg2->Add(graph_x2_0_lesq);
  mg2->Add(graph_x2_0_proj);
  mg2->Draw("AL");
  mg2->GetXaxis()->SetTitle("t");
  mg2->GetYaxis()->SetTitle("x2");
  mg2->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2->AddEntry(graph_x2_2, "x2 (index 2)", "l");
  legend2->AddEntry(graph_x2_0_lesq, "x2 (index 0, least squares)", "l");
  legend2->AddEntry(graph_x2_0_proj, "x2 (index 0, projection)", "l");
  legend2->Draw();

  canvas->cd(3);
  gPad->SetGrid();
  auto mg3 = new TMultiGraph();
  mg3->SetTitle("State 3");
  mg3->Add(graph_x3_2);
  mg3->Add(graph_x3_0_lesq);
  mg3->Add(graph_x3_0_proj);
  mg3->Draw("AL");
  mg3->GetXaxis()->SetTitle("t");
  mg3->GetYaxis()->SetTitle("x3");
  mg3->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend3 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend3->AddEntry(graph_x3_2, "x3 (index 2)", "l");
  legend3->AddEntry(graph_x3_0_lesq, "x3 (index 0, least squares)", "l");
  legend3->AddEntry(graph_x3_0_proj, "x3 (index 0, projection)", "l");
  legend3->Draw();

  canvas->cd(4);
  gPad->SetGrid();
  auto mg4 = new TMultiGraph();
  mg4->SetTitle("State 4");
  mg4->Add(graph_x4_2);
  mg4->Add(graph_x4_0_lesq);
  mg4->Add(graph_x4_0_proj);
  mg4->Draw("AL");
  mg4->GetXaxis()->SetTitle("t");
  mg4->GetYaxis()->SetTitle("x4");
  mg4->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend4 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend4->AddEntry(graph_x4_2, "x4 (index 2)", "l");
  legend4->AddEntry(graph_x4_0_lesq, "x4 (index 0, least squares)", "l");
  legend4->AddEntry(graph_x4_0_proj, "x4 (index 0, projection)", "l");
  legend4->Draw();

  canvas->cd(5);
  gPad->SetGrid();
  auto mg5 = new TMultiGraph();
  mg5->SetTitle("Invariant manifold (index 0)");
  mg5->Add(graph_h1_2);
  mg5->Add(graph_h2_2);
  mg5->Add(graph_h1_0_lesq);
  mg5->Add(graph_h2_0_lesq);
  mg5->Add(graph_h1_0_proj);
  mg5->Add(graph_h2_0_proj);
  mg5->Draw("AL");
  mg5->GetXaxis()->SetTitle("t");
  mg5->GetYaxis()->SetTitle("h1, h2, h3 (-)");
  mg5->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend5 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend5->AddEntry(graph_h1_2, "h1 (index 2)", "l");
  legend5->AddEntry(graph_h2_2, "h2 (index 2)", "l");
  legend5->AddEntry(graph_h1_0_lesq, "h1 (l.s. index 0)", "l");
  legend5->AddEntry(graph_h2_0_lesq, "h2 (l.s. index 0)", "l");
  legend5->AddEntry(graph_h1_0_proj, "h1 (proj. index 0)", "l");
  legend5->AddEntry(graph_h2_0_proj, "h2 (proj. index 0)", "l");
  legend5->Draw();

  canvas->Update();
  app.Run();
#endif

  return 0;
}
