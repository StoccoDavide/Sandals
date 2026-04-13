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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR_3, INTEGRATOR_0)                  \
  PROBLEM##Problem<Real, PROBLEM##index3<Real>, INTEGRATOR_3<Real, 10, 0>> \
      problem_index_3;                                                     \
  PROBLEM##Problem<Real, PROBLEM##index0<Real>, INTEGRATOR_0<Real, 10, 3>> \
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
  PROBLEM_INIT(PendulumOCP, RadauIIA5, RK4)

  // Set verbose mode
  constexpr bool verbose{true};
  problem_index_3.verbose_mode(verbose);
  problem_index_0.verbose_mode(verbose);
  problem_index_3.integrator()->verbose_mode(false);
  problem_index_0.integrator()->verbose_mode(false);

  // Set reverse mode
  constexpr bool reverse{false};
  problem_index_3.integrator()->reverse_mode(reverse);
  problem_index_0.integrator()->reverse_mode(reverse);

  // Set solver tolerance
  problem_index_3.tolerance(1.0e-9);
  problem_index_0.tolerance(1.0e-9);

  // Set solver maximum number of iterations
  problem_index_3.max_iterations(100);
  problem_index_0.max_iterations(100);

  // Set solution parameters
  constexpr Integer num_subintervals{1};
  problem_index_3.subintervals(num_subintervals);
  problem_index_0.subintervals(num_subintervals);

  // Set time mesh
  constexpr Integer num_points{50};
  Eigen::Vector<Real, Eigen::Dynamic> time(
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(
          num_points,
          problem_index_3.time_start(),
          problem_index_3.time_end()));

  // Set initial guess
  Eigen::Matrix<Real, 10, Eigen::Dynamic> guess(problem_index_3.guess(time));

  // Solve the problems with shooting
  Solution<Real, 10, 0> sol_index_3(time.size());
  try {
    std::cout << "Solving problem with index 3..." << std::endl;
    problem_index_3.sigma(1.0);
    problem_index_0.lambda(1.0e-3);
    problem_index_3.integrator()->projection_mode(true);
    problem_index_3.multiple_shooting(time, guess);
    sol_index_3 = problem_index_3.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 3: " << e.what()
              << std::endl;
  }

  Solution<Real, 10, 3> sol_index_0_lesq(time.size());
  try {
    std::cout << "Solving problem with index 0 (least squares)..." << std::endl;
    sol_index_0_lesq.x = guess;
    std::vector<Real> lambda_vec{1.0e-3, 1.0e-6, 1.0e-6, 1.0e-9, 0.0};
    std::vector<Real> sigma_vec{1.0, 10.0, 100.0, 10.0, 10.0};
    for (Integer i{0}; i < 4; ++i) {
      problem_index_0.sigma(sigma_vec[i]);
      problem_index_0.lambda(lambda_vec[i]);
      problem_index_0.integrator()->projection_mode(false);
      problem_index_0.multiple_shooting(time, sol_index_0_lesq.x);
      sol_index_0_lesq = problem_index_0.solution();
    }
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 0 (least squares): "
              << e.what() << std::endl;
  }

  Solution<Real, 10, 3> sol_index_0_proj(time.size());
  try {
    std::cout << "Solving problem with index 0 (projection)..." << std::endl;
    problem_index_0.sigma(1.0);
    problem_index_0.lambda(1.0e-3);
    problem_index_0.integrator()->projection_mode(true);
    problem_index_0.multiple_shooting(time, guess);
    sol_index_0_proj = problem_index_0.solution();
  } catch (const std::exception &e) {
    std::cerr << "Error solving problem with index 0 (projection): " << e.what()
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

  TGraph *graph_x_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(0));
  TGraph *graph_y_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(1));
  TGraph *graph_u_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(2));
  TGraph *graph_v_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(3));
  TGraph *graph_l_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(4));
  TGraph *graph_l1_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(5));
  TGraph *graph_l2_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(6));
  TGraph *graph_l3_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(7));
  TGraph *graph_l4_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(8));
  TGraph *graph_l5_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_x(9));
  TGraph *graph_h1_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_h(0));
  TGraph *graph_h2_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_h(1));
  TGraph *graph_h3_0_lesq =
      to_TGraph(sol_index_0_lesq.t, sol_index_0_lesq.eigen_h(2));

  TGraph *graph_x_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(0));
  TGraph *graph_y_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(1));
  TGraph *graph_u_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(2));
  TGraph *graph_v_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(3));
  TGraph *graph_l_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(4));
  TGraph *graph_l1_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(5));
  TGraph *graph_l2_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(6));
  TGraph *graph_l3_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(7));
  TGraph *graph_l4_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(8));
  TGraph *graph_l5_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_x(9));
  TGraph *graph_h1_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_h(0));
  TGraph *graph_h2_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_h(1));
  TGraph *graph_h3_0_proj =
      to_TGraph(sol_index_0_proj.t, sol_index_0_proj.eigen_h(2));

  graph_x_3->SetLineColor(colors[0]);
  graph_y_3->SetLineColor(colors[0]);
  graph_u_3->SetLineColor(colors[0]);
  graph_v_3->SetLineColor(colors[0]);
  graph_l_3->SetLineColor(colors[0]);
  graph_l1_3->SetLineColor(colors[0]);
  graph_l2_3->SetLineColor(colors[0]);
  graph_l3_3->SetLineColor(colors[0]);
  graph_l4_3->SetLineColor(colors[0]);
  graph_l5_3->SetLineColor(colors[0]);
  graph_x_3->SetLineWidth(2.0);
  graph_y_3->SetLineWidth(2.0);
  graph_u_3->SetLineWidth(2.0);
  graph_v_3->SetLineWidth(2.0);
  graph_l_3->SetLineWidth(2.0);
  graph_l1_3->SetLineWidth(2.0);
  graph_l2_3->SetLineWidth(2.0);
  graph_l3_3->SetLineWidth(2.0);
  graph_l4_3->SetLineWidth(2.0);
  graph_l5_3->SetLineWidth(2.0);

  graph_x_0_lesq->SetLineColor(colors[1]);
  graph_y_0_lesq->SetLineColor(colors[1]);
  graph_u_0_lesq->SetLineColor(colors[1]);
  graph_v_0_lesq->SetLineColor(colors[1]);
  graph_l_0_lesq->SetLineColor(colors[1]);
  graph_l1_0_lesq->SetLineColor(colors[1]);
  graph_l2_0_lesq->SetLineColor(colors[1]);
  graph_l3_0_lesq->SetLineColor(colors[1]);
  graph_l4_0_lesq->SetLineColor(colors[1]);
  graph_l5_0_lesq->SetLineColor(colors[1]);
  graph_x_0_lesq->SetLineWidth(2.0);
  graph_y_0_lesq->SetLineWidth(2.0);
  graph_u_0_lesq->SetLineWidth(2.0);
  graph_v_0_lesq->SetLineWidth(2.0);
  graph_l_0_lesq->SetLineWidth(2.0);
  graph_l1_0_lesq->SetLineWidth(2.0);
  graph_l2_0_lesq->SetLineWidth(2.0);
  graph_l3_0_lesq->SetLineWidth(2.0);
  graph_l4_0_lesq->SetLineWidth(2.0);
  graph_l5_0_lesq->SetLineWidth(2.0);

  graph_x_0_lesq->SetLineColor(colors[2]);
  graph_y_0_lesq->SetLineColor(colors[2]);
  graph_u_0_lesq->SetLineColor(colors[2]);
  graph_v_0_lesq->SetLineColor(colors[2]);
  graph_l_0_lesq->SetLineColor(colors[2]);
  graph_l1_0_lesq->SetLineColor(colors[2]);
  graph_l2_0_lesq->SetLineColor(colors[2]);
  graph_l3_0_lesq->SetLineColor(colors[2]);
  graph_l4_0_lesq->SetLineColor(colors[2]);
  graph_l5_0_lesq->SetLineColor(colors[2]);
  graph_x_0_lesq->SetLineWidth(2.0);
  graph_y_0_lesq->SetLineWidth(2.0);
  graph_u_0_lesq->SetLineWidth(2.0);
  graph_v_0_lesq->SetLineWidth(2.0);
  graph_l_0_lesq->SetLineWidth(2.0);
  graph_l1_0_lesq->SetLineWidth(2.0);
  graph_l2_0_lesq->SetLineWidth(2.0);
  graph_l3_0_lesq->SetLineWidth(2.0);
  graph_l4_0_lesq->SetLineWidth(2.0);
  graph_l5_0_lesq->SetLineWidth(2.0);

  canvas->Divide(2, 3);

  canvas->cd(1);
  gPad->SetGrid();
  auto mg1 = new TMultiGraph();
  mg1->SetTitle("Cartesian coordinates");
  mg1->Add(graph_x_3);
  mg1->Add(graph_y_3);
  mg1->Add(graph_x_0_lesq);
  mg1->Add(graph_y_0_lesq);
  mg1->Add(graph_x_0_proj);
  mg1->Add(graph_y_0_proj);
  mg1->Draw("AL");
  mg1->GetXaxis()->SetTitle("t");
  mg1->GetYaxis()->SetTitle("x, y (m)");
  mg1->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->AddEntry(graph_x_3, "x (index 3)", "l");
  legend1->AddEntry(graph_y_3, "y (index 3)", "l");
  legend1->AddEntry(graph_x_0_lesq, "x (l.s. index 0)", "l");
  legend1->AddEntry(graph_y_0_lesq, "y (l.s. index 0)", "l");
  legend1->AddEntry(graph_x_0_proj, "x (proj. index 0)", "l");
  legend1->AddEntry(graph_y_0_proj, "y (proj. index 0)", "l");
  legend1->Draw();

  canvas->cd(2);
  gPad->SetGrid();
  auto mg2 = new TMultiGraph();
  mg2->SetTitle("Velocities");
  mg2->Add(graph_u_3);
  mg2->Add(graph_v_3);
  mg2->Add(graph_u_0_lesq);
  mg2->Add(graph_v_0_lesq);
  mg2->Add(graph_u_0_proj);
  mg2->Add(graph_v_0_proj);
  mg2->Draw("AL");
  mg2->GetXaxis()->SetTitle("t");
  mg2->GetYaxis()->SetTitle("u, v (m/s)");
  mg2->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2->AddEntry(graph_u_3, "u (index 3)", "l");
  legend2->AddEntry(graph_v_3, "v (index 3)", "l");
  legend2->AddEntry(graph_u_0_lesq, "u (l.s. index 0)", "l");
  legend2->AddEntry(graph_v_0_lesq, "v (l.s. index 0)", "l");
  legend2->AddEntry(graph_u_0_proj, "u (proj. index 0)", "l");
  legend2->AddEntry(graph_v_0_proj, "v (proj. index 0)", "l");
  legend2->Draw();

  canvas->cd(3);
  gPad->SetGrid();
  auto mg3 = new TMultiGraph();
  mg3->SetTitle("Lagrange multiplier (DAE)");
  mg3->Add(graph_l_3);
  mg3->Add(graph_l_0_lesq);
  mg3->Add(graph_l_0_proj);
  mg3->Draw("AL");
  mg3->GetXaxis()->SetTitle("t");
  mg3->GetYaxis()->SetTitle("lambda (N)");
  mg3->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend3 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend3->AddEntry(graph_l_3, "l (index 3)", "l");
  legend3->AddEntry(graph_l_0_lesq, "l (l.s. index 0)", "l");
  legend3->AddEntry(graph_l_0_proj, "l (proj. index 0)", "l");
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
  mg4->Add(graph_l1_0_lesq);
  mg4->Add(graph_l2_0_lesq);
  mg4->Add(graph_l3_0_lesq);
  mg4->Add(graph_l4_0_lesq);
  mg4->Add(graph_l5_0_lesq);
  mg4->Add(graph_l1_0_proj);
  mg4->Add(graph_l2_0_proj);
  mg4->Add(graph_l3_0_proj);
  mg4->Add(graph_l4_0_proj);
  mg4->Add(graph_l5_0_proj);
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
  legend4->AddEntry(graph_l1_0_lesq, "l1 (l.s. index 0)", "l");
  legend4->AddEntry(graph_l2_0_lesq, "l2 (l.s. index 0)", "l");
  legend4->AddEntry(graph_l3_0_lesq, "l3 (l.s. index 0)", "l");
  legend4->AddEntry(graph_l4_0_lesq, "l4 (l.s. index 0)", "l");
  legend4->AddEntry(graph_l5_0_lesq, "l5 (l.s. index 0)", "l");
  legend4->AddEntry(graph_l1_0_proj, "l1 (proj. index 0)", "l");
  legend4->AddEntry(graph_l2_0_proj, "l2 (proj. index 0)", "l");
  legend4->AddEntry(graph_l3_0_proj, "l3 (proj. index 0)", "l");
  legend4->AddEntry(graph_l4_0_proj, "l4 (proj. index 0)", "l");
  legend4->AddEntry(graph_l5_0_proj, "l5 (proj. index 0)", "l");
  legend4->Draw();

  canvas->cd(5);
  gPad->SetGrid();
  auto mg5 = new TMultiGraph();
  mg5->SetTitle("Invariant manifold (index 0)");
  mg5->Add(graph_h1_0_lesq);
  mg5->Add(graph_h2_0_lesq);
  mg5->Add(graph_h3_0_lesq);
  mg5->Add(graph_h1_0_proj);
  mg5->Add(graph_h2_0_proj);
  mg5->Add(graph_h3_0_proj);
  mg5->Draw("AL");
  mg5->GetXaxis()->SetTitle("t");
  mg5->GetYaxis()->SetTitle("h1, h2, h3 (-)");
  mg5->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff());
  auto legend5 = new TLegend(0.7, 0.8, 0.9, 0.9);
  legend5->AddEntry(graph_h1_0_lesq, "h1 (l.s. index 0)", "l");
  legend5->AddEntry(graph_h2_0_lesq, "h2 (l.s. index 0)", "l");
  legend5->AddEntry(graph_h3_0_lesq, "h3 (l.s. index 0)", "l");
  legend5->AddEntry(graph_h1_0_proj, "h1 (proj. index 0)", "l");
  legend5->AddEntry(graph_h2_0_proj, "h2 (proj. index 0)", "l");
  legend5->AddEntry(graph_h3_0_proj, "h3 (proj. index 0)", "l");
  legend5->Draw();

  canvas->Update();
  app.Run();
#endif

  return 0;
}
