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
#include "Sandals/Problem.hh"
#include "Sandals/RungeKutta/RK4.hh" // ERK
#include "Sandals/RungeKutta/LobattoIIIA2.hh" // DIRK
#include "Sandals/RungeKutta/RadauIIA5.hh" // IRK

#include "BasicExplicit.hh"
#include "BasicImplicit.hh"
#include "BasicSemiExplicit.hh"

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

int main(int argc, char **argv) {

#ifdef SANDALS_ENABLE_PLOTTING
  TApplication app("plot", &argc, argv);
#else
  // Prevent compiler warnings
  (void)argc;
  (void)argv;
#endif

  static constexpr long num_points{10};

  Eigen::VectorXd time(Eigen::VectorXd::LinSpaced(num_points, 0.0, 1.0));

  std::cout << "\n=== Basic Explicit Problem ===\n" << std::endl;

  BasicExplicitProblem problem_explicit(std::make_shared<LobattoIIIA2<Real, 2>>());
  problem_explicit.verbose_mode(true);
  problem_explicit.integrator()->verbose_mode(false);
  {
    std::vector<Eigen::Vector<Real, 2>> guess(num_points, Eigen::Vector<Real, 2>::Zero());
    problem_explicit.multiple_shooting(time, guess);
  }

  std::cout << "\n=== Basic Implicit Problem ===\n" << std::endl;

  BasicImplicitProblem problem_implicit(std::make_shared<LobattoIIIA2<Real, 2>>());
  problem_implicit.verbose_mode(true);
  problem_implicit.integrator()->verbose_mode(false);
  {
    std::vector<Eigen::Vector<Real, 2>> guess(num_points, Eigen::Vector<Real, 2>::Zero());
    problem_implicit.multiple_shooting(time, guess);
  }

  std::cout << "\n=== Basic Semi-Explicit Problem ===\n" << std::endl;

  BasicSemiExplicitProblem problem_semiexplicit(std::make_shared<LobattoIIIA2<Real, 2>>());
  problem_semiexplicit.verbose_mode(true);
  problem_semiexplicit.integrator()->verbose_mode(false);
  {
    std::vector<Eigen::Vector<Real, 2>> guess(num_points, Eigen::Vector<Real, 2>::Zero());
    problem_semiexplicit.multiple_shooting(time, guess);
  }

  #ifdef SANDALS_ENABLE_PLOTTING
    auto esol = problem_explicit.solution();
    auto isol = problem_implicit.solution();
    auto ssol = problem_semiexplicit.solution();
    auto asol = problem_semiexplicit.analytical_solution(esol->t);

    auto colors = matlab_lines_colormap();

    TCanvas *canvas = new TCanvas("canvas", "Solution Comparison", 1200, 400);
    canvas->Divide(3, 1);

    canvas->cd(1);
    TGraph *graph_ex = to_TGraph(esol->t, esol->x.row(0));
    TGraph *graph_ey = to_TGraph(esol->t, esol->x.row(1));
    TGraph *graph_ax = to_TGraph(esol->t, asol.row(0).transpose().eval());
    TGraph *graph_ay = to_TGraph(esol->t, asol.row(1).transpose().eval());
    graph_ex->SetTitle("Explicit vs Analytical");
    graph_ex->SetLineColor(colors[0]); graph_ex->Draw("AL");
    graph_ey->SetLineColor(colors[1]); graph_ey->Draw("L SAME");
    graph_ax->SetMarkerColor(colors[0]); graph_ax->SetMarkerStyle(8); graph_ax->Draw("P SAME");
    graph_ay->SetMarkerColor(colors[1]); graph_ay->SetMarkerStyle(8); graph_ay->Draw("P SAME");
    graph_ex->GetXaxis()->SetTitle("t (s)");
    graph_ex->GetYaxis()->SetTitle("x, y (-)");
    graph_ex->GetXaxis()->SetLimits(time(0), time(Eigen::last));
    graph_ex->GetYaxis()->SetRangeUser(-1.25, 0.25);
    TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg1->AddEntry(graph_ex, "x (explicit)", "l");
    leg1->AddEntry(graph_ey, "y (explicit)", "l");
    leg1->AddEntry(graph_ax, "x (analytical)", "p");
    leg1->AddEntry(graph_ay, "y (analytical)", "p");
    leg1->Draw();

    canvas->cd(2);
    TGraph *graph_ix = to_TGraph(isol->t, isol->x.row(0));
    TGraph *graph_iy = to_TGraph(isol->t, isol->x.row(1));
    graph_ix->SetTitle("Implicit vs Analytical");
    graph_ix->SetLineColor(colors[0]); graph_ix->SetLineStyle(9); graph_ix->Draw("AL");
    graph_iy->SetLineColor(colors[1]); graph_iy->SetLineStyle(9); graph_iy->Draw("L SAME");
    graph_ax->SetMarkerColor(colors[0]); graph_ax->Draw("P SAME");
    graph_ay->SetMarkerColor(colors[1]); graph_ay->Draw("P SAME");
    graph_ix->GetXaxis()->SetTitle("t (s)");
    graph_ix->GetYaxis()->SetTitle("x, y (-)");
    graph_ix->GetXaxis()->SetLimits(time(0), time(Eigen::last));
    graph_ix->GetYaxis()->SetRangeUser(-1.25, 0.25);
    TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->AddEntry(graph_ix, "x (implicit)", "l");
    leg2->AddEntry(graph_iy, "y (implicit)", "l");
    leg2->AddEntry(graph_ax, "x (analytical)", "p");
    leg2->AddEntry(graph_ay, "y (analytical)", "p");
    leg2->Draw();

    canvas->cd(3);
    TGraph *graph_sx = to_TGraph(ssol->t, ssol->x.row(0));
    TGraph *graph_sy = to_TGraph(ssol->t, ssol->x.row(1));
    graph_sx->SetTitle("Semi-Explicit vs Analytical");
    graph_sx->SetLineColor(colors[0]); graph_sx->SetLineStyle(7); graph_sx->Draw("AL");
    graph_sy->SetLineColor(colors[1]); graph_sy->SetLineStyle(7); graph_sy->Draw("L SAME");
    graph_ax->SetMarkerColor(colors[0]); graph_ax->Draw("P SAME");
    graph_ay->SetMarkerColor(colors[1]); graph_ay->Draw("P SAME");
    graph_sx->GetXaxis()->SetTitle("t (s)");
    graph_sx->GetYaxis()->SetTitle("x, y (-)");
    graph_sx->GetXaxis()->SetLimits(time(0), time(Eigen::last));
    graph_sx->GetYaxis()->SetRangeUser(-1.25, 0.25);
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

