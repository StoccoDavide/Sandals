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
#include "Sandals/RungeKutta/RK4.hh"
#include "Sandals/RungeKutta/Fehlberg45.hh"

#include "OscillatorImplicit.hh"
#include "OscillatorExplicit.hh"
#include "OscillatorSemiExplicit.hh"

#include "SimplePendulumImplicit.hh"
#include "SimplePendulumExplicit.hh"
#include "SimplePendulumSemiExplicit.hh"

#include "SinCosImplicit.hh"
#include "SinCosExplicit.hh"
#include "SinCosSemiExplicit.hh"

#include "SinImplicit.hh"
#include "SinExplicit.hh"
#include "SinSemiExplicit.hh"

#include "ThreeBodyImplicit.hh"
#include "ThreeBodyExplicit.hh"
#include "ThreeBodySemiExplicit.hh"

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

  Fehlberg45<Real, 2, 0>::Time time = Eigen::VectorXd::LinSpaced(1000, 0.0, 10.0);

  SinCosImplicit system_implicit;
  Fehlberg45<Real, 2, 0> integrator_implicit(std::make_shared<SinCosImplicit<Real>>(system_implicit));
  integrator_implicit.enable_projection();
  integrator_implicit.disable_reverse_mode();
  Solution<Real, 2, 0> solution_implicit;
  Solution<Real, 2, 0> solution_implicit_adaptive;
  integrator_implicit.solve(time, system_implicit.ics(), solution_implicit);
  integrator_implicit.adaptive_solve(time, system_implicit.ics(), solution_implicit_adaptive);

#ifdef SANDALS_ENABLE_PLOTTING
  TCanvas *c1a = new TCanvas("c1a", "Oscillator Implicit", 800, 600);
  TGraph *g1a_x0 = to_TGraph(solution_implicit.t, solution_implicit.x.row(0));
  TGraph *g1a_x1 = to_TGraph(solution_implicit.t, solution_implicit.x.row(1));
  TGraph *g1a_x0a = to_TGraph(solution_implicit_adaptive.t, solution_implicit_adaptive.x.row(0));
  TGraph *g1a_x1a = to_TGraph(solution_implicit_adaptive.t, solution_implicit_adaptive.x.row(1));

  auto colors = matlab_lines_colormap();
  g1a_x0->SetLineColor(colors[0]);  g1a_x0->Draw("AL");
  g1a_x1->SetLineColor(colors[0]);  g1a_x1->Draw("L SAME");
  g1a_x0a->SetLineColor(colors[1]); g1a_x0a->SetLineStyle(9); g1a_x0a->Draw("L SAME");
  g1a_x1a->SetLineColor(colors[1]); g1a_x1a->SetLineStyle(9); g1a_x1a->Draw("L SAME");
  g1a_x0->SetTitle("Sine Implicit");
  g1a_x0->GetXaxis()->SetTitle("t (s)");
  g1a_x0->GetYaxis()->SetTitle("x (-)");
  c1a->SetTitle("Sine Implicit");

  TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(g1a_x0,  "Implicit x0 (fixed)", "l"); // "l" = line, "p" = point
  leg->AddEntry(g1a_x1,  "Implicit x1 (fixed)", "l");
  leg->AddEntry(g1a_x0a, "Implicit x0 (adaptive)", "l");
  leg->AddEntry(g1a_x1a, "Implicit x1 (adaptive)", "l");
  leg->Draw();
  c1a->Update();
#endif

  SinCosExplicit system_explicit;
  Fehlberg45<Real, 2, 0> integrator_explicit(std::make_shared<SinCosExplicit<Real>>(system_explicit));
  integrator_explicit.enable_projection();
  integrator_explicit.disable_reverse_mode();
  Solution<Real, 2, 0> solution_explicit;
  Solution<Real, 2, 0> solution_explicit_adaptive;
  integrator_explicit.solve(time, system_explicit.ics(), solution_explicit);
  integrator_explicit.adaptive_solve(time, system_explicit.ics(), solution_explicit_adaptive);

#ifdef SANDALS_ENABLE_PLOTTING
  TCanvas *c1b = new TCanvas("c1b", "Sine Explicit", 800, 600);
  TGraph *g1b_x0 = to_TGraph(solution_explicit.t, solution_explicit.x.row(0));
  TGraph *g1b_x1 = to_TGraph(solution_explicit.t, solution_explicit.x.row(1));
  TGraph *g1b_x0a = to_TGraph(solution_explicit_adaptive.t, solution_explicit_adaptive.x.row(0));
  TGraph *g1b_x1a = to_TGraph(solution_explicit_adaptive.t, solution_explicit_adaptive.x.row(1));

  g1b_x0->SetLineColor(colors[0]);  g1b_x0->Draw("AL");
  g1b_x1->SetLineColor(colors[0]);  g1b_x1->Draw("L SAME");
  g1b_x0a->SetLineColor(colors[1]); g1b_x0a->SetLineStyle(9); g1b_x0a->Draw("L SAME");
  g1b_x1a->SetLineColor(colors[1]); g1b_x1a->SetLineStyle(9); g1b_x1a->Draw("L SAME");
  g1b_x0->SetTitle("Sine Explicit");
  g1b_x0->GetXaxis()->SetTitle("t (s)");
  g1b_x0->GetYaxis()->SetTitle("x (-)");
  c1b->SetTitle("Sine Explicit");

  leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(g1b_x0, "Explicit x0 (fixed)", "l");
  leg->AddEntry(g1b_x1, "Explicit x1 (fixed)", "l");
  leg->AddEntry(g1b_x0a, "Explicit x0 (adaptive)", "l");
  leg->AddEntry(g1b_x1a, "Explicit x1 (adaptive)", "l");
  leg->Draw();
  c1b->Update();
#endif

  SinCosSemiExplicit system_semiexplicit;
  Fehlberg45<Real, 2, 0> integrator_semiexplicit(std::make_shared<SinCosSemiExplicit<Real>>(system_semiexplicit));
  integrator_semiexplicit.enable_projection();
  integrator_semiexplicit.disable_reverse_mode();
  Solution<Real, 2, 0> solution_semiexplicit;
  Solution<Real, 2, 0> solution_semiexplicit_adaptive;
  integrator_semiexplicit.solve(time, system_semiexplicit.ics(), solution_semiexplicit);
  integrator_semiexplicit.adaptive_solve(time, system_semiexplicit.ics(), solution_semiexplicit_adaptive);

#ifdef SANDALS_ENABLE_PLOTTING
  TCanvas *c1c = new TCanvas("c1c", "Oscillator SemiExplicit", 800, 600);
  TGraph *g1c_x0 = to_TGraph(solution_semiexplicit.t, solution_semiexplicit.x.row(0));
  TGraph *g1c_x1 = to_TGraph(solution_semiexplicit.t, solution_semiexplicit.x.row(1));
  TGraph *g1c_x0a = to_TGraph(solution_semiexplicit_adaptive.t, solution_semiexplicit_adaptive.x.row(0));
  TGraph *g1c_x1a = to_TGraph(solution_semiexplicit_adaptive.t, solution_semiexplicit_adaptive.x.row(1));

  g1c_x0->SetLineColor(colors[0]);  g1c_x0->Draw("AL");
  g1c_x1->SetLineColor(colors[0]);  g1c_x1->Draw("L SAME");
  g1c_x0a->SetLineColor(colors[1]); g1c_x0a->SetLineStyle(9); g1c_x0a->Draw("L SAME");
  g1c_x1a->SetLineColor(colors[1]); g1c_x1a->SetLineStyle(9); g1c_x1a->Draw("L SAME");
  g1c_x0->SetTitle("Oscillator SemiExplicit");
  g1c_x0->GetXaxis()->SetTitle("t (s)");
  g1c_x0->GetYaxis()->SetTitle("x (-)");
  c1c->SetTitle("Oscillator SemiExplicit");

  leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(g1c_x0, "SemiExplicit x0 (fixed)", "l");
  leg->AddEntry(g1c_x1, "SemiExplicit x1 (fixed)", "l");
  leg->AddEntry(g1c_x0a, "SemiExplicit x0 (adaptive)", "l");
  leg->AddEntry(g1c_x1a, "SemiExplicit x1 (adaptive)", "l");
  leg->Draw();
  c1c->Update();
#endif

#ifdef SANDALS_ENABLE_PLOTTING
  app.Run();
#endif

  return 0;
}

