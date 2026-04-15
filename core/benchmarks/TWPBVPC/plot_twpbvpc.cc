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
#include "Sandals/RungeKutta/GaussLegendre4.hh"
#include "TWPBVPC.hh"

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
#define PROBLEM_INIT(PROBLEM, INTEGRATOR)                                     \
  constexpr Integer D{PROBLEM##Explicit<Real>::equations()};                  \
  PROBLEM##Problem<Real, PROBLEM##Explicit<Real>, INTEGRATOR<Real, D, 0>>     \
      PROBLEM##_explicit;                                                     \
  PROBLEM##Problem<Real, PROBLEM##Implicit<Real>, INTEGRATOR<Real, D, 0>>     \
      PROBLEM##_implicit;                                                     \
  PROBLEM##Problem<Real, PROBLEM##SemiExplicit<Real>, INTEGRATOR<Real, D, 0>> \
      PROBLEM##_semiexplicit;
#endif

#ifndef PROBLEM_SOLVE
#define PROBLEM_SOLVE(PROBLEM)                                \
  /* Set stiffness parameter */                               \
  constexpr Real lambda{1.0e-3};                              \
  static_cast<PROBLEM##Explicit<Real> *>(                     \
      PROBLEM##_explicit.integrator()->system())              \
      ->lambda(lambda);                                       \
  static_cast<PROBLEM##Implicit<Real> *>(                     \
      PROBLEM##_implicit.integrator()->system())              \
      ->lambda(lambda);                                       \
  static_cast<PROBLEM##SemiExplicit<Real> *>(                 \
      PROBLEM##_semiexplicit.integrator()->system())          \
      ->lambda(lambda);                                       \
  /* Set verbose mode */                                      \
  constexpr bool verbose{false};                              \
  PROBLEM##_explicit.verbose_mode(verbose);                   \
  PROBLEM##_implicit.verbose_mode(verbose);                   \
  PROBLEM##_semiexplicit.verbose_mode(verbose);               \
  PROBLEM##_explicit.integrator()->verbose_mode(false);       \
  PROBLEM##_implicit.integrator()->verbose_mode(false);       \
  PROBLEM##_semiexplicit.integrator()->verbose_mode(false);   \
  /* Set reverse mode */                                      \
  constexpr bool reverse{true};                               \
  PROBLEM##_explicit.integrator()->reverse_mode(reverse);     \
  PROBLEM##_implicit.integrator()->reverse_mode(reverse);     \
  PROBLEM##_semiexplicit.integrator()->reverse_mode(reverse); \
  /* Set solver parameters */                                 \
  PROBLEM##_explicit.sigma(1.0);                              \
  PROBLEM##_implicit.sigma(1.0);                              \
  PROBLEM##_semiexplicit.sigma(1.0);                          \
  /* Set solver tolerance */                                  \
  PROBLEM##_explicit.tolerance(1.0e-14);                      \
  PROBLEM##_implicit.tolerance(1.0e-14);                      \
  PROBLEM##_semiexplicit.tolerance(1.0e-14);                  \
  /* Set solver maximum number of iterations */               \
  PROBLEM##_explicit.max_iterations(100);                     \
  PROBLEM##_implicit.max_iterations(100);                     \
  PROBLEM##_semiexplicit.max_iterations(100);                 \
  /* Set solution parameters */                               \
  constexpr Integer num_subintervals{1};                      \
  PROBLEM##_explicit.subintervals(num_subintervals);          \
  PROBLEM##_implicit.subintervals(num_subintervals);          \
  PROBLEM##_semiexplicit.subintervals(num_subintervals);      \
  /* Set time mesh */                                         \
  constexpr Integer num_points{100};                          \
  Eigen::Vector<Real, Eigen::Dynamic> time(                   \
      Eigen::Vector<Real, Eigen::Dynamic>::LinSpaced(         \
          num_points,                                         \
          PROBLEM##_explicit.time_start(),                    \
          PROBLEM##_explicit.time_end()));                    \
  /* Set initial guess */                                     \
  Eigen::Matrix<Real, D, Eigen::Dynamic> guess(               \
      PROBLEM##_explicit.guess(time));                        \
  /* Solve the problems*/                                     \
  PROBLEM##_explicit.multiple_shooting(time, guess);          \
  PROBLEM##_implicit.multiple_shooting(time, guess);          \
  PROBLEM##_semiexplicit.multiple_shooting(time, guess);
#endif

#ifndef PROBLEM_PLOT
#define PROBLEM_PLOT(PROBLEM)                                       \
  auto sol_e = PROBLEM##_explicit.solution();                       \
  auto sol_i = PROBLEM##_implicit.solution();                       \
  auto sol_s = PROBLEM##_semiexplicit.solution();                   \
  auto sol_a = PROBLEM##_explicit.exact_solution(                   \
      reverse ? sol_e.t.reverse().eval() : sol_e.t);                \
  auto colors     = matlab_lines_colormap();                        \
  TCanvas *canvas = new TCanvas(#PROBLEM, #PROBLEM, 1200, 400);     \
  canvas->Divide(3, 1);                                             \
  canvas->cd(1);                                                    \
  TGraph *graph_e = to_TGraph(sol_e.t, sol_e.eigen_x(0));           \
  TGraph *graph_a = to_TGraph(sol_e.t, sol_a);                      \
  graph_e->SetTitle("Explicit vs Exact");                           \
  graph_e->SetLineColor(colors[0]);                                 \
  graph_e->Draw("AL");                                              \
  graph_a->SetLineColor(colors[1]);                                 \
  graph_a->Draw("L");                                               \
  graph_e->GetXaxis()->SetTitle("t");                               \
  graph_e->GetYaxis()->SetTitle("x");                               \
  graph_e->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff()); \
  TLegend *legend_1 = new TLegend(0.7, 0.8, 0.9, 0.9);              \
  legend_1->AddEntry(graph_e, "Numerical", "l");                    \
  legend_1->AddEntry(graph_a, "Exact", "l");                        \
  legend_1->Draw();                                                 \
  canvas->cd(2);                                                    \
  TGraph *graph_i = to_TGraph(sol_i.t, sol_i.eigen_x(0));           \
  graph_i->SetTitle("Implicit vs Exact");                           \
  graph_i->SetLineColor(colors[0]);                                 \
  graph_i->Draw("AL");                                              \
  graph_a->SetLineColor(colors[1]);                                 \
  graph_a->Draw("L");                                               \
  graph_i->GetXaxis()->SetTitle("t");                               \
  graph_i->GetYaxis()->SetTitle("x");                               \
  graph_i->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff()); \
  TLegend *legend_2 = new TLegend(0.7, 0.8, 0.9, 0.9);              \
  legend_2->AddEntry(graph_i, "Numerical", "l");                    \
  legend_2->AddEntry(graph_a, "Exact", "l");                        \
  legend_2->Draw();                                                 \
  canvas->cd(3);                                                    \
  TGraph *graph_s = to_TGraph(sol_s.t, sol_s.eigen_x(0));           \
  graph_s->SetTitle("Semi-Explicit vs Exact");                      \
  graph_s->SetLineColor(colors[0]);                                 \
  graph_s->Draw("AL");                                              \
  graph_a->SetLineColor(colors[1]);                                 \
  graph_a->Draw("L");                                               \
  graph_s->GetXaxis()->SetTitle("t");                               \
  graph_s->GetYaxis()->SetTitle("x");                               \
  graph_s->GetXaxis()->SetLimits(time.minCoeff(), time.maxCoeff()); \
  TLegend *legend_3 = new TLegend(0.7, 0.8, 0.9, 0.9);              \
  legend_3->AddEntry(graph_s, "Numerical", "l");                    \
  legend_3->AddEntry(graph_a, "Exact", "l");                        \
  legend_3->Draw();                                                 \
  canvas->Update();
#endif

#ifndef GENERATE_PLOT
#define GENERATE_PLOT(PROBLEM, INTEGRATOR)                           \
  try {                                                              \
    PROBLEM_INIT(PROBLEM, INTEGRATOR)                                \
    PROBLEM_SOLVE(PROBLEM)                                           \
    PROBLEM_PLOT(PROBLEM)                                            \
  } catch (const std::exception &e) {                                \
    std::cerr << "Error in " #PROBLEM ": " << e.what() << std::endl; \
  }
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
  GENERATE_PLOT(BVPT1, GaussLegendre4)
  GENERATE_PLOT(BVPT2, GaussLegendre4)
  GENERATE_PLOT(BVPT3, GaussLegendre4)
  GENERATE_PLOT(BVPT4, GaussLegendre4)
  GENERATE_PLOT(BVPT5, GaussLegendre4)
  GENERATE_PLOT(BVPT6, GaussLegendre4)
  GENERATE_PLOT(BVPT7, GaussLegendre4)
  GENERATE_PLOT(BVPT8, GaussLegendre4)
  GENERATE_PLOT(BVPT9, GaussLegendre4)
  GENERATE_PLOT(BVPT10, GaussLegendre4)
  GENERATE_PLOT(BVPT11, GaussLegendre4)
  GENERATE_PLOT(BVPT12, GaussLegendre4)
  GENERATE_PLOT(BVPT13, GaussLegendre4)
  GENERATE_PLOT(BVPT14, GaussLegendre4)
  GENERATE_PLOT(BVPT15, GaussLegendre4)
  GENERATE_PLOT(BVPT16, GaussLegendre4)
  GENERATE_PLOT(BVPT17, GaussLegendre4)
  GENERATE_PLOT(BVPT18, GaussLegendre4)
  GENERATE_PLOT(BVPT19, GaussLegendre4)
  GENERATE_PLOT(BVPT20, GaussLegendre4)
  GENERATE_PLOT(BVPT21, GaussLegendre4)
  GENERATE_PLOT(BVPT22, GaussLegendre4)
  GENERATE_PLOT(BVPT23, GaussLegendre4)
  GENERATE_PLOT(BVPT24, GaussLegendre4)
  GENERATE_PLOT(BVPT25, GaussLegendre4)
  GENERATE_PLOT(BVPT26, GaussLegendre4)
  GENERATE_PLOT(BVPT27, GaussLegendre4)
  GENERATE_PLOT(BVPT28, GaussLegendre4)
  GENERATE_PLOT(BVPT29, GaussLegendre4)
  GENERATE_PLOT(BVPT30, GaussLegendre4)
  GENERATE_PLOT(BVPT31, GaussLegendre4)
  GENERATE_PLOT(BVPT32, GaussLegendre4)
  GENERATE_PLOT(BVPT33, GaussLegendre4)

#ifdef SANDALS_ENABLE_PLOTTING
  app.Run();
#endif

  return 0;
}
