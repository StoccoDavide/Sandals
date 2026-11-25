/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_GERDTS_SHOOTING_HH
#define TESTS_GERDTS_SHOOTING_HH

#include "Sandals/Ipopt.hh"
#include "Sandals/MultipleShooting.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"

#include "DAEProblemsAD.hh"

// Test utilities for the Gerdts examples: a problem transcribed by multiple shooting and solved by
// IPOPT with the tolerance 1e-12, and the Richardson extrapolation of order-two results.
namespace TestOCP
{
  using namespace Sandals;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;

  // Quiet IPOPT with the default tolerance (1e-12); SANDALS_IPOPT_DEBUG=1 prints the iterations.
  inline IpoptOptions options()
  {
    IpoptOptions opts;
    opts.integers["print_level"] = std::getenv("SANDALS_IPOPT_DEBUG") ? 5 : 0;
    opts.strings["sb"] = "yes";
    opts.integers["max_iter"] = 500;
    return opts;
  }

  // A problem transcribed by multiple shooting with the Radau IIA (order 5) method, and solved by IPOPT.
  template <Integer N, Integer U, Integer M>
  struct Shooting
  {
    using FunctionU = std::function<Eigen::Vector<Real, U>(Real)>;

    MultipleShooting<Real, 3, N, U, M> nlp;
    Real J{0.0};
    VectorX multipliers;

    template <typename Problem>
    Shooting(std::shared_ptr<Problem> problem, Integer const intervals, Integer const substeps = 1)
      : nlp(problem, RadauIIA5Tableau<Real>(), intervals, substeps) {}

    // Start from a forward simulation (zero controls by default), and solve
    bool solve(Eigen::Vector<Real, N> const & x_ini, FunctionU const & u_guess = nullptr, IpoptOptions const & opts = options())
    {
      auto const zero{[] (Real) {return Eigen::Vector<Real, U>::Zero().eval();}};
      return nlp.simulate_guess(x_ini, u_guess ? u_guess : zero) && resolve(opts);
    }

    // Solve again from the current solution (e.g., after changing the bounds)
    bool resolve(IpoptOptions const & opts = options())
    {return ipopt_solve(nlp, opts, J, nullptr, &multipliers) == Ipopt::Solve_Succeeded;}

    Eigen::Vector<Real, N> x(Integer const k) const {return nlp.x(nlp.guess(), k);}
    Eigen::Vector<Real, U> u(Integer const k) const {return nlp.u(nlp.guess(), k);}
    Real t(Integer const k) const {return nlp.t_mesh()(k);}
    Integer intervals() const {return nlp.intervals();}

    // The piecewise constant controls of the solution, as a function of time (to warm-start another
    // transcription on the same mesh)
    FunctionU controls() const
    {
      return [this] (Real const t) {
        Integer k{0};
        while (k + 1 < this->intervals() && this->t(k + 1) <= t) {++k;}
        return this->u(k);
      };
    }

    // Multiplier of the boundary condition i
    Real boundary_multiplier(Integer const i) const {return multipliers(nlp.boundary_row() + i);}
  };

  // Richardson extrapolation of second-order approximations on meshes of n and 2n intervals.
  inline Real richardson(Real const coarse, Real const fine) {return fine + (fine - coarse)/3.0;}

} // namespace TestOCP

#endif // TESTS_GERDTS_SHOOTING_HH
