/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_NONLINEAR_SOLVER_HXX
#define SANDALS_NONLINEAR_SOLVER_HXX

namespace Sandals
{

  /*\
   |   _   _             _ _                       ____        _
   |  | \ | | ___  _ __ | (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __
   |  |  \| |/ _ \| '_ \| | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__|
   |  | |\  | (_) | | | | | | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |
   |  |_| \_|\___/|_| |_|_|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|
   |
  \*/

  //! \brief Nonlinear solver base class container.
  //! Nonlinear solver base class container for a problem of the form \f$
  //! \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$ where \f$ \mathbf{F}: \mathbb{R}^n
  //! \rightarrow \mathbb{R}^n \f$.
  template <unsigned N>
  class NonlinearSolver
  {
  public:
    using vec = Eigen::Vector<real, N>;                  //!< Templetized vector type.
    using mat = Eigen::Matrix<real, N, N>;               //!< Templetized matrix type.
    using fun = std::function<void(vec const &, vec &)>; //!< Non-linear function type.
    using jac = std::function<void(vec const &, mat &)>; //!< Jacobian function type.

  protected:
    fun      m_fun;                      //!< Non-linear function.
    jac      m_jac;                      //!< Jacobian function.
    real     m_tolerance{EPSILON_LOW};  //!< Solver tolerance.
    real     m_alpha{0.8};               //!< Relaxation factor.
    unsigned m_max_iterations{100};      //!< Maximum number of allowed algorithm iterations.
    unsigned m_max_fun_evaluations{100}; //!< Maximum number of allowed function evaluations.
    unsigned m_max_jac_evaluations{100}; //!< Maximum number of allowed Jacobian evaluations.
    unsigned m_max_relaxations{10};      //!< Maximum number of allowed algorithm relaxations.
    bool     m_verbose{false};           //!< Verbose mode boolean.
    bool     m_converged{false};         //!< Convergence boolean.
    real     m_residuals{0.0};           //!< Function residuals.
    unsigned m_iterations{0};            //!< Algorithm iterations number.
    unsigned m_fun_evaluations{0};       //!< Function evaluations number.
    unsigned m_jac_evaluations{0};       //!< Jacobian evaluations number.
    unsigned m_relaxations{0};           //!< Algorithm relaxations number.

  public:
    //! Class constructor for a nonlinear solver.
    NonlinearSolver(void) {}

    //!< Get the tolerance \f$ \epsilon \f$.
    //! \return The tolerance \f$ \epsilon \f$.
    real tolerance(void) const {return this->m_tolerance;}

    //! Set the tolerance \f$ \epsilon \f$ for which the nonlinear solver
    //! stops, i.e., \f$ \left\| \mathbf{F}(\mathbf{x}) \right\| < \epsilon \f$.
    //! \param t_tolerance The tolerance \f$ \epsilon \f$.
    void
    tolerance(real t_tolerance)
    {
      SANDALS_ASSERT(
        !std::isnan(t_tolerance) && std::isfinite(t_tolerance) && t_tolerance > real(0.0),
        "Sandals::NonlinearSolver::tolerance(...): invalid input detected.");
      this->m_tolerance = t_tolerance;
    }

    //!< Get the number of maximum allowed iterations.
    //! \return The number of maximum allowed iterations.
    unsigned max_iterations(void) const {return this->max_iterations;}

    //! Set the number of maximum allowed iterations.
    //! \param t_max_evaluations The number of maximum allowed iterations.
    void
    max_iterations(unsigned t_max_iterations)
    {
      SANDALS_ASSERT(
        !std::isnan(t_max_iterations) && std::isfinite(t_max_iterations),
        "Sandals::NonlinearSolver::max_iterations(...): invalid input detected.");
      this->m_max_iterations = t_max_iterations;
    }

    //!< Get the number of maximum allowed function evaluations.
    //! \return The number of maximum allowed function evaluations.
    unsigned max_function_evaluations(void) const {return this->m_max_fun_evaluations;}

    //! Set the number of maximum allowed function evaluations.
    //! \param t_max_fun_evaluations The number of maximum allowed function evaluations.
    void
    max_function_evaluations(unsigned t_max_fun_evaluations)
    {
      SANDALS_ASSERT(
        !std::isnan(t_max_fun_evaluations) && std::isfinite(t_max_fun_evaluations),
        "Sandals::NonlinearSolver::max_function_evaluations(...): invalid input detected.");

      this->m_max_fun_evaluations = t_max_fun_evaluations;
    }

    //!< Get the number of maximum allowed Jacobian evaluations.
    //! \return The number of maximum allowed Jacobian evaluations.
    unsigned max_jacobian_evaluations(void) const {return this->m_max_jac_evaluations;}

    //! Set the number of maximum allowed Jacobian evaluations.
    //! \param t_max_jac_evaluations The number of maximum allowed Jacobian evaluations.
    void
    max_jacobian_evaluations(unsigned t_max_jac_evaluations)
    {
      SANDALS_ASSERT(
        !std::isnan(t_max_jac_evaluations) && std::isfinite(t_max_jac_evaluations),
        "Sandals::NonlinearSolver::max_jacobian_evaluations(...): invalid input detected.");

      this->m_max_jac_evaluations = t_max_jac_evaluations;
    }

    //!< Get the number of maximum allowed algorithm relaxations.
    //! \return The number of maximum allowed algorithm relaxations.
    unsigned max_relaxations(void) const {return this->max_relaxations;}

    //! Set the number of maximum allowed algorithm relaxations.
    //! \param t_max_relaxations The number of maximum allowed algorithm relaxations.
    void
    max_relaxations(unsigned t_max_relaxations)
    {
      SANDALS_ASSERT(
        !std::isnan(t_max_relaxations) && std::isfinite(t_max_relaxations),
        "Sandals::NonlinearSolver::max_relaxations(...): invalid input detected.");
      this->m_max_relaxations = t_max_relaxations;
    }

    //! Get relaxation factor \f$\alpha\f$.
    //! \return The relaxation factor \f$ \alpha \f$.
    real alpha(void) const {return this->m_alpha;}

    //! Set relaxation factor \f$ \alpha \f$.
    //! \param t_alpha The relaxation factor \f$ \alpha \f$.
    void
    alpha(real t_alpha)
    {
      SANDALS_ASSERT(
        !std::isnan(t_alpha) && std::isfinite(t_alpha) &&
        real(0.0) <= t_alpha && t_alpha <= real(1.0),
        "Sandals::NonlinearSolver::alpha(...): invalid input detected.");
      this->m_alpha = t_alpha;
    }

    //! Enable solver's verbose mode.
    void enable_verbose(void) {this->m_verbose = true;}

    //! Disable solver's verbose mode.
    void disable_verbose(void) {this->m_verbose = false;}

    //! Get the number of algorithm iterations.
    unsigned iterations(void) const {return this->m_iterations;}

    //! Get the number of function evaluations.
    //! \return The number of function evaluations.
    unsigned function_evaluations(void) const {return this->m_fun_evaluations;}

    //! Get the number of function evaluations.
    //! \return The number of function evaluations.
    unsigned jacobian_evaluations(void) const {return this->m_jac_evaluations;}

    //! Get the number of algorithm relaxations.
    //! \return The number of algorithm relaxations.
    unsigned relaxations(void) const {return this->m_relaxations;}

    //! Get the number of function evaluations.
    //! \return The number of function evaluations.
    real residuals(void) const {return this->m_residuals;}

    //! Get the convergence boolean value.
    //! \return The convergence boolean value.
    bool converged(void) const {return this->m_converged;}

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    //! \param t_fun The non-linear function pointer.
    //! \param t_jac The jacobian of the non-linear function pointer.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    bool solve(fun t_fun, jac t_jac, vec const &x_ini, vec &x_sol)
    {
      this->m_fun = t_fun;
      this->m_jac = t_jac;
      return this->solve(x_ini, x_sol);
    }

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$
    //! with dumping factor \f$ \alpha \f$.
    //! \param t_fun The non-linear function pointer.
    //! \param t_jac The jacobian of the non-linear function pointer.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    bool solve_dumped(fun t_fun, jac t_jac, vec const &x_ini, vec &x_sol)
    {
      this->m_fun = t_fun;
      this->m_jac = t_jac;
      return this->solve_dumped(x_ini, x_sol);
    }

    //! Get the nonlinear solver name.
    //! \return The nonlinear solver name.
    virtual std::string name(void) const = 0;

  protected:
    //! Reset solver internal counters and variables.
    void reset(void)
    {
      this->m_iterations      = unsigned(0);
      this->m_fun_evaluations = unsigned(0);
      this->m_jac_evaluations = unsigned(0);
      this->m_relaxations     = unsigned(0);
      this->m_residuals       = real(0.0);
      this->m_converged       = false;
    }

    //! Evaluate function \f$ \mathbf{F}(\mathbf{x}) \f$.
    //! \param x Input point.
    //! \param F Output value.
    void evaluate_function(vec const &x, vec &F)
    {
      ++this->m_fun_evaluations;
      this->m_fun(x, F);
    }

    //! Evaluate jacobian \f$ \mathbf{J}(\mathbf{x}) \f$.
    //! \param x Input point.
    //! \param J Output value.
    void evaluate_jacobian(vec const &x, mat &J)
    {
      ++this->m_jac_evaluations;
      this->m_jac(x, J);
    }

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    //! \param t_fun The function pointer.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    virtual bool solve(vec const &x_ini, vec &x_sol) = 0;

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$
    //! with dumping factor \f$ \alpha \f$.
    //! \param t_fun The function pointer.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    virtual bool solve_dumped(vec const &x_ini, vec &x_sol) = 0;

  }; // class NonlinearSolver

} // namespace Sandals

#endif // SANDALS_NONLINEAR_SOLVER_HXX
