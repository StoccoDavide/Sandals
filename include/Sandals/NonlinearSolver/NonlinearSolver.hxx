/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

  /**
  * \brief Class container for the nonlinear solver.
  *
  * \includedoc docs/markdown/NonlinearSolver.md
  *
  * \tparam N The dimension of the nonlinear system of equations.
  */
  template <Size N>
  class NonlinearSolver
  {
  public:
    using Vector   = Eigen::Vector<Real, N>;                        /**< Templetized vector type. */
    using Matrix   = Eigen::Matrix<Real, N, N>;                     /**< Templetized matrix type. */
    using Function = std::function<void(Vector const &, Vector &)>; /**< Nonlinear function type. */
    using Jacobian = std::function<void(Vector const &, Matrix &)>; /**< Jacobian function type. */

  protected:
    Function m_function;                      /**< Nonlinear function \f$ \mathbf{F} \f$. */
    Jacobian m_jacobian;                      /**< Jacobian function \f$ \mathbf{JF}_{\mathbf{x}} \f$. */
    Real     m_tolerance{EPSILON_HIGH};       /**< Solver tolerance \f$ \epsilon \f$. */
    Real     m_alpha{0.8};                    /**< Relaxation factor \f$ \alpha \f$. */
    Size     m_max_iterations{100};           /**< Maximum number of allowed algorithm iterations. */
    Size     m_max_function_evaluations{100}; /**< Maximum number of allowed function evaluations. */
    Size     m_max_jacobian_evaluations{100}; /**< Maximum number of allowed Jacobian evaluations. */
    Size     m_max_relaxations{10};           /**< Maximum number of allowed algorithm relaxations. */
    bool     m_verbose{false};                /**< Verbose mode boolean flag. */
    bool     m_converged{false};              /**< Convergence boolean flag. */
    Real     m_residuals{0.0};                /**< Function residuals. */
    Size     m_iterations{0};                 /**< Algorithm iterations number. */
    Size     m_function_evaluations{0};       /**< Function evaluations number. */
    Size     m_jacobian_evaluations{0};       /**< Jacobian evaluations number. */
    Size     m_relaxations{0};                /**< Algorithm relaxations number. */

  public:
    /**
    * Class constructor for the nonlinear solver.
    */
    NonlinearSolver() {}

    /**
    * Get the tolerance \f$ \epsilon \f$.
    * \return The tolerance \f$ \epsilon \f$.
    */
    Real tolerance() const {return this->m_tolerance;}

    /**
    * Set the tolerance \f$ \epsilon \f$ for which the nonlinear solver
    * stops, i.e., \f$ \left\| \mathbf{F}(\mathbf{x}) \right\| < \epsilon \f$.
    * \param[in] t_tolerance The tolerance \f$ \epsilon \f$.
    */
    void tolerance(Real t_tolerance) {
      SANDALS_ASSERT(
        !std::isnan(t_tolerance) && std::isfinite(t_tolerance) && t_tolerance > Real(0.0),
        "Sandals::NonlinearSolver::tolerance(...): invalid input detected.");
      this->m_tolerance = t_tolerance;
    }

    /**
    * Get the number of maximum allowed iterations.
    * \return The number of maximum allowed iterations.
    */
    Size max_iterations() const {return this->max_iterations;}

    /**
    * Set the number of maximum allowed iterations.
    * \param[in] t_max_iterations The number of maximum allowed iterations.
    */
    void max_iterations(Size t_max_iterations) {
      SANDALS_ASSERT(!std::isnan(t_max_iterations) && std::isfinite(t_max_iterations),
        "Sandals::NonlinearSolver::max_iterations(...): invalid input detected.");
      this->m_max_iterations = t_max_iterations;
    }

    /**
    * Get the number of maximum allowed function evaluations.
    * \return The number of maximum allowed function evaluations.
    */
    Size max_function_evaluations() const {return this->m_max_function_evaluations;}

    /**
    * Set the number of maximum allowed function evaluations.
    * \param[in] t_max_fun_evaluations The number of maximum allowed function evaluations.
    */
    void max_function_evaluations(Size t_max_fun_evaluations)
    {
      SANDALS_ASSERT(!std::isnan(t_max_fun_evaluations) && std::isfinite(t_max_fun_evaluations),
        "Sandals::NonlinearSolver::max_function_evaluations(...): invalid input detected.");
      this->m_max_function_evaluations = t_max_fun_evaluations;
    }

    /**
    * Get the number of maximum allowed Jacobian evaluations.
    * \return The number of maximum allowed Jacobian evaluations.
    */
    Size max_jacobian_evaluations() const {return this->m_max_jacobian_evaluations;}

    /**
    * Set the number of maximum allowed Jacobian evaluations.
    * \param[in] t_max_jac_evaluations The number of maximum allowed Jacobian evaluations.
    */
    void max_jacobian_evaluations(Size t_max_jac_evaluations)
    {
      SANDALS_ASSERT(!std::isnan(t_max_jac_evaluations) && std::isfinite(t_max_jac_evaluations),
        "Sandals::NonlinearSolver::max_jacobian_evaluations(...): invalid input detected.");
      this->m_max_jacobian_evaluations = t_max_jac_evaluations;
    }

    /**
    * Get the number of maximum allowed relaxations.
    * \return The number of maximum allowed relaxations.
    */
    Size max_relaxations() const {return this->max_relaxations;}

    /**
    * Set the number of maximum allowed relaxations.
    * \param[in] t_max_relaxations The number of maximum allowed relaxations.
    */
    void max_relaxations(Size t_max_relaxations)
    {
      SANDALS_ASSERT(!std::isnan(t_max_relaxations) && std::isfinite(t_max_relaxations),
        "Sandals::NonlinearSolver::max_relaxations(...): invalid input detected.");
      this->m_max_relaxations = t_max_relaxations;
    }

    /**
    * Get relaxation factor \f$ \alpha \f$.
    * \return The relaxation factor \f$ \alpha \f$.
    */
    Real alpha() const {return this->m_alpha;}

    /**
    * Set relaxation factor \f$ \alpha \f$.
    * \param[in] t_alpha The relaxation factor \f$ \alpha \f$.
    */
    void alpha(Real t_alpha)
    {
      SANDALS_ASSERT(!std::isnan(t_alpha) && std::isfinite(t_alpha) && Real(0.0) <= t_alpha && t_alpha <= Real(1.0),
        "Sandals::NonlinearSolver::alpha(...): invalid input detected.");
      this->m_alpha = t_alpha;
    }

    /**
    * Enable solver's verbose mode.
    */
    void enable_verbose_mode() {this->m_verbose = true;}

    /**
    * Disable solver's verbose mode.
    */
    void disable_verbose_mode() {this->m_verbose = false;}

    /**
    * Set solver's verbose mode.
    * \param[in] t_verbose The verbose mode boolean flag.
    */
    void verbose_mode(bool t_verbose) {this->m_verbose = t_verbose;}

    /**
    * Get the number of algorithm iterations.
    */
    Size iterations() const {return this->m_iterations;}

    /**
    * Get the number of function evaluations.
    * \return The number of function evaluations.
    */
    Size function_evaluations() const {return this->m_function_evaluations;}

    /**
    * Get the number of function evaluations.
    * \return The number of function evaluations.
    */
    Size jacobian_evaluations() const {return this->m_jacobian_evaluations;}

    /**
    * Get the number of algorithm relaxations.
    * \return The number of algorithm relaxations.
    */
    Size relaxations() const {return this->m_relaxations;}

    /**
    * Get the number of function evaluations.
    * \return The number of function evaluations.
    */
    Real residuals() const {return this->m_residuals;}

    /**
    * Get the convergence boolean flag.
    * \return The convergence boolean flag.
    */
    bool converged() const {return this->m_converged;}

    /**
    * Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    * \param[in] t_function The nonlinear function pointer.
    * \param[in] t_jacobian The Jacobian of the nonlinear function pointer.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    bool solve(Function t_function, Jacobian t_jacobian, Vector const &x_ini, Vector &x_sol)
    {
      this->m_function = t_function;
      this->m_jacobian = t_jacobian;
      return this->solve(x_ini, x_sol);
    }

    /**
    * Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$ with damping.
    * \param[in] t_function The nonlinear function pointer.
    * \param[in] t_jacobian The Jacobian of the nonlinear function pointer.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    bool solve_damped(Function t_function, Jacobian t_jacobian, Vector const &x_ini, Vector &x_sol)
    {
      this->m_function = t_function;
      this->m_jacobian = t_jacobian;
      return this->solve_damped(x_ini, x_sol);
    }

    /**
    * Get the nonlinear solver name.
    * \return The nonlinear solver name.
    */
    virtual std::string name() const = 0;

  protected:
    /**
    * Reset solver internal counters and variables.
    */
    void reset()
    {
      this->m_iterations           = Size(0);
      this->m_function_evaluations = Size(0);
      this->m_jacobian_evaluations = Size(0);
      this->m_relaxations          = Size(0);
      this->m_residuals            = Real(0.0);
      this->m_converged            = false;
    }

    /**
    * Evaluate function \f$ \mathbf{F}(\mathbf{x}) \f$.
    * \param[in] x Input point.
    * \param[out] function Output value.
    */
    void evaluate_function(Vector const &x, Vector &function)
    {
      ++this->m_function_evaluations;
      this->m_function(x, function);
    }

    /**
    * Evaluate Jacobian \f$ \mathbf{J}(\mathbf{x}) \f$.
    * \param[in] x Input point.
    * \param[out] jacobian Output value.
    */
    void evaluate_jacobian(Vector const &x, Matrix &jacobian)
    {
      ++this->m_jacobian_evaluations;
      this->m_jacobian(x, jacobian);
    }

    /**
    * Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    virtual bool solve(Vector const &x_ini, Vector &x_sol) = 0;

    /**
    * Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$ with damping.
    * \param[in] x_ini The initialization point.
    * \param[out] x_sol The solution point.
    */
    virtual bool solve_damped(Vector const &x_ini, Vector &x_sol) = 0;

  }; // class NonlinearSolver

} // namespace Sandals

#endif // SANDALS_NONLINEAR_SOLVER_HXX
