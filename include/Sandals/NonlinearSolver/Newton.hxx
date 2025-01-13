/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_NEWTON_HXX
#define SANDALS_NEWTON_HXX

namespace Sandals
{

  /*\
   |   _   _               _
   |  | \ | | _____      _| |_ ___  _ __
   |  |  \| |/ _ \ \ /\ / / __/ _ \| '_ \
   |  | |\  |  __/\ V  V /| || (_) | | | |
   |  |_| \_|\___| \_/\_/  \__\___/|_| |_|
   |
  \*/

  //! \brief Class container for the (damped) Newton's method with affine invariant step.
  //!
  //! \includedoc docs/markdown/Newton.md
  template <Size N>
  class Newton : public NonlinearSolver<N>
  {
  public:
    using Vector   = typename NonlinearSolver<N>::Vector;   //!< Templetized vector type.
    using Matrix   = typename NonlinearSolver<N>::Matrix;   //!< Templetized matrix type.
    using Function = typename NonlinearSolver<N>::Function; //!< Nonlinear function type.
    using Jacobian = typename NonlinearSolver<N>::Jacobian; //!< Jacobian function type.
    using NonlinearSolver<N>::solve;
    using NonlinearSolver<N>::solve_damped;

  private:
    Eigen::FullPivLU<Matrix> m_lu; //!< LU decomposition.

  public:
    //! Class constructor for the Newton solver.
    Newton() {}

    //! Get the Newton solver name.
    //! \return The Newton solver name.
    std::string name() const override {return "Newton";}

    //! Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    //! \param[in] x_ini The initialization point.
    //! \param[out] x_sol The solution point.
    bool solve(Vector const &x_ini, Vector &x_sol) override
    {
      // Reset internal variables
      this->reset();

      // Initialize variables
      Real residuals, step_norm;
      Vector x_old, x_new, function, step;
      Matrix jacobian;

      // Set initial iteration
      x_old = x_ini;
      this->evaluate_function(x_old, function);
      this->evaluate_jacobian(x_old, jacobian);

      // Algorithm iterations
      Real tolerance_residuals{this->m_tolerance};
      Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
      for (this->m_iterations = Size(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
      {

        // Calculate step
        this->m_lu.compute(jacobian);
        SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Newton::solve(...): singular Jacobian detected.");
        step = -this->m_lu.solve(function);

        // Check convergence
        residuals = function.norm();
        step_norm = step.norm();
        if (residuals < tolerance_residuals || step_norm < tolerance_step_norm) {
          this->m_converged = true;
          break;
        }

        // Update point
        x_new = x_old + step;
        this->evaluate_function(x_new, function);
        this->evaluate_jacobian(x_new, jacobian);

        // Update internal variables
        x_old = x_new;
      }

      // Convergence data
      x_sol             = x_old;
      this->m_residuals = residuals;
      return this->m_converged;
    }

    //! Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$ with damping.
    //! \param[in] x_ini The initialization point.
    //! \param[out] x_sol The solution point.
    bool solve_damped(Vector const &x_ini, Vector &x_sol) override
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      Real residuals_old, residuals_new, step_norm_old, step_norm_new, tau;
      Vector x_old, x_new, function_old, function_new, step_old, step_new;
      Matrix jacobian;

      // Set initial iteration
      x_old = x_ini;
      this->evaluate_function(x_old, function_old);
      this->evaluate_jacobian(x_old, jacobian);

      // Algorithm iterations
      Real tolerance_residuals{this->m_tolerance};
      Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
      for (this->m_iterations = Size(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
      {

        // Calculate step
        this->m_lu.compute(jacobian);
        SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Newton::solve_damped(...): singular Jacobian detected.");
        step_old = -this->m_lu.solve(function_old);

        // Check convergence
        residuals_old = function_old.norm();
        step_norm_old = step_old.norm();
        if (residuals_old < tolerance_residuals || step_norm_old < tolerance_step_norm) {
          this->m_converged = true;
          break;
        }

        // Relax the iteration process
        tau = Real(1.0);
        for (this->m_relaxations = Size(0); this->m_relaxations < this->m_max_relaxations; ++this->m_relaxations)
        {
          // Update point
          step_new = tau * step_old;
          x_new = x_old + step_new;
          this->evaluate_function(x_new, function_new);

          // Check relaxation
          residuals_new = function_new.norm();
          step_norm_new = step_new.norm();
          if (residuals_new < residuals_old || step_norm_new < (Real(1.0)-tau/Real(2.0))*step_norm_old) {
            this->evaluate_jacobian(x_new, jacobian);
            break;
          } else {
            tau *= this->m_alpha;
          }
        }

        // Update internal variables
        x_old        = x_new;
        function_old = function_new;
        step_old     = step_new;
      }

      // Convergence data
      x_sol             = x_old;
      this->m_residuals = residuals_old;
      return this->m_converged;
    }

  }; // class Newton

} // namespace Sandals

#endif // SANDALS_NEWTON_HXX
