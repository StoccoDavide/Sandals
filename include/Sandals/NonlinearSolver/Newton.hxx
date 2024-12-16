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

#ifndef SANDALS_NEWTON_HXX
#define SANDALS_NEWTON_HXX

#include "Sandals/NonlinearSolver/NonlinearSolver.hxx"

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

  //! \brief Class container for a dumped Newton's method with affine invariant step.
  //! Class container for a dumped Newton's method with affine invariant step.
  //! Given a zeros of a vectorial function problem of the form \f$ \mathbf{F}
  //! (\mathbf{x}) = \mathbf{0} \f$, where \f$ \mathbf{F}: \mathbb{R}^n \rightarrow
  //! \mathbb{R}^n \f$, then the Newton's method is defined as:
  //!
  //! \f[
  //!   \mathbf{JF}_{\mathbf{x}}(\mathbf{x}_k) \mathbf{h} = -\mathbf{F}(\mathbf{x}_k).
  //! \f]
  //!
  //! The dumped step is defined as:
  //!
  //! \f[
  //!   \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{h}
  //! \f]
  //!
  //! where \f$ \alpha_k \f$ is a dumping coefficient that satisfies:
  //!
  //! \f[
  //! \left\| \mathbf{JF}_{\mathbf{x}}(\mathbf{x}_k)^{-1} \mathbf{F}(\mathbf{x}_{k+1}) \right\|
  //! \leq \left(1 - \dfrac{\alpha_k}{2}\right) \left\| \mathbf{JF}_{\mathbf{x}}(\mathbf{x}_k)^{-1}
  //! \mathbf{F}(\mathbf{x}_k) \right\| = \left(1 - \dfrac{\alpha_k}{2} \right)
  //! \left\| \mathbf{h}  \right\|.
  //! \f]
  //!
  //! For more details on the Newton's method with affine invariant step refer to:
  //! https://www.zib.de/deuflhard/research/algorithm/ainewton.en.html.
  //! Newton affine solver class container.
  template <unsigned N>
  class Newton : public NonlinearSolver<N>
  {
  public:
    using vec = typename NonlinearSolver<N>::vec; //!< Templetized vector type.
    using mat = typename NonlinearSolver<N>::mat; //!< Templetized matrix type.
    using fun = typename NonlinearSolver<N>::fun; //!< Non-linear function type.
    using jac = typename NonlinearSolver<N>::jac; //!< Jacobian function type.
    using NonlinearSolver<N>::solve;
    using NonlinearSolver<N>::solve_dumped;

    //! Class constructor for a quasi-Newton solver.
    Newton(void) {}

    //! Get the quasi-Newton solver name.
    //! \return The quasi-Newton solver name.
    std::string name(void) const override {return "Newton";}

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    bool solve(vec const &x_ini, vec &x_sol) override
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      mat J0, J1;
      vec X0, F0, D0, X1, F1, D1;
      real F0_norm = real(0.0);
      real D0_norm = real(0.0);

      // Set initial iteration
      X0 = x_ini;
      this->evaluate_function(X0, F0);
      this->evaluate_jacobian(X0, J0);

      // Algorithm iterations
      real tolerance_F_norm = this->m_tolerance;
      real tolerance_D_norm = tolerance_F_norm * tolerance_F_norm;
      this->m_converged = false;
      for (this->m_iterations = unsigned(1);
          this->m_iterations < this->m_max_iterations;
           ++this->m_iterations) {

        // Calculate step
        this->step(F0, J0, D0);

        // Check convergence
        F0_norm = F0.norm();
        D0_norm = D0.norm();
        if (F0_norm < tolerance_F_norm || D0_norm < tolerance_D_norm) {
          this->m_converged = true;
          break;
        }

        // Update point
        X1 = X0 + D0;
        this->evaluate_function(X1, F1);
        this->evaluate_jacobian(X1, J1);

        // Update internal variables
        X0 = X1;
        F0 = F1;
        D0 = D1;
        J0 = J1;
      }

      // Convergence data
      x_sol             = X0;
      this->m_residuals = F0_norm;
      return this->m_converged;
    }

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$
    //! with dumping factor \f$ \alpha \f$.
    //! \param t_fun The function pointer.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    bool solve_dumped(vec const &x_ini, vec &x_sol) override
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      mat J0, J1;
      vec X0, F0, D0, X1, F1, D1;
      real F0_norm = real(0.0);
      real D0_norm = real(0.0);
      real F1_norm = real(0.0);
      real D1_norm = real(0.0);

      // Set initial iteration
      X0 = x_ini;
      this->evaluate_function(X0, F0);
      this->evaluate_jacobian(X0, J0);

      // Algorithm iterations
      real tolerance_F_norm = this->m_tolerance;
      real tolerance_D_norm = this->m_tolerance * this->m_tolerance;
      this->m_converged     = false;
      real tau              = real(1.0);
      for (this->m_iterations = unsigned(1);
           this->m_iterations < this->m_max_iterations;
           ++this->m_iterations) {
        // Calculate step
        this->step(F0, J0, D0);

        // Check convergence
        F0_norm = F0.norm();
        D0_norm = D0.norm();
        if (F0_norm < tolerance_F_norm || D0_norm < tolerance_D_norm) {
          this->m_converged = true;
          break;
        }

        // Relax the iteration process
        tau = real(1.0);
        for (this->m_relaxations = unsigned(0);
             this->m_relaxations < this->m_max_relaxations;
             ++this->m_relaxations) {
          // Update point
          D1 = tau * D0;
          X1 = X0 + D1;
          this->evaluate_function(X1, F1);
          this->evaluate_jacobian(X1, J1);

          // Check relaxation
          F1_norm = F1.norm();
          D1_norm = D1.norm();
          if (F1_norm < F0_norm || D1_norm < (real(1.0)-tau/real(2.0))*D0_norm) {
            break;
          } else {
            tau *= this->m_alpha;
          }
        }

        // Update internal variables
        X0 = X1;
        F0 = F1;
        D0 = D1;
        J0 = J1;
      }

      // Convergence data
      x_sol             = X0;
      this->m_residuals = F0_norm;
      return this->m_converged;
    }

    //! Calculate the step.
    //! \param F Function.
    //! \param J Jacobian approximation.
    //! \param D Step.
    void step(vec const &F, mat const &J, vec &D) const {D = -J.inverse() * F;}

  }; // class Newton

} // namespace Sandals

#endif // SANDALS_NEWTON_HXX
