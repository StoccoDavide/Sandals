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

#ifndef SANDALS_BROYDEN_HXX
#define SANDALS_BROYDEN_HXX

namespace Sandals
{

  /*\
   |   ____                      _
   |  | __ ) _ __ ___  _   _  __| | ___ _ __
   |  |  _ \| '__/ _ \| | | |/ _` |/ _ \ '_ \
   |  | |_) | | | (_) | |_| | (_| |  __/ | | |
   |  |____/|_|  \___/ \__, |\__,_|\___|_| |_|
   |                   |___/
  \*/

  //! Broyden solver class container.
  template <unsigned N>
  class Broyden : public NonlinearSolver<N>
  {
  public:
    using vec = typename NonlinearSolver<N>::vec; //!< Templetized vector type.
    using mat = typename NonlinearSolver<N>::mat; //!< Templetized matrix type.
    using fun = typename NonlinearSolver<N>::fun; //!< Non-linear function type.
    using jac = typename NonlinearSolver<N>::jac; //!< Jacobian function type.
    using NonlinearSolver<N>::solve;
    using NonlinearSolver<N>::solve_dumped;

    //! Class constructor for a Broyden solver.
    Broyden(void) {}

    //! Get the Broyden solver name.
    //! \return The Broyden solver name.
    std::string name(void) const override {return "Broyden";}

    //! Solve non-linear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    //! \param t_fun The function pointer.
    //! \param x_ini The initialization point.
    //! \param x_sol The solution point.
    bool solve(vec const &x_ini, vec &x_sol) override
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      mat J0, J1;
      vec X0, F0, D0, DX0, DF0, X1, F1, D1, DX1, DF1;
      real F0_norm = real(0.0);
      real D0_norm = real(0.0);

      // Set initial iteration
      X0 = x_ini;
      this->evaluate_function(X0, F0);
      this->evaluate_jacobian(X0, J0);

      // Algorithm iterations
      real tolerance_F_norm = this->m_tolerance;
      real tolerance_D_norm = this->m_tolerance * this->m_tolerance;
      this->m_converged     = false;
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

        // Update jacobian approximation
        DX1 = X1 - X0;
        DF1 = F1 - F0;
        this->update(
          DX0, DF0, J0, // Old step data
          DX1, DF1, J1  // New step data
        );

        // Update internal variables
        X0  = X1;
        F0  = F1;
        DX0 = DX1;
        DF0 = DF1;
        D0  = D1;
        J0  = J1;
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
      vec X0, F0, D0, DX0, DF0, X1, F1, D1, DX1, DF1;
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

          // Check relaxation
          F1_norm = F1.norm();
          D1_norm = D1.norm();
          if (F1_norm < F0_norm || D1_norm < (real(1.0)-tau/real(2.0))*D0_norm) {
            break;
          } else {
            tau *= this->m_alpha;
          }
        }

        // Update jacobian approximation
        DX1 = X1 - X0;
        DF1 = F1 - F0;
        this->update(
          DX0, DF0, J0, // Old step data
          DX1, DF1, J1  // New step data
        );

        // Update internal variables
        X0  = X1;
        F0  = F1;
        DX0 = DX1;
        DF0 = DF1;
        D0  = D1;
        J0  = J1;
      }

      // Convergence data
      x_sol             = X0;
      this->m_residuals = F0_norm;
      return this->m_converged;
    }

    //! Jacobian approximation update rule.
    //! \param DX0 Old point difference.
    //! \param DF0 Old function value difference.
    //! \param J0  Old jacobian approximation.
    //! \param DX1 New point difference.
    //! \param DF1 New function value difference.
    //! \param J1  New jacobian approximation.
    void update(
      vec const &DX0, vec const &DF0, mat const &J0,
      vec const &DX1, vec const &DF1, mat       &J1
    ) {
      vec J0DF1(J0 * DF1);
      real DX1J0DF1 = DX1.transpose() * J0DF1;
      real DF1DF1   = DF1.transpose() * DF1;
      // Selection criteria
      // |DX1'*DX0|/|DX1'*J0*DF1| < |DF1'*DF0|/(DF1'*DF1)
      if (std::abs(DX1.transpose() * DX0) / std::abs(DX1J0DF1) <
          std::abs(DF1.transpose() * DF0) / DF1DF1 ||
          this->iterations() < unsigned(2)) {
        // Broyden's Good solver
        // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C = J0'*DX1;
        vec C(J0.transpose() * DX1);
        J1 = J0 - (J0DF1 - DX1) / (C.transpose() * DF1) * C.transpose();
      } else {
        // Broyden's Bad solver
        // J1 = J0 - (J0*DF1-DX1)/(C'*DF1)*C', where C = DF1;
        J1 = J0 - (J0DF1 - DX1) / (DF1DF1)*DF1.transpose();
      }
    }

    //! Calculate the step.
    //! \param F Function.
    //! \param J Jacobian approximation.
    //! \param D Step.
    void step(vec const &F, mat const &J, vec &D) const {D = -J * F;}

  }; // class Broyden

} // namespace Sandals

#endif // SANDALS_BROYDEN_HXX
