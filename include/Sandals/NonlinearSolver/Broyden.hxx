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

  //! \brief Class container for the (damped) Broyden's methods: \em good, \em bad, and \em combined.
  //!
  //! Class container for the (damped) Broyden's  methods: \em good, \em bad, and \em combined. The
  //! \em combined Broyden's method is a combination of the Broyden's Good and Broyden's Bad solvers.
  //! Given a zeros of a vectorial function problem of the form \f$ \mathbf{F}(\mathbf{x}) =
  //! \mathbf{0} \f$, where \f$ \mathbf{F}: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n} \f$, the
  //! generic Broyden's method is defined as
  //!
  //! \f[ \mathbf{H}_k(\mathbf{x}_k) \mathbf{h} = -\mathbf{F}(\mathbf{x}_k) \text{,} \f]
  //!
  //! where \f$ \mathbf{H} \f$ is the (inverse) Jacobian approximation. The advancing
  //! step is then defined as
  //!
  //! \f[ \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{h} \text{,} \f]
  //!
  //! What distinguishes the Broyden's \em combined method from the generic Broyden's method is the
  //! update of the Jacobian approximation. The Broyden's \em combined method uses the following
  //! update rule, which is a combination of the Broyden's Good and Broyden's Bad solvers
  //!
  //! \f[ \left\{\begin{array}{ll} \mathbf{H}_{k+1}^{-1} = \mathbf{H}_k^{-1} - \displaystyle\frac{
  //! \mathbf{H}_k^{-1} \Delta\mathbf{F}_k - \Delta\mathbf{x}_k}{\mathbf{C}_g \Delta\mathbf{F}_k} \mathbf{C}_g
  //! & \left\| \displaystyle\frac{\Delta\mathbf{x}_k^\top \Delta\mathbf{x}_{k-1}}{\Delta\mathbf{x}_k^\top
  //! \mathbf{H}_k^{-1} \Delta\mathbf{F}_k} \right\| < \left\| \displaystyle\frac{\Delta\mathbf{F}_k^\top
  //! \Delta\mathbf{F}_{k-1}}{\Delta\mathbf{F}_k^\top\Delta\mathbf{F}_k} \right\| \\
  //! \mathbf{H}_{k+1}^{-1} = \mathbf{H}_k^{-1} - \displaystyle\frac{\mathbf{H}_k^{-1} \Delta\mathbf{F}_k
  //! -\Delta\mathbf{x}_k}{\mathbf{C}_b \Delta\mathbf{F}_k} \mathbf{C}_b & \text{otherwise}
  //! \end{array}\right. \text{,} \f]
  //!
  //! with \f$ \mathbf{C}_g = \Delta\mathbf{x}_k^\top \mathbf{H}_k^{-1} \f$, \f$ \mathbf{C}_b =
  //! \Delta\mathbf{F}_k^\top \f$, \f$ \Delta\mathbf{x}_k = \mathbf{x}_{k+1} - \mathbf{x}_k \f$,
  //! and \f$ \Delta\mathbf{F}_k = \mathbf{F}(\mathbf{x}_{k+1}) - \mathbf{F}(\mathbf{x}_k) \f$.
  //! For more details on the Broyden's \em combined method refer to the reference: <em> "Sobre dois
  //! métodos de Broyden" </em> by J. M. Martínez and J. M. Martínez, Matemática Aplicada e Computacional,
  //! IV Congresso Nacional de Matemática Aplicada e Computacional, Rio de Janeiro, Brasil, setembro de 1981.
  //! \note The implemented Broyden's class for the \em combined method can be used as a \em good or
  //! \em bad solver by setting the appropriate parameters.
  template <Size N>
  class Broyden : public NonlinearSolver<N>
  {
  public:
    using Type = enum class Type : Size {GOOD=0, BAD=1, COMBINED=2}; //!< Broyden solver type.
    using Vector   = typename NonlinearSolver<N>::Vector;   //!< Templetized vector type.
    using Matrix   = typename NonlinearSolver<N>::Matrix;   //!< Templetized matrix type.
    using Function = typename NonlinearSolver<N>::Function; //!< nonlinear function type.
    using Jacobian = typename NonlinearSolver<N>::Jacobian; //!< Jacobian function type.
    using NonlinearSolver<N>::solve;
    using NonlinearSolver<N>::solve_damped;

  private:
    Type m_type{Type::COMBINED}; //!< ODE system type.

  public:
    //! Class constructor for the Broyden solver.
    //!
    Broyden() {}

    //! Get the Broyden solver name.
    //! \return The Broyden solver name.
    std::string name() const override {return "Broyden";}

    //! Get the enumeration type of the Broyden solver.
    //! \return The Broyden solver enumeration type.
    Type type() const {return this->m_type;}

    //! Set the enumeration type of the Broyden solver.
    //! \param[in] t_type The Broyden solver enumeration type.
    void type(Type t_type) {this->m_type = t_type;}

    //! Enable the \em good Broyden solver.
    //!
    void enable_good_mode() {this->m_type = Type::GOOD;}

    //! Enable the \em bad Broyden solver.
    //!
    void enable_bad_mode() {this->m_type = Type::BAD;}

    //! Enable the \em combined Broyden solver.
    //!
    void enable_combined_mode() {this->m_type = Type::COMBINED;}

    //! Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$.
    //! \param[in] x_ini The initialization point.
    //! \param[out] x_sol The solution point.
    bool solve(Vector const &x_ini, Vector &x_sol) override
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      Real residuals, step_norm;
      Vector x_old, x_new, function_old, function_new, step_old, step_new, delta_x_old, delta_x_new,
        delta_function_old, delta_function_new;
      Matrix jacobian_old, jacobian_new;

      // Set initial iteration
      x_old = x_ini;
      this->evaluate_function(x_old, function_old);
      this->evaluate_jacobian(x_old, jacobian_old);

      // Algorithm iterations
      Real tolerance_residuals{this->m_tolerance};
      Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
      for (this->m_iterations = Size(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
      {

        // Calculate step
        step_old = -jacobian_old * function_old;

        // Check convergence
        residuals = function_old.norm();
        step_norm = step_old.norm();
        if (residuals < tolerance_residuals || step_norm < tolerance_step_norm) {
          this->m_converged = true;
          break;
        }

        // Update point
        x_new = x_old + step_old;
        this->evaluate_function(x_new, function_new);

        // Update jacobian approximation
        delta_x_new = x_new - x_old;
        delta_function_new = function_new - function_old;
        this->update(
          delta_x_old, delta_function_old, jacobian_old, // Old step data
          delta_x_new, delta_function_new, jacobian_new  // New step data
        );

        // Update internal variables
        x_old              = x_new;
        function_old       = function_new;
        delta_x_old        = delta_x_new;
        delta_function_old = delta_function_new;
        step_old           = step_new;
        jacobian_old       = jacobian_new;
      }

      // Convergence data
      x_sol             = x_old;
      this->m_residuals = residuals;
      return this->m_converged;
    }

    //! Solve nonlinear system of equations \f$ \mathbf{F}(\mathbf{x}) = \mathbf{0} \f$
    //! with damping factor \f$ \alpha \f$.
    //! \param[in] x_ini The initialization point.
    //! \param[out] x_sol The solution point.
    bool solve_damped(Vector const &x_ini, Vector &x_sol) override
    {
      // Setup internal variables
      this->reset();

      // Initialize variables
      Real residuals_old, residuals_new, step_norm_old, step_norm_new, tau;
      Vector x_old, x_new, function_old, function_new, step_old, step_new, delta_x_old, delta_x_new,
        delta_function_old, delta_function_new;
      Matrix jacobian_old, jacobian_new;

      // Set initial iteration
      x_old = x_ini;
      this->evaluate_function(x_old, function_old);
      this->evaluate_jacobian(x_old, jacobian_old);

      // Algorithm iterations
      Real tolerance_residuals{this->m_tolerance};
      Real tolerance_step_norm{this->m_tolerance * this->m_tolerance};
      for (this->m_iterations = Size(1); this->m_iterations < this->m_max_iterations; ++this->m_iterations)
      {

        // Calculate step
        step_old = -jacobian_old * function_old;

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
            break;
          } else {
            tau *= this->m_alpha;
          }
        }

        // Update jacobian approximation
        delta_x_new = x_new - x_old;
        delta_function_new = function_new - function_old;
        this->update(
          delta_x_old, delta_function_old, jacobian_old, // Old step data
          delta_x_new, delta_function_new, jacobian_new  // New step data
        );

        // Update internal variables
        x_old              = x_new;
        function_old       = function_new;
        delta_x_old        = delta_x_new;
        delta_function_old = delta_function_new;
        step_old           = step_new;
        jacobian_old       = jacobian_new;
      }

      // Convergence data
      x_sol             = x_old;
      this->m_residuals = residuals_old;
      return this->m_converged;
    }

    //! Jacobian approximation update rule for the Broyden's method.
    //! \param[in] delta_x_old Old difference between points.
    //! \param[in] delta_function_old Old difference between function values.
    //! \param[in] jacobian_old Old jacobian approximation.
    //! \param[in] delta_x_new New difference between points.
    //! \param[in] delta_function_new New difference between function values.
    //! \param[out] jacobian_new New jacobian approximation.
    void update(
      Vector const &delta_x_old, Vector const &delta_function_old, Matrix const &jacobian_old,
      Vector const &delta_x_new, Vector const &delta_function_new, Matrix       &jacobian_new
    ) {
      Vector tmp_1(jacobian_old * delta_function_new);
      Real tmp_2{delta_function_new.squaredNorm()};
      // Selection criteria
      // |(dx_new'*dx_old) / (dx_new'*J_old*dF_new)| < |(dF_new'*dF_old) / (dF_new'*dF_new)|
      if (this->m_type == Type::COMBINED || this->m_type == Type::GOOD || this->iterations() < Size(2) ||
          std::abs(delta_x_new.transpose() * delta_x_old) / std::abs(delta_x_new.transpose() * tmp_1)
          < std::abs(delta_function_new.transpose() * delta_function_old) / tmp_2) {
        // Broyden's Good solver
        // J_new = J_old - (J_old*dF_new-dx_new) / (C'*dF_new)*C', where C = J_old'*dx_new;
        Vector C_g(jacobian_old.transpose() * delta_x_new);
        jacobian_new = jacobian_old - (tmp_1 - delta_x_new) / (C_g.transpose() * delta_function_new) * C_g.transpose();
      } else {
        // Broyden's Bad solver
        // J_new = J_old - (J_old*dF_old-dx_new) / (C'*dF_old)*C', where C = dF_old;
        jacobian_new = jacobian_old - (tmp_1 - delta_x_new) / tmp_2 * delta_function_old.transpose();
      }
    }

  }; // class Broyden

} // namespace Sandals

#endif // SANDALS_BROYDEN_HXX
