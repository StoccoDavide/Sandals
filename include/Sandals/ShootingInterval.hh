/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_SHOOTING_INTERVAL_HH
#define SANDALS_SHOOTING_INTERVAL_HH

#include <Sandals.hh>
#include <Sandals/Tableau.hh>
#include <Sandals/OptimalControl.hh>

namespace Sandals
{

  /*\
   |   ____  _                 _   _             ___       _                       _
   |  / ___|| |__   ___   ___ | |_(_)_ __   __ _|_ _|_ __ | |_ ___ _ ____   ____ _| |
   |  \___ \| '_ \ / _ \ / _ \| __| | '_ \ / _` || || '_ \| __/ _ \ '__\ \ / / _` | |
   |   ___) | | | | (_) | (_) | |_| | | | | (_| || || | | | ||  __/ |   \ V / (_| | |
   |  |____/|_| |_|\___/ \___/ \__|_|_| |_|\__, |___|_| |_|\__\___|_|    \_/ \__,_|_|
   |                                       |___/
  \*/

  /**
  * \brief Integration of one shooting interval, with first- and second-order sensitivities.
  *
  * The interval \f$ [t_0, t_1] \f$ is integrated with \f$ s \f$ equal substeps of a Runge-Kutta
  * method, from the states \f$ \mathbf{x} \f$ and with the constant controls \f$ \mathbf{u} \f$. The
  * running cost is integrated together with the states, as the extra state \f$ q \f$ with \f$ q^\prime
  * = L(\mathbf{x}, \mathbf{u}, t) \f$ and \f$ q(t_0) = 0 \f$, so the integrated end point is
  * \f$ \mathbf{y} = [\mathbf{X}; Q] \f$ (the "augmented states"). Each substep solves the stage
  * equations, with \f$ \mathbf{K}_i = h\,\mathbf{y}^\prime_i \f$,
  *
  * \f[
  * \mathbf{G}_i(\mathbf{K}) = \mathbf{F}_a\left(\mathbf{y} + \textstyle\sum_j a_{ij}\mathbf{K}_j,
  * \mathbf{K}_i/h, \mathbf{u}, t + c_i h\right) = \mathbf{0}, \quad \mathbf{y}^+ = \mathbf{y} +
  * \textstyle\sum_i b_i\mathbf{K}_i \text{,}
  * \f]
  *
  * with Newton's method. The derivatives of the end point with respect to the parameters
  * \f$ \mathbf{p} = [\mathbf{x}; \mathbf{u}] \f$ follow from the implicit function theorem:
  * - first order (forward, through the substeps): \f$ \partial_{\mathbf{p}}\mathbf{K} =
  *   -(\partial_{\mathbf{K}}\mathbf{G})^{-1}\partial_{\mathbf{p}}\mathbf{G} \f$;
  * - second order (Hessian of \f$ \mathbf{w}^\top\mathbf{y} \f$, backward through the substeps): with
  *   the stage adjoints \f$ \boldsymbol{\nu} = -(\partial_{\mathbf{K}}\mathbf{G})^{-\top}(\mathbf{b}
  *   \otimes \boldsymbol{\lambda}) \f$, each substep contributes \f$ \sum_i \mathbf{S}_i^\top
  *   \nabla^2(\boldsymbol{\nu}_i^\top\mathbf{F}_a)\,\mathbf{S}_i \f$, where \f$ \mathbf{S}_i \f$ maps
  *   the parameters to the arguments of the \f$ i \f$-th stage, and the adjoint \f$ \boldsymbol{\lambda}
  *   \f$ is propagated backward by the transposed substep Jacobian.
  * \tparam Real The scalar number type.
  * \tparam S The number of stages of the Runge-Kutta method.
  * \tparam N The number of states.
  * \tparam U The number of controls.
  * \tparam M The number of invariants of the system.
  */
  template <typename Real, Integer S, Integer N, Integer U, Integer M = 0>
  class ShootingInterval
  {
  public:
    static constexpr Integer A{N + 1};   /**< Augmented states [x; q]. */
    static constexpr Integer P{N + U};   /**< Parameters [x; u] of the interval. */
    static constexpr Integer Q{N};       /**< Index of the running cost integral q in the augmented states. */

    using Problem   = OptimalControlProblem<Real, N, U, M>; /**< Optimal control problem type. */
    using VectorF   = Eigen::Vector<Real, N>;      /**< States vector type. */
    using VectorU   = Eigen::Vector<Real, U>;      /**< Controls vector type. */
    using VectorA   = Eigen::Vector<Real, A>;      /**< Augmented states vector type. */
    using MatrixAP  = Eigen::Matrix<Real, A, P>;   /**< Sensitivity of the augmented states to the parameters. */
    using MatrixPP  = Eigen::Matrix<Real, P, P>;   /**< Hessian in the parameters. */

  private:
    using MatrixAA  = Eigen::Matrix<Real, A, A>;
    using MatrixAU  = Eigen::Matrix<Real, A, U>;
    using MatrixK   = Eigen::Matrix<Real, A, S>;           // stage variables, one column per stage
    using VectorK   = Eigen::Vector<Real, A*S>;            // stage variables, stacked
    using MatrixKK  = Eigen::Matrix<Real, A*S, A*S>;       // stage Jacobian dG/dK
    using MatrixKL  = Eigen::Matrix<Real, A*S, A+U>;       // stage sensitivity to the substep parameters [y; u]
    using MatrixW   = Eigen::Matrix<Real, 2*A+U, 2*A+U>;   // Hessian in the stage arguments (y, y', u)
    using MatrixWL  = Eigen::Matrix<Real, 2*A+U, A+U>;     // map from [y; u] to the stage arguments
    using MatrixLL  = Eigen::Matrix<Real, A+U, A+U>;       // Hessian in the substep parameters [y; u]

    /** Data of one substep, stored by the forward pass for the derivatives. */
    struct Substep
    {
      Real    t{0.0};                       // initial time
      VectorA y;                            // initial augmented states
      MatrixK K;                            // stage variables at convergence
      Eigen::PartialPivLU<MatrixKK> lu;     // factorization of dG/dK at convergence
      MatrixKL K_l;                         // dK/d[y; u]
      MatrixAA phi_y;                       // d y+/d y
      MatrixAP Y;                           // d y/d p (sensitivity of the initial states)
    };

    Problem const &   m_problem;   /**< Optimal control problem (system and running cost). */
    Tableau<Real, S>  m_tableau;   /**< Butcher tableau. */
    Real              m_t_ini;     /**< Initial time of the interval. */
    Real              m_h;         /**< Substep size. */
    Integer           m_substeps;  /**< Number of substeps. */

    // Results of the last integration, valid for the parameters m_p
    Eigen::Vector<Real, P> m_p;
    bool                   m_valid{false};
    std::vector<Substep>   m_steps;
    VectorA                m_y_end;
    MatrixAP               m_Y_end;

    Integer m_max_newton_iterations{50}; /**< Maximum number of Newton iterations per substep. */
    bool    m_method_checked{false};     /**< True once the method has been checked against the system. */

  public:
    /**
    * Class constructor.
    * \param[in] t_problem The optimal control problem.
    * \param[in] t_tableau The Butcher tableau of the Runge-Kutta method.
    * \param[in] t_t_ini The initial time of the interval.
    * \param[in] t_t_end The final time of the interval.
    * \param[in] t_substeps The number of equal substeps.
    */
    ShootingInterval(Problem const & t_problem, Tableau<Real, S> const & t_tableau, Real const t_t_ini,
      Real const t_t_end, Integer const t_substeps)
      : m_problem(t_problem), m_tableau(t_tableau), m_t_ini(t_t_ini),
        m_h((t_t_end - t_t_ini) / t_substeps), m_substeps(t_substeps)
    {
      SANDALS_ASSERT(t_substeps >= 1 && t_t_end > t_t_ini, "Sandals::ShootingInterval: invalid interval.");
    }

    /** \return The number of substeps. */
    Integer substeps() const {return this->m_substeps;}

    /**
    * Integrate the interval (the result is cached for the given states and controls).
    * \param[in] x The initial states.
    * \param[in] u The controls.
    * \return True if every substep converged, false otherwise.
    */
    bool integrate(VectorF const & x, VectorU const & u)
    {
      Eigen::Vector<Real, P> p;
      p << x, u;
      if (this->m_valid && p == this->m_p) {return true;}

      if (!this->m_method_checked) {
        this->check_method(x, u);
        this->m_method_checked = true;
      }

      this->m_valid = false;
      this->m_p = p;
      this->m_steps.resize(this->m_substeps);

      VectorA y;
      y << x, Real(0.0); // the running cost integral starts from zero
      MatrixAP Y(MatrixAP::Zero());
      Y.template topLeftCorner<N, N>().setIdentity();

      for (Integer j{0}; j < this->m_substeps; ++j) {
        Substep & step{this->m_steps[j]};
        step.t = this->m_t_ini + j*this->m_h;
        step.y = y;
        step.Y = Y;
        if (!this->solve_stages(step)) {return false;}

        // Advance the states and the first-order sensitivity: Y+ = phi_y Y + [0, phi_u]
        y += step.K * this->m_tableau.b;
        MatrixAU const phi_u{this->substep_jacobian_u(step)};
        Y = step.phi_y * Y;
        Y.template rightCols<U>() += phi_u;
      }

      this->m_y_end = y;
      this->m_Y_end = Y;
      this->m_valid = true;
      return true;
    }

    /** \return The integrated end point \f$ [\mathbf{X}; Q] \f$. */
    VectorA const & end() const {return this->m_y_end;}

    /** \return The sensitivity \f$ \partial[\mathbf{X}; Q]/\partial[\mathbf{x}; \mathbf{u}] \f$. */
    MatrixAP const & sensitivity() const {return this->m_Y_end;}

    /**
    * Compute the Hessian of \f$ \mathbf{w}^\top[\mathbf{X}; Q] \f$ with respect to \f$ [\mathbf{x};
    * \mathbf{u}] \f$ (the interval must have been integrated).
    * \param[in] w The weights of the augmented end point.
    * \return The (symmetric) Hessian.
    */
    MatrixPP hessian(VectorA const & w) const
    {
      SANDALS_ASSERT(this->m_valid, "Sandals::ShootingInterval::hessian(...): integrate first.");
      MatrixPP H(MatrixPP::Zero());
      VectorA lambda{w}; // adjoint of the states at the end of the current substep

      for (Integer j{this->m_substeps - 1}; j >= 0; --j) {
        Substep const & step{this->m_steps[j]};

        // Stage adjoints: (dG/dK)^T nu = -(b (x) lambda)
        VectorK rhs;
        for (Integer i{0}; i < S; ++i) {rhs.template segment<A>(i*A) = -this->m_tableau.b(i) * lambda;}
        VectorK const nu{step.lu.transpose().solve(rhs)};

        // Hessian of the substep in its parameters [y; u], chained to [x; u] through Q = d[y; u]/d[x; u]
        MatrixLL H_step(MatrixLL::Zero());
        for (Integer i{0}; i < S; ++i) {
          MatrixWL const S_i{this->stage_map(step, i)};
          MatrixW const W_i{this->stage_hessian(step, i, nu.template segment<A>(i*A))};
          H_step.noalias() += S_i.transpose() * W_i * S_i;
        }
        Eigen::Matrix<Real, A+U, P> Q_j(Eigen::Matrix<Real, A+U, P>::Zero());
        Q_j.template topRows<A>() = step.Y;
        Q_j.template bottomRightCorner<U, U>().setIdentity();
        H.noalias() += Q_j.transpose() * H_step * Q_j;

        // Adjoint at the beginning of the substep
        lambda = step.phi_y.transpose() * lambda;
      }
      return Real(0.5)*(H + H.transpose());
    }

  private:
    /*\
     |  Augmented model F_a([x; q], [x'; q'], u, t) = [F(x, x', u, t); q' - L(x, u, t)]
    \*/

    VectorA residual(VectorA const & y, VectorA const & y_dot, VectorU const & u, Real const t) const
    {
      VectorF const x{y.template head<N>()}, x_dot{y_dot.template head<N>()};
      VectorA r;
      r << this->m_problem.system()->F(x, x_dot, u, t), y_dot(Q) - this->m_problem.lagrange(x, u, t);
      return r;
    }

    /** Jacobians of the augmented model with respect to y, y', and u. */
    void jacobians(VectorA const & y, VectorA const & y_dot, VectorU const & u, Real const t,
      MatrixAA & J_y, MatrixAA & J_y_dot, MatrixAU & J_u) const
    {
      auto const & sys{this->m_problem.system()};
      VectorF const x{y.template head<N>()}, x_dot{y_dot.template head<N>()};
      Eigen::Vector<Real, P> const dL{this->m_problem.lagrange_gradient(x, u, t)};

      J_y.setZero();
      J_y.template topLeftCorner<N, N>() = sys->JF_x(x, x_dot, u, t);
      J_y.template block<1, N>(Q, 0)     = -dL.template head<N>().transpose();

      J_y_dot.setZero();
      J_y_dot.template topLeftCorner<N, N>() = sys->JF_x_dot(x, x_dot, u, t);
      J_y_dot(Q, Q) = Real(1.0);

      J_u.template topRows<N>() = sys->JF_u(x, x_dot, u, t);
      J_u.template bottomRows<1>() = -dL.template tail<U>().transpose();
    }

    /** Hessian of the weighted augmented model nu^T F_a in the stage arguments (y, y', u). */
    MatrixW stage_hessian(Substep const & step, Integer const i, VectorA const & nu) const
    {
      VectorA y_i, y_dot_i;
      VectorU u;
      Real t_i;
      this->stage_arguments(step, i, y_i, y_dot_i, u, t_i);
      VectorF const x{y_i.template head<N>()}, x_dot{y_dot_i.template head<N>()};
      VectorF const nu_F{nu.template head<N>()};
      auto const & sys{this->m_problem.system()};

      // Blocks in (x, x', u): the q components enter linearly, so their rows and columns are zero
      constexpr Integer X{0}, X_DOT{A}, CONTROLS{2*A};
      MatrixW W(MatrixW::Zero());
      W.template block<N, N>(X, X)                = sys->HF_x_x(x, x_dot, u, t_i, nu_F);
      W.template block<N, N>(X, X_DOT)            = sys->HF_x_x_dot(x, x_dot, u, t_i, nu_F);
      W.template block<N, N>(X_DOT, X_DOT)        = sys->HF_x_dot_x_dot(x, x_dot, u, t_i, nu_F);
      W.template block<N, U>(X, CONTROLS)         = sys->HF_x_u(x, x_dot, u, t_i, nu_F);
      W.template block<N, U>(X_DOT, CONTROLS)     = sys->HF_x_dot_u(x, x_dot, u, t_i, nu_F);
      W.template block<U, U>(CONTROLS, CONTROLS)  = sys->HF_u_u(x, x_dot, u, t_i, nu_F);

      // Running cost row: -nu_q L(x, u, t)
      Eigen::Matrix<Real, P, P> const HL{this->m_problem.lagrange_hessian(x, u, t_i)};
      W.template block<N, N>(X, X)               -= nu(Q) * HL.template topLeftCorner<N, N>();
      W.template block<N, U>(X, CONTROLS)        -= nu(Q) * HL.template topRightCorner<N, U>();
      W.template block<U, U>(CONTROLS, CONTROLS) -= nu(Q) * HL.template bottomRightCorner<U, U>();

      // Symmetric completion
      W.template block<N, N>(X_DOT, X)        = W.template block<N, N>(X, X_DOT).transpose();
      W.template block<U, N>(CONTROLS, X)     = W.template block<N, U>(X, CONTROLS).transpose();
      W.template block<U, N>(CONTROLS, X_DOT) = W.template block<N, U>(X_DOT, CONTROLS).transpose();
      return W;
    }

    /*\
     |  Stages
    \*/

    /** Arguments (y_i, y'_i, u, t_i) of the i-th stage of a substep. */
    void stage_arguments(Substep const & step, Integer const i, VectorA & y_i, VectorA & y_dot_i,
      VectorU & u, Real & t_i) const
    {
      y_i     = step.y + step.K * this->m_tableau.A.row(i).transpose();
      y_dot_i = step.K.col(i) / this->m_h;
      u       = this->m_p.template tail<U>();
      t_i     = step.t + this->m_tableau.c(i) * this->m_h;
    }

    /** Map from the substep parameters [y; u] to the arguments (y_i, y'_i, u) of the i-th stage. */
    MatrixWL stage_map(Substep const & step, Integer const i) const
    {
      MatrixWL S_i(MatrixWL::Zero());
      S_i.template topLeftCorner<A, A>().setIdentity();
      for (Integer j{0}; j < S; ++j) {
        S_i.template topRows<A>() += this->m_tableau.A(i, j) * step.K_l.template middleRows<A>(j*A);
      }
      S_i.template middleRows<A>(A) = step.K_l.template middleRows<A>(i*A) / this->m_h;
      S_i.template bottomRightCorner<U, U>().setIdentity();
      return S_i;
    }

    /** Stage residuals G(K) and Jacobians dG/dK and dG/d[y; u]. */
    void stage_system(Substep const & step, VectorK & G, MatrixKK & G_K, MatrixKL & G_l) const
    {
      G_K.setZero();
      for (Integer i{0}; i < S; ++i) {
        VectorA y_i, y_dot_i;
        VectorU u_i;
        Real t_i;
        this->stage_arguments(step, i, y_i, y_dot_i, u_i, t_i);

        MatrixAA J_y, J_y_dot;
        MatrixAU J_u;
        this->jacobians(y_i, y_dot_i, u_i, t_i, J_y, J_y_dot, J_u);

        G.template segment<A>(i*A) = this->residual(y_i, y_dot_i, u_i, t_i);
        for (Integer j{0}; j < S; ++j) {
          G_K.template block<A, A>(i*A, j*A) = this->m_tableau.A(i, j) * J_y;
        }
        G_K.template block<A, A>(i*A, i*A) += J_y_dot / this->m_h;
        G_l.template block<A, A>(i*A, 0) = J_y;
        G_l.template block<A, U>(i*A, A) = J_u;
      }
    }

    /**
    * Only stiffly accurate methods (whose last stage is the end of the step, as Radau IIA) integrate
    * systems with algebraic states, i.e., with a singular \f$ \partial\mathbf{F}/\partial\mathbf{x}' \f$
    * (DAEs of index one or higher): the other methods (explicit, Gauss-Legendre, ...) do not give
    * consistent end points, and are allowed for index-0 systems (ODEs, possibly with invariants) only.
    */
    void check_method(VectorF const & x, VectorU const & u) const
    {
      Tableau<Real, S> const & T{this->m_tableau};
      bool const stiffly_accurate{(T.A.row(S-1).transpose() - T.b).template lpNorm<Eigen::Infinity>() <= Real(1.0e-14)};
      if (stiffly_accurate) {return;}
      using MatrixNN = Eigen::Matrix<Real, N, N>;
      MatrixNN const J{this->m_problem.system()->JF_x_dot(x, VectorF::Zero(), u, this->m_t_ini)};
      bool const ode{Eigen::FullPivLU<MatrixNN>(J).isInvertible()};
      SANDALS_ASSERT(ode,
        "Sandals::ShootingInterval: the system has algebraic states (singular dF/dx'), the method " << T.name <<
        " is not stiffly accurate: use a stiffly accurate method (e.g., Radau IIA) for DAEs.");
    }

    /** Solve the stage equations of a substep with Newton's method, and store its derivatives. */
    bool solve_stages(Substep & step)
    {
      VectorK G;
      MatrixKK G_K;
      MatrixKL G_l;
      step.K.setZero();

      // Converged when the residual vanishes, or when the Newton step is negligible with respect to
      // the stage states (the residual alone cannot tell: its rounding floor grows with the magnitude
      // of the terms of the model, and with the index). The refinement step below then takes the
      // stages to machine precision.
      bool converged{false};
      for (Integer iteration{0}; iteration < this->m_max_newton_iterations; ++iteration) {
        this->stage_system(step, G, G_K, G_l);
        if (!G.allFinite()) {return false;}
        if (G.template lpNorm<Eigen::Infinity>() <= Real(1.0e-13) * (Real(1.0) + step.K.template lpNorm<Eigen::Infinity>())) {
          converged = true;
          break;
        }
        step.lu.compute(G_K);
        VectorK const dK{step.lu.solve(-G)};
        if (!dK.allFinite()) {return false;}
        Eigen::Map<VectorK>(step.K.data()) += dK;
        if (this->negligible(step, dK)) {
          this->stage_system(step, G, G_K, G_l); // at the converged stages, for the refinement step
          converged = G.allFinite();
          break;
        }
      }
      if (!converged) {return false;}

      // One more Newton step from the converged stages takes them to machine precision, so that the
      // end point is a smooth function of the parameters (no solver noise in its differences)
      step.lu.compute(G_K);
      Eigen::Map<VectorK>(step.K.data()) += step.lu.solve(-G);
      this->stage_system(step, G, G_K, G_l);
      if (!G.allFinite()) {return false;}

      // Derivatives at the solution
      step.lu.compute(G_K);
      step.K_l = -step.lu.solve(G_l);
      if (!step.K_l.allFinite()) {return false;}
      step.phi_y = MatrixAA::Identity();
      for (Integer i{0}; i < S; ++i) {
        step.phi_y += this->m_tableau.b(i) * step.K_l.template block<A, A>(i*A, 0);
      }
      return true;
    }

    /** True if the Newton step is negligible (1e-10 relative) with respect to every stage state. */
    bool negligible(Substep const & step, VectorK const & dK) const
    {
      Eigen::Map<MatrixK const> const step_K(dK.data());
      for (Integer i{0}; i < S; ++i) {
        auto const scale{Real(1.0) + step.y.array().abs() + step.K.col(i).array().abs()};
        if ((step_K.col(i).array().abs() > Real(1.0e-10) * scale).any()) {return false;}
      }
      return true;
    }

    /** Sensitivity of the substep end point to the controls, phi_u = sum_i b_i dK_i/du. */
    MatrixAU substep_jacobian_u(Substep const & step) const
    {
      MatrixAU phi_u(MatrixAU::Zero());
      for (Integer i{0}; i < S; ++i) {
        phi_u += this->m_tableau.b(i) * step.K_l.template block<A, U>(i*A, A);
      }
      return phi_u;
    }
  };

} // namespace Sandals

#endif // SANDALS_SHOOTING_INTERVAL_HH
