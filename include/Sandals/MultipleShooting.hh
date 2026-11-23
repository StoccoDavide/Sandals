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

#ifndef SANDALS_MULTIPLE_SHOOTING_HH
#define SANDALS_MULTIPLE_SHOOTING_HH

#include <Sandals.hh>
#include <Sandals/ShootingInterval.hh>

namespace Sandals
{

  /*\
   |   __  __       _ _   _       _      ____  _                 _   _
   |  |  \/  |_   _| | |_(_)_ __ | | ___/ ___|| |__   ___   ___ | |_(_)_ __   __ _
   |  | |\/| | | | | | __| | '_ \| |/ _ \___ \| '_ \ / _ \ / _ \| __| | '_ \ / _` |
   |  | |  | | |_| | | |_| | |_) | |  __/___) | | | | (_) | (_) | |_| | | | | (_| |
   |  |_|  |_|\__,_|_|\__|_| .__/|_|\___|____/|_| |_|\___/ \___/ \__|_|_| |_|\__, |
   |                       |_|                                               |___/
  \*/

  /**
  * \brief Direct multiple shooting transcription of an optimal control problem of a DAE with invariants.
  *
  * On the fixed mesh \f$ t_0 < \dots < t_m \f$, the variables are the states at the nodes
  * \f$ \mathbf{x}_k \f$ and the piecewise-constant controls \f$ \mathbf{u}_k \f$. Each interval is
  * integrated by a Runge-Kutta method with a fixed number of substeps (more where the dynamics is
  * stiff), giving the end point \f$ \mathbf{X}_{k+1} \f$ and the running cost integral \f$ Q_k \f$,
  * with exact first and second derivatives (see Sandals::ShootingInterval).
  *
  * <b>Invariants.</b> At every node the rectangular system \f$ [\mathbf{x}_k - \mathbf{X}_k = \mathbf{0};
  * \mathbf{h}(\mathbf{x}_k) = \mathbf{0}] \f$ (\f$ N + M \f$ rows, \f$ N \f$ unknowns) is solved in the
  * least-squares sense,
  *
  * \f[
  * \min_{\mathbf{x}_k} \frac{1}{2}\|\mathbf{x}_k - \mathbf{X}_k\|^2 + \frac{1}{2\varepsilon}\|\mathbf{h}(
  * \mathbf{x}_k)\|^2 \text{,}
  * \f]
  *
  * through its optimality conditions in augmented form (the stable alternative to the normal
  * equations, as a QR factorization would solve it), with the multipliers \f$ \boldsymbol{\mu}_k \f$
  * as variables:
  *
  * \f[
  * \mathbf{x}_k - \mathbf{X}_k + \mathbf{Jh}_{\mathbf{x}}(\mathbf{x}_k)^\top\boldsymbol{\mu}_k =
  * \mathbf{0}, \quad \mathbf{h}(\mathbf{x}_k) - \varepsilon\boldsymbol{\mu}_k = \mathbf{0} \text{.}
  * \f]
  *
  * The weight \f$ \varepsilon \geq 0 \f$ balances continuity and invariants: \f$ \varepsilon = 0 \f$
  * makes the invariants exact (the minimum-norm correction of the continuity), \f$ \varepsilon = 1 \f$
  * weights them equally. The system stays square, so the constraints are independent even when the
  * integrator preserves the invariants. The incoming state of the first node, \f$ \mathbf{X}_0 \f$, is a
  * variable on which the boundary conditions act, so that slightly inconsistent boundary data is
  * absorbed by the same least-squares correction. Without invariants (\f$ M = 0 \f$) the conditions
  * reduce to the continuity \f$ \mathbf{x}_k = \mathbf{X}_k \f$, and the boundary conditions act on
  * \f$ \mathbf{x}_0 \f$.
  *
  * Layout of the variables and of the constraints:
  * - variables: \f$ [\mathbf{x}_0, \mathbf{u}_0, \dots, \mathbf{x}_{m-1}, \mathbf{u}_{m-1}, \mathbf{x}_m,
  *   \mathbf{X}_0, \boldsymbol{\mu}_0, \dots, \boldsymbol{\mu}_m] \f$ (the last ones only if \f$ M > 0 \f$);
  * - constraints: continuity at the nodes \f$ 1, \dots, m \f$, continuity at the node \f$ 0 \f$ (if
  *   \f$ M > 0 \f$), invariants at the nodes \f$ 0, \dots, m \f$, boundary conditions on \f$
  *   (\mathbf{X}_0, \mathbf{x}_m) \f$, and path constraints at the nodes (at the last node with the
  *   controls of the last interval).
  *
  * \note DAEs: the node values of the algebraic states only shift the stage variables of stiffly
  * accurate methods, so their initial value must be fixed by the boundary conditions.
  * \tparam Real The scalar number type.
  * \tparam S The number of stages of the Runge-Kutta method.
  * \tparam N The number of states.
  * \tparam U The number of controls.
  * \tparam M The number of invariants.
  */
  template <typename Real, Integer S, Integer N, Integer U, Integer M = 0>
  class MultipleShooting
  {
  public:
    using Problem    = OptimalControlProblem<Real, N, U, M>; /**< Optimal control problem type. */
    using ProblemPtr = std::shared_ptr<Problem>;            /**< Optimal control problem pointer type. */
    using Interval   = ShootingInterval<Real, S, N, U, M>;  /**< Shooting interval type. */
    using VectorF    = Eigen::Vector<Real, N>;              /**< States vector type. */
    using VectorU    = Eigen::Vector<Real, U>;              /**< Controls vector type. */
    using VectorH    = Eigen::Vector<Real, M>;              /**< Invariants vector type. */
    using MatrixF    = Eigen::Matrix<Real, N, N>;           /**< States matrix type. */
    using VectorX    = Eigen::Vector<Real, Eigen::Dynamic>; /**< Dynamic vector type. */
    using MatrixX    = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; /**< Dynamic matrix type. */
    using MatrixXN   = Eigen::Matrix<Real, N, Eigen::Dynamic>; /**< States on the nodes type. */
    using MatrixXU   = Eigen::Matrix<Real, U, Eigen::Dynamic>; /**< Controls on the intervals type. */
    using FunctionX  = std::function<VectorF(Real const)>;  /**< States guess function type. */
    using FunctionU  = std::function<VectorU(Real const)>;  /**< Controls guess function type. */

    static constexpr Integer P{N + U}; /**< Variables [x_k; u_k] of an interval. */

  private:
    using Index = std::pair<Integer, Integer>;

    ProblemPtr            m_problem;   /**< Optimal control problem. */
    VectorX               m_t_mesh;    /**< Shooting mesh. */
    Integer               m_m{0};      /**< Number of intervals. */
    Integer               m_nb{0};     /**< Number of boundary conditions. */
    Integer               m_ng{0};     /**< Number of path constraints. */
    Real                  m_epsilon;   /**< Weight of the invariants in the node least squares. */
    mutable std::vector<Interval> m_intervals; /**< Shooting intervals (with their cached integration). */
    VectorX               m_guess;     /**< Starting point. */

    std::vector<Integer>     m_jac_rows, m_jac_cols;   /**< Jacobian sparsity pattern. */
    std::vector<Integer>     m_hess_rows, m_hess_cols; /**< Hessian (lower triangle) sparsity pattern. */
    std::map<Index, Integer> m_hess_map;               /**< Hessian (row, column) to nonzero index. */

  public:
    /**
    * Class constructor.
    * \param[in] t_problem The optimal control problem.
    * \param[in] t_tableau The Butcher tableau of the Runge-Kutta method.
    * \param[in] t_t_mesh The shooting mesh, from the initial to the final time of the problem.
    * \param[in] t_substeps The number of substeps of each interval.
    * \param[in] t_epsilon The weight \f$ \varepsilon \geq 0 \f$ of the invariants (0: exact invariants).
    */
    MultipleShooting(ProblemPtr t_problem, Tableau<Real, S> const & t_tableau, VectorX const & t_t_mesh,
      std::vector<Integer> const & t_substeps, Real const t_epsilon = 0.0)
      : m_problem(t_problem), m_t_mesh(t_t_mesh), m_epsilon(t_epsilon)
    {
      #define CMD "Sandals::MultipleShooting::MultipleShooting(...): "
      this->m_m = static_cast<Integer>(t_t_mesh.size()) - 1;
      SANDALS_ASSERT(this->m_m >= 1, CMD "at least two mesh points are required.");
      SANDALS_ASSERT(static_cast<Integer>(t_substeps.size()) == this->m_m, CMD "one number of substeps per interval is required.");
      SANDALS_ASSERT(t_epsilon >= 0.0, CMD "the weight of the invariants must be non-negative.");
      SANDALS_ASSERT(std::abs(t_t_mesh(0) - t_problem->t_ini()) < 1.0e-12 &&
        std::abs(t_t_mesh(this->m_m) - t_problem->t_end()) < 1.0e-12, CMD "the mesh must span the problem time interval.");
      #undef CMD
      this->m_nb = t_problem->boundary_conditions_number();
      this->m_ng = t_problem->path_constraints_number();
      for (Integer k{0}; k < this->m_m; ++k) {
        this->m_intervals.emplace_back(*t_problem, t_tableau, t_t_mesh(k), t_t_mesh(k+1), t_substeps[k]);
      }
      this->build_jacobian_structure();
      this->build_hessian_structure();
      this->m_guess = VectorX::Zero(this->variables_number());
    }

    /**
    * Class constructor with the same number of substeps on every interval.
    * \param[in] t_problem The optimal control problem.
    * \param[in] t_tableau The Butcher tableau of the Runge-Kutta method.
    * \param[in] t_intervals The number of equal intervals.
    * \param[in] t_substeps The number of substeps of each interval.
    * \param[in] t_epsilon The weight \f$ \varepsilon \geq 0 \f$ of the invariants (0: exact invariants).
    */
    MultipleShooting(ProblemPtr t_problem, Tableau<Real, S> const & t_tableau, Integer const t_intervals,
      Integer const t_substeps = 1, Real const t_epsilon = 0.0)
      : MultipleShooting(t_problem, t_tableau, VectorX::LinSpaced(t_intervals + 1, t_problem->t_ini(), t_problem->t_end()),
        std::vector<Integer>(t_intervals, t_substeps), t_epsilon)
    {}

    /*\
     |  Layout
    \*/

    /** \return The number of intervals \f$ m \f$. */
    Integer intervals() const {return this->m_m;}

    /** \return The shooting mesh. */
    VectorX const & t_mesh() const {return this->m_t_mesh;}

    /** \return The weight \f$ \varepsilon \f$ of the invariants. */
    Real epsilon() const {return this->m_epsilon;}

    /** \return The number of NLP variables. */
    Integer variables_number() const
    {
      Integer const nodes{this->m_m*P + N};
      return M > 0 ? nodes + N + (this->m_m + 1)*M : nodes;
    }

    /** \return The number of NLP constraints. */
    Integer constraints_number() const {return this->path_row(this->m_m + 1);}

    /** \param[in] k Node index. \return The index of \f$ \mathbf{x}_k \f$. */
    Integer x_index(Integer const k) const {return k*P;}

    /** \param[in] k Interval index. \return The index of \f$ \mathbf{u}_k \f$. */
    Integer u_index(Integer const k) const {return k*P + N;}

    /** \return The index of the incoming state of the first node \f$ \mathbf{X}_0 \f$ (if \f$ M > 0 \f$). */
    Integer X0_index() const {return this->m_m*P + N;}

    /** \param[in] k Node index. \return The index of \f$ \boldsymbol{\mu}_k \f$ (if \f$ M > 0 \f$). */
    Integer mu_index(Integer const k) const {return this->X0_index() + N + k*M;}

    /** \param[in] k Node index (\f$ k \geq 1 \f$, or \f$ k = 0 \f$ if \f$ M > 0 \f$). \return The first continuity row of the node. */
    Integer continuity_row(Integer const k) const {return k > 0 ? (k - 1)*N : this->m_m*N;}

    /** \param[in] k Node index. \return The first invariants row of the node (if \f$ M > 0 \f$). */
    Integer invariants_row(Integer const k) const {return this->m_m*N + (M > 0 ? N : 0) + k*M;}

    /** \return The first boundary conditions row. */
    Integer boundary_row() const {return this->invariants_row(this->m_m + 1);}

    /** \param[in] k Node index. \return The first path constraints row of the node. */
    Integer path_row(Integer const k) const {return this->boundary_row() + this->m_nb + k*this->m_ng;}

    /** \return The index of the states on which the initial boundary conditions act. */
    Integer boundary_start_index() const {return M > 0 ? this->X0_index() : this->x_index(0);}

    /** \param[in] z NLP variables. \param[in] k Node index. \return The states \f$ \mathbf{x}_k \f$. */
    VectorF x(VectorX const & z, Integer const k) const {return z.template segment<N>(this->x_index(k));}

    /** \param[in] z NLP variables. \param[in] k Interval index. \return The controls \f$ \mathbf{u}_k \f$. */
    VectorU u(VectorX const & z, Integer const k) const {return z.template segment<U>(this->u_index(k));}

    /** \param[in] z NLP variables. \param[in] k Node index. \return The multipliers \f$ \boldsymbol{\mu}_k \f$. */
    VectorH mu(VectorX const & z, Integer const k) const {return z.template segment<M>(this->mu_index(k));}

    /** \param[in] z NLP variables. \return The states on the nodes. */
    MatrixXN x_mesh(VectorX const & z) const
    {
      MatrixXN X(N, this->m_m + 1);
      for (Integer k{0}; k <= this->m_m; ++k) {X.col(k) = this->x(z, k);}
      return X;
    }

    /** \param[in] z NLP variables. \return The controls on the intervals. */
    MatrixXU u_mesh(VectorX const & z) const
    {
      MatrixXU Uk(U, this->m_m);
      for (Integer k{0}; k < this->m_m; ++k) {Uk.col(k) = this->u(z, k);}
      return Uk;
    }

    /*\
     |  Initial guess
    \*/

    /** \return The starting point. */
    VectorX const & guess() const {return this->m_guess;}

    /** \param[in] t_guess The starting point. */
    void guess(VectorX const & t_guess)
    {
      SANDALS_ASSERT(t_guess.size() == this->variables_number(), "Sandals::MultipleShooting::guess(...): wrong size.");
      this->m_guess = t_guess;
    }

    /**
    * Set the starting point from states and controls functions of time (the multipliers are zero).
    * \param[in] x_guess The states guess \f$ \mathbf{x}(t) \f$.
    * \param[in] u_guess The controls guess \f$ \mathbf{u}(t) \f$ (at the middle of the intervals).
    */
    void guess(FunctionX const & x_guess, FunctionU const & u_guess)
    {
      VectorX & z{this->m_guess};
      z.setZero();
      for (Integer k{0}; k <= this->m_m; ++k) {z.template segment<N>(this->x_index(k)) = x_guess(this->m_t_mesh(k));}
      for (Integer k{0}; k < this->m_m; ++k) {z.template segment<U>(this->u_index(k)) = u_guess(this->middle(k));}
      if constexpr (M > 0) {z.template segment<N>(this->X0_index()) = this->x(z, 0);}
    }

    /**
    * Set the starting point by integrating the system forward from the given initial states, with the
    * given controls. Each node is the least-squares correction of the incoming integrated state, so the
    * starting point satisfies the continuity and invariants constraints.
    * \param[in] x_ini The initial states (the incoming state of the first node).
    * \param[in] u_guess The controls guess \f$ \mathbf{u}(t) \f$ (at the middle of the intervals).
    * \return True if the integration and the corrections succeeded, false otherwise.
    */
    bool simulate_guess(VectorF const & x_ini, FunctionU const & u_guess)
    {
      VectorX & z{this->m_guess};
      z.setZero();
      VectorF X_k{x_ini};
      if constexpr (M > 0) {z.template segment<N>(this->X0_index()) = x_ini;}
      for (Integer k{0}; k <= this->m_m; ++k) {
        VectorF x_k;
        VectorH mu_k;
        if (!this->node_least_squares(X_k, this->m_t_mesh(k), x_k, mu_k)) {return false;}
        z.template segment<N>(this->x_index(k)) = x_k;
        if constexpr (M > 0) {z.template segment<M>(this->mu_index(k)) = mu_k;}
        if (k == this->m_m) {break;}

        VectorU const u_k{u_guess(this->middle(k))};
        z.template segment<U>(this->u_index(k)) = u_k;
        if (!this->m_intervals[k].integrate(x_k, u_k)) {return false;}
        X_k = this->m_intervals[k].end().template head<N>();
      }
      return true;
    }

    /*\
     |  Nonlinear program
    \*/

    /** Bounds of the variables and of the constraints (infinite bounds as infinity). */
    void bounds(VectorX & z_lower, VectorX & z_upper, VectorX & c_lower, VectorX & c_upper) const
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      z_lower = VectorX::Constant(this->variables_number(), -inf);
      z_upper = VectorX::Constant(this->variables_number(), inf);
      for (Integer k{0}; k <= this->m_m; ++k) {
        z_lower.template segment<N>(this->x_index(k)) = this->m_problem->x_lower();
        z_upper.template segment<N>(this->x_index(k)) = this->m_problem->x_upper();
      }
      for (Integer k{0}; k < this->m_m; ++k) {
        z_lower.template segment<U>(this->u_index(k)) = this->m_problem->u_lower();
        z_upper.template segment<U>(this->u_index(k)) = this->m_problem->u_upper();
      }
      c_lower = VectorX::Zero(this->constraints_number());
      c_upper = VectorX::Zero(this->constraints_number());
      if (this->m_nb > 0) {
        c_lower.segment(this->boundary_row(), this->m_nb) = this->m_problem->boundary_conditions_lower();
        c_upper.segment(this->boundary_row(), this->m_nb) = this->m_problem->boundary_conditions_upper();
      }
      for (Integer k{0}; k <= this->m_m && this->m_ng > 0; ++k) {
        c_lower.segment(this->path_row(k), this->m_ng) = this->m_problem->path_constraints_lower();
        c_upper.segment(this->path_row(k), this->m_ng) = this->m_problem->path_constraints_upper();
      }
    }

    /** \param[in] z NLP variables. \return The objective: Mayer term plus the integrated running cost. */
    Real objective(VectorX const & z) const
    {
      Real cost{this->m_problem->mayer(this->x(z, 0), this->x(z, this->m_m))};
      for (Integer k{0}; k < this->m_m; ++k) {cost += this->interval(z, k).end()(Interval::Q);}
      return cost;
    }

    /** \param[in] z NLP variables. \return The gradient of the objective. */
    VectorX objective_gradient(VectorX const & z) const
    {
      VectorX gradient(VectorX::Zero(this->variables_number()));
      Eigen::Vector<Real, 2*N> const g_mayer{this->m_problem->mayer_gradient(this->x(z, 0), this->x(z, this->m_m))};
      gradient.template segment<N>(this->x_index(0)) += g_mayer.template head<N>();
      gradient.template segment<N>(this->x_index(this->m_m)) += g_mayer.template tail<N>();
      for (Integer k{0}; k < this->m_m; ++k) {
        gradient.template segment<P>(this->x_index(k)) += this->interval(z, k).sensitivity().row(Interval::Q).transpose();
      }
      return gradient;
    }

    /** \param[in] z NLP variables. \return The constraints. */
    VectorX constraints(VectorX const & z) const
    {
      VectorX c(this->constraints_number());
      for (Integer k{1}; k <= this->m_m; ++k) {
        VectorF const X_k{this->interval(z, k-1).end().template head<N>()};
        c.template segment<N>(this->continuity_row(k)) = this->continuity(z, k, X_k);
      }
      if constexpr (M > 0) {
        c.template segment<N>(this->continuity_row(0)) = this->continuity(z, 0, z.template segment<N>(this->X0_index()));
        for (Integer k{0}; k <= this->m_m; ++k) {
          c.template segment<M>(this->invariants_row(k)) = this->system().h(this->x(z, k), this->m_t_mesh(k)) - this->m_epsilon*this->mu(z, k);
        }
      }
      if (this->m_nb > 0) {
        c.segment(this->boundary_row(), this->m_nb) = this->m_problem->boundary_conditions(this->boundary_start(z), this->x(z, this->m_m));
      }
      for (Integer k{0}; k <= this->m_m && this->m_ng > 0; ++k) {
        c.segment(this->path_row(k), this->m_ng) = this->m_problem->path_constraints(this->x(z, k), this->u(z, this->path_interval(k)), this->m_t_mesh(k));
      }
      return c;
    }

    /** \return The number of nonzeros of the constraints Jacobian. */
    Integer jacobian_nnz() const {return static_cast<Integer>(this->m_jac_rows.size());}

    /** \return The Jacobian sparsity rows. */
    std::vector<Integer> const & jacobian_rows() const {return this->m_jac_rows;}

    /** \return The Jacobian sparsity columns. */
    std::vector<Integer> const & jacobian_cols() const {return this->m_jac_cols;}

    /** \param[in] z NLP variables. \return The Jacobian nonzeros, in the order of the sparsity pattern. */
    VectorX jacobian_values(VectorX const & z) const
    {
      std::vector<Integer> rows, cols;
      std::vector<Real> values;
      this->template jacobian_traversal<true>(&z, rows, cols, values);
      return Eigen::Map<VectorX>(values.data(), values.size());
    }

    /** \return The number of nonzeros of the lower triangle of the Lagrangian Hessian. */
    Integer hessian_nnz() const {return static_cast<Integer>(this->m_hess_rows.size());}

    /** \return The Hessian sparsity rows. */
    std::vector<Integer> const & hessian_rows() const {return this->m_hess_rows;}

    /** \return The Hessian sparsity columns. */
    std::vector<Integer> const & hessian_cols() const {return this->m_hess_cols;}

    /**
    * Evaluate the lower triangle of the Hessian of the Lagrangian \f$ \sigma f + \boldsymbol{\lambda}^\top
    * \mathbf{c} \f$, in the order of the sparsity pattern.
    * \param[in] z NLP variables.
    * \param[in] sigma The objective factor.
    * \param[in] lambda The constraints multipliers.
    * \return The Hessian nonzeros.
    */
    VectorX hessian_values(VectorX const & z, Real const sigma, VectorX const & lambda) const
    {
      VectorX values(VectorX::Zero(this->hessian_nnz()));

      // Intervals: the end point X_{k+1} enters the continuity at the node k+1 with a minus sign, and
      // the running cost integral enters the objective
      for (Integer k{0}; k < this->m_m; ++k) {
        typename Interval::VectorA w;
        w << -lambda.template segment<N>(this->continuity_row(k+1)), sigma;
        this->accumulate(this->interval_indices(k), this->interval(z, k).hessian(w), values);
      }

      // Node least squares: x_k - X_k + Jh(x_k)^T mu_k and h(x_k) - epsilon mu_k
      if constexpr (M > 0) {
        for (Integer k{0}; k <= this->m_m; ++k) {
          this->accumulate(this->node_indices(k), this->node_hessian(z, k, lambda), values);
        }
      }

      // Mayer term on (x_0, x_m), boundary conditions on (X_0, x_m) (or (x_0, x_m) without invariants)
      this->accumulate(this->pair_indices(this->x_index(0)), sigma*this->m_problem->mayer_hessian(this->x(z, 0), this->x(z, this->m_m)), values);
      if (this->m_nb > 0) {
        VectorX const lambda_b{lambda.segment(this->boundary_row(), this->m_nb)};
        this->accumulate(this->pair_indices(this->boundary_start_index()),
          this->m_problem->boundary_conditions_hessian(this->boundary_start(z), this->x(z, this->m_m), lambda_b), values);
      }

      // Path constraints at the nodes
      for (Integer k{0}; k <= this->m_m && this->m_ng > 0; ++k) {
        VectorX const lambda_g{lambda.segment(this->path_row(k), this->m_ng)};
        this->accumulate(this->path_indices(k), this->m_problem->path_constraints_hessian(
          this->x(z, k), this->u(z, this->path_interval(k)), this->m_t_mesh(k), lambda_g), values);
      }
      return values;
    }

  private:
    /*\
     |  Helpers
    \*/

    auto const & system() const {return *this->m_problem->system();}

    Real middle(Integer const k) const {return (this->m_t_mesh(k) + this->m_t_mesh(k+1)) / Real(2.0);}

    /** The interval at which the path constraints of the node k are evaluated (the last one at the end). */
    Integer path_interval(Integer const k) const {return std::min(k, this->m_m - 1);}

    VectorF boundary_start(VectorX const & z) const {return z.template segment<N>(this->boundary_start_index());}

    /** The interval k, integrated from the current nodes and controls (an evaluation error otherwise). */
    Interval const & interval(VectorX const & z, Integer const k) const
    {
      if (!this->m_intervals[k].integrate(this->x(z, k), this->u(z, k))) {
        SANDALS_ERROR("Sandals::MultipleShooting: integration failed on interval " << k << ".");
      }
      return this->m_intervals[k];
    }

    /** Continuity (least-squares stationarity) of the node k: x_k - X_k + Jh(x_k)^T mu_k. */
    VectorF continuity(VectorX const & z, Integer const k, VectorF const & X_k) const
    {
      VectorF c{this->x(z, k) - X_k};
      if constexpr (M > 0) {c += this->system().Jh_x(this->x(z, k), this->m_t_mesh(k)).transpose() * this->mu(z, k);}
      return c;
    }

    /**
    * Solve the node least squares for given incoming states X, with Newton's method on its augmented
    * optimality conditions [x - X + Jh^T mu; h - epsilon mu] = 0 (used to build the starting point).
    */
    bool node_least_squares(VectorF const & X, Real const t, VectorF & x, VectorH & mu) const
    {
      x = X;
      mu.setZero();
      if constexpr (M > 0) {
        using VectorL = Eigen::Vector<Real, N+M>;
        using MatrixL = Eigen::Matrix<Real, N+M, N+M>;
        Real const tolerance{Real(1.0e-14) * std::max(Real(1.0), X.norm())};
        for (Integer iteration{0}; iteration < 50; ++iteration) {
          Eigen::Matrix<Real, M, N> const Jh{this->system().Jh_x(x, t)};
          VectorL residual;
          residual << x - X + Jh.transpose()*mu, this->system().h(x, t) - this->m_epsilon*mu;
          if (residual.template lpNorm<Eigen::Infinity>() <= tolerance) {return true;}
          MatrixL J;
          J << MatrixF::Identity() + this->system().Hh_x_x(x, t, mu), Jh.transpose(),
               Jh, -this->m_epsilon*Eigen::Matrix<Real, M, M>::Identity();
          VectorL const step{J.fullPivLu().solve(-residual)};
          if (!step.allFinite()) {return false;}
          x  += step.template head<N>();
          mu += step.template tail<M>();
        }
        return false;
      }
      return true;
    }

    /**
    * Hessian of the node least-squares terms lambda_c^T Jh(x_k)^T mu_k + lambda_h^T h(x_k) in
    * (x_k, mu_k), with lambda_c and lambda_h the multipliers of the continuity and invariants rows:
    * - xx: Th(mu, lambda_c) + Hh(lambda_h), the third derivatives of the invariants and their Hessian;
    * - x mu_j: Hh(e_j) lambda_c, the Hessian of the j-th invariant times lambda_c.
    */
    Eigen::Matrix<Real, N+M, N+M> node_hessian(VectorX const & z, Integer const k, VectorX const & lambda) const
    {
      VectorF const x_k{this->x(z, k)};
      Real const t_k{this->m_t_mesh(k)};
      VectorF const lambda_c{lambda.template segment<N>(this->continuity_row(k))};
      VectorH const lambda_h{lambda.template segment<M>(this->invariants_row(k))};

      Eigen::Matrix<Real, N+M, N+M> H(Eigen::Matrix<Real, N+M, N+M>::Zero());
      H.template topLeftCorner<N, N>() = this->system().Th_x_x(x_k, t_k, this->mu(z, k), lambda_c)
        + this->system().Hh_x_x(x_k, t_k, lambda_h);
      for (Integer j{0}; j < M; ++j) {
        VectorF const column{this->system().Hh_x_x(x_k, t_k, VectorH::Unit(j)) * lambda_c};
        H.template block<N, 1>(0, N + j) = column;
        H.template block<1, N>(N + j, 0) = column.transpose();
      }
      return H;
    }

    /*\
     |  Sparsity patterns
    \*/

    /** Indices of the variables [x_k; u_k] of the interval k. */
    std::vector<Integer> interval_indices(Integer const k) const {return this->range(this->x_index(k), P);}

    /** Indices of the variables [x_k; mu_k] of the node k. */
    std::vector<Integer> node_indices(Integer const k) const
    {
      std::vector<Integer> idx{this->range(this->x_index(k), N)};
      for (Integer a{0}; a < M; ++a) {idx.push_back(this->mu_index(k) + a);}
      return idx;
    }

    /** Indices of the states [x_start; x_m], with x_start the states at the given index. */
    std::vector<Integer> pair_indices(Integer const start) const
    {
      std::vector<Integer> idx{this->range(start, N)};
      for (Integer a{0}; a < N; ++a) {idx.push_back(this->x_index(this->m_m) + a);}
      return idx;
    }

    /** Indices of the arguments [x_k; u] of the path constraints of the node k. */
    std::vector<Integer> path_indices(Integer const k) const
    {
      std::vector<Integer> idx{this->range(this->x_index(k), N)};
      for (Integer a{0}; a < U; ++a) {idx.push_back(this->u_index(this->path_interval(k)) + a);}
      return idx;
    }

    static std::vector<Integer> range(Integer const first, Integer const size)
    {
      std::vector<Integer> idx(size);
      for (Integer a{0}; a < size; ++a) {idx[a] = first + a;}
      return idx;
    }

    /** Traverse the constraints Jacobian in a fixed order (pattern only, or pattern and values). */
    template <bool Values>
    void jacobian_traversal(VectorX const * z, std::vector<Integer> & rows, std::vector<Integer> & cols,
      std::vector<Real> & values) const
    {
      auto push_block = [&] (Integer const row, std::vector<Integer> const & columns, MatrixX const & block) {
        for (Integer r{0}; r < block.rows(); ++r) {
          for (Integer c{0}; c < block.cols(); ++c) {
            rows.push_back(row + r);
            cols.push_back(columns[c]);
            if constexpr (Values) {values.push_back(block(r, c));}
          }
        }
      };
      auto value = [&] (auto const & compute, Integer const r, Integer const c) {
        if constexpr (Values) {return MatrixX(compute());} else {(void) compute; return MatrixX(MatrixX::Zero(r, c));}
      };

      // Continuity at the nodes 1..m: -dX_k/d[x_{k-1}; u_{k-1}], and d/dx_k, d/dmu_k of the correction
      for (Integer k{1}; k <= this->m_m; ++k) {
        push_block(this->continuity_row(k), this->interval_indices(k-1),
          value([&] {return MatrixX(-this->interval(*z, k-1).sensitivity().template topRows<N>());}, N, P));
        this->push_node_correction(push_block, value, z, k);
      }

      // Continuity at the node 0 (incoming state X_0) and invariants at every node
      if constexpr (M > 0) {
        push_block(this->continuity_row(0), this->range(this->X0_index(), N), value([] {return MatrixX(-MatrixF::Identity());}, N, N));
        this->push_node_correction(push_block, value, z, 0);
        for (Integer k{0}; k <= this->m_m; ++k) {
          push_block(this->invariants_row(k), this->range(this->x_index(k), N),
            value([&] {return MatrixX(this->system().Jh_x(this->x(*z, k), this->m_t_mesh(k)));}, M, N));
          push_block(this->invariants_row(k), this->range(this->mu_index(k), M),
            value([&] {return MatrixX(-this->m_epsilon*MatrixX::Identity(M, M));}, M, M));
        }
      }

      // Boundary conditions on (x_start, x_m)
      if (this->m_nb > 0) {
        push_block(this->boundary_row(), this->pair_indices(this->boundary_start_index()),
          value([&] {return this->m_problem->boundary_conditions_jacobian(this->boundary_start(*z), this->x(*z, this->m_m));}, this->m_nb, 2*N));
      }

      // Path constraints at the nodes
      for (Integer k{0}; k <= this->m_m && this->m_ng > 0; ++k) {
        push_block(this->path_row(k), this->path_indices(k), value([&] {return this->m_problem->path_constraints_jacobian(
          this->x(*z, k), this->u(*z, this->path_interval(k)), this->m_t_mesh(k));}, this->m_ng, N+U));
      }
    }

    /** Derivatives of x_k + Jh(x_k)^T mu_k with respect to x_k and mu_k (continuity rows of the node k). */
    template <typename PushBlock, typename Value>
    void push_node_correction(PushBlock & push_block, Value & value, VectorX const * z, Integer const k) const
    {
      push_block(this->continuity_row(k), this->range(this->x_index(k), N), value([&] {
        MatrixF d{MatrixF::Identity()};
        if constexpr (M > 0) {d += this->system().Hh_x_x(this->x(*z, k), this->m_t_mesh(k), this->mu(*z, k));}
        return MatrixX(d);
      }, N, N));
      if constexpr (M > 0) {
        push_block(this->continuity_row(k), this->range(this->mu_index(k), M),
          value([&] {return MatrixX(this->system().Jh_x(this->x(*z, k), this->m_t_mesh(k)).transpose());}, N, M));
      }
    }

    void build_jacobian_structure()
    {
      std::vector<Real> unused;
      this->template jacobian_traversal<false>(nullptr, this->m_jac_rows, this->m_jac_cols, unused);
    }

    /** Add the lower triangle of a dense block over the given (distinct) indices to the Hessian pattern. */
    void add_hessian_block(std::vector<Integer> const & idx)
    {
      for (Integer a : idx) {
        for (Integer b : idx) {
          if (a < b || this->m_hess_map.count({a, b})) {continue;}
          this->m_hess_map.emplace(Index{a, b}, static_cast<Integer>(this->m_hess_rows.size()));
          this->m_hess_rows.push_back(a);
          this->m_hess_cols.push_back(b);
        }
      }
    }

    void build_hessian_structure()
    {
      for (Integer k{0}; k < this->m_m; ++k) {this->add_hessian_block(this->interval_indices(k));}
      if constexpr (M > 0) {
        for (Integer k{0}; k <= this->m_m; ++k) {this->add_hessian_block(this->node_indices(k));}
      }
      this->add_hessian_block(this->pair_indices(this->x_index(0)));
      if (this->m_nb > 0) {this->add_hessian_block(this->pair_indices(this->boundary_start_index()));}
      for (Integer k{0}; k <= this->m_m && this->m_ng > 0; ++k) {this->add_hessian_block(this->path_indices(k));}
    }

    /** Accumulate the lower triangle of a dense symmetric block over the given indices. */
    template <typename Matrix>
    void accumulate(std::vector<Integer> const & idx, Matrix const & H, VectorX & values) const
    {
      Integer const size{static_cast<Integer>(idx.size())};
      for (Integer a{0}; a < size; ++a) {
        for (Integer b{0}; b < size; ++b) {
          if (idx[a] >= idx[b]) {values(this->m_hess_map.at(Index{idx[a], idx[b]})) += H(a, b);}
        }
      }
    }
  };

} // namespace Sandals

#endif // SANDALS_MULTIPLE_SHOOTING_HH
