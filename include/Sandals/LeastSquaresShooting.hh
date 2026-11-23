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

#ifndef SANDALS_LEAST_SQUARES_SHOOTING_HH
#define SANDALS_LEAST_SQUARES_SHOOTING_HH

#include <Sandals.hh>
#include <Sandals/Problem.hh>

namespace Sandals
{

  /*\
   |   _                   _   ____                                  ____  _                 _   _
   |  | |    ___  __ _ ___| |_/ ___|  __ _ _   _  __ _ _ __ ___  ___/ ___|| |__   ___   ___ | |_(_)_ __   __ _
   |  | |   / _ \/ _` / __| __\___ \ / _` | | | |/ _` | '__/ _ \/ __\___ \| '_ \ / _ \ / _ \| __| | '_ \ / _` |
   |  | |__|  __/ (_| \__ \ |_ ___) | (_| | |_| | (_| | | |  __/\__ \___) | | | | (_) | (_) | |_| | | | | (_| |
   |  |_____\___|\__,_|___/\__|____/ \__, |\__,_|\__,_|_|  \___||___/____/|_| |_|\___/ \___/ \__|_|_| |_|\__, |
   |                                    |_|                                                              |___/
  \*/

  /**
  * \brief Least-squares multiple shooting for ODE/DAE boundary value problems, as a nonlinear program.
  *
  * Given the shooting nodes \f$ \mathbf{x}_0, \dots, \mathbf{x}_m \f$ on the mesh \f$ t_0 < \dots < t_m
  * \f$, each interval is integrated with \f$ s \f$ fixed substeps of the Runge-Kutta method of the
  * problem, \f$ \mathbf{X}_k = \boldsymbol{\Phi}(\mathbf{x}_k) \f$, and the (possibly rectangular)
  * residual system is
  *
  * \f[
  * \mathbf{g}(\mathbf{x}) = \begin{bmatrix}
  *   \boldsymbol{\Phi}(\mathbf{x}_k) - \mathbf{x}_{k+1} & k = 0, \dots, m-1 \\
  *   \mathbf{h}(\mathbf{x}_k, t_k) & k = 0, \dots, m \\
  *   \mathbf{b}(\mathbf{x}_0, \mathbf{x}_m)
  * \end{bmatrix} \text{.}
  * \f]
  *
  * The invariants \f$ \mathbf{h} \f$ are handled in one of two ways:
  * - Invariants::LEAST_SQUARES: they are rows of the residual system, which is then overdetermined
  *   and solved in the least-squares sense;
  * - Invariants::PROJECTION: the integrator projects on them after every substep (the sensitivities
  *   include the projection), and they are not rows of the residual system.
  *
  * Instead of applying Gauss-Newton to \f$ \frac{1}{2}\|\mathbf{g}(\mathbf{x})\|^2 \f$, the residuals
  * \f$ \mathbf{r} \f$ are made variables of the nonlinear program
  *
  * \f[
  * \min_{\mathbf{x}, \mathbf{r}} \frac{1}{2}\|\mathbf{r}\|^2 \quad \text{s.t.} \quad
  * \mathbf{g}(\mathbf{x}) - \mathbf{r} = \mathbf{0} \text{,}
  * \f]
  *
  * whose constraint Jacobian \f$ [\partial_{\mathbf{x}}\mathbf{g}, -\mathbf{I}] \f$ always has full
  * rank, however overdetermined \f$ \mathbf{g} \f$ is. The second derivatives of \f$ \mathbf{g} \f$
  * through the integrator are not available, so the Hessian of the Lagrangian is the Gauss-Newton one,
  * identity on \f$ \mathbf{r} \f$ and zero on \f$ \mathbf{x} \f$. With it, the Newton system of an
  * interior-point solver is the augmented Gauss-Newton system
  *
  * \f[
  * \begin{bmatrix} \delta\mathbf{I} & \mathbf{0} & \partial_{\mathbf{x}}\mathbf{g}^\top \\
  *   \mathbf{0} & \mathbf{I} & -\mathbf{I} \\ \partial_{\mathbf{x}}\mathbf{g} & -\mathbf{I} &
  *   \mathbf{0} \end{bmatrix} \text{,}
  * \f]
  *
  * where the solver's inertia regularization \f$ \delta \f$ plays the role of a Tikhonov parameter and
  * its line search globalizes the iteration. The class is independent of the solver (see
  * Sandals::IpoptAdapter for IPOPT).
  *
  * \note The integrator of the problem is set to fixed steps (no adaptive stepping), and its
  * projection is enabled or disabled according to the treatment of the invariants.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The number of invariants.
  * \tparam Integrator The Runge-Kutta integrator type.
  */
  template <typename Real, Integer N, Integer M, typename Integrator>
  class LeastSquaresShooting
  {
  public:
    using ProblemType = Problem<Real, N, M, Integrator>; /**< Boundary value problem type. */
    using ProblemPtr  = std::shared_ptr<ProblemType>; /**< Boundary value problem shared pointer type. */
    using VectorF     = Eigen::Vector<Real, N>; /**< States vector type. */
    using MatrixF     = Eigen::Matrix<Real, N, N>; /**< States matrix type. */
    using VectorH     = Eigen::Vector<Real, M>; /**< Invariants vector type. */
    using MatrixJH    = Eigen::Matrix<Real, M, N>; /**< Invariants Jacobian type. */
    using VectorX     = Eigen::Vector<Real, Eigen::Dynamic>; /**< Dynamic vector type. */
    using MatrixXN    = Eigen::Matrix<Real, N, Eigen::Dynamic>; /**< States on the nodes type. */
    using FunctionX   = std::function<VectorF(Real const)>; /**< States guess function type. */

    /** Treatment of the invariants. */
    enum class Invariants {LEAST_SQUARES, PROJECTION};

  private:
    ProblemPtr m_problem;      /**< Boundary value problem (system, integrator, boundary conditions). */
    VectorX    m_t_mesh;       /**< Shooting mesh. */
    Integer    m_m{0};         /**< Number of shooting intervals. */
    Integer    m_substeps{1};  /**< Number of integration substeps per interval. */
    Invariants m_invariants{Invariants::LEAST_SQUARES}; /**< Treatment of the invariants. */
    VectorX    m_guess;        /**< Starting point [x; r]. */

    std::vector<Integer> m_jac_rows; /**< Jacobian sparsity rows. */
    std::vector<Integer> m_jac_cols; /**< Jacobian sparsity columns. */
    std::vector<Integer> m_hess_idx; /**< Hessian sparsity rows and columns (diagonal of the residuals). */

    // Cache of the integration of the intervals, valid for the nodes in m_cache_nodes
    mutable VectorX              m_cache_nodes; /**< Nodes of the cached integration. */
    mutable std::vector<VectorF> m_cache_end;   /**< Integrated (and possibly projected) end states. */
    mutable std::vector<MatrixF> m_cache_jac;   /**< Sensitivities of the end states to the nodes. */

  public:
    /**
    * Class constructor.
    * \param[in] t_problem The boundary value problem.
    * \param[in] t_t_mesh The shooting mesh.
    * \param[in] t_substeps The number of fixed integration substeps per interval.
    * \param[in] t_invariants The treatment of the invariants.
    */
    LeastSquaresShooting(ProblemPtr t_problem, VectorX const & t_t_mesh, Integer const t_substeps = 1,
      Invariants const t_invariants = Invariants::LEAST_SQUARES)
      : m_problem(t_problem), m_t_mesh(t_t_mesh), m_substeps(t_substeps), m_invariants(t_invariants)
    {
      #define CMD "Sandals::LeastSquaresShooting::LeastSquaresShooting(...): "
      SANDALS_ASSERT(this->m_t_mesh.size() >= 2, CMD "at least two mesh points are required.");
      SANDALS_ASSERT(this->m_substeps >= 1, CMD "at least one substep is required.");
      #undef CMD
      this->m_m = static_cast<Integer>(this->m_t_mesh.size()) - 1;
      this->m_problem->integrator()->disable_adaptive_mode();
      this->m_problem->integrator()->projection(this->projected());
      this->build_jacobian_structure();
      this->m_hess_idx.resize(this->residuals_number());
      for (Integer i{0}; i < this->residuals_number(); ++i) {this->m_hess_idx[i] = this->r_index(i);}
      this->m_guess = VectorX::Zero(this->variables_number());
    }

    /** \return True if the invariants are enforced by projection. */
    bool projected() const {return M > 0 && this->m_invariants == Invariants::PROJECTION;}

    /** \return True if the invariants are rows of the residual system. */
    bool least_squares_invariants() const {return M > 0 && this->m_invariants == Invariants::LEAST_SQUARES;}

    /** \return The shooting mesh. */
    VectorX const & t_mesh() const {return this->m_t_mesh;}

    /** \return The number of shooting intervals \f$ m \f$. */
    Integer intervals() const {return this->m_m;}

    /** \return The number of residuals (rows of the residual system). */
    Integer residuals_number() const
    {
      Integer const invariants_rows{this->least_squares_invariants() ? (this->m_m + 1)*M : 0};
      return this->m_m*N + invariants_rows + N;
    }

    /** \return The number of NLP variables (nodes and residuals). */
    Integer variables_number() const {return (this->m_m + 1)*N + this->residuals_number();}

    /** \return The number of NLP constraints (one per residual). */
    Integer constraints_number() const {return this->residuals_number();}

    /** \param[in] k Node index. \return The index of the node \f$ \mathbf{x}_k \f$. */
    Integer x_index(Integer const k) const {return k*N;}

    /** \param[in] i Residual index. \return The index of the residual variable \f$ r_i \f$. */
    Integer r_index(Integer const i) const {return (this->m_m + 1)*N + i;}

    /** \param[in] k Interval index. \return The first row of the continuity residuals. */
    Integer continuity_row(Integer const k) const {return k*N;}

    /** \param[in] k Node index. \return The first row of the invariants residuals. */
    Integer invariants_row(Integer const k) const {return this->m_m*N + k*M;}

    /** \return The first row of the boundary conditions residuals. */
    Integer boundary_row() const {return this->residuals_number() - N;}

    /** \param[in] z NLP variables. \param[in] k Node index. \return The node \f$ \mathbf{x}_k \f$. */
    VectorF x(VectorX const & z, Integer const k) const {return z.template segment<N>(this->x_index(k));}

    /** \param[in] z NLP variables. \return The residual variables \f$ \mathbf{r} \f$. */
    VectorX r(VectorX const & z) const {return z.tail(this->residuals_number());}

    /** \param[in] z NLP variables. \return The states on the nodes (\f$ N \times (m+1) \f$). */
    MatrixXN x_mesh(VectorX const & z) const
    {
      MatrixXN X(N, this->m_m + 1);
      for (Integer k{0}; k <= this->m_m; ++k) {X.col(k) = this->x(z, k);}
      return X;
    }

    /*\
     |  Initial guess
    \*/

    /** \return The starting point. */
    VectorX const & guess() const {return this->m_guess;}

    /** \param[in] t_guess The starting point. */
    void guess(VectorX const & t_guess)
    {
      SANDALS_ASSERT(t_guess.size() == this->variables_number(), "Sandals::LeastSquaresShooting::guess(...): wrong size.");
      this->m_guess = t_guess;
    }

    /**
    * Set the starting point from a states guess; the residual variables are set to the residuals of
    * the guess, so that the starting point is feasible.
    * \param[in] x_guess The states guess \f$ \mathbf{x}(t) \f$.
    */
    void guess(FunctionX const & x_guess)
    {
      VectorX & z{this->m_guess};
      for (Integer k{0}; k <= this->m_m; ++k) {z.template segment<N>(this->x_index(k)) = x_guess(this->m_t_mesh(k));}
      z.tail(this->residuals_number()) = this->residuals(z);
    }

    /*\
     |  Nonlinear program
    \*/

    /**
    * Get the bounds of the NLP variables and constraints: all the variables are free and all the
    * constraints are equalities.
    */
    void bounds(VectorX & z_lower, VectorX & z_upper, VectorX & c_lower, VectorX & c_upper) const
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      z_lower = VectorX::Constant(this->variables_number(), -inf);
      z_upper = VectorX::Constant(this->variables_number(), inf);
      c_lower = VectorX::Zero(this->constraints_number());
      c_upper = VectorX::Zero(this->constraints_number());
    }

    /** \param[in] z NLP variables. \return The objective \f$ \frac{1}{2}\|\mathbf{r}\|^2 \f$. */
    Real objective(VectorX const & z) const {return Real(0.5)*this->r(z).squaredNorm();}

    /** \param[in] z NLP variables. \return The objective gradient \f$ [\mathbf{0}; \mathbf{r}] \f$. */
    VectorX objective_gradient(VectorX const & z) const
    {
      VectorX gradient(VectorX::Zero(this->variables_number()));
      gradient.tail(this->residuals_number()) = this->r(z);
      return gradient;
    }

    /** \param[in] z NLP variables. \return The constraints \f$ \mathbf{g}(\mathbf{x}) - \mathbf{r} \f$. */
    VectorX constraints(VectorX const & z) const {return this->residuals(z) - this->r(z);}

    /**
    * Evaluate the residual system \f$ \mathbf{g}(\mathbf{x}) \f$.
    * \param[in] z NLP variables (only the nodes are used).
    * \return The residuals.
    */
    VectorX residuals(VectorX const & z) const
    {
      this->integrate(z);
      VectorX g(this->residuals_number());
      for (Integer k{0}; k < this->m_m; ++k) {
        g.template segment<N>(this->continuity_row(k)) = this->m_cache_end[k] - this->x(z, k+1);
      }
      if (this->least_squares_invariants()) {
        auto const & sys{this->m_problem->integrator()->system()};
        for (Integer k{0}; k <= this->m_m; ++k) {
          g.template segment<M>(this->invariants_row(k)) = sys->h(this->x(z, k), this->m_t_mesh(k));
        }
      }
      g.template segment<N>(this->boundary_row()) = this->m_problem->b(this->x(z, 0), this->x(z, this->m_m));
      return g;
    }

    /** \return The number of nonzeros of the constraints Jacobian. */
    Integer jacobian_nnz() const {return static_cast<Integer>(this->m_jac_rows.size());}

    /** \return The Jacobian sparsity rows. */
    std::vector<Integer> const & jacobian_rows() const {return this->m_jac_rows;}

    /** \return The Jacobian sparsity columns. */
    std::vector<Integer> const & jacobian_cols() const {return this->m_jac_cols;}

    /**
    * Evaluate the nonzeros of the constraints Jacobian \f$ [\partial_{\mathbf{x}}\mathbf{g}, -\mathbf{I}]
    * \f$, in the order given by jacobian_rows() and jacobian_cols().
    * \param[in] z NLP variables.
    * \return The Jacobian nonzeros.
    */
    VectorX jacobian_values(VectorX const & z) const
    {
      std::vector<Integer> rows, cols;
      std::vector<Real> values;
      this->template jacobian_traversal<true>(&z, rows, cols, values);
      return Eigen::Map<VectorX>(values.data(), values.size());
    }

    /** \return The number of nonzeros of the (Gauss-Newton) Hessian, the diagonal of the residuals. */
    Integer hessian_nnz() const {return this->residuals_number();}

    /** \return The Hessian sparsity rows. */
    std::vector<Integer> const & hessian_rows() const {return this->m_hess_idx;}

    /** \return The Hessian sparsity columns. */
    std::vector<Integer> const & hessian_cols() const {return this->m_hess_idx;}

    /**
    * Evaluate the Gauss-Newton Hessian of the Lagrangian, \f$ \sigma\mathbf{I} \f$ on the residuals
    * (exact for the objective, the curvature of the constraints is neglected).
    */
    VectorX hessian_values(VectorX const & /*z*/, Real const sigma, VectorX const & /*lambda*/) const
    {
      return VectorX::Constant(this->residuals_number(), sigma);
    }

  private:
    /**
    * Integrate every interval from its node with fixed substeps, propagating the sensitivities (and
    * projecting after every substep, if requested). The result is cached for the given nodes.
    */
    void integrate(VectorX const & z) const
    {
      VectorX const nodes{z.head((this->m_m + 1)*N)};
      if (this->m_cache_end.size() == static_cast<std::size_t>(this->m_m) && nodes == this->m_cache_nodes) {return;}

      auto const & integrator{this->m_problem->integrator()};
      this->m_cache_end.resize(this->m_m);
      this->m_cache_jac.resize(this->m_m);
      Solution<Real, N, M> solution;
      for (Integer k{0}; k < this->m_m; ++k) {
        VectorX const t_substeps{VectorX::LinSpaced(this->m_substeps + 1, this->m_t_mesh(k), this->m_t_mesh(k+1))};
        if (!integrator->solve(t_substeps, this->x(z, k), solution, this->m_cache_jac[k])) {
          this->m_cache_end.clear(); // invalidate the cache
          SANDALS_ERROR("Sandals::LeastSquaresShooting::integrate(...): integration failed on interval " << k << ".");
        }
        this->m_cache_end[k] = solution.x.col(solution.size() - 1);
      }
      this->m_cache_nodes = nodes;
    }

    /**
    * Traverse the constraints Jacobian in a fixed order, collecting the sparsity pattern and,
    * optionally, the values.
    */
    template <bool Values>
    void jacobian_traversal(VectorX const * z, std::vector<Integer> & rows, std::vector<Integer> & cols,
      std::vector<Real> & values) const
    {
      auto push = [&] (Integer const r, Integer const c, Real const v) {
        rows.push_back(r); cols.push_back(c);
        if constexpr (Values) {values.push_back(v);} else {(void)v;}
      };
      auto push_block = [&] (Integer const row, Integer const col, auto const & block) {
        for (Integer r{0}; r < block.rows(); ++r) {
          for (Integer c{0}; c < block.cols(); ++c) {push(row + r, col + c, block(r, c));}
        }
      };
      if constexpr (Values) {this->integrate(*z);}

      // Continuity: sensitivity of the end state on x_k, -I on x_{k+1}
      for (Integer k{0}; k < this->m_m; ++k) {
        MatrixF const A{Values ? this->m_cache_jac[k] : MatrixF::Zero()};
        push_block(this->continuity_row(k), this->x_index(k), A);
        for (Integer r{0}; r < N; ++r) {push(this->continuity_row(k) + r, this->x_index(k+1) + r, Real(-1.0));}
      }

      // Invariants on every node
      if (this->least_squares_invariants()) {
        auto const & sys{this->m_problem->integrator()->system()};
        for (Integer k{0}; k <= this->m_m; ++k) {
          MatrixJH const Jh{Values ? MatrixJH(sys->Jh_x(this->x(*z, k), this->m_t_mesh(k))) : MatrixJH::Zero()};
          push_block(this->invariants_row(k), this->x_index(k), Jh);
        }
      }

      // Boundary conditions on x_0 and x_m
      MatrixF Jb_ini(MatrixF::Zero()), Jb_end(MatrixF::Zero());
      if constexpr (Values) {
        VectorF const x_ini{this->x(*z, 0)}, x_end{this->x(*z, this->m_m)};
        Jb_ini = this->m_problem->Jb_x_ini(x_ini, x_end);
        Jb_end = this->m_problem->Jb_x_end(x_ini, x_end);
      }
      push_block(this->boundary_row(), this->x_index(0), Jb_ini);
      push_block(this->boundary_row(), this->x_index(this->m_m), Jb_end);

      // Residual variables: -I
      for (Integer i{0}; i < this->residuals_number(); ++i) {push(i, this->r_index(i), Real(-1.0));}
    }

    /** Build the Jacobian sparsity pattern. */
    void build_jacobian_structure()
    {
      std::vector<Real> unused;
      this->m_jac_rows.clear();
      this->m_jac_cols.clear();
      this->template jacobian_traversal<false>(nullptr, this->m_jac_rows, this->m_jac_cols, unused);
    }

  }; // class LeastSquaresShooting

} // namespace Sandals

#endif // SANDALS_LEAST_SQUARES_SHOOTING_HH
