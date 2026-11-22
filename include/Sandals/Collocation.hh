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

#ifndef SANDALS_COLLOCATION_HH
#define SANDALS_COLLOCATION_HH

#include <Sandals.hh>
#include <Sandals/Tableau.hh>
#include <Sandals/RungeKutta.hh>
#include <Sandals/System/Implicit.hh>
#include <Sandals/OptimalControl.hh>

namespace Sandals
{

  /*\
   |    ____      _ _                 _   _
   |   / ___|___ | | | ___   ___ __ _| |_(_) ___  _ __
   |  | |   / _ \| | |/ _ \ / __/ _` | __| |/ _ \| '_ \
   |  | |__| (_) | | | (_) | (_| (_| | |_| | (_) | | | |
   |   \____\___/|_|_|\___/ \___\__,_|\__|_|\___/|_| |_|
   |
  \*/

  /**
  * \brief Direct collocation transcription of an optimal control problem with a Runge-Kutta method.
  *
  * The optimal control problem is transcribed into a sparse nonlinear program (NLP) on the mesh
  * \f$ t_0 < t_1 < \dots < t_n \f$, with step \f$ h_k = t_{k+1} - t_k \f$. Each interval holds the
  * states \f$ \mathbf{x}_k \f$, the stage variables \f$ \mathbf{K}_k = [\mathbf{K}_{k,1}, \dots,
  * \mathbf{K}_{k,S}] \f$ (with \f$ \mathbf{K}_{k,i} = h_k \mathbf{x}^{\prime}_{k,i} \f$, as in the
  * Runge-Kutta integrators), and piecewise-constant controls \f$ \mathbf{u}_k \f$. The final states
  * \f$ \mathbf{x}_n \f$ close the variables vector
  *
  * \f[
  * \mathbf{z} = [\mathbf{x}_0, \mathbf{K}_0, \mathbf{u}_0, \dots, \mathbf{x}_{n-1}, \mathbf{K}_{n-1},
  * \mathbf{u}_{n-1}, \mathbf{x}_n] \text{.}
  * \f]
  *
  * The constraints are, for each interval, the stage equations and the continuity conditions
  *
  * \f[
  * \mathbf{F}\left(\mathbf{x}_k + \textstyle\sum_j a_{ij}\mathbf{K}_{k,j}, \mathbf{K}_{k,i}/h_k,
  * \mathbf{u}_k, t_k + c_i h_k\right) = \mathbf{0}, \quad \mathbf{x}_k + \textstyle\sum_i b_i
  * \mathbf{K}_{k,i} - \mathbf{x}_{k+1} = \mathbf{0} \text{,}
  * \f]
  *
  * followed by the boundary conditions and the path constraints at every mesh node (at the last node
  * the controls of the last interval are used). The Lagrange term is integrated with the quadrature
  * of the tableau at the stage nodes. The first derivatives (gradient and Jacobian) and the second
  * derivatives (Hessian of the Lagrangian, lower triangle) are exact up to the derivatives provided
  * by the system and the problem. The sparsity patterns are fixed and computed once.
  *
  * \note For DAEs the initial states must be consistent: the boundary conditions must also fix the
  * initial value of the algebraic states (e.g., by imposing the algebraic equations at \f$ t_0 \f$).
  * With stiffly accurate methods (e.g., RadauIIA) the algebraic part of \f$ \mathbf{x}_0 \f$ only
  * enters through the stage values, so leaving it free makes the NLP singular.
  * \tparam Real The scalar number type.
  * \tparam S The number of stages of the Runge-Kutta method.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam U The number of controls.
  */
  template <typename Real, Integer S, Integer N, Integer U>
  class Collocation
  {
  public:
    using Problem    = OptimalControlProblem<Real, N, U>; /**< Optimal control problem type. */
    using ProblemPtr = std::shared_ptr<Problem>; /**< Optimal control problem shared pointer type. */
    using VectorF    = Eigen::Vector<Real, N>; /**< States vector type. */
    using VectorU    = Eigen::Vector<Real, U>; /**< Controls vector type. */
    using MatrixF    = Eigen::Matrix<Real, N, N>; /**< States matrix type. */
    using MatrixK    = Eigen::Matrix<Real, N, S>; /**< Stage variables matrix type. */
    using VectorX    = Eigen::Vector<Real, Eigen::Dynamic>; /**< Dynamic vector type. */
    using MatrixX    = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; /**< Dynamic matrix type. */
    using MatrixXN   = Eigen::Matrix<Real, N, Eigen::Dynamic>; /**< States trajectory type. */
    using MatrixXU   = Eigen::Matrix<Real, U, Eigen::Dynamic>; /**< Controls trajectory type. */
    using FunctionX  = std::function<VectorF(Real const)>; /**< States guess function type. */
    using FunctionU  = std::function<VectorU(Real const)>; /**< Controls guess function type. */

    static constexpr Integer D{N + N*S + U}; /**< Number of variables per interval. */
    static constexpr Integer C{N*S + N};     /**< Number of dynamics constraints per interval. */

  private:
    using Index = std::pair<Integer, Integer>;

    ProblemPtr       m_problem;  /**< Optimal control problem. */
    Tableau<Real, S> m_tableau;  /**< Butcher tableau of the collocation method. */
    VectorX          m_t_mesh;   /**< Mesh of the independent variable. */
    Integer          m_n{0};     /**< Number of intervals. */
    Integer          m_nb{0};    /**< Number of boundary conditions. */
    Integer          m_ng{0};    /**< Number of path constraints. */
    VectorX          m_guess;    /**< Initial guess of the NLP variables. */

    std::vector<Integer>    m_jac_rows;  /**< Jacobian sparsity rows. */
    std::vector<Integer>    m_jac_cols;  /**< Jacobian sparsity columns. */
    std::vector<Integer>    m_hess_rows; /**< Hessian (lower triangle) sparsity rows. */
    std::vector<Integer>    m_hess_cols; /**< Hessian (lower triangle) sparsity columns. */
    std::map<Index, Integer> m_hess_map; /**< Hessian (row, column) to nonzero index map. */

  public:
    /**
    * Class constructor for the collocation transcription on a given mesh.
    * \param[in] t_problem The optimal control problem.
    * \param[in] t_tableau The Butcher tableau of the collocation method.
    * \param[in] t_t_mesh The mesh, from the initial to the final time of the problem.
    */
    Collocation(ProblemPtr t_problem, Tableau<Real, S> const & t_tableau, VectorX const & t_t_mesh)
      : m_problem(t_problem), m_tableau(t_tableau), m_t_mesh(t_t_mesh)
    {
      #define CMD "Sandals::Collocation::Collocation(...): "
      SANDALS_ASSERT(this->m_t_mesh.size() >= 2, CMD "at least two mesh points are required.");
      SANDALS_ASSERT(std::abs(this->m_t_mesh(0) - t_problem->t_ini()) < 1.0e-12 &&
        std::abs(this->m_t_mesh(this->m_t_mesh.size()-1) - t_problem->t_end()) < 1.0e-12,
        CMD "the mesh must span the problem time interval.");
      for (Integer k{0}; k < this->m_t_mesh.size()-1; ++k) {
        SANDALS_ASSERT(this->m_t_mesh(k+1) > this->m_t_mesh(k), CMD "the mesh must be increasing.");
      }
      #undef CMD
      this->m_n  = static_cast<Integer>(this->m_t_mesh.size()) - 1;
      this->m_nb = t_problem->boundary_conditions_number();
      this->m_ng = t_problem->path_constraints_number();
      this->build_jacobian_structure();
      this->build_hessian_structure();
      this->m_guess = VectorX::Zero(this->variables_number());
    }

    /**
    * Class constructor for the collocation transcription on a uniform mesh.
    * \param[in] t_problem The optimal control problem.
    * \param[in] t_tableau The Butcher tableau of the collocation method.
    * \param[in] t_intervals The number of mesh intervals.
    */
    Collocation(ProblemPtr t_problem, Tableau<Real, S> const & t_tableau, Integer const t_intervals)
      : Collocation(t_problem, t_tableau, VectorX::LinSpaced(t_intervals + 1, t_problem->t_ini(), t_problem->t_end()))
    {}

    /** \return The optimal control problem. */
    ProblemPtr const & problem() const {return this->m_problem;}

    /** \return The Butcher tableau. */
    Tableau<Real, S> const & tableau() const {return this->m_tableau;}

    /** \return The mesh. */
    VectorX const & t_mesh() const {return this->m_t_mesh;}

    /** \return The number of mesh intervals \f$ n \f$. */
    Integer intervals() const {return this->m_n;}

    /** \return The number of NLP variables. */
    Integer variables_number() const {return this->m_n*D + N;}

    /** \return The number of NLP constraints. */
    Integer constraints_number() const {return this->m_n*C + this->m_nb + (this->m_n + 1)*this->m_ng;}

    /** \return The number of nonzeros of the constraints Jacobian. */
    Integer jacobian_nnz() const {return static_cast<Integer>(this->m_jac_rows.size());}

    /** \return The number of nonzeros of the lower triangle of the Lagrangian Hessian. */
    Integer hessian_nnz() const {return static_cast<Integer>(this->m_hess_rows.size());}

    /** \return The Jacobian sparsity rows. */
    std::vector<Integer> const & jacobian_rows() const {return this->m_jac_rows;}

    /** \return The Jacobian sparsity columns. */
    std::vector<Integer> const & jacobian_cols() const {return this->m_jac_cols;}

    /** \return The Hessian (lower triangle) sparsity rows. */
    std::vector<Integer> const & hessian_rows() const {return this->m_hess_rows;}

    /** \return The Hessian (lower triangle) sparsity columns. */
    std::vector<Integer> const & hessian_cols() const {return this->m_hess_cols;}

    /** \param[in] k Mesh node index (\f$ 0 \leq k \leq n \f$). \return The index of \f$ \mathbf{x}_k \f$. */
    Integer x_index(Integer const k) const {return k*D;}

    /** \param[in] k Interval index. \param[in] i Stage index. \return The index of \f$ \mathbf{K}_{k,i} \f$. */
    Integer K_index(Integer const k, Integer const i) const {return k*D + N + i*N;}

    /** \param[in] k Interval index. \return The index of \f$ \mathbf{u}_k \f$. */
    Integer u_index(Integer const k) const {return k*D + N + N*S;}

    /** \param[in] k Interval index. \param[in] i Stage index. \return The first row of the stage equations. */
    Integer stage_row(Integer const k, Integer const i) const {return k*C + i*N;}

    /** \param[in] k Interval index. \return The first row of the continuity conditions. */
    Integer continuity_row(Integer const k) const {return k*C + S*N;}

    /** \return The first row of the boundary conditions. */
    Integer boundary_row() const {return this->m_n*C;}

    /** \param[in] k Mesh node index. \return The first row of the path constraints at the node. */
    Integer path_row(Integer const k) const {return this->m_n*C + this->m_nb + k*this->m_ng;}

    /** \param[in] z NLP variables. \param[in] k Mesh node index. \return The states \f$ \mathbf{x}_k \f$. */
    VectorF x(VectorX const & z, Integer const k) const {return z.template segment<N>(this->x_index(k));}

    /** \param[in] z NLP variables. \param[in] k Interval index. \return The controls \f$ \mathbf{u}_k \f$. */
    VectorU u(VectorX const & z, Integer const k) const {return z.template segment<U>(this->u_index(k));}

    /** \param[in] z NLP variables. \param[in] k Interval index. \return The stage variables \f$ \mathbf{K}_k \f$. */
    MatrixK K(VectorX const & z, Integer const k) const
    {return Eigen::Map<MatrixK const>(z.data() + this->K_index(k, 0));}

    /** \param[in] z NLP variables. \return The states on the mesh nodes (\f$ N \times (n+1) \f$). */
    MatrixXN x_mesh(VectorX const & z) const
    {
      MatrixXN X(N, this->m_n + 1);
      for (Integer k{0}; k <= this->m_n; ++k) {X.col(k) = this->x(z, k);}
      return X;
    }

    /** \param[in] z NLP variables. \return The piecewise-constant controls (\f$ U \times n \f$). */
    MatrixXU u_mesh(VectorX const & z) const
    {
      MatrixXU Uk(U, this->m_n);
      for (Integer k{0}; k < this->m_n; ++k) {Uk.col(k) = this->u(z, k);}
      return Uk;
    }

    /*\
     |  Initial guess
    \*/

    /** \return The initial guess of the NLP variables. */
    VectorX const & guess() const {return this->m_guess;}

    /** \param[in] t_guess The initial guess of the NLP variables. */
    void guess(VectorX const & t_guess)
    {
      SANDALS_ASSERT(t_guess.size() == this->variables_number(), "Sandals::Collocation::guess(...): wrong size.");
      this->m_guess = t_guess;
    }

    /**
    * Set the initial guess from states and controls functions of time. The stage variables are
    * initialized with the states increment over each interval.
    * \param[in] x_guess The states guess \f$ \mathbf{x}(t) \f$.
    * \param[in] u_guess The controls guess \f$ \mathbf{u}(t) \f$ (sampled at the intervals midpoints).
    */
    void guess(FunctionX const & x_guess, FunctionU const & u_guess)
    {
      VectorX & z{this->m_guess};
      for (Integer k{0}; k <= this->m_n; ++k) {z.template segment<N>(this->x_index(k)) = x_guess(this->m_t_mesh(k));}
      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        VectorF const dx{x_guess(t_k + h_k) - x_guess(t_k)};
        for (Integer i{0}; i < S; ++i) {z.template segment<N>(this->K_index(k, i)) = dx;}
        z.template segment<U>(this->u_index(k)) = u_guess(t_k + h_k/Real(2.0));
      }
    }

    /**
    * Set the initial guess by integrating the system forward with the same Runge-Kutta method, from
    * the given initial states and with the given (piecewise-constant) controls. The resulting guess
    * satisfies the dynamics constraints.
    * \param[in] x_ini The initial states \f$ \mathbf{x}_0 \f$ (should be consistent for DAEs).
    * \param[in] u_guess The controls guess \f$ \mathbf{u}(t) \f$ (sampled at the intervals midpoints).
    * \return True if the integration succeeded, false otherwise.
    */
    bool simulate_guess(VectorF const & x_ini, FunctionU const & u_guess)
    {
      auto sys{this->m_problem->system()};
      VectorU u_k;
      auto adapter = std::make_shared<ImplicitWrapper<Real, N, 0>>("CollocationGuess",
        [&] (VectorF const & x, VectorF const & x_dot, Real const t) {return sys->F(x, x_dot, u_k, t);},
        [&] (VectorF const & x, VectorF const & x_dot, Real const t) {return sys->JF_x(x, x_dot, u_k, t);},
        [&] (VectorF const & x, VectorF const & x_dot, Real const t) {return sys->JF_x_dot(x, x_dot, u_k, t);});
      RungeKutta<Real, S, N, 0> rk(this->m_tableau, adapter);
      rk.disable_adaptive_mode();
      rk.disable_projection();

      VectorX & z{this->m_guess};
      VectorF x_k{x_ini}, x_new;
      MatrixK K;
      Real h_new;
      z.template segment<N>(this->x_index(0)) = x_k;
      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        u_k = u_guess(t_k + h_k/Real(2.0));
        if (!rk.step(x_k, t_k, h_k, x_new, h_new, K)) {return false;}
        Eigen::Map<MatrixK>(z.data() + this->K_index(k, 0)) = K;
        z.template segment<U>(this->u_index(k)) = u_k;
        z.template segment<N>(this->x_index(k+1)) = x_new;
        x_k = x_new;
      }
      return true;
    }

    /*\
     |  Bounds
    \*/

    /**
    * Get the bounds of the NLP variables and constraints (infinite bounds are kept as infinity).
    * \param[out] z_lower Variables lower bounds.
    * \param[out] z_upper Variables upper bounds.
    * \param[out] c_lower Constraints lower bounds.
    * \param[out] c_upper Constraints upper bounds.
    */
    void bounds(VectorX & z_lower, VectorX & z_upper, VectorX & c_lower, VectorX & c_upper) const
    {
      Real const inf{std::numeric_limits<Real>::infinity()};
      z_lower = VectorX::Constant(this->variables_number(), -inf);
      z_upper = VectorX::Constant(this->variables_number(), inf);
      for (Integer k{0}; k <= this->m_n; ++k) {
        z_lower.template segment<N>(this->x_index(k)) = this->m_problem->x_lower();
        z_upper.template segment<N>(this->x_index(k)) = this->m_problem->x_upper();
      }
      for (Integer k{0}; k < this->m_n; ++k) {
        z_lower.template segment<U>(this->u_index(k)) = this->m_problem->u_lower();
        z_upper.template segment<U>(this->u_index(k)) = this->m_problem->u_upper();
      }
      c_lower = VectorX::Zero(this->constraints_number());
      c_upper = VectorX::Zero(this->constraints_number());
      if (this->m_ng > 0) {
        VectorX const g_lower{this->m_problem->path_constraints_lower()};
        VectorX const g_upper{this->m_problem->path_constraints_upper()};
        for (Integer k{0}; k <= this->m_n; ++k) {
          c_lower.segment(this->path_row(k), this->m_ng) = g_lower;
          c_upper.segment(this->path_row(k), this->m_ng) = g_upper;
        }
      }
    }

    /*\
     |  Objective
    \*/

    /**
    * Evaluate the objective function.
    * \param[in] z NLP variables.
    * \return The objective value.
    */
    Real objective(VectorX const & z) const
    {
      Real obj{this->m_problem->mayer(this->x(z, 0), this->x(z, this->m_n))};
      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        VectorF const x_k{this->x(z, k)};
        VectorU const u_k{this->u(z, k)};
        MatrixK const K_k{this->K(z, k)};
        for (Integer i{0}; i < S; ++i) {
          VectorF const x_i{x_k + K_k * this->m_tableau.A.row(i).transpose()};
          obj += h_k * this->m_tableau.b(i) * this->m_problem->lagrange(x_i, u_k, t_k + this->m_tableau.c(i)*h_k);
        }
      }
      return obj;
    }

    /**
    * Evaluate the gradient of the objective function.
    * \param[in] z NLP variables.
    * \return The gradient (dense).
    */
    VectorX objective_gradient(VectorX const & z) const
    {
      VectorX grad(VectorX::Zero(this->variables_number()));
      auto const g_phi{this->m_problem->mayer_gradient(this->x(z, 0), this->x(z, this->m_n))};
      grad.template segment<N>(this->x_index(0))       += g_phi.template head<N>();
      grad.template segment<N>(this->x_index(this->m_n)) += g_phi.template tail<N>();
      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        VectorF const x_k{this->x(z, k)};
        VectorU const u_k{this->u(z, k)};
        MatrixK const K_k{this->K(z, k)};
        for (Integer i{0}; i < S; ++i) {
          VectorF const x_i{x_k + K_k * this->m_tableau.A.row(i).transpose()};
          Real const w{h_k * this->m_tableau.b(i)};
          auto const g_L{this->m_problem->lagrange_gradient(x_i, u_k, t_k + this->m_tableau.c(i)*h_k)};
          VectorF const g_x{w * g_L.template head<N>()};
          grad.template segment<N>(this->x_index(k)) += g_x;
          for (Integer j{0}; j < S; ++j) {grad.template segment<N>(this->K_index(k, j)) += this->m_tableau.A(i, j) * g_x;}
          grad.template segment<U>(this->u_index(k)) += w * g_L.template tail<U>();
        }
      }
      return grad;
    }

    /*\
     |  Constraints
    \*/

    /**
    * Evaluate the constraints.
    * \param[in] z NLP variables.
    * \return The constraints values.
    */
    VectorX constraints(VectorX const & z) const
    {
      VectorX c(this->constraints_number());
      auto const & sys{this->m_problem->system()};
      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        VectorF const x_k{this->x(z, k)};
        VectorU const u_k{this->u(z, k)};
        MatrixK const K_k{this->K(z, k)};
        for (Integer i{0}; i < S; ++i) {
          VectorF const x_i{x_k + K_k * this->m_tableau.A.row(i).transpose()};
          c.template segment<N>(this->stage_row(k, i)) =
            sys->F(x_i, K_k.col(i)/h_k, u_k, t_k + this->m_tableau.c(i)*h_k);
        }
        c.template segment<N>(this->continuity_row(k)) = x_k + K_k * this->m_tableau.b - this->x(z, k+1);
      }
      if (this->m_nb > 0) {
        c.segment(this->boundary_row(), this->m_nb) =
          this->m_problem->boundary_conditions(this->x(z, 0), this->x(z, this->m_n));
      }
      if (this->m_ng > 0) {
        for (Integer k{0}; k <= this->m_n; ++k) {
          c.segment(this->path_row(k), this->m_ng) =
            this->m_problem->path_constraints(this->x(z, k), this->u(z, std::min(k, this->m_n-1)), this->m_t_mesh(k));
        }
      }
      return c;
    }

    /**
    * Evaluate the nonzeros of the constraints Jacobian, in the order given by jacobian_rows() and
    * jacobian_cols().
    * \param[in] z NLP variables.
    * \return The Jacobian nonzeros.
    */
    VectorX jacobian_values(VectorX const & z) const
    {
      std::vector<Integer> rows, cols;
      std::vector<Real> vals;
      this->template jacobian_traversal<true>(&z, rows, cols, vals);
      return Eigen::Map<VectorX>(vals.data(), vals.size());
    }

    /**
    * Evaluate the constraints Jacobian as a dense matrix (for testing and debugging).
    * \param[in] z NLP variables.
    * \return The dense Jacobian.
    */
    MatrixX jacobian_dense(VectorX const & z) const
    {
      MatrixX J(MatrixX::Zero(this->constraints_number(), this->variables_number()));
      VectorX const v{this->jacobian_values(z)};
      for (Integer e{0}; e < this->jacobian_nnz(); ++e) {J(this->m_jac_rows[e], this->m_jac_cols[e]) += v(e);}
      return J;
    }

    /*\
     |  Hessian of the Lagrangian
    \*/

    /**
    * Evaluate the nonzeros of the lower triangle of the Hessian of the Lagrangian \f$ \sigma f(
    * \mathbf{z}) + \boldsymbol{\lambda}^\top\mathbf{c}(\mathbf{z}) \f$, in the order given by
    * hessian_rows() and hessian_cols().
    * \param[in] z NLP variables.
    * \param[in] sigma Objective factor \f$ \sigma \f$.
    * \param[in] lambda Constraints multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The Hessian nonzeros.
    */
    VectorX hessian_values(VectorX const & z, Real const sigma, VectorX const & lambda) const
    {
      using MatrixW = Eigen::Matrix<Real, 2*N+U, 2*N+U>;
      using MatrixPW = Eigen::Matrix<Real, 2*N+U, D>;
      using MatrixPL = Eigen::Matrix<Real, N+U, D>;
      using MatrixD = Eigen::Matrix<Real, D, D>;

      VectorX vals(VectorX::Zero(this->hessian_nnz()));
      auto const & sys{this->m_problem->system()};
      std::vector<Integer> idx_interval(D);

      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        VectorF const x_k{this->x(z, k)};
        VectorU const u_k{this->u(z, k)};
        MatrixK const K_k{this->K(z, k)};
        MatrixD H(MatrixD::Zero());

        for (Integer i{0}; i < S; ++i) {
          Real const t_i{t_k + this->m_tableau.c(i)*h_k};
          VectorF const x_i{x_k + K_k * this->m_tableau.A.row(i).transpose()};
          VectorF const x_dot_i{K_k.col(i)/h_k};

          // Map from the interval variables to the stage arguments (x, x', u)
          MatrixPW P(MatrixPW::Zero());
          P.template block<N, N>(0, 0).setIdentity();
          for (Integer j{0}; j < S; ++j) {P.template block<N, N>(0, N + j*N) = this->m_tableau.A(i, j) * MatrixF::Identity();}
          P.template block<N, N>(N, N + i*N) = MatrixF::Identity() / h_k;
          P.template block<U, U>(2*N, N + N*S).setIdentity();

          // Stage equations contribution
          VectorF const l_i{lambda.template segment<N>(this->stage_row(k, i))};
          MatrixW W;
          W.template block<N, N>(0, 0) = sys->HF_x_x(x_i, x_dot_i, u_k, t_i, l_i);
          W.template block<N, N>(0, N) = sys->HF_x_x_dot(x_i, x_dot_i, u_k, t_i, l_i);
          W.template block<N, N>(N, 0) = W.template block<N, N>(0, N).transpose();
          W.template block<N, N>(N, N) = sys->HF_x_dot_x_dot(x_i, x_dot_i, u_k, t_i, l_i);
          W.template block<N, U>(0, 2*N) = sys->HF_x_u(x_i, x_dot_i, u_k, t_i, l_i);
          W.template block<N, U>(N, 2*N) = sys->HF_x_dot_u(x_i, x_dot_i, u_k, t_i, l_i);
          W.template block<U, N>(2*N, 0) = W.template block<N, U>(0, 2*N).transpose();
          W.template block<U, N>(2*N, N) = W.template block<N, U>(N, 2*N).transpose();
          W.template block<U, U>(2*N, 2*N) = sys->HF_u_u(x_i, x_dot_i, u_k, t_i, l_i);
          H.noalias() += P.transpose() * W * P;

          // Lagrange term contribution (variables x, u)
          if (sigma != Real(0.0)) {
            MatrixPL Q;
            Q.template topRows<N>() = P.template topRows<N>();
            Q.template bottomRows<U>() = P.template bottomRows<U>();
            H.noalias() += (sigma * h_k * this->m_tableau.b(i)) * Q.transpose() *
              this->m_problem->lagrange_hessian(x_i, u_k, t_i) * Q;
          }
        }

        // Path constraints contribution at the node k (variables x_k, u_k)
        if (this->m_ng > 0) {
          MatrixPL R(MatrixPL::Zero());
          R.template block<N, N>(0, 0).setIdentity();
          R.template block<U, U>(N, N + N*S).setIdentity();
          H.noalias() += R.transpose() * this->m_problem->path_constraints_hessian(x_k, u_k, t_k,
            lambda.segment(this->path_row(k), this->m_ng)) * R;
        }

        for (Integer a{0}; a < D; ++a) {idx_interval[a] = this->x_index(k) + a;}
        this->accumulate(idx_interval, H, vals);
      }

      // Path constraints at the last node (variables x_n, u_{n-1})
      if (this->m_ng > 0) {
        Integer const n{this->m_n};
        std::vector<Integer> idx(this->path_indices(n));
        this->accumulate(idx, this->m_problem->path_constraints_hessian(this->x(z, n), this->u(z, n-1),
          this->m_t_mesh(n), lambda.segment(this->path_row(n), this->m_ng)), vals);
      }

      // Mayer term and boundary conditions (variables x_0, x_n)
      {
        VectorF const x_0{this->x(z, 0)}, x_n{this->x(z, this->m_n)};
        Eigen::Matrix<Real, 2*N, 2*N> B(sigma * this->m_problem->mayer_hessian(x_0, x_n));
        if (this->m_nb > 0) {
          B += this->m_problem->boundary_conditions_hessian(x_0, x_n, lambda.segment(this->boundary_row(), this->m_nb));
        }
        this->accumulate(this->boundary_indices(), B, vals);
      }

      return vals;
    }

    /**
    * Evaluate the Hessian of the Lagrangian as a dense symmetric matrix (for testing and debugging).
    * \param[in] z NLP variables.
    * \param[in] sigma Objective factor \f$ \sigma \f$.
    * \param[in] lambda Constraints multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The dense Hessian.
    */
    MatrixX hessian_dense(VectorX const & z, Real const sigma, VectorX const & lambda) const
    {
      MatrixX H(MatrixX::Zero(this->variables_number(), this->variables_number()));
      VectorX const v{this->hessian_values(z, sigma, lambda)};
      for (Integer e{0}; e < this->hessian_nnz(); ++e) {
        H(this->m_hess_rows[e], this->m_hess_cols[e]) += v(e);
        if (this->m_hess_rows[e] != this->m_hess_cols[e]) {H(this->m_hess_cols[e], this->m_hess_rows[e]) += v(e);}
      }
      return H;
    }

  private:
    /**
    * Traverse the constraints Jacobian in a fixed order, collecting the sparsity pattern and,
    * optionally, the values.
    */
    template <bool Values>
    void jacobian_traversal(VectorX const * z, std::vector<Integer> & rows, std::vector<Integer> & cols,
      std::vector<Real> & vals) const
    {
      auto const & sys{this->m_problem->system()};
      auto push = [&] (Integer const r, Integer const c, Real const v) {
        rows.push_back(r); cols.push_back(c);
        if constexpr (Values) {vals.push_back(v);} else {(void)v;}
      };

      for (Integer k{0}; k < this->m_n; ++k) {
        Real const t_k{this->m_t_mesh(k)}, h_k{this->m_t_mesh(k+1) - t_k};
        VectorF x_k; VectorU u_k; MatrixK K_k;
        if constexpr (Values) {x_k = this->x(*z, k); u_k = this->u(*z, k); K_k = this->K(*z, k);}

        // Stage equations (dense N x D blocks)
        for (Integer i{0}; i < S; ++i) {
          Eigen::Matrix<Real, N, D> B(Eigen::Matrix<Real, N, D>::Zero());
          if constexpr (Values) {
            Real const t_i{t_k + this->m_tableau.c(i)*h_k};
            VectorF const x_i{x_k + K_k * this->m_tableau.A.row(i).transpose()};
            VectorF const x_dot_i{K_k.col(i)/h_k};
            MatrixF const JF_x{sys->JF_x(x_i, x_dot_i, u_k, t_i)};
            B.template block<N, N>(0, 0) = JF_x;
            for (Integer j{0}; j < S; ++j) {B.template block<N, N>(0, N + j*N) = this->m_tableau.A(i, j) * JF_x;}
            B.template block<N, N>(0, N + i*N) += sys->JF_x_dot(x_i, x_dot_i, u_k, t_i) / h_k;
            B.template block<N, U>(0, N + N*S) = sys->JF_u(x_i, x_dot_i, u_k, t_i);
          }
          for (Integer r{0}; r < N; ++r) {
            for (Integer c{0}; c < D; ++c) {push(this->stage_row(k, i) + r, this->x_index(k) + c, B(r, c));}
          }
        }

        // Continuity conditions (diagonal blocks)
        for (Integer r{0}; r < N; ++r) {
          Integer const row{this->continuity_row(k) + r};
          push(row, this->x_index(k) + r, Real(1.0));
          for (Integer i{0}; i < S; ++i) {push(row, this->K_index(k, i) + r, this->m_tableau.b(i));}
          push(row, this->x_index(k+1) + r, Real(-1.0));
        }
      }

      // Boundary conditions (dense blocks on x_0 and x_n)
      if (this->m_nb > 0) {
        MatrixX Jb(MatrixX::Zero(this->m_nb, 2*N));
        if constexpr (Values) {Jb = this->m_problem->boundary_conditions_jacobian(this->x(*z, 0), this->x(*z, this->m_n));}
        for (Integer r{0}; r < this->m_nb; ++r) {
          for (Integer c{0}; c < N; ++c) {push(this->boundary_row() + r, this->x_index(0) + c, Jb(r, c));}
          for (Integer c{0}; c < N; ++c) {push(this->boundary_row() + r, this->x_index(this->m_n) + c, Jb(r, N + c));}
        }
      }

      // Path constraints (dense blocks on x_k and u_k)
      if (this->m_ng > 0) {
        for (Integer k{0}; k <= this->m_n; ++k) {
          MatrixX Jg(MatrixX::Zero(this->m_ng, N+U));
          if constexpr (Values) {
            Jg = this->m_problem->path_constraints_jacobian(this->x(*z, k), this->u(*z, std::min(k, this->m_n-1)), this->m_t_mesh(k));
          }
          std::vector<Integer> const idx(this->path_indices(k));
          for (Integer r{0}; r < this->m_ng; ++r) {
            for (Integer c{0}; c < N+U; ++c) {push(this->path_row(k) + r, idx[c], Jg(r, c));}
          }
        }
      }
    }

    /** Build the Jacobian sparsity pattern. */
    void build_jacobian_structure()
    {
      std::vector<Real> dummy;
      this->m_jac_rows.clear();
      this->m_jac_cols.clear();
      this->template jacobian_traversal<false>(nullptr, this->m_jac_rows, this->m_jac_cols, dummy);
    }

    /** \return The NLP indices of the path constraints arguments [x_k, u_k] at node k. */
    std::vector<Integer> path_indices(Integer const k) const
    {
      std::vector<Integer> idx(N+U);
      for (Integer a{0}; a < N; ++a) {idx[a] = this->x_index(k) + a;}
      for (Integer a{0}; a < U; ++a) {idx[N + a] = this->u_index(std::min(k, this->m_n-1)) + a;}
      return idx;
    }

    /** \return The NLP indices of the boundary arguments [x_0, x_n]. */
    std::vector<Integer> boundary_indices() const
    {
      std::vector<Integer> idx(2*N);
      for (Integer a{0}; a < N; ++a) {idx[a] = this->x_index(0) + a; idx[N + a] = this->x_index(this->m_n) + a;}
      return idx;
    }

    /** Add the lower triangle of a dense block over the given (distinct) indices to the pattern. */
    void add_hessian_block(std::vector<Integer> const & idx)
    {
      for (Integer a : idx) {
        for (Integer b : idx) {
          if (a < b) {continue;}
          Index const key{a, b};
          if (this->m_hess_map.find(key) == this->m_hess_map.end()) {
            this->m_hess_map.emplace(key, static_cast<Integer>(this->m_hess_rows.size()));
            this->m_hess_rows.push_back(a);
            this->m_hess_cols.push_back(b);
          }
        }
      }
    }

    /** Build the Hessian sparsity pattern (lower triangle, duplicates merged). */
    void build_hessian_structure()
    {
      this->m_hess_rows.clear();
      this->m_hess_cols.clear();
      this->m_hess_map.clear();
      std::vector<Integer> idx(D);
      for (Integer k{0}; k < this->m_n; ++k) {
        for (Integer a{0}; a < D; ++a) {idx[a] = this->x_index(k) + a;}
        this->add_hessian_block(idx);
      }
      if (this->m_ng > 0) {this->add_hessian_block(this->path_indices(this->m_n));}
      this->add_hessian_block(this->boundary_indices());
    }

    /** Accumulate the lower triangle of a dense symmetric block over the given indices. */
    template <typename Matrix>
    void accumulate(std::vector<Integer> const & idx, Matrix const & H, VectorX & vals) const
    {
      Integer const size{static_cast<Integer>(idx.size())};
      for (Integer a{0}; a < size; ++a) {
        for (Integer b{0}; b < size; ++b) {
          if (idx[a] < idx[b]) {continue;}
          vals(this->m_hess_map.at(Index{idx[a], idx[b]})) += H(a, b);
        }
      }
    }

  }; // class Collocation

} // namespace Sandals

#endif // SANDALS_COLLOCATION_HH
