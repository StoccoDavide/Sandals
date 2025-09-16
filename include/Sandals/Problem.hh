/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_PROBLEM_HH
#define SANDALS_PROBLEM_HH

#include <Sandals/RungeKutta.hh>

namespace Sandals
{

  /*\
   |.  ____            _     _
   |. |  _ \ _ __ ___ | |__ | | ___ _ __ ___
   |. | |_) | '__/ _ \| '_ \| |/ _ \ '_ ` _ \
   |. |  __/| | | (_) | |_) | |  __/ | | | | |
   |. |_|   |_|  \___/|_.__/|_|\___|_| |_| |_|
   |.
  \*/

  /**
  * \brief Class container for an boundary value problem (BVP) of ODEs/DAEs.
  *
  * Class container for an boundary value problem (BVP) of ordinary differential equations (ODEs) or
  * differential algebraic equations (DAEs).
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the implicit ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M, Integer S>
  class Problem
  {
  public:
    using TimeSpan = std::pair<Real, Real>; /**< Pair of initial and final time. */
    using SystemPtr = std::shared_ptr<Implicit<Real, N, M>>; /**< Shared pointer to an ODE/DAE system. */
    using RungeKuttaPtr = std::shared_ptr<RungeKutta<Real, S, N, M>>; /**< Shared pointer to a Runge-Kutta method. */
    using SolutionPtr = std::shared_ptr<Solution<Real, N, M>>; /**< Shared pointer to a solution. */

    using typename RungeKutta<Real, S, N, M>::VectorX;
    using typename RungeKutta<Real, S, N, M>::MatrixJX;
    using typename Implicit<Real, N, M>::VectorF;
    using typename Implicit<Real, N, M>::MatrixJF;

  private:
    std::string   m_name{"(undefined)"}; /**< Name of the problem. */
    SystemPtr     m_system; /**< ODE/DAE system. */
    RungeKuttaPtr m_rungekutta; /**< Runge-Kutta method. */
    SolutionPtr   m_solution; /**< Solution of the problem. */

  public:
    /**
    * Class constructor for the boundary value problem (BVP).
    * \param[in] t_name The name of the problem.
    * \param[in] t_system The ODE/DAE system shared pointer.
    * \param[in] t_rungekutta The Runge-Kutta method shared pointer.
    */
    Problem(std::string t_name, SystemPtr t_system, RungeKuttaPtr t_rungekutta)
      : m_system(t_system), m_rungekutta(t_rungekutta), m_name(t_name)
    {
      this->m_solution = std::make_shared<Solution<Real, N, M>>();
    }

    /**
    * Class destructor.
    */
    virtual ~Problem() {}

    /**
    * Get the name of the problem.
    * \return The name of the problem.
    */
    std::string name() const {return this->m_name;}

    /**
    * Set the name of the problem.
    * \param[in] t_name The name of the problem.
    */
    void name(std::string const t_name) {this->m_name = t_name;}

    /**
    * Get the ODE/DAE system pointer.
    * \return The ODE/DAE system pointer.
    */
    SystemPtr system() {return this->m_system;}

    /**
    * Set the ODE/DAE system pointer.
    * \param[in] t_system The ODE/DAE system pointer.
    */
    void system(SystemPtr t_system) {this->m_system = t_system;}

    /**
    * Get the Runge-Kutta method pointer.
    * \return The Runge-Kutta method pointer.
    */
    RungeKuttaPtr rungekutta() {return this->m_rungekutta;}

    /**
    * Set the Runge-Kutta method pointer.
    * \param[in] t_rungekutta The Runge-Kutta method pointer.
    */
    void rungekutta(RungeKuttaPtr t_rungekutta) {this->m_rungekutta = t_rungekutta;}

    /**
    * Evaluate the boundary conditions function of the problem \f$ \mathbf{b}(\mathbf{x}_{\text{ini}},
    * \mathbf{x}_{\text{end}}) \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_{\text{ini}} \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_{\text{end}} \f$.
    * \return The residual of the boundary conditions.
    */
    virtual VectorF b(VectorF const & x_ini, VectorF const & x_end) const = 0;

    /**
    * Evaluate the Jacobian of the boundary conditions function \f$ \mathbf{b}(\mathbf{x}_{\text{ini}},
    * \mathbf{x}_{\text{end}}) \f$ with respect to the initial states \f$ \mathbf{x}_{\text{ini}} \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_{\text{ini}} \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_{\text{end}} \f$.
    * \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}_{\text{ini}}}(\mathbf{x}_{\text{ini}},
    * \mathbf{x}_{\text{end}}) \f$.
    */
    virtual MatrixJF Jb_x_ini(VectorF const & x_ini, VectorF const & x_end) const = 0;

    /**
    * Evaluate the Jacobian of the boundary conditions function \f$ \mathbf{b}(\mathbf{x}_{\text{ini}},
    * \mathbf{x}_{\text{end}}) \f$ with respect to the final states \f$ \mathbf{x}_{\text{end}} \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_{\text{ini}} \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_{\text{end}} \f$.
    * \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}_{\text{end}}}(\mathbf{x}_{\text{ini}},
    * \mathbf{x}_{\text{end}}) \f$.
    */
    virtual MatrixJF Jb_x_end(VectorF const & x_ini, VectorF const & x_end) const = 0;

    /**
    * Solve the boundary value problem (BVP) using the Runge-Kutta method uning a single shooting method.
    * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t} \f$.
    * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    * \param[out] sol The solution of the system over the mesh of independent variable.
    * \return True if the system is successfully solved, false otherwise.
    * \warning Do not use the solution for internal backtracking, as the step callback may directly
    * modify the solution.
    */
    bool single_shooting(VectorX const & t_mesh, VectorF const & ics) const
    {
      using Eigen::last;
      using F = Eigen::Matrix<Real, 2*N, 1>;
      using JF = Eigen::Matrix<Real, 2*N, 2*N>;

      #define CMD "Sandals::BVP::single_shooting(...): "

      // Temporary variables
      F x_sol, dx, fun;
      JF jac;

      // Initialize the guess and the solution
      x_sol << ics, ics;
      this->m_solution.clear();

      // Solve the boundary value problem using a linearized Newton method
      Integer max_iter{100};
      MatrixJX Jx; Jx.setIdentity();
      for (Integer i{0}; i < max_iter; ++i) {
        if (!this->m_rungekutta.solve(t_mesh, x_sol.template head<N>(), this->m_solution, Jx)) {
          SANDALS_ERROR(CMD "failed to integrate the system " << this->m_name << ".");
          return false;
        }

        // Evaluate the boundary conditions and its Jacobian
        fun.template head<N>() = this->m_solution.x.col(last) - x_sol.template tail<N>();
        fun.template tail<N>() = -this->b(x_sol, this->m_solution.x.col(last));
        jac.template block<N, N>(0, 0) = -Jx;
        jac.template block<N, N>(0, N) = JF::Identity();
        jac.template block<N, N>(N, 0) = this->Jb_x_ini(x_sol, this->m_solution.x.col(last));
        jac.template block<N, N>(N, N) = this->Jb_x_end(x_sol, this->m_solution.x.col(last));

        // Solve the linear system
        Eigen::FullPivLU<JF> lu(jac);
        SANDALS_ASSERT(lu.rank() == 2*N, CMD "singular Jacobian detected.");
        dx = -lu.solve(fun);

        // Update the guess
        x_sol += dx;
        if (dx.norm() < 0.01) {return true;}

      }
      return false;

      #undef CMD
    }

  }; // class BVP

} // namespace Sandals

#endif // SANDALS_PROBLEM_HH
