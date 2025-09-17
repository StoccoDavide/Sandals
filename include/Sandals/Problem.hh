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
  template <typename Real, Integer N, Integer M, typename Integrator>
  class Problem
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */
    const Real SQRT_EPSILON{std::sqrt(EPSILON)}; /**< Square root of machine epsilon epsilon static constant value. */ \

    using SystemPtr = typename Implicit<Real, N, M>::Pointer; /**< Shared pointer to an ODE/DAE system. */
    using IntegratorPtr = std::shared_ptr<Integrator>; /**< Shared pointer to a Runge-Kutta method. */
    using SolutionPtr = std::shared_ptr<Solution<Real, N, M>>; /**< Shared pointer to a solution. */

    using VectorX = typename Integrator::VectorX;
    using MatrixJX = typename Integrator::MatrixJX;
    using VectorF = typename Implicit<Real, N, M>::VectorF;
    using MatrixJF = typename Implicit<Real, N, M>::MatrixJF;

  private:
    std::string   m_name{"(undefined)"}; /**< Name of the problem. */
    SystemPtr     m_system; /**< ODE/DAE system. */
    IntegratorPtr m_integrator; /**< Runge-Kutta method. */
    SolutionPtr   m_solution; /**< Solution of the problem. */

    bool    m_verbose{false};                   /**< Verbose mode boolean. */
    Real    m_tolerance{EPSILON_HIGH}; /**< Tolerance for the solution. */
    Integer m_max_iterations{100};     /**< Maximum number of iterations. */

  public:
    /**
    * Class constructor for the boundary value problem (BVP).
    * \param[in] t_name The name of the problem.
    * \param[in] t_system The ODE/DAE system shared pointer.
    * \param[in] t_integrator The integrator shared pointer.
    */
    Problem(std::string t_name, SystemPtr t_system, IntegratorPtr t_integrator)
      : m_name(t_name), m_integrator(t_integrator)
    {
      this->m_integrator->system(t_system);
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
    SystemPtr system() {return this->m_integrator->system();}

    /**
    * Set the ODE/DAE system pointer.
    * \param[in] t_system The ODE/DAE system pointer.
    */
    void system(SystemPtr t_system) {this->m_integrator->system(t_system);}

    /**
    * Get the integrator pointer.
    * \return The integrator pointer.
    */
    IntegratorPtr integrator() {return this->m_integrator;}

    /**
    * Set the integrator pointer.
    * \param[in] t_integrator The integrator pointer.
    */
    void integrator(IntegratorPtr t_integrator) {this->m_integrator = t_integrator;}

    /**
    * Get the solution pointer.
    * \return The solution pointer.
    */
    SolutionPtr const solution() const {return this->m_solution;}

    /**
    * Get the verbose mode.
    * \return The verbose mode.
    */
    bool verbose_mode() {return this->m_verbose;}

    /**
    * Set the verbose mode.
    * \param[in] t_verbose The verbose mode.
    */
    void verbose_mode(bool t_verbose) {
      this->m_verbose = t_verbose;
      this->m_integrator->verbose_mode(t_verbose);
    }

    /**
    * Enable the verbose mode.
    */
    void enable_verbose_mode() {this->verbose_mode(true);}

    /**
    * Disable the verbose mode.
    */
    void disable_verbose_mode() {this->verbose_mode(false);}

    /**
    * Get the solution tolerance.
    * \return The solution tolerance.
    */
    Real tolerance() {return this->m_tolerance;}

    /**
    * Set the solution tolerance.
    * \param[in] t_tolerance The solution tolerance.
    */
    void tolerance(Real const t_tolerance)
    {this->m_tolerance = t_tolerance;}

    /**
    * Get the maximum number of iterations.
    * \return The maximum number of iterations.
    */
    Integer & max_iterations() {return this->m_max_iterations;}

    /**
    * Set the maximum number of iterations.
    * \param[in] t_max_iterations The maximum number of iterations.
    */
    void max_iterations(Integer const t_max_iterations)
    {this->m_max_iterations = t_max_iterations;}

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
    bool single_shooting(VectorX const & t_mesh, VectorF const & ics)
    {
      using ShootingF = Eigen::Matrix<Real, 2*N, 1>;
      using ShootingJF = Eigen::Matrix<Real, 2*N, 2*N>;

      #define CMD "Sandals::Problem::single_shooting(...): "

      // Temporary variables
      VectorF bcs, x_ini, x_end;
      ShootingF x_sol, x_step, b;
      ShootingJF A(ShootingJF::Zero());
      A.template block<N, N>(0, N) = MatrixJX::Identity();

      // Initialize the guess and the solution
      x_sol << ics, ics;
      this->m_solution->clear();

      // Solve the boundary value problem using a linearized Newton method
      MatrixJX Jx(MatrixJX::Identity());
      Eigen::FullPivLU<ShootingJF> lu;
      for (Integer i{0}; i < this->m_max_iterations; ++i) {

        /* Single shooting method
                  [A]           {x}   =          {f}
         /   -Jx       I     \ /    \   / x_end - x_sol_end \
         |                   | | dx | = |                   |
         \ Jb_x_ini Jb_x_end / \    /   \        -b         /
        */

        // Integrate the system and propagate the Jacobian
        if (!this->m_integrator->solve(t_mesh, x_sol.template head<N>(), *this->m_solution, Jx)) {
          SANDALS_ERROR(CMD "failed to integrate the system " << this->m_integrator->name() << ".");
          return false;
        }

        // Retrieve the initial and final states
        x_ini = x_sol.template head<N>();
        x_end = this->m_solution->x.col(this->m_solution->t.size() - 1);

        // Evaluate the residual of the boundary conditions
        bcs = this->b(x_ini, x_end);

        // Print the iteration info
        if (this->m_verbose) {
          std::cout << "Iteration " << i << ": |b| = " << bcs.norm() << std::endl
                    << "  x(" << t_mesh.template head<1>() << ") = " << x_ini.transpose() << std::endl
                    << "  x(" << t_mesh.template tail<1>() << ") = " << x_end.transpose() << std::endl;
        }

        // Check if the solution is found (i.e., if the boundary conditions are satisfied)
        if (bcs.norm() < this->m_tolerance) {return true;}

        // Build the linear system
        b.template head<N>() = x_end - x_sol.template tail<N>();
        b.template tail<N>() = -bcs;
        A.template block<N, N>(0, 0) = -Jx;
        A.template block<N, N>(N, 0) = this->Jb_x_ini(x_ini, x_end);
        A.template block<N, N>(N, N) = this->Jb_x_end(x_ini, x_end);

        // Compute the solution of the linear system
        lu.compute(A);
        SANDALS_ASSERT(lu.rank() == 2*N, CMD "singular Jacobian detected.");
        x_step = lu.solve(b);

        // Update the solution
        x_sol += x_step;
      }

      // If the loop completes without returning, indicate failure
      if (this->m_verbose) {SANDALS_WARNING(CMD "maximum number of iterations reached.");}
      return false;

      #undef CMD
    }

  }; // class Problem

} // namespace Sandals

#endif // SANDALS_PROBLEM_HH
