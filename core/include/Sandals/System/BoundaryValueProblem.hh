/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_BOUNDARY_VALUE_PROBLEM_HH
#define SANDALS_BOUNDARY_VALUE_PROBLEM_HH

// Sparse linear algebra solvers
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

// Sandals Runge-Kutta integrator
#include <Sandals/RungeKutta.hh>

namespace Sandals {

  /*\
   |   ______     ______
   |  | __ ) \   / /  _ \
   |  |  _ \\ \ / /| |_) |
   |  | |_) |\ V / |  __/
   |  |____/  \_/  |_|
   |
  \*/

  /**
   * \brief Class container for an boundary value problem (BVP) of ODEs/DAEs.
   *
   * Class container for an boundary value problem (BVP) of ordinary
   * differential equations (ODEs) or differential algebraic equations (DAEs).
   * \tparam Real The scalar number type.
   * \tparam N The dimension of the implicit ODE/DAE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M, typename Integrator>
  class BoundaryValueProblem {
   public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */
    const Real SQRT_EPSILON{
      std::sqrt(EPSILON)}; /**< Square root of machine epsilon epsilon static
                              constant value. */

    using ShootingChoice = enum class ShootingChoice : Integer {
      SINGLE   = 0,
      MULTIPLE = 1
    }; /**< Shooting method choice. */
    using ContinuationChoice = enum class ContinuationChoice : Integer {
      NONE   = 0,
      SIMPLE = 1
    }; /**< Continuation method choice. */

    using System =
        Implicit<Real, N, M>; /**< Unique pointer to an ODE/DAE system. */
    using SystemPtr =
        typename Implicit<Real, N, M>::Pointer;        /**< Unique pointer to an
                                                          ODE/DAE system. */
    using IntegratorPtr = std::unique_ptr<Integrator>; /**< Unique pointer to a
                                                          Runge-Kutta method. */
    using SolutionPtr =
        std::unique_ptr<Solution<Real, N, M>>;         /**< Unique pointer to a
                                                          solution. */

    using VectorX = Eigen::Vector<Real, Eigen::Dynamic>; /**< Dynamic vector of
                                                            real numbers. */
    using MatrixX =
        Eigen::Matrix<Real, N, Eigen::Dynamic>; /**< Dynamic matrix of real
                                                   numbers. */
    using MatrixJX = typename Integrator::MatrixJX;
    using VectorF  = typename Implicit<Real, N, M>::VectorF;
    using MatrixJF = typename Implicit<Real, N, M>::MatrixJF;
    using VectorH  = typename Implicit<Real, N, M>::VectorH;
    using MatrixJH = typename Implicit<Real, N, M>::MatrixJH;

   private:
    std::string m_name{"(undefined name)"}; /**< Name of the BVP. */
    IntegratorPtr m_integrator;             /**< Runge-Kutta method. */
    SolutionPtr m_solution;                 /**< Solution of the BVP. */

    bool m_verbose{false};                  /**< Verbose mode boolean. */
    Real m_tolerance{std::sqrt(EPSILON)};   /**< Tolerance for the solution. */
    Real m_sigma{0.01}; /**< Invariants manifold weight for the multiple
                           shooting method. */
    Integer m_max_iterations{100}; /**< Maximum number of iterations. */
    Integer m_subintervals{
      1}; /**< Number of subintervals for the shooting methods. */

   public:
    /**
     * Class constructor for the boundary value problem.
     * \param[in] t_name The name of the BVP.
     * \param[in] t_system The ODE/DAE system unique pointer.
     * \param[in] t_integrator The integrator unique pointer.
     */
    BoundaryValueProblem(std::string t_name,
                         SystemPtr t_system,
                         IntegratorPtr t_integrator)
        : m_name(t_name), m_integrator(std::move(t_integrator)) {
      this->m_integrator->system(std::move(t_system));
      this->m_solution = std::make_unique<Solution<Real, N, M>>();
    }

    /**
     * Class destructor.
     */
    virtual ~BoundaryValueProblem() {}

    /**
     * Get the name of the BVP.
     * \return The name of the BVP.
     */
    std::string name() const {
      return this->m_name;
    }

    /**
     * Set the name of the BVP.
     * \param[in] t_name The name of the BVP.
     */
    void name(const std::string t_name) {
      this->m_name = t_name;
    }

    /**
     * Get the ODE/DAE system pointer.
     * \return The ODE/DAE system pointer.
     */
    System *system() {
      return this->m_integrator->system();
    }

    /**
     * Get the ODE/DAE system const pointer.
     * \return The ODE/DAE system const pointer.
     */
    const System *system() const {
      return this->m_integrator->system();
    }

    /**
     * Set the ODE/DAE system pointer.
     * \param[in] t_system The ODE/DAE system pointer.
     */
    void system(SystemPtr t_system) {
      this->m_integrator->system(t_system);
    }

    /**
     * Get the integrator pointer.
     * \return The integrator pointer.
     */
    Integrator *integrator() {
      return this->m_integrator.get();
    }

    /**
     * Get the integrator const pointer.
     * \return The integrator const pointer.
     */
    const Integrator *integrator() const {
      return this->m_integrator.get();
    }

    /**
     * Set the integrator pointer.
     * \param[in] t_integrator The integrator pointer.
     */
    void integrator(IntegratorPtr t_integrator) {
      this->m_integrator = std::move(t_integrator);
    }

    /**
     * Get the solution pointer.
     * \return The solution pointer.
     */
    Solution<Real, N, M> &solution() {
      return *this->m_solution;
    }

    /**
     * Get the solution const pointer.
     * \return The solution const pointer.
     */
    const Solution<Real, N, M> &solution() const {
      return *this->m_solution;
    }

    /**
     * Get the verbose mode.
     * \return The verbose mode.
     */
    bool verbose_mode() {
      return this->m_verbose;
    }

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
    void enable_verbose_mode() {
      this->verbose_mode(true);
    }

    /**
     * Disable the verbose mode.
     */
    void disable_verbose_mode() {
      this->verbose_mode(false);
    }

    /**
     * Get the solution tolerance.
     * \return The solution tolerance.
     */
    Real tolerance() {
      return this->m_tolerance;
    }

    /**
     * Set the solution tolerance.
     * \param[in] t_tolerance The solution tolerance.
     */
    void tolerance(const Real t_tolerance) {
      this->m_tolerance = t_tolerance;
    }

    /**
     * Get the weight \f$ \sigma \f$ of the invariants manifold for the multiple
     * shooting method.
     * \return The weight \f$ \sigma \f$.
     */
    Real sigma() {
      return this->m_sigma;
    }

    /**
     * Set the weight \f$ \sigma \f$ of the invariants manifold for the multiple
     * \param[in] t_sigma The weight \f$ \sigma \f$ of the invariants manifold
     * for the multiple shooting method.
     */
    void sigma(const Real t_sigma) {
      this->m_sigma = t_sigma;
    }

    /**
     * Get the maximum number of iterations.
     * \return The maximum number of iterations.
     */
    Integer &max_iterations() {
      return this->m_max_iterations;
    }

    /**
     * Set the maximum number of iterations.
     * \param[in] t_max_iterations The maximum number of iterations.
     */
    void max_iterations(const Integer t_max_iterations) {
      this->m_max_iterations = t_max_iterations;
    }

    /**
     * Get the number of subintervals for the shooting methods.
     * \return The number of subintervals for the shooting methods.
     */
    Integer &subintervals() {
      return this->m_subintervals;
    }

    /**
     * Set the number of subintervals for the shooting methods.
     * \param[in] t_subintervals The number of subintervals for the shooting
     * methods.
     */
    void subintervals(const Integer t_subintervals) {
#define CMD "Sandals::BoundaryValueProblem::subintervals(...): "

      SANDALS_ASSERT(t_subintervals > 0,
                     CMD "number of subintervals must be positive.");
      this->m_subintervals = t_subintervals;

#undef CMD
    }

    /**
     * Evaluate the boundary conditions function of the problem \f$
     * \mathbf{b}(\mathbf{x}_{\text{ini}},
     * \mathbf{x}_{\text{end}}) \f$.
     * \param[in] x_ini Initial states \f$ \mathbf{x}_{\text{ini}} \f$.
     * \param[in] x_end Final states \f$ \mathbf{x}_{\text{end}} \f$.
     * \return The residual of the boundary conditions.
     */
    virtual VectorF b(const VectorF &x_ini, const VectorF &x_end) const = 0;

    /**
     * Evaluate the Jacobian of the boundary conditions function \f$
     * \mathbf{b}(\mathbf{x}_{\text{ini}},
     * \mathbf{x}_{\text{end}}) \f$ with respect to the initial states \f$
     * \mathbf{x}_{\text{ini}} \f$.
     * \param[in] x_ini Initial states \f$ \mathbf{x}_{\text{ini}} \f$.
     * \param[in] x_end Final states \f$ \mathbf{x}_{\text{end}} \f$.
     * \return The Jacobian \f$
     * \mathbf{Jb}_{\mathbf{x}_{\text{ini}}}(\mathbf{x}_{\text{ini}},
     * \mathbf{x}_{\text{end}}) \f$.
     */
    virtual MatrixJF Jb_x_ini(const VectorF &x_ini,
                              const VectorF &x_end) const = 0;

    /**
     * Evaluate the Jacobian of the boundary conditions function \f$
     * \mathbf{b}(\mathbf{x}_{\text{ini}},
     * \mathbf{x}_{\text{end}}) \f$ with respect to the final states \f$
     * \mathbf{x}_{\text{end}} \f$.
     * \param[in] x_ini Initial states \f$ \mathbf{x}_{\text{ini}} \f$.
     * \param[in] x_end Final states \f$ \mathbf{x}_{\text{end}} \f$.
     * \return The Jacobian \f$
     * \mathbf{Jb}_{\mathbf{x}_{\text{end}}}(\mathbf{x}_{\text{ini}},
     * \mathbf{x}_{\text{end}}) \f$.
     */
    virtual MatrixJF Jb_x_end(const VectorF &x_ini,
                              const VectorF &x_end) const = 0;

    /**
     * Solve the BVP using the single shooting method.
     * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t} \f$.
     * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
     * \param[out] sol The solution of the system over the mesh of independent
     * variable.
     * \return True if the system is successfully solved, false otherwise.
     * \warning Single shooting method is not recommended for stiff problems,
     * as it may fail to converge. Nonetheless, invariants preservation is not
     * implemented for the single shooting method.
     */
    bool single_shooting(const VectorX &t_mesh, const VectorF &ics) {
#define CMD "Sandals::BoundaryValueProblem::single_shooting(...): "

      using VectorShooting = Eigen::Vector<Real, 2 * N>;
      using MatrixShooting = Eigen::Matrix<Real, 2 * N, 2 * N>;

      // Disable adaptive time-stepping
      this->m_integrator->disable_adaptive_mode();

      // Temporary variables
      VectorF b, x_ini, x_end;

      // Shooting method scheme
      VectorShooting x_sol, x_step, b_sys;
      MatrixShooting A_sys;

      // Temporary variables
      const Integer local_intervals{std::max(1, this->m_subintervals)};
      const Integer num_intervals{static_cast<Integer>(t_mesh.size()) - 1};
      Integer idx_x_ini{0}, idx_x_end{num_intervals};
      if (this->m_integrator->reverse_mode()) {
        idx_x_ini = num_intervals;
        idx_x_end = 0;
      }

      SANDALS_ASSERT_WARNING(local_intervals == 1,
                             CMD
                             "please avoid using multiple subintervals with "
                             "single shooting method, just resize the "
                             "mesh accordingly.");

      // Prepare the linear system for the Newton step
      for (Integer i{0}; i < N; ++i) {
        A_sys.insert(i, i) = -1.0;
      }

      // Initialize the guess and the solution
      x_sol << ics, ics;
      this->m_solution->clear();

      // Solve the boundary value problem using a linearized Newton method
      MatrixJX Jx(MatrixJX::Identity());
      Eigen::FullPivHouseholderQR<MatrixShooting> qr;
      for (Integer iter{0}; iter < this->m_max_iterations; ++iter) {
        /* Single shooting method scheme
                [A_sys]          {x}  =       {b_sys}
         /   -Jx       I     \ /    \   / x_end - x_sol_end \
         |                   | | dx | = |                   |
         \ Jb_x_ini Jb_x_end / \    /   \        -b         /
        */

        // Integrate the system and propagate the Jacobian
        if (!this->m_integrator->solve(t_mesh,
                                       x_sol.template head<N>(),
                                       *this->m_solution,
                                       Jx)) {
          SANDALS_ERROR(CMD "failed to integrate the system "
                        << this->m_integrator->name() << ".");
          return false;
        }

        // Retrieve the initial and final states
        x_ini = this->m_solution->x.col(idx_x_ini);
        x_end = this->m_solution->x.col(idx_x_end);

        // Evaluate the residual of the boundary conditions
        b = this->b(x_ini, x_end);

        // Print the iteration info
        if (this->m_verbose) {
          std::cout << "Iteration " << iter << ": |b| = " << b.norm()
                    << std::endl
                    << "  x(" << t_mesh(idx_x_ini)
                    << ") = " << x_ini.transpose() << std::endl
                    << "  x(" << t_mesh(idx_x_end)
                    << ") = " << x_end.transpose() << std::endl;
        }

        // Check convergence
        if (b.norm() < this->m_tolerance) {
          return true;
        }

        // Build the linear system
        A_sys.template block<N, N>(0, 0) = -Jx;
        if (this->m_integrator->reverse_mode()) {
          A_sys.template block<N, N>(N, N) = this->Jb_x_ini(x_ini, x_end);
          A_sys.template block<N, N>(N, 0) = this->Jb_x_end(x_ini, x_end);
        } else {
          A_sys.template block<N, N>(N, 0) = this->Jb_x_ini(x_ini, x_end);
          A_sys.template block<N, N>(N, N) = this->Jb_x_end(x_ini, x_end);
        }
        b_sys.template head<N>() = x_end - x_sol.template tail<N>();
        b_sys.template tail<N>() = -b;

        // Compute the solution of the linear system
        A_sys.makeCompressed();
        qr.compute(A_sys);
        SANDALS_ASSERT(qr.rank() == qr.cols(),
                       CMD "singular Jacobian detected.");
        x_sol += qr.solve(b_sys);
      }

      // If the loop completes without returning, indicate failure
      if (this->m_verbose) {
        SANDALS_WARNING(CMD "maximum number of iterations reached.");
      }
      return false;

#undef CMD
    }

    /**
     * Solve the BVP using the multiple shooting method.
     * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t} \f$.
     * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
     * \param[in] x_guess Initial guess for the states at the mesh nodes.
     * \return True if the system is successfully solved, false otherwise.
     */
    bool multiple_shooting(const VectorX &t_mesh, const MatrixX &x_guess) {
#define CMD "Sandals::BoundaryValueProblem::multiple_shooting(...): "

      using VectorShooting = Eigen::Vector<Real, Eigen::Dynamic>;
      using MatrixShooting = Eigen::SparseMatrix<Real>;

      // Check if the mesh and the guess have compatible sizes
      SANDALS_ASSERT(t_mesh.size() > 1,
                     CMD "expected at least two time mesh points.");
      SANDALS_ASSERT(x_guess.cols() == t_mesh.size(),
                     CMD
                     "incompatible sizes between time mesh and states guess.");

      // Disable adaptive time-stepping
      this->m_integrator->disable_adaptive_mode();

      // Temporary variables
      VectorF x_ini, x_end;
      const Integer num_intervals{static_cast<Integer>(t_mesh.size()) - 1};
      const Integer local_intervals{std::max(1, this->m_subintervals)};
      const Integer c_size{num_intervals * N};
      const Integer x_size{(num_intervals + 1) * N};
      const Integer h_size{
        this->m_integrator->projection_mode() ? 0 : (num_intervals + 1) * M};
      const Real sqrt_sigma{std::sqrt(this->m_sigma)};
      Integer idx_x_ini{0}, idx_x_end{num_intervals};
      if (this->m_integrator->reverse_mode()) {
        idx_x_ini = num_intervals;
        idx_x_end = 0;
      }

      // Prepare the linear system for the Newton step
      VectorShooting b_sys(c_size + h_size + N);
      MatrixShooting A_sys(x_size + h_size, x_size);
      A_sys.reserve(c_size * (N + 1) + h_size * N + 2 * N * N);

      // Initialize the solution
      this->m_solution->clear();
      this->m_solution->resize(t_mesh.size());
      this->m_solution->t = t_mesh;

      // Initial guess for the states as the guess at the mesh nodes
      MatrixX x_sol(x_guess), delta_x_sol(x_guess);

      // Solve the boundary value problem using a linearized Newton method
      VectorX t_local_mesh;
      Solution<Real, N, M> local_sol;
      MatrixJX Jx;
      MatrixJH Jh;
      Eigen::SparseQR<MatrixShooting, Eigen::COLAMDOrdering<Integer>> qr;
      std::vector<Eigen::Triplet<Real>> triplets;
      triplets.reserve(c_size * (N + 1) + h_size * N + 2 * N * N);
      for (Integer iter{0}; iter < this->m_max_iterations; ++iter) {
        /* Multiple shooting method scheme
                [A_sys]                   {x}  =         {b_sys}
         /   -Jx_1     I              \           /  x_ini_2 - x_sol_1  \
         |       .        .           |           |          :          |
         |          .        .        | /  :  \ = |          :          |
         |             -Jx_n     I    | |  :  |   | x_ini_n+1 - x_sol_n |
         |  σ½*Jh_0    0              | | dx  |   |      -σ^½*h_0       |
         |       .        .           | |  :  |   |          :          |
         |          .        .        | \  :  /   |          :          |
         |         σ½*Jh_m       0    |           |     -σ^½*h_m        |
         \ Jb_x_ini          Jb_x_end /           \         -b          /
        */

        // Reset the linear system
        triplets.clear();
        A_sys.setZero();
        b_sys.setZero();

        // Integrate each interval with the local mesh
        for (Integer k{0}; k < num_intervals; ++k) {
          t_local_mesh =
              VectorX::LinSpaced(local_intervals + 1, t_mesh(k), t_mesh(k + 1));
          Jx.setIdentity();
          if (!this->m_integrator->template solve<true>(t_local_mesh,
                                                        x_sol.col(k),
                                                        local_sol,
                                                        Jx)) {
            SANDALS_ERROR(CMD "failed to integrate interval " << k << ".");
            return false;
          }

          // Store the local solution
          if (k == 0) {
            this->m_solution->x.col(0) = local_sol.x.col(0);
            this->m_solution->h.col(0) = local_sol.h.col(0);
          }
          this->m_solution->x.col(k + 1) = local_sol.x.col(local_intervals);
          this->m_solution->h.col(k + 1) = local_sol.h.col(local_intervals);

          // Jacobian propagation for the current interval
          for (Integer i{0}; i < N; ++i) {
            triplets.emplace_back(k * N + i, (k + 1) * N + i, 1.0);
            for (Integer j{0}; j < N; ++j) {
              triplets.emplace_back(k * N + i, k * N + j, -Jx(i, j));
            }
          }

          // Continuity residual for interior nodes
          b_sys.template segment<N>(k * N) =
              this->m_solution->x.col(k + 1) - x_sol.col(k + 1);

          // Compute the Jacobian of contraints manifold
          if constexpr (M > 0) {
            if (!this->m_integrator->projection_mode()) {
              Jh = this->m_integrator->system()->Jh_x(
                  this->m_solution->x.col(k + 1),
                  t_mesh(k + 1));
              // Insert the Jacobian in the linear system
              for (Integer i{0}; i < M; ++i) {
                for (Integer j{0}; j < N; ++j) {
                  triplets.emplace_back(c_size + k * M + i,
                                        k * N + j,
                                        sqrt_sigma * Jh(i, j));
                }
              }

              // Residual for the constraints manifold
              b_sys.template segment<M>(c_size + k * M) =
                  -sqrt_sigma * this->m_solution->h.col(k + 1);
            }
          }
        }

        // Update the boundary condition states
        x_ini = this->m_solution->x.col(idx_x_ini);
        x_end = this->m_solution->x.col(idx_x_end);

        // Boundary condition residuals
        b_sys.template tail<N>() = -this->b(x_ini, x_end);

        // Print the iteration info
        if (this->m_verbose) {
          std::cout << "Iteration " << iter << ": |b| = " << b_sys.norm()
                    << std::endl
                    << "  x(" << t_mesh(idx_x_ini)
                    << ") = " << x_ini.transpose() << std::endl
                    << "  x(" << t_mesh(idx_x_end)
                    << ") = " << x_end.transpose() << std::endl;
        }

        // Check convergence
        if (b_sys.norm() < this->m_tolerance) {
          return true;
        }

        // Update the boundary condition Jacobian blocks
        MatrixJX Jb_x_ini(this->Jb_x_ini(x_ini, x_end));
        MatrixJX Jb_x_end(this->Jb_x_end(x_ini, x_end));
        for (Integer i{0}; i < N; ++i) {
          for (Integer j{0}; j < N; ++j) {
            triplets.emplace_back(c_size + h_size + i,
                                  idx_x_ini * N + j,
                                  Jb_x_ini(i, j));
            triplets.emplace_back(c_size + h_size + i,
                                  idx_x_end * N + j,
                                  Jb_x_end(i, j));
          }
        }

        // Solve the linear system
        A_sys.setFromTriplets(triplets.begin(), triplets.end());
        A_sys.makeCompressed();
        qr.compute(A_sys);
        SANDALS_ASSERT(qr.info() == Eigen::Success,
                       CMD "failed to factorize the linear system.");

        // Update the solution
        delta_x_sol = qr.solve(b_sys).reshaped(N, num_intervals + 1);
        if (!delta_x_sol.allFinite()) {
          SANDALS_ERROR(CMD "invalid solution of the linear system.");
          return false;
        } else if (delta_x_sol.norm() < this->m_tolerance * this->m_tolerance) {
          SANDALS_WARNING(CMD "small update step, possible convergence.");
          return false;
        }
        x_sol += delta_x_sol;
      }

      // If the loop completes without returning, indicate failure
      if (this->m_verbose) {
        SANDALS_WARNING(CMD "maximum number of iterations reached.");
      }
      return false;

#undef CMD
    }

    /**
     * Solve the boundary value problem using a shooting method.
     * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t}
     * \f$.
     * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
     * \param[in] x_guess Initial guess for the states at the mesh nodes (only
     * for multiple shooting).
     * \tparam ShootingType Type of shooting method to use (single or
     * multiple).
     */
    template <ShootingChoice ShootingType = ShootingChoice::MULTIPLE>
    bool solve(const VectorX &t_mesh, const VectorF &ics, MatrixX &x_guess) {
#define CMD "Sandals::BoundaryValueProblem::solve(...): "

      if constexpr (ShootingType == ShootingChoice::SINGLE) {
        return this->single_shooting(t_mesh, ics);
      } else if constexpr (ShootingType == ShootingChoice::MULTIPLE) {
        return this->multiple_shooting(t_mesh, x_guess);
      } else {
        SANDALS_ERROR(CMD "unknown shooting method.");
        return false;
      }

#undef CMD
    }

  };  // BoundaryValueProblem

}  // namespace Sandals

#endif  // SANDALS_BOUNDARY_VALUE_PROBLEM_HH
