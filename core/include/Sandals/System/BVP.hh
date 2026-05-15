/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_BVP_HH
#define SANDALS_BVP_HH

// Sparse linear algebra solvers
#include <Eigen/SparseQR>

// Sandals Runge-Kutta integrator
#include <Sandals/RungeKutta.hh>

namespace Sandals {

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
  class BVP {
   public:
    SANDALS_BASIC_CONSTANTS(Real)
    const Real SQRT_EPSILON{std::sqrt(EPSILON)};

    using ShootingChoice = enum class ShootingChoice : Integer {
      SINGLE   = 0,
      MULTIPLE = 1
    }; /**< Shooting method choice. */

    using SolutionChoice = enum class SolutionChoice : Integer {
      GN_STANDARD  = 0,
      GN_AUGMENTED = 1,
      GN_KKT       = 2
    }; /**< Solution method choice. */

    using System        = Implicit<Real, N, M>;
    using SystemPtr     = typename Implicit<Real, N, M>::Pointer;
    using IntegratorPtr = std::unique_ptr<Integrator>;
    using SolutionPtr   = std::unique_ptr<Solution<Real, N, M>>;

    using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>;
    using MatrixX  = Eigen::Matrix<Real, N, Eigen::Dynamic>;
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
    Real m_sigma{1.0};  /**< Invariants manifold weight for the multiple
     shooting method. */
    Real m_lambda{0.0}; /**< Augmentation parameter for the least squares
                           multiple shooting method. */
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
    BVP(std::string t_name, SystemPtr t_system, IntegratorPtr t_integrator)
        : m_name(t_name), m_integrator(std::move(t_integrator)) {
      this->m_integrator->system(std::move(t_system));
      this->m_solution = std::make_unique<Solution<Real, N, M>>();
    }

    /**
     * Class destructor.
     */
    virtual ~BVP() {}

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
     * Get the augmentation parameter \f$ \lambda \f$ for the least squares
     * multiple shooting method.
     * \return The augmentation parameter \f$ \lambda \f$.
     */
    Real lambda() {
      return this->m_lambda;
    }

    /**
     * Set the augmentation parameter \f$ \lambda \f$ for the least squares
     * multiple shooting method.
     * \param[in] t_lambda The augmentation parameter \f$ \lambda \f$ for the
     * least squares multiple shooting method.
     */
    void lambda(const Real t_lambda) {
      this->m_lambda = t_lambda;
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
#define CMD "Sandals::BVP::subintervals(...): "

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
#define CMD "Sandals::BVP::single_shooting(...): "

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
      MatrixJF Jx(MatrixJF::Identity());
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
                       CMD "singular linear system in BVP solution.");
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
     * \tparam SolutionMethod The solution method choice.
     * \param[in] t_mesh Independent variable mesh.
     * \param[in] x_guess Initial guess for the states at the mesh nodes.
     * \return True if the system is successfully solved, false otherwise.
     */
    template <SolutionChoice SolutionMethod = SolutionChoice::GN_KKT>
    bool multiple_shooting_damped(const VectorX &t_mesh,
                                  const MatrixX &x_guess) {
#define CMD "Sandals::BVP::multiple_shooting(...): "

      using VectorShooting = Eigen::Vector<Real, Eigen::Dynamic>;
      using MatrixShooting = Eigen::SparseMatrix<Real>;

      // Check input dimensions
      SANDALS_ASSERT(t_mesh.size() > 1,
                     CMD "expected at least two mesh points.");
      SANDALS_ASSERT(x_guess.cols() == t_mesh.size(),
                     CMD "incompatible mesh and guess sizes.");

      // Disable adaptive time-stepping
      this->m_integrator->disable_adaptive_mode();

      // Problem dimensions
      const Integer num_intervals{static_cast<Integer>(t_mesh.size()) - 1};
      const Integer local_intervals{std::max(1, this->m_subintervals)};
      const Integer c_size{num_intervals * N};
      const Integer x_size{(num_intervals + 1) * N};
      const Integer h_size{
        this->m_integrator->projection_mode() ? 0 : (num_intervals + 1) * M};
      const Integer lsq_rows{c_size + h_size};
      const Integer lsq_cols{x_size};
      const Integer A_lsq_nnz{(c_size + h_size + 2 * N) * N + c_size};
      const Real sqrt_sigma{std::sqrt(this->m_sigma)};

      Integer idx_x_ini{0};
      Integer idx_x_end{num_intervals};

      if (this->m_integrator->reverse_mode()) {
        idx_x_ini = num_intervals;
        idx_x_end = 0;
      }

      // Least-squares operator
      MatrixShooting A_lsq(lsq_rows, lsq_cols);
      A_lsq.reserve(A_lsq_nnz);

      VectorShooting b_lsq(lsq_rows);

      std::vector<Eigen::Triplet<Real>> triplets_lsq;
      triplets_lsq.reserve(A_lsq_nnz);

      // Boundary constraint matrix
      MatrixX Bmat(MatrixX::Zero(N, x_size));

      // Augmented KKT system
      MatrixShooting A_aug;
      VectorShooting b_aug;
      std::vector<Eigen::Triplet<Real>> triplets_aug;

      // Normal-equation KKT system
      MatrixShooting A_kkt;
      VectorShooting b_kkt;
      std::vector<Eigen::Triplet<Real>> triplets_kkt;

      // Allocate augmented KKT system
      if constexpr (SolutionMethod == SolutionChoice::GN_AUGMENTED) {
        // Symmetric augmented KKT system
        // / λI  Aᵀ  Bᵀ \
        // |  A  -I   0 |
        // \  B   0   0 /

        const Integer aug_size{x_size + lsq_rows + N};

        A_aug.resize(aug_size, aug_size);
        A_aug.reserve(2 * A_lsq_nnz + 4 * N * N + aug_size);

        b_aug.resize(aug_size);

        triplets_aug.reserve(2 * A_lsq_nnz + 4 * N * N + aug_size);
      }

      // Allocate normal-equation KKT system
      if constexpr (SolutionMethod == SolutionChoice::GN_KKT) {
        // Symmetric normal-equation KKT system
        // / AᵀA  Bᵀ \
        // \  B    0 /

        const Integer kkt_size{x_size + N};

        A_kkt.resize(kkt_size, kkt_size);
        A_kkt.reserve(4 * A_lsq_nnz + 4 * N * N);

        b_kkt.resize(kkt_size);

        triplets_kkt.reserve(4 * A_lsq_nnz + 4 * N * N);
      }

      // Initialize solution container
      this->m_solution->clear();
      this->m_solution->resize(t_mesh.size());
      this->m_solution->t = t_mesh;

      // State vectors
      MatrixX x_sol(x_guess);
      MatrixX delta_x_sol(MatrixX::Zero(N, num_intervals + 1));
      MatrixX x_trial(MatrixX::Zero(N, num_intervals + 1));

      // Temporary variables
      VectorF x_ini, x_end;
      VectorX t_local_mesh;

      Solution<Real, N, M> local_sol;

      MatrixJF Jx;
      MatrixJF Jb_x_ini;
      MatrixJF Jb_x_end;

      MatrixJH Jh;

      // Solvers
      Eigen::SparseQR<MatrixShooting, Eigen::COLAMDOrdering<Integer>> qr;
      // Eigen::SparseQR<MatrixShooting, Eigen::COLAMDOrdering<Integer>> ldlt;
      Eigen::SimplicialLDLT<MatrixShooting> ldlt;

      // Merit function
      auto residual_norm = [&](const MatrixX &x_eval) -> Real {
        Real norm2{0};

        Solution<Real, N, M> sol_eval;
        MatrixJF Jtmp;

        for (Integer k{0}; k < num_intervals; ++k) {
          VectorX t_eval =
              VectorX::LinSpaced(local_intervals + 1, t_mesh(k), t_mesh(k + 1));

          Jtmp.setIdentity();

          if (!this->m_integrator->template solve<true>(t_eval,
                                                        x_eval.col(k),
                                                        sol_eval,
                                                        Jtmp)) {
            return std::numeric_limits<Real>::infinity();
          }

          VectorF continuity =
              sol_eval.x.col(local_intervals) - x_eval.col(k + 1);

          norm2 += continuity.squaredNorm();

          if constexpr (M > 0) {
            if (!this->m_integrator->projection_mode()) {
              norm2 += this->m_sigma *
                       this->m_integrator->system()
                           ->h(sol_eval.x.col(local_intervals), t_mesh(k + 1))
                           .squaredNorm();

              if (k == 0) {
                norm2 += this->m_sigma * this->m_integrator->system()
                                             ->h(x_eval.col(0), t_mesh(0))
                                             .squaredNorm();
              }
            }
          }
        }

        norm2 +=
            this->b(x_eval.col(idx_x_ini), x_eval.col(idx_x_end)).squaredNorm();

        return 0.5 * norm2;
      };

      // Initial residual norm
      Real residual{std::numeric_limits<Real>::infinity()};

      // Newton iterations
      for (Integer iter{0}; iter < this->m_max_iterations; ++iter) {
        // Reset least-squares system
        triplets_lsq.clear();
        b_lsq.setZero();
        Bmat.setZero();

        // Assemble least-squares operator
        // A_lsq * dx = b_lsq
        for (Integer k{0}; k < num_intervals; ++k) {
          // Local integration mesh
          t_local_mesh =
              VectorX::LinSpaced(local_intervals + 1, t_mesh(k), t_mesh(k + 1));

          // Reset state transition Jacobian
          Jx.setIdentity();

          // Integrate interval
          if (!this->m_integrator->template solve<true>(t_local_mesh,
                                                        x_sol.col(k),
                                                        local_sol,
                                                        Jx)) {
            SANDALS_ERROR(CMD "failed to integrate interval " << k << ".");
            return false;
          }

          // Store local solution
          if (k == 0) {
            this->m_solution->x.col(0) = local_sol.x.col(0);
            this->m_solution->h.col(0) = local_sol.h.col(0);
          }

          this->m_solution->x.col(k + 1) = local_sol.x.col(local_intervals);
          this->m_solution->h.col(k + 1) = local_sol.h.col(local_intervals);

          // Continuity equations
          for (Integer i{0}; i < N; ++i) {
            // +I block
            triplets_lsq.emplace_back(k * N + i, (k + 1) * N + i, 1.0);

            // -Jx block
            for (Integer j{0}; j < N; ++j) {
              if (std::abs(Jx(i, j)) > EPSILON) {
                triplets_lsq.emplace_back(k * N + i, k * N + j, -Jx(i, j));
              }
            }
          }

          // Continuity residual
          b_lsq.template segment<N>(k * N) =
              this->m_solution->x.col(k + 1) - x_sol.col(k + 1);

          // Constraint manifold equations
          if constexpr (M > 0) {
            if (!this->m_integrator->projection_mode()) {
              // Initial node
              if (k == 0) {
                Jh =
                    sqrt_sigma *
                    this->m_integrator->system()->Jh_x(x_sol.col(0), t_mesh(0));

                for (Integer i{0}; i < M; ++i) {
                  for (Integer j{0}; j < N; ++j) {
                    if (std::abs(Jh(i, j)) > EPSILON) {
                      triplets_lsq.emplace_back(c_size + i, j, Jh(i, j));
                    }
                  }
                }

                b_lsq.template segment<M>(c_size) =
                    -sqrt_sigma *
                    this->m_integrator->system()->h(x_sol.col(0), t_mesh(0));
              }

              // Interior/final nodes
              Jh = sqrt_sigma *
                   this->m_integrator->system()->Jh_x(x_sol.col(k + 1),
                                                      t_mesh(k + 1));

              for (Integer i{0}; i < M; ++i) {
                for (Integer j{0}; j < N; ++j) {
                  if (std::abs(Jh(i, j)) > EPSILON) {
                    triplets_lsq.emplace_back(c_size + (k + 1) * M + i,
                                              (k + 1) * N + j,
                                              Jh(i, j));
                  }
                }
              }

              // Constraint residual
              b_lsq.template segment<M>(c_size + (k + 1) * M) =
                  -sqrt_sigma *
                  this->m_integrator->system()->h(x_sol.col(k + 1),
                                                  t_mesh(k + 1));
            }
          }
        }

        // Boundary states
        x_ini = x_sol.col(idx_x_ini);
        x_end = x_sol.col(idx_x_end);

        // Boundary Jacobians
        Jb_x_ini = this->Jb_x_ini(x_ini, x_end);
        Jb_x_end = this->Jb_x_end(x_ini, x_end);

        // Assemble boundary matrix
        for (Integer i{0}; i < N; ++i) {
          for (Integer j{0}; j < N; ++j) {
            Bmat(i, idx_x_ini * N + j) += Jb_x_ini(i, j);
            Bmat(i, idx_x_end * N + j) += Jb_x_end(i, j);
          }
        }

        // Build least-squares operator
        A_lsq.setFromTriplets(triplets_lsq.begin(), triplets_lsq.end());
        A_lsq.makeCompressed();

        // Residual norm
        residual = std::sqrt(2.0 * residual_norm(x_sol));

        // Iteration information
        if (this->m_verbose) {
          std::cout << "Newton iteration " << iter << ": |r| = " << residual
                    << std::endl
                    << "  x(" << t_mesh(idx_x_ini)
                    << ") = " << x_ini.transpose() << std::endl
                    << "  x(" << t_mesh(idx_x_end)
                    << ") = " << x_end.transpose() << std::endl;
        }

        // Convergence check
        if (residual < this->m_tolerance) {
          return true;
        }

        // Standard least-squares solve
        if constexpr (SolutionMethod == SolutionChoice::GN_STANDARD) {
          MatrixShooting A_full(lsq_rows + N, x_size);

          std::vector<Eigen::Triplet<Real>> triplets_full;
          triplets_full.reserve(A_lsq.nonZeros() + 2 * N * N);

          for (Integer k{0}; k < A_lsq.outerSize(); ++k) {
            for (typename MatrixShooting::InnerIterator it(A_lsq, k); it;
                 ++it) {
              triplets_full.emplace_back(it.row(), it.col(), it.value());
            }
          }

          for (Integer i{0}; i < N; ++i) {
            for (Integer j{0}; j < x_size; ++j) {
              if (std::abs(Bmat(i, j)) > EPSILON) {
                triplets_full.emplace_back(lsq_rows + i, j, Bmat(i, j));
              }
            }
          }

          A_full.setFromTriplets(triplets_full.begin(), triplets_full.end());

          VectorShooting b_full(lsq_rows + N);

          b_full.head(lsq_rows) = b_lsq;
          b_full.tail(N)        = -this->b(x_ini, x_end);

          qr.compute(A_full);

          SANDALS_ASSERT(qr.info() == Eigen::Success,
                         CMD "failed to factorize least-squares system.");

          delta_x_sol = qr.solve(b_full).reshaped(N, num_intervals + 1);
        }

        // Symmetric augmented KKT solve
        else if constexpr (SolutionMethod == SolutionChoice::GN_AUGMENTED) {
          triplets_aug.clear();
          b_aug.setZero();

          // Unknown ordering: (dx, z, lambda)^T
          const Integer dx_offset{0};
          const Integer z_offset{x_size};
          const Integer lambda_offset{x_size + lsq_rows};

          // A and Aᵀ blocks
          for (Integer k{0}; k < A_lsq.outerSize(); ++k) {
            for (typename MatrixShooting::InnerIterator it(A_lsq, k); it;
                 ++it) {
              // Aᵀ block
              triplets_aug.emplace_back(dx_offset + it.col(),
                                        z_offset + it.row(),
                                        it.value());
              // A block
              triplets_aug.emplace_back(z_offset + it.row(),
                                        dx_offset + it.col(),
                                        it.value());
            }
          }

          // -I block
          for (Integer i{0}; i < lsq_rows; ++i) {
            triplets_aug.emplace_back(z_offset + i, z_offset + i, -1.0);
          }

          // B and Bᵀ blocks
          for (Integer i{0}; i < N; ++i) {
            for (Integer j{0}; j < x_size; ++j) {
              if (std::abs(Bmat(i, j)) > EPSILON) {
                // Bᵀ block
                triplets_aug.emplace_back(dx_offset + j,
                                          lambda_offset + i,
                                          Bmat(i, j));
                // B block
                triplets_aug.emplace_back(lambda_offset + i,
                                          dx_offset + j,
                                          Bmat(i, j));
              }
            }
          }

          // RHS: (0, r, -b)^T
          b_aug.segment(z_offset, lsq_rows) = b_lsq;
          b_aug.tail(N)                     = -this->b(x_ini, x_end);

          // Build system
          A_aug.setFromTriplets(triplets_aug.begin(), triplets_aug.end());
          A_aug.makeCompressed();

          // Factorize
          ldlt.compute(A_aug);

          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD "failed to factorize augmented KKT system.");

          // Solve
          VectorShooting sol_aug(ldlt.solve(b_aug));

          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD "failed to solve augmented KKT system.");

          // Extract Newton step
          delta_x_sol =
              sol_aug.segment(dx_offset, x_size).reshaped(N, num_intervals + 1);
        }

        // Normal-equation KKT solve
        else if constexpr (SolutionMethod == SolutionChoice::GN_KKT) {
          triplets_kkt.clear();
          b_kkt.setZero();

          // Normal equations
          MatrixShooting AtA(A_lsq.transpose() * A_lsq);

          // AᵀA block
          for (Integer k{0}; k < AtA.outerSize(); ++k) {
            for (typename MatrixShooting::InnerIterator it(AtA, k); it; ++it) {
              triplets_kkt.emplace_back(it.row(), it.col(), it.value());
            }
          }

          // B and Bᵀ blocks
          for (Integer i{0}; i < N; ++i) {
            for (Integer j{0}; j < x_size; ++j) {
              if (std::abs(Bmat(i, j)) > EPSILON) {
                triplets_kkt.emplace_back(j, x_size + i, Bmat(i, j));

                triplets_kkt.emplace_back(x_size + i, j, Bmat(i, j));
              }
            }
          }

          // RHS: (Aᵀr - λr_bar, -b)^T
          b_kkt.head(x_size) =
              A_lsq.transpose() * b_lsq;  // - this->m_lambda * r_bar;
          b_kkt.tail(N) = -this->b(x_ini, x_end);

          // Build KKT matrix
          A_kkt.setFromTriplets(triplets_kkt.begin(), triplets_kkt.end());
          A_kkt.makeCompressed();

          // Factorize
          ldlt.compute(A_kkt);

          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD "failed to factorize KKT system.");

          // Solve
          VectorShooting sol_kkt(ldlt.solve(b_kkt));

          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD "failed to solve KKT system.");

          // Extract Newton step
          delta_x_sol = sol_kkt.head(x_size).reshaped(N, num_intervals + 1);
        }

        // Invalid step
        if (!delta_x_sol.allFinite()) {
          SANDALS_ERROR(CMD "invalid Newton update.");
          return false;
        }

        // Small step
        const Real step_norm{delta_x_sol.norm()};
        const Real x_norm{x_sol.norm()};
        constexpr Real stagnation_factor{10.0};
        if (step_norm < this->m_tolerance * (1.0 + x_norm)) {
          // Accept if residual already reasonably small
          if (residual < stagnation_factor * this->m_tolerance) {
            if (this->m_verbose) {
              std::cout << "Converged by small step criterion." << std::endl;
            }
            return true;
          }
          SANDALS_WARNING(CMD "Newton stagnation detected.");
          return false;
        }

        // Backtracking line search
        Real alpha{1.0};
        constexpr Real alpha_min{1.0e-6};
        constexpr Real contraction{0.5};
        constexpr Real armijo{1.0e-4};

        const Real merit_old{residual_norm(x_sol)};
        const VectorShooting gradient{-A_lsq.transpose() * b_lsq};
        const Real gradient_norm{gradient.norm()};

        if (residual < this->m_tolerance &&
            gradient_norm < 10.0 * this->m_tolerance) {
          if (this->m_verbose) {
            std::cout << "Converged: residual + gradient criteria satisfied.\n";
          }
          return true;
        }

        const VectorShooting delta_vec{
          Eigen::Map<VectorShooting>(delta_x_sol.data(), x_size)};
        const Real directional_derivative{gradient.dot(delta_vec)};

        // if (directional_derivative >= 0.0) {
        //   SANDALS_WARNING(CMD "non-descent direction detected.");
        //   return false;
        // }

        bool accepted{false};
        Integer iter_ls{0};
        while (alpha > alpha_min) {
          x_trial = x_sol + alpha * delta_x_sol;

          const Real merit_new{residual_norm(x_trial)};

          // Iteration information
          if (this->m_verbose) {
            std::cout << "  Line search iteration " << iter_ls
                      << ": (alpha = " << alpha << "), merit = " << merit_new
                      << std::endl;
          }

          const Real armijo_rhs{merit_old +
                                armijo * alpha * directional_derivative};

          if (std::isfinite(merit_new) && merit_new <= armijo_rhs) {
            if (this->m_verbose) {
              std::cout << "  Accepted step with alpha = " << alpha
                        << std::endl;
            }
            // Accept the trial step
            accepted = true;
            break;
          }
          alpha *= contraction;
          ++iter_ls;
        }

        if (!accepted) {
          SANDALS_WARNING(CMD "line search failed.");
          return false;
        }

        // Update solution
        x_sol = x_trial;
      }

      // Maximum iterations reached
      if (this->m_verbose) {
        SANDALS_WARNING(CMD "maximum iterations reached.");
      }
      return false;

#undef CMD
    }

    /**
     * Solve the BVP using the multiple shooting method.
     * \tparam SolutionMethod The solution method choice.
     * \param[in] t_mesh Independent variable mesh.
     * \param[in] x_guess Initial guess for the states at the mesh nodes.
     * \return True if the system is successfully solved, false otherwise.
     */
    template <SolutionChoice SolutionMethod = SolutionChoice::GN_AUGMENTED>
    bool multiple_shooting(const VectorX &t_mesh, const MatrixX &x_guess) {
#define CMD "Sandals::BVP::multiple_shooting(...): "

      using VectorShooting = Eigen::Vector<Real, Eigen::Dynamic>;
      using MatrixShooting = Eigen::SparseMatrix<Real>;

      // Check input dimensions
      SANDALS_ASSERT(t_mesh.size() > 1,
                     CMD "expected at least two mesh points.");
      SANDALS_ASSERT(x_guess.cols() == t_mesh.size(),
                     CMD "incompatible mesh and guess sizes.");

      // Disable adaptive time-stepping
      this->m_integrator->disable_adaptive_mode();

      // Problem dimensions
      const Integer num_intervals{static_cast<Integer>(t_mesh.size()) - 1};
      const Integer local_intervals{std::max(1, this->m_subintervals)};
      const Integer c_size{num_intervals * N};
      const Integer x_size{(num_intervals + 1) * N};
      const Integer h_size{
        this->m_integrator->projection_mode() ? 0 : (num_intervals + 1) * M};
      const Integer lsq_rows{c_size + h_size + N};
      const Integer lsq_cols{x_size};
      const Integer A_lsq_nnz{(c_size + h_size + 2 * N) * N + c_size};
      const Real sqrt_sigma{std::sqrt(this->m_sigma)};
      Integer idx_x_ini{0};
      Integer idx_x_end{num_intervals};
      if (this->m_integrator->reverse_mode()) {
        idx_x_ini = num_intervals;
        idx_x_end = 0;
      }

      constexpr Real epsilon_x{1.0e-4};
      constexpr Real epsilon_mu{1.0e-8};

      // Least-squares operator
      MatrixShooting A_lsq(lsq_rows, lsq_cols);
      A_lsq.reserve(A_lsq_nnz);
      VectorShooting b_lsq(lsq_rows);
      std::vector<Eigen::Triplet<Real>> triplets_lsq;
      triplets_lsq.reserve(A_lsq_nnz);

      // Augmented KKT system
      MatrixShooting A_aug;
      VectorShooting b_aug;
      std::vector<Eigen::Triplet<Real>> triplets_aug;

      // Normal-equation KKT system
      MatrixShooting A_kkt;
      VectorShooting b_kkt;
      std::vector<Eigen::Triplet<Real>> triplets_kkt;

      // Allocate augmented system
      if constexpr (SolutionMethod == SolutionChoice::GN_AUGMENTED) {
        // Symmetric augmented KKT system
        // / λI  Aᵀ   Bᵀ \
        // | A   -I    0 |
        // \ B    0  -λI /
        const Integer aug_size{x_size + lsq_rows + N};
        A_aug.resize(aug_size, aug_size);
        A_aug.reserve(2 * A_lsq_nnz + 4 * N * N + aug_size + N);
        b_aug.resize(aug_size);
        triplets_aug.reserve(2 * A_lsq_nnz + 4 * N * N + aug_size + N);
      }

      // Allocate normal-equation KKT system
      if constexpr (SolutionMethod == SolutionChoice::GN_KKT) {
        // Symmetric normal-equation KKT system
        // / AᵀA  Bᵀ \
        // \  B    0 /
        const Integer kkt_size{x_size + N};
        A_kkt.resize(kkt_size, kkt_size);
        A_kkt.reserve(4 * A_lsq_nnz + 4 * N * N);
        b_kkt.resize(kkt_size);
        triplets_kkt.reserve(4 * A_lsq_nnz + 4 * N * N);
      }

      // Initialize the solution
      this->m_solution->clear();
      this->m_solution->resize(t_mesh.size());
      this->m_solution->t = t_mesh;

      // State vectors
      MatrixX x_sol(x_guess);
      MatrixX delta_x_sol(MatrixX::Zero(N, num_intervals + 1));

      // Temporary variables
      VectorF x_ini, x_end;
      VectorX t_local_mesh;
      Solution<Real, N, M> local_sol;
      MatrixJF Jx, Jb_x_ini, Jb_x_end;
      MatrixJH Jh;

      // Solvers
      Eigen::SparseQR<MatrixShooting, Eigen::COLAMDOrdering<Integer>> qr;
      // Eigen::SimplicialLDLT<MatrixShooting> ldlt;
      Eigen::SparseQR<MatrixShooting, Eigen::COLAMDOrdering<Integer>> ldlt;

      // Newton iterations
      for (Integer iter{0}; iter < this->m_max_iterations; ++iter) {
        // Reset least-squares system
        triplets_lsq.clear();
        b_lsq.setZero();

        // Assemble the rectangular least-squares operator
        // A_lsq * dx = b_lsq
        for (Integer k{0}; k < num_intervals; ++k) {
          // Local integration mesh
          t_local_mesh =
              VectorX::LinSpaced(local_intervals + 1, t_mesh(k), t_mesh(k + 1));

          // Reset state transition Jacobian
          Jx.setIdentity();

          // Integrate interval
          if (!this->m_integrator->template solve<true>(t_local_mesh,
                                                        x_sol.col(k),
                                                        local_sol,
                                                        Jx)) {
            SANDALS_ERROR(CMD "failed to integrate interval " << k << ".");
            return false;
          }

          // Store local solution
          if (k == 0) {
            this->m_solution->x.col(0) = local_sol.x.col(0);
            this->m_solution->h.col(0) = local_sol.h.col(0);
          }
          this->m_solution->x.col(k + 1) = local_sol.x.col(local_intervals);
          this->m_solution->h.col(k + 1) = local_sol.h.col(local_intervals);

          // Continuity equations: x_{k+1} - Φ_k(x_k)
          for (Integer i{0}; i < N; ++i) {
            // +I block
            triplets_lsq.emplace_back(k * N + i, (k + 1) * N + i, 1.0);

            // -Jx block
            for (Integer j{0}; j < N; ++j) {
              if (std::abs(Jx(i, j)) > EPSILON) {
                triplets_lsq.emplace_back(k * N + i, k * N + j, -Jx(i, j));
              }
            }
          }

          // Continuity residual
          b_lsq.template segment<N>(k * N) =
              this->m_solution->x.col(k + 1) - x_sol.col(k + 1);

          // Constraint manifold equations
          if constexpr (M > 0) {
            if (!this->m_integrator->projection_mode()) {
              // Initial node
              if (k == 0) {
                Jh =
                    sqrt_sigma *
                    this->m_integrator->system()->Jh_x(x_sol.col(0), t_mesh(0));
                for (Integer i{0}; i < M; ++i) {
                  for (Integer j{0}; j < N; ++j) {
                    if (std::abs(Jh(i, j)) > EPSILON) {
                      triplets_lsq.emplace_back(c_size + i, j, Jh(i, j));
                    }
                  }
                }
                b_lsq.template segment<M>(c_size) =
                    -sqrt_sigma *
                    this->m_integrator->system()->h(x_sol.col(0), t_mesh(0));
              }

              // Interior/final nodes
              Jh = sqrt_sigma *
                   this->m_integrator->system()->Jh_x(x_sol.col(k + 1),
                                                      t_mesh(k + 1));

              for (Integer i{0}; i < M; ++i) {
                for (Integer j{0}; j < N; ++j) {
                  if (std::abs(Jh(i, j)) > EPSILON) {
                    triplets_lsq.emplace_back(c_size + (k + 1) * M + i,
                                              (k + 1) * N + j,
                                              Jh(i, j));
                  }
                }
              }

              // Constraint residual
              b_lsq.template segment<M>(c_size + (k + 1) * M) =
                  -sqrt_sigma *
                  this->m_integrator->system()->h(x_sol.col(k + 1),
                                                  t_mesh(k + 1));
            }
          }
        }

        // Boundary states
        x_ini = x_sol.col(idx_x_ini);
        x_end = x_sol.col(idx_x_end);

        // Boundary Jacobians
        Jb_x_ini = this->Jb_x_ini(x_ini, x_end);
        Jb_x_end = this->Jb_x_end(x_ini, x_end);

        // Boundary residual
        b_lsq.template tail<N>() = -this->b(x_ini, x_end);

        // Boundary equations
        // B * dx = -b
        for (Integer i{0}; i < N; ++i) {
          for (Integer j{0}; j < N; ++j) {
            if (std::abs(Jb_x_ini(i, j)) > EPSILON) {
              triplets_lsq.emplace_back(c_size + h_size + i,
                                        idx_x_ini * N + j,
                                        Jb_x_ini(i, j));
            }
            if (std::abs(Jb_x_end(i, j)) > EPSILON) {
              triplets_lsq.emplace_back(c_size + h_size + i,
                                        idx_x_end * N + j,
                                        Jb_x_end(i, j));
            }
          }
        }

        // Build least-squares operator
        A_lsq.setFromTriplets(triplets_lsq.begin(), triplets_lsq.end());
        A_lsq.makeCompressed();

        // Iteration information
        if (this->m_verbose) {
          std::cout << "Iteration " << iter
                    << ": |b| = " << this->b(x_ini, x_end).norm()
                    << ", |r| = " << (A_lsq.transpose() * b_lsq).norm()
                    << std::endl
                    << "  x(" << t_mesh(idx_x_ini)
                    << ") = " << x_ini.transpose() << std::endl
                    << "  x(" << t_mesh(idx_x_end)
                    << ") = " << x_end.transpose() << std::endl;
        }

        // Convergence check
        if (b_lsq.norm() < this->m_tolerance) {
          return true;
        }

        // Rectangular least-squares solve
        if constexpr (SolutionMethod == SolutionChoice::GN_STANDARD) {
          qr.compute(A_lsq);

          SANDALS_ASSERT(qr.info() == Eigen::Success,
                         CMD
                         "failed to factorize "
                         "least-squares system.");

          delta_x_sol = qr.solve(b_lsq).reshaped(N, num_intervals + 1);
        }

        // Symmetric augmented KKT solve
        else if constexpr (SolutionMethod == SolutionChoice::GN_AUGMENTED) {
          triplets_aug.clear();
          b_aug.setZero();

          // Unknown ordering: (dx, z, λ)^T
          const Integer dx_offset{0};
          const Integer z_offset{x_size};
          const Integer lambda_offset{x_size + lsq_rows};

          // A and Aᵀ blocks
          for (Integer k{0}; k < A_lsq.outerSize(); ++k) {
            for (typename MatrixShooting ::InnerIterator it(A_lsq, k); it;
                 ++it) {
              // Aᵀ block
              triplets_aug.emplace_back(dx_offset + it.col(),
                                        z_offset + it.row(),
                                        it.value());
              // A block
              triplets_aug.emplace_back(z_offset + it.row(),
                                        dx_offset + it.col(),
                                        it.value());
            }
          }

          // λI block
          for (Integer i{0}; i < z_offset; ++i) {
            triplets_aug.emplace_back(i, i, epsilon_x);
          }

          // -I block
          for (Integer i{0}; i < lsq_rows; ++i) {
            triplets_aug.emplace_back(z_offset + i, z_offset + i, -1.0);
          }

          // -λI block
          for (Integer i{0}; i < N; ++i) {
            triplets_aug.emplace_back(lambda_offset + i,
                                      lambda_offset + i,
                                      -epsilon_mu);
          }

          // B and Bᵀ blocks
          for (Integer i{0}; i < N; ++i) {
            for (Integer j{0}; j < N; ++j) {
              if (std::abs(Jb_x_ini(i, j)) > EPSILON) {
                const Integer col{idx_x_ini * N + j};
                // Bᵀ block
                triplets_aug.emplace_back(dx_offset + col,
                                          lambda_offset + i,
                                          Jb_x_ini(i, j));
                // B block
                triplets_aug.emplace_back(lambda_offset + i,
                                          dx_offset + col,
                                          Jb_x_ini(i, j));
              }
              if (std::abs(Jb_x_end(i, j)) > EPSILON) {
                const Integer col{idx_x_end * N + j};
                // Bᵀ block
                triplets_aug.emplace_back(dx_offset + col,
                                          lambda_offset + i,
                                          Jb_x_end(i, j));
                // B block
                triplets_aug.emplace_back(lambda_offset + i,
                                          dx_offset + col,
                                          Jb_x_end(i, j));
              }
            }
          }

          // Right-hand side: (0, -r, -b)^T
          // b_aug.head(z_offset) = A_lsq.transpose() * b_lsq;
          b_aug.head(z_offset).setZero();
          b_aug.segment(z_offset, lsq_rows) = b_lsq;
          b_aug.tail(N)                     = -this->b(x_ini, x_end);

          // Build augmented system
          A_aug.setFromTriplets(triplets_aug.begin(), triplets_aug.end());
          A_aug.makeCompressed();

          // Factorize
          ldlt.compute(A_aug);
          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD
                         "failed to factorize "
                         "augmented KKT system.");

          // Solve
          VectorShooting sol_aug(ldlt.solve(b_aug));
          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD
                         "failed to solve "
                         "augmented KKT system.");

          // Extract dx block
          delta_x_sol =
              sol_aug.segment(dx_offset, x_size).reshaped(N, num_intervals + 1);
        }

        // Normal-equation KKT solve
        else if constexpr (SolutionMethod == SolutionChoice::GN_KKT) {
          triplets_kkt.clear();
          b_kkt.setZero();

          // Normal equations
          MatrixShooting AtA(A_lsq.transpose() * A_lsq);

          // AᵀA block
          for (Integer k{0}; k < AtA.outerSize(); ++k) {
            for (typename MatrixShooting ::InnerIterator it(AtA, k); it; ++it) {
              triplets_kkt.emplace_back(it.row(), it.col(), it.value());
            }
          }

          // B and Bᵀ blocks
          for (Integer i{0}; i < N; ++i) {
            for (Integer j{0}; j < N; ++j) {
              if (std::abs(Jb_x_ini(i, j)) > EPSILON) {
                triplets_kkt.emplace_back(j, x_size + i, Jb_x_ini(i, j));
                triplets_kkt.emplace_back(x_size + i, j, Jb_x_ini(i, j));
              }
              if (std::abs(Jb_x_end(i, j)) > EPSILON) {
                triplets_kkt.emplace_back(idx_x_end * N + j,
                                          x_size + i,
                                          Jb_x_end(i, j));
                triplets_kkt.emplace_back(x_size + i,
                                          idx_x_end * N + j,
                                          Jb_x_end(i, j));
              }
            }
          }

          // Right-hand side: (Aᵀr, -b)^T
          b_kkt.head(x_size) = A_lsq.transpose() * b_lsq;
          b_kkt.tail(N)      = -this->b(x_ini, x_end);

          // Build KKT matrix
          A_kkt.setFromTriplets(triplets_kkt.begin(), triplets_kkt.end());
          A_kkt.makeCompressed();

          // Factorize
          ldlt.compute(A_kkt);
          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD
                         "failed to factorize "
                         "KKT system.");

          // Solve
          VectorShooting sol_kkt(ldlt.solve(b_kkt));
          SANDALS_ASSERT(ldlt.info() == Eigen::Success,
                         CMD
                         "failed to solve "
                         "KKT system.");

          // Extract dx block
          delta_x_sol = sol_kkt.head(x_size).reshaped(N, num_intervals + 1);
        }

        // Check update validity
        if (!delta_x_sol.allFinite()) {
          SANDALS_ERROR(CMD "invalid Newton update.");
          return false;
        }
        // Small update
        if (delta_x_sol.norm() < this->m_tolerance) {
          SANDALS_WARNING(CMD "small update step.");
          return false;
        }

        // Newton update
        x_sol += delta_x_sol;
      }

      // Maximum iterations reached
      if (this->m_verbose) {
        SANDALS_WARNING(CMD "maximum iterations reached.");
      }

      return false;

#undef CMD
    }

    /**
     * Solve the boundary value problem using a shooting method.
     * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t}
     * \f$.
     * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
     * \param[in] x_guess Initial guess for the states at the mesh nodes
     * (only for multiple shooting).
     * \tparam ShootingType Type of shooting method to use (single or
     * multiple).
     */
    template <ShootingChoice ShootingType = ShootingChoice::MULTIPLE>
    bool solve(const VectorX &t_mesh, const VectorF &ics, MatrixX &x_guess) {
#define CMD "Sandals::BVP::solve(...): "

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

  };  // BVP

}  // namespace Sandals

#endif  // SANDALS_BVP_HH
