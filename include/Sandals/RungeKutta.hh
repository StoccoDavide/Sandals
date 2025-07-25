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

#ifndef SANDALS_RUNGEKUTTA_HH
#define SANDALS_RUNGEKUTTA_HH

#include <optional>

#include <Sandals.hh>

#include <Sandals/System/Implicit.hh>
#include <Sandals/System/Explicit.hh>
#include <Sandals/System/SemiExplicit.hh>
#include <Sandals/System/Linear.hh>

#include <Sandals/Tableau.hh>
#include <Sandals/Solution.hh>

namespace Sandals {

  /*\
   |   ____                         _  __     _   _
   |  |  _ \ _   _ _ __   __ _  ___| |/ /   _| |_| |_ __ _
   |  | |_) | | | | '_ \ / _` |/ _ \ ' / | | | __| __/ _` |
   |  |  _ <| |_| | | | | (_| |  __/ . \ |_| | |_| || (_| |
   |  |_| \_\\__,_|_| |_|\__, |\___|_|\_\__,_|\__|\__\__,_|
   |                     |___/
  \*/

  /**
  * \brief Class container for the generic \em implicit, \em explicit, and <em> diagonally implicit
  * </em> Runge-Kutta methods.
  *
  * \includedoc docs/markdown/RungeKutta.md
  *
  * \tparam Real The scalar number type.
  * \tparam S The number of stages of the Runge-Kutta method.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer S, Integer N, Integer M>
  class RungeKutta
  {
    using VectorX = Eigen::Vector<Real, Eigen::Dynamic>; /**< \f$ N \times 1 \f$ vector of Real number type (column vector). */
    using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; /**< \f$ N \times N \f$ matrix of Real number type. */
    using VectorK = Eigen::Vector<Real, N*S>; /**< Templetized vector type. */
    using MatrixK = Eigen::Matrix<Real, N, S>; /**< Templetized matrix type. */
    using MatrixJK = Eigen::Matrix<Real, N*N, S>; /**< Templetized matrix type. */
    using MatrixJ = Eigen::Matrix<Real, N*S, N*S>; /**< Templetized matrix type. */
    using VectorP = Eigen::Matrix<Real, N+M, 1>; /**< Templetized vector type. */
    using MatrixP = Eigen::Matrix<Real, N+M, N+M>; /**< Templetized matrix type. */
    using NewtonX = Optimist::RootFinder::Newton<Real, N>; /**< Templetized Newton solver for ERK and DIRK methods. */
    using NewtonK = Optimist::RootFinder::Newton<Real, N*S>; /**< Templetized Newton solver for IRK methods. */
    using VectorS = typename Tableau<Real, S>::Vector; /**< Templetized vector type. */
    using MatrixS = typename Tableau<Real, S>::Matrix; /**< Templetized matrix type. */
    using VectorN = typename Implicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixN = typename Implicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using VectorM = typename Implicit<Real, N, M>::VectorH; /**< Templetized vector type. */
    using MatrixM = typename Implicit<Real, N, M>::MatrixJH; /**< Templetized matrix type. */
    using FunctionSC = std::function<void(Integer const, VectorX const &, Real const)>; /**< Step callback function type. */

  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */
    const Real SQRT_EPSILON{std::sqrt(EPSILON)}; /**< Square root of machine epsilon epsilon static constant value. */ \

    using System = typename Implicit<Real, N, M>::Pointer; /**< Shared pointer to an implicit ODE/DAE system. */
    using Type = typename Tableau<Real, S>::Type; /**< Runge-Kutta type enumeration. */
    using Time = Eigen::Vector<Real, Eigen::Dynamic>; /**< Templetized vector type for the independent variable (or time). */

  private:
    mutable NewtonX m_newtonX;                     /**< Newton solver for ERK and DIRK methods. */
    mutable NewtonK m_newtonK;                     /**< Newton solver for IRK methods. */
    mutable Eigen::FullPivLU<MatrixP> m_lu;        /**< LU decomposition for the projection matrix. */

    Tableau<Real, S> m_tableau;                    /**< Butcher tableau of the Runge-Kutta method. */
    System     m_system;                           /**< ODE/DAE system object pointer. */
    Real       m_absolute_tolerance{1e-6};         /**< Absolute tolerance for adaptive step \f$ \epsilon_{\text{abs}} \f$. */
    Real       m_relative_tolerance{1e-3};         /**< Relative tolerance for adaptive step \f$ \epsilon_{\text{rel}} \f$. */
    Real       m_safety_factor{0.9};               /**< Safety factor for adaptive step \f$ f \f$. */
    Real       m_min_safety_factor{0.1};           /**< Minimum safety factor for adaptive step \f$ f_{\max} \f$. */
    Real       m_max_safety_factor{10.0};          /**< Maximum safety factor for adaptive step \f$ f_{\min} \f$. */
    Real       m_min_step{EPSILON_HIGH};           /**< Minimum step for advancing \f$ h_{\min} \f$. */
    Integer    m_max_substeps{5};                  /**< Maximum number of substeps. */
    bool       m_adaptive{true};                   /**< Adaptive step mode boolean. */
    bool       m_verbose{false};                   /**< Verbose mode boolean. */
    bool       m_reverse{false};                   /**< Time reverse mode boolean. */
    FunctionSC m_step_callback{nullptr};           /**< Step callback function. */

    Real    m_projection_tolerance{EPSILON_HIGH}; /**< Projection tolerance \f$ \epsilon_{\text{proj}} \f$. */
    Integer m_max_projection_iterations{5};       /**< Maximum number of projection steps. */
    bool    m_projection{true};                   /**< Adaptive step mode boolean. */

  public:
    /**
    * Copy constructor for the timer.
    */
    RungeKutta(const RungeKutta &) = delete;

    /**
    * Assignment operator for the timer.
    */
    RungeKutta & operator=(RungeKutta const & ) = delete;

    /**
    * Class constructor for the Runge-Kutta method.
    * \param[in] t_tableau The Tableau reference.
    */
    RungeKutta(Tableau<Real, S> const & t_tableau)
      : m_tableau(t_tableau) {
      this->verbose_mode(this->m_verbose);
    }

    /**
    * Class constructor for the Runge-Kutta method.
    * \param[in] t_tableau The Tableau reference.
    * \param[in] t_system The ODE/DAE system shared pointer.
    */
    RungeKutta(Tableau<Real, S> const & t_tableau, System t_system)
      : m_tableau(t_tableau), m_system(t_system) {
      this->verbose_mode(this->m_verbose);
    }

    /**
    * Get the enumeration type of the Runge-Kutta method.
    * \return The enumeration type of the Runge-Kutta method.
    */
    Type type() const {return this->m_tableau.type;}

    /**
    * Check if the method is an explicit Runge-Kutta (ERK) method.
    * \return True if the method is an explicit Runge-Kutta method, false otherwise.
    */
    bool is_erk() const {return this->m_tableau.type == Type::ERK;}

    /**
    * Check if the method is an implicit Runge-Kutta (IRK) method.
    * \return True if the method is an implicit Runge-Kutta method, false otherwise.
    */
    bool is_irk() const {return this->m_tableau.type == Type::IRK;}

    /**
    * Check if the method is a diagonally implicit Runge-Kutta (DIRK) method.
    * \return True if the method is a diagonally implicit Runge-Kutta method, false otherwise.
    */
    bool is_dirk() const {return this->m_tableau.type == Type::DIRK;}

    /**
    * Get the Butcher Tableau reference.
    * \return The Tableau reference.
    */
    Tableau<Real, S> & tableau() {return this->m_tableau;}

    /**
    * Get the Butcher Tableau const reference.
    * \return The Tableau const reference.
    */
    Tableau<Real, S> const & tableau() const {return this->m_tableau;}

    /**
    * Get the stages \f$ s \f$ number of the Runge-Kutta method.
    * \return The stages \f$ s \f$ number of the Runge-Kutta method.
    */
    Integer stages() const {return S;}

    /**
    * Get the name of the Runge-Kutta method.
    * \return The name of the Runge-Kutta method.
    */
    std::string name() const {return this->m_tableau.name;}

    /**
    * Get the order \f$ p \f$ of the Runge-Kutta method.
    * \return The order \f$ p \f$ of the Runge-Kutta method.
    */
    Integer order() const {return this->m_tableau.order;}

    /**
    * Check if the Runge-Kutta method is embedded.
    * \return True if the Runge-Kutta method is embedded, false otherwise.
    */
    bool is_embedded() const {return this->m_tableau.is_embedded;}

    /**
    * Get the Butcher tableau matrix \f$ \mathbf{A} \f$.
    * \return The Butcher tableau matrix \f$ \mathbf{A} \f$.
    */
    MatrixS A() const {return this->m_tableau.A;}

    /**
    * Get the Butcher tableau weights vector \f$ \mathbf{b} \f$.
    * \return The Butcher tableau weights vector \f$ \mathbf{b} \f$.
    */
    VectorS b() const {return this->m_tableau.b;}

    /**
    * Get the Butcher tableau embedded weights vector \f$ \hat{\mathbf{b}} \f$.
    * \return The Butcher tableau embedded weights vector \f$ \hat{\mathbf{b}} \f$.
    */
    VectorS b_embedded() const {return this->m_tableau.b_e;}

    /**
    * Get the Butcher tableau nodes vector \f$ \mathbf{c} \f$.
    * \return The Butcher tableau nodes vector \f$ \mathbf{c} \f$.
    */
    VectorS c() const {return this->m_tableau.c;}

    /**
    * Get the ODE/DAE system pointer.
    * \return The ODE/DAE system pointer.
    */
    System system() {return this->m_system;}

    /**
    * Set the ODE/DAE system pointer.
    * \param[in] t_system The ODE/DAE system pointer.
    */
    void system(System t_system) {this->m_system = t_system;}

    /**
    * Set the implicit ODE/DAE system with lambda functions.
    * \param[in] F The ODE/DAE system function.
    * \param[in] JF_x The Jacobian of the ODE/DAE system function with respect to the states.
    * \param[in] JF_x_dot The Jacobian of the ODE/DAE system function with respect to the states derivative.
    * \param[in] h The system's invariants.
    * \param[in] Jh_x The Jacobian of the system's invariants with respect to the states.
    * \param[in] in_domain The in-domain function.
    */
    void implicit_system(
      typename ImplicitWrapper<Real, N, M>::FunctionF F,
      typename ImplicitWrapper<Real, N, M>::FunctionJF JF_x,
      typename ImplicitWrapper<Real, N, M>::FunctionJF JF_x_dot,
      typename ImplicitWrapper<Real, N, M>::FunctionH h = ImplicitWrapper<Real, N, M>::DefaultH,
      typename ImplicitWrapper<Real, N, M>::FunctionJH Jh_x = ImplicitWrapper<Real, N, M>::DefaultJH,
      typename ImplicitWrapper<Real, N, M>::FunctionID in_domain = ImplicitWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<ImplicitWrapper<Real, N, M>>(F, JF_x, JF_x_dot, h, Jh_x, in_domain);
    }

    /**
    * Set the implicit ODE/DAE system with lambda functions.
    * \param[in] name The name of the implicit ODE/DAE system.
    * \param[in] F The ODE/DAE system function.
    * \param[in] JF_x The Jacobian of the ODE/DAE system function with respect to the states.
    * \param[in] JF_x_dot The Jacobian of the ODE/DAE system function with respect to the states derivative.
    * \param[in] h The system's invariants.
    * \param[in] Jh_x The Jacobian of the system's invariants with respect to the states.
    * \param[in] in_domain The in-domain function.
    */
    void implicit_system(
      std::string name,
      typename ImplicitWrapper<Real, N, M>::FunctionF F,
      typename ImplicitWrapper<Real, N, M>::FunctionJF JF_x,
      typename ImplicitWrapper<Real, N, M>::FunctionJF JF_x_dot,
      typename ImplicitWrapper<Real, N, M>::FunctionH h = ImplicitWrapper<Real, N, M>::DefaultH,
      typename ImplicitWrapper<Real, N, M>::FunctionJH Jh_x = ImplicitWrapper<Real, N, M>::DefaultJH,
      typename ImplicitWrapper<Real, N, M>::FunctionID in_domain = ImplicitWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<ImplicitWrapper<Real, N, M>>(name, F, JF_x, JF_x_dot, h, Jh_x, in_domain);
    }

    /**
    * Set the explicit ODE/DAE system with lambda functions.
    * \param[in] f The explicit ODE system function.
    * \param[in] Jf_x The Jacobian of the explicit ODE system function with respect to the states.
    * \param[in] h The system's invariants.
    * \param[in] Jh_x The Jacobian of the system's invariants with respect to the states.
    * \param[in] in_domain The in-domain function.
    */
    void explicit_system(
      typename ExplicitWrapper<Real, N, M>::FunctionF f,
      typename ExplicitWrapper<Real, N, M>::FunctionJF Jf_x,
      typename ExplicitWrapper<Real, N, M>::FunctionH h = ExplicitWrapper<Real, N, M>::DefaultH,
      typename ExplicitWrapper<Real, N, M>::FunctionJH Jh_x = ExplicitWrapper<Real, N, M>::DefaultJH,
      typename ExplicitWrapper<Real, N, M>::FunctionID in_domain = ExplicitWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<ExplicitWrapper<Real, N, M>>(f, Jf_x, h, Jh_x, in_domain);
    }

    /**
    * Set the explicit ODE/DAE system with lambda functions.
    * \param[in] name The name of the explicit ODE/DAE system.
    * \param[in] f The explicit ODE system function.
    * \param[in] Jf_x The Jacobian of the explicit ODE system function with respect to the states.
    * \param[in] h The system's invariants.
    * \param[in] Jh_x The Jacobian of the system's invariants with respect to the states.
    * \param[in] in_domain The in-domain function.
    */
    void explicit_system(
      std::string name,
      typename ExplicitWrapper<Real, N, M>::FunctionF f,
      typename ExplicitWrapper<Real, N, M>::FunctionJF Jf_x,
      typename ExplicitWrapper<Real, N, M>::FunctionH h = ExplicitWrapper<Real, N, M>::DefaultH,
      typename ExplicitWrapper<Real, N, M>::FunctionJH Jh_x = ExplicitWrapper<Real, N, M>::DefaultJH,
      typename ExplicitWrapper<Real, N, M>::FunctionID in_domain = ExplicitWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<ExplicitWrapper<Real, N, M>>(name, f, Jf_x, h, Jh_x, in_domain);
    }

    /**
    * Set the linear ODE/DAE system with lambda functions.
    * \param[in] E The mass matrix function.
    * \param[in] A The system matrix function.
    * \param[in] b The system vector function.
    * \param[in] h The invariants function.
    * \param[in] Jh_x The Jacobian of the invariants function.
    * \param[in] in_domain The in-domain function.
    */
    void linear_system(
      typename LinearWrapper<Real, N, M>::FunctionE E,
      typename LinearWrapper<Real, N, M>::FunctionA A,
      typename LinearWrapper<Real, N, M>::FunctionB b,
      typename LinearWrapper<Real, N, M>::FunctionH h = LinearWrapper<Real, N, M>::DefaultH,
      typename LinearWrapper<Real, N, M>::FunctionJH Jh_x = LinearWrapper<Real, N, M>::DefaultJH,
      typename LinearWrapper<Real, N, M>::FunctionID in_domain = LinearWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<LinearWrapper<Real, N, M>>(E, A, b, h, Jh_x, in_domain);
    }

    /**
    * Set the linear ODE/DAE system with lambda functions.
    * \param[in] name The name of the linear ODE/DAE system.
    * \param[in] E The mass matrix function.
    * \param[in] A The system matrix function.
    * \param[in] b The system vector function.
    * \param[in] h The invariants function.
    * \param[in] Jh_x The Jacobian of the invariants function.
    * \param[in] in_domain The in-domain function.
    */
    void linear_system(
      std::string name,
      typename LinearWrapper<Real, N, M>::FunctionE E,
      typename LinearWrapper<Real, N, M>::FunctionA A,
      typename LinearWrapper<Real, N, M>::FunctionB b,
      typename LinearWrapper<Real, N, M>::FunctionH h = LinearWrapper<Real, N, M>::DefaultH,
      typename LinearWrapper<Real, N, M>::FunctionJH Jh_x = LinearWrapper<Real, N, M>::DefaultJH,
      typename LinearWrapper<Real, N, M>::FunctionID in_domain = LinearWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<LinearWrapper<Real, N, M>>(name, E, A, b, h, Jh_x, in_domain);
    }

    /**
    * Set the semi-explicit ODE/DAE system with lambda functions.
    * \param[in] A The function for the mass matrix.
    * \param[in] TA_x The function for the mass matrix.
    * \param[in] b The function for the right-hand-side.
    * \param[in] Jb_x The function for the right-hand-side Jacobian.
    * \param[in] h The invariants function.
    * \param[in] Jh_x The Jacobian of the invariants function.
    * \param[in] in_domain The in-domain function.
    */
    void semi_explicit_system(
      typename SemiExplicitWrapper<Real, N, M>::FunctionA A,
      typename SemiExplicitWrapper<Real, N, M>::FunctionTA TA_x,
      typename SemiExplicitWrapper<Real, N, M>::FunctionB b,
      typename SemiExplicitWrapper<Real, N, M>::FunctionJB Jb_x,
      typename SemiExplicitWrapper<Real, N, M>::FunctionH h = SemiExplicitWrapper<Real, N, M>::DefaultH,
      typename SemiExplicitWrapper<Real, N, M>::FunctionJH Jh_x = SemiExplicitWrapper<Real, N, M>::DefaultJH,
      typename SemiExplicitWrapper<Real, N, M>::FunctionID in_domain = SemiExplicitWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<SemiExplicitWrapper<Real, N, M>>(A, TA_x, b, Jb_x, h, Jh_x, in_domain);
    }

    /**
    * Set the semi-explicit ODE/DAE system with lambda functions.
    * \param[in] name The name of the semi-explicit ODE/DAE system.
    * \param[in] A The function for the mass matrix.
    * \param[in] TA_x The function for the mass matrix.
    * \param[in] b The function for the right-hand-side.
    * \param[in] Jb_x The function for the right-hand-side Jacobian.
    * \param[in] h The invariants function.
    * \param[in] Jh_x The Jacobian of the invariants function.
    * \param[in] in_domain The in-domain function.
    */
    void semi_explicit_system(
      std::string name,
      typename SemiExplicitWrapper<Real, N, M>::FunctionA A,
      typename SemiExplicitWrapper<Real, N, M>::FunctionTA TA_x,
      typename SemiExplicitWrapper<Real, N, M>::FunctionB b,
      typename SemiExplicitWrapper<Real, N, M>::FunctionJB Jb_x,
      typename SemiExplicitWrapper<Real, N, M>::FunctionH h = SemiExplicitWrapper<Real, N, M>::DefaultH,
      typename SemiExplicitWrapper<Real, N, M>::FunctionJH Jh_x = SemiExplicitWrapper<Real, N, M>::DefaultJH,
      typename SemiExplicitWrapper<Real, N, M>::FunctionID in_domain = SemiExplicitWrapper<Real, N, M>::DefaultID
    ) {
      this->m_system = std::make_shared<SemiExplicitWrapper<Real, N, M>>(name, A, TA_x, b, Jb_x, h, Jh_x, in_domain);
    }

    /**
    * Check if the ODE/DAE system pointer is set.
    * \return True if the ODE/DAE system pointer is set, false otherwise.
    */
    bool has_system() {return this->m_system != nullptr;}

    /**
    * Get the adaptive step absolute tolerance \f$ \epsilon_{\text{abs}} \f$.
    * \return The adaptive step absolute tolerance \f$ \epsilon_{\text{abs}} \f$.
    */
    Real absolute_tolerance() {return this->m_absolute_tolerance;}

    /**
    * Get the adaptive step absolute tolerance reference \f$ \epsilon_{\text{abs}} \f$.
    * \param[in] t_absolute_tolerance The adaptive step absolute tolerance \f$ \epsilon_{\text{abs}} \f$.
    */
    void absolute_tolerance(Real const t_absolute_tolerance) {this->m_absolute_tolerance = t_absolute_tolerance;}

    /**
    * Get the adaptive step relative tolerance \f$ \epsilon_{\text{rel}} \f$.
    * \return The adaptive step relative tolerance \f$ \epsilon_{\text{rel}} \f$.
    */
    Real relative_tolerance() {return this->m_relative_tolerance;}

    /**
    * Get the adaptive step relative tolerance reference \f$ \epsilon_{\text{rel}} \f$.
    * \param[in] t_relative_tolerance The adaptive step relative tolerance \f$ \epsilon_{\text{rel}} \f$.
    */
    void relative_tolerance(Real const t_relative_tolerance) {this->m_relative_tolerance = t_relative_tolerance;}

    /**
    * Get the safety factor for adaptive step \f$ f \f$.
    * \return The safety factor for adaptive step \f$ f \f$.
    */
    Real & safety_factor() {return this->m_safety_factor;}

    /**
    * Set safety factor for adaptive step \f$ f \f$.
    * \param[in] t_safety_factor The safety factor for adaptive step \f$ f \f$.
    */
    void safety_factor(Real const t_safety_factor) {this->m_safety_factor = t_safety_factor;}

    /**
    * Get the minimum safety factor for adaptive step \f$ f_{\min} \f$.
    * \return The minimum safety factor for adaptive step \f$ f_{\min} \f$.
    */
    Real & min_safety_factor() {return this->m_min_safety_factor;}

    /**
    * Set the minimum safety factor for adaptive step \f$ f_{\min} \f$.
    * \param[in] t_min_safety_factor The minimum safety factor for adaptive step \f$ f_{\min} \f$.
    */
    void min_safety_factor(Real const t_min_safety_factor) {this->m_min_safety_factor = t_min_safety_factor;}

    /**
    * Get the maximum safety factor for adaptive step \f$ f_{\max} \f$.
    * \return The maximum safety factor for adaptive step \f$ f_{\max} \f$.
    */
    Real & max_safety_factor() {return this->m_max_safety_factor;}

    /**
    * Set the maximum safety factor for adaptive step \f$ f_{\max} \f$.
    * \param[in] t_max_safety_factor The maximum safety factor for adaptive step \f$ f_{\max} \f$.
    */
    void max_safety_factor(Real const t_max_safety_factor) {this->m_max_safety_factor = t_max_safety_factor;}

    /**
    * Get the minimum step for advancing \f$ h_{\min} \f$.
    * \return The minimum step for advancing \f$ h_{\min} \f$.
    */
    Real & min_step() {return this->m_min_step;}

    /**
    * Set the minimum step for advancing \f$ h_{\min} \f$.
    * \param[in] t_min_step The minimum step for advancing \f$ h_{\min} \f$.
    */
    void min_step(Real const t_min_step) {this->m_min_step = t_min_step;}

    /**
    * Get the maximum number of substeps.
    * \return The maximum number of substeps.
    */
    Integer & max_substeps() {return this->m_max_substeps;}

    /**
    * Set the maximum number of substeps.
    * \param[in] t_max_substeps The maximum number of substeps.
    */
    void max_substeps(Integer const t_max_substeps) {this->m_max_substeps = t_max_substeps;}

    /**
    * Get the adaptive step mode.
    * \return The adaptive step mode.
    */
    bool adaptive_mode() {return this->m_adaptive;}

    /**
    * Set the adaptive step mode.
    * \param[in] t_adaptive The adaptive step mode.
    */
    void adaptive(bool t_adaptive) {this->m_adaptive = t_adaptive;}

    /**
    * Enable the adaptive step mode.
    */
    void enable_adaptive_mode() {this->m_adaptive = true;}

    /**
    * Disable the adaptive step mode.
    */
    void disable_adaptive_mode() {this->m_adaptive = false;}

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
      this->m_newtonX.verbose_mode(t_verbose);
      this->m_newtonK.verbose_mode(t_verbose);
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
    * Get the time reverse mode.
    * \return The time reverse mode.
    */
    bool reverse_mode() {return this->m_reverse;}

    /**
    * Set the time reverse mode.
    * \param[in] t_reverse The time reverse mode.
    */
    void reverse(bool t_reverse) {this->m_reverse = t_reverse;}

    /**
    * Enable the time reverse mode.
    */
    void enable_reverse_mode() {this->m_reverse = true;}

    /**
    * Disable the time reverse mode.
    */
    void disable_reverse_mode() {this->m_reverse = false;}

    /**
    * Get the step callback function.
    * \return The step callback function.
    */
    FunctionSC step_callback() {return this->m_step_callback;}

    /**
    * Set the step callback function.
    * \param[in] t_step_callback The step callback function.
    */
    void step_callback(FunctionSC const & t_step_callback) {this->m_step_callback = t_step_callback;}

    /**
    * Get the projection tolerance.
    * \return The projection tolerance.
    */
    Real projection_tolerance() {return this->m_projection_tolerance;}

    /**
    * Set the projection tolerance.
    * \param[in] t_projection_tolerance The projection tolerance.
    */
    void projection_tolerance(Real const t_projection_tolerance)
      {this->m_projection_tolerance = t_projection_tolerance;}

    /**
    * Get the maximum number of projection iterations.
    * \return The maximum number of projection iterations.
    */
    Integer & max_projection_iterations() {return this->m_max_projection_iterations;}

    /**
    * Set the maximum number of projection iterations.
    * \param[in] t_max_projection_iterations The maximum number of projection iterations.
    */
    void max_projection_iterations(Integer const t_max_projection_iterations)
      {this->m_max_projection_iterations = t_max_projection_iterations;}

    /**
    * Get projection mode.
    * \return The projection mode.
    */
    bool projection() {return this->m_projection;}

    /**
    * Set the projection mode.
    * \param[in] t_projection The projection mode.
    */
    void projection(bool t_projection) {this->m_projection = t_projection;}

    /**
    * Enable the projection mode.
    */
    void enable_projection() {this->m_projection = true;}

    /**
    * Disable the projection mode.
    */
    void disable_projection() {this->m_projection = false;}

    /**
    * Estimate the optimal step size for the next advancing step according to the error control method.
    * The error control method used is based on the local truncation error, which is computed as
    *
    * \f[
    * \varepsilon_{\text{t}} = \max\left(\left| \mathbf{x} - \hat{\mathbf{x}} \right|\right)
    * \f]
    *
    * where \f$ \mathbf{x} \f$ is the approximation of the states computed with the higher order method,
    * and \f$ \hat{\mathbf{x}} \f$ is the approximation of the states computed with the lower order method.
    * The desired error is computed as
    *
    * \f[
    * \varepsilon_{\text{d}} = \epsilon_{\text{abs}} + \epsilon_{\text{rel}}
    *   \max\left(\left| \mathbf{x} \right|, \left| \hat{\mathbf{x}} \right|\right)
    * \f]
    *
    * where \f$ \epsilon_{\text{abs}} \f$ is the absolute tolerance and \f$ \epsilon_{\text{rel}} \f$
    * is the relative tolerance. The optimal step size is then estimated as
    *
    * \f[
    * h_{\text{opt}} = h_k \min\left( f_{\max}, \max\left( f_{\min}, f \left(
    *   \displaystyle\frac{\varepsilon_{\text{d}}}{\varepsilon_{\text{t}}} \right)^{\frac{1}{q+1}}
    * \right)\right)
    * \f]
    *
    * where \f$ f \f$ is the safety factor, \f$ f_{\max} \f$ is the maximum safety factor,
    * \f$ f_{\min} \f$ is the minimum safety factor, and \f$ q \f$ is the order of the embedded method.
    *
    * \param[in] x States approximation \f$ \mathbf{x}_{k+1} \f$.
    * \param[in] x_e Embedded method's states approximation \f$ \hat{\mathbf{x}}_{k+1} \f$.
    * \param[in] h_k Actual advancing step \f$ h_k \f$.
    * \return The suggested step for the next integration step \f$ h_{k+1}^\star \f$.
    */
    Real estimate_step(VectorN const & x, VectorN const & x_e, Real const h_k) const
    {
      Real desired_error{this->m_absolute_tolerance + this->m_relative_tolerance *
        std::max(x.array().abs().maxCoeff(), x_e.array().abs().maxCoeff())};
      Real truncation_error{(x - x_e).array().abs().maxCoeff()};
      return h_k * std::min(this->m_max_safety_factor, std::max(this->m_min_safety_factor,
        this->m_safety_factor * std::pow(desired_error/truncation_error,
        1.0/std::max(this->m_tableau.order, this->m_tableau.order_e))));
    }

    /**
    * Print the Runge-Kutta method information to the output stream.
    * \return The Runge-Kutta method information as a string.
    */
    std::string info() const
    {
      std::ostringstream os;
      os
        << "Runge-Kutta method:\t" << this->name() << std::endl
        << "\t- order:\t" << this->order() << std::endl
        << "\t- stages:\t" << this->stages() << std::endl
        << "\t- type:\t";
      switch (this->type()) {
        case Type::ERK: os << "explicit"; break;
        case Type::IRK: os << "implicit"; break;
        case Type::DIRK: os << "diagonally implicit"; break;
      }
      os
        << std::endl
        << "\t- embedded:\t" << this->is_embedded() << std::endl;
      if (this->has_system()) {
        os << "\t- system:\t" << this->m_system->name() << std::endl;
      } else {
        os << "\t- system:\t" << "none" << std::endl;
      }
      return os.str();
    }

    /**
    * Print the Runge-Kutta method information on a stream.
    * \param[in,out] os Output stream.
    */
    void info(std::ostream &os) {os << this->info();}

    /*\
     |   _____ ____  _  __
     |  | ____|  _ \| |/ /
     |  |  _| | |_) | ' /
     |  | |___|  _ <| . \
     |  |_____|_| \_\_|\_\
     |
    \*/

    /**
    * Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    * + h_k \f$ using an explicit Runge-Kutta method, given an explicit system of the form \f$
    * \mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t) \f$. If the method is embedded, the step for
    * the next advancing step \f$ h_{k+1}^\star \f$ is also computed according to the error control
    * method.
    * \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    * \param[in] t_old Independent variable (or time) \f$ t_k \f$ at the \f$ k \f$-th step.
    * \param[in] h_old Advancing step \f$ h_k \f$ at the \f$ k \f$-th step.
    * \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    * \param[out] h_new The suggested step \f$ h_{k+1}^\star \f$ for the next advancing step.
    * \param[out] K The \f$ \mathbf{K} \f$ variables of the Runge-Kutta method.
    * \return True if the step is successfully computed, false otherwise.
    */
    bool erk_explicit_step(VectorN const & x_old, Real const t_old, Real const h_old, VectorN & x_new,
      Real & h_new, MatrixK & K) const
    {
      using Eigen::all;
      using Eigen::seqN;

      // Compute the K variables in the case of an explicit method and explicit system
      VectorN x_node;
      for (Integer i{0}; i < S; ++i) {
        x_node = x_old + K(all, seqN(0, i)) * this->m_tableau.A(i, seqN(0, i)).transpose();
        if (!this->m_reverse) {
          K.col(i) = h_old * static_cast<Explicit<Real, N, M> const *>(this->m_system.get())->f(x_node, t_old + h_old*this->m_tableau.c(i));
        } else {
          K.col(i) = h_old * static_cast<Explicit<Real, N, M> const *>(this->m_system.get())->f_reverse(x_node, t_old + h_old*this->m_tableau.c(i));
        }
      }
      if (!K.allFinite()) {return false;}

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next step
      if (this->m_adaptive && this->m_tableau.is_embedded) {
        VectorN x_emb = x_old + K * this->m_tableau.b_e;
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    /**
    * Compute the residual of system to be solved, which is given by the values of the system
    *
    * \f[
    * \begin{array}{l}
    *   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^{n-1} a_{nj}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    * \end{array} \text{.}
    * \f]
    *
    * where \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$.
    * \note If \f$ h \to 0 \f$, then the first guess to solve the nonlinear system is given by \f$
    * \tilde{\mathbf{K}} = \mathbf{0} \f$.
    * \param[in] s Stage index \f$ s \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] h Advancing step \f$ h \f$.
    * \param[in] K Variables \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ of the system to be solved.
    * \param[out] fun The residual of system to be solved.
    */
    void erk_implicit_function(Integer const s, VectorN const & x, Real const t, Real const h,
      MatrixK const & K, VectorN & fun) const
    {
      using Eigen::all;
      using Eigen::seqN;
      VectorN x_node(x + K(all, seqN(0, s)) * this->m_tableau.A(s, seqN(0, s)).transpose());
      if (!this->m_reverse) {
        fun = this->m_system->F(x_node, K.col(s)/h, t + h * this->m_tableau.c(s));
      } else {
        fun = this->m_system->F_reverse(x_node, K.col(s)/h, t + h * this->m_tableau.c(s));
      }
    }

    /**
    * Compute the Jacobian of the system of equations
    *
    * \f[
    * \begin{array}{l}
    *   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^{n-1} a_{nj}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    * \end{array} \text{.}
    * \f]
    *
    * with respect to the \f$ \tilde{\mathbf{K}}_n = h \mathbf{K}_n \f$ variables. The Jacobian
    * matrix of the \f$ n \f$-th equation with respect to the \f$ \tilde{\mathbf{K}}_n \f$ variables
    * is given by
    *
    * \f[
    *   \frac{\partial\mathbf{F}_n}{\partial\tilde{\mathbf{K}}_n}\left(\mathbf{x} + \displaystyle
    *   \sum_{j=1}^{n-1} a_{nj}\tilde{\mathbf{K}}_j, \displaystyle\frac{\tilde{\mathbf{K}}_n}{h},
    *   t + h c_n \right) = a_{nj} \mathbf{JF}_x + \displaystyle\frac{1}{h}\begin{cases}
    *   \mathbf{JF}_{x^{\prime}} & n = j \\
    *   0 & n \neq j
    * \end{cases} \text{,}
    * \f]
    *
    * for \f$ j = 1, 2, \ldots, s \f$.
    * \param[in] s Stage index \f$ s \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] h Advancing step \f$ h \f$.
    * \param[in] K Variables \f$ h \mathbf{K} \f$ of the system to be solved.
    * \param[out] jac The Jacobian of system to be solved.
    */
    void erk_implicit_jacobian(Integer const s, VectorN const & x, Real const t, Real const h,
      MatrixK const & K, MatrixN & jac) const
    {
      using Eigen::all;
      using Eigen::seqN;
      VectorN x_node(x + K(all, seqN(0, s)) * this->m_tableau.A(s, seqN(0, s)).transpose());
      if (!this->m_reverse) {
        jac = this->m_system->JF_x_dot(x_node, K.col(s)/h, t + h * this->m_tableau.c(s)) / h;
      } else {
        jac = this->m_system->JF_x_dot_reverse(x_node, K.col(s)/h, t + h * this->m_tableau.c(s)) / h;
      }
    }

    /**
    * Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    * + h_k \f$ using an explicit Runge-Kutta method, given an explicit system of the form \f$
    * \mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t) \f$. If the method is embedded, the step for
    * the next advancing step \f$ h_{k+1}^\star \f$ is also computed according to the error control
    * method.
    * \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    * \param[in] t_old Independent variable (or time) \f$ t_k \f$ at the \f$ k \f$-th step.
    * \param[in] h_old Advancing step \f$ h_k \f$ at the \f$ k \f$-th step.
    * \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    * \param[out] h_new The suggested step \f$ h_{k+1}^\star \f$ for the next advancing step.
    * \param[out] K The \f$ \mathbf{K} \f$ variables of the Runge-Kutta method.
    * \return True if the step is successfully computed, false otherwise.
    */
    bool erk_implicit_step(VectorN const & x_old, Real const t_old, Real const h_old, VectorN & x_new,
      Real & h_new, MatrixK & K) const
    {
      VectorN K_sol;
      VectorN K_ini(VectorN::Zero());

      // Check if the solver converged
      for (Integer s{0}; s < S; ++s) {
        if (this->m_newtonX.solve(
            [this, s, &K, &x_old, t_old, h_old](VectorN const & K_fun, VectorN & fun)
              {K.col(s) = K_fun; this->erk_implicit_function(s, x_old, t_old, h_old, K, fun);},
            [this, s, &K, &x_old, t_old, h_old](VectorN const & K_jac, MatrixN & jac)
              {K.col(s) = K_jac; this->erk_implicit_jacobian(s, x_old, t_old, h_old, K, jac);},
            K_ini, K_sol)) {
          K.col(s) = K_sol;
        } else {
          return false;
        }
      }

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next step
      if (this->m_adaptive && this->m_tableau.is_embedded) {
        VectorN x_emb(x_old + K * this->m_tableau.b_e);
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    /**
    * Compute the Jacobian of the variables \f$ \mathbf{K} \f$ with respect to the states
    * \f$ \mathbf{x} \f$ as
    *
    * \f[
    *   \frac{\partial\mathbf{K}}{\partial\mathbf{x}} =
    *   \begin{bmatrix}
    *     \mathbf{I} & \mathbf{0} & \cdots & \mathbf{0} \\
    *


    & \mathbf{0} & \cdots & \mathbf{0} \\
    *     \mathbf{0} & \mathbf{I} & \cdots & \mathbf{0} \\
    *     \vdots & \vdots & \ddots & \vdots \\
    *     \mathbf{0} & \mathbf{0} & \cdots & \mathbf{I}
    *   \end{bmatrix} \text{,}
    * \f]
    *
    * where \f$ \mathbf{I} \f$ is the identity matrix of size \f$ N \times N \f$.
    * \param[out] jac The Jacobian of the variables \f$ \mathbf{K} \f$ with respect to the states
    * \f$ \mathbf{x} \f$.
    */
   // void irk_jacobian_K(MatrixK const & K, MatrixJK & jac) const
   // {
//
   // }

    /*\
     |   ___ ____  _  __
     |  |_ _|  _ \| |/ /
     |   | || |_) | ' /
     |   | ||  _ <| . \
     |  |___|_| \_\_|\_\
     |
    \*/

    /**
    * Compute the residual of system to be solved, which is given by the values of the system
    *
    * \f[
    * \begin{cases}
    *   \mathbf{F}_1 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{1j}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_1}{h}, t + h c_1\right) \\
    *   \mathbf{F}_2 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{2j}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_2}{h}, t + h c_2\right) \\[-0.5em]
    *   \vdots \\[-0.5em]
    *   \mathbf{F}_s \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{sj}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_s}{h}, t + h c_s\right)
    * \end{cases}
    * \f]
    *
    * where \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$.
    * \note If \f$ h \to 0 \f$, then the first guess to solve the nonlinear system is given by \f$
    * \tilde{\mathbf{K}} = \mathbf{0} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] h Advancing step \f$ h \f$.
    * \param[in] K Variables \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ of the system to be solved.
    * \param[out] fun The residual of system to be solved.
    */
    void irk_function(VectorN const & x, Real const t, Real const h, VectorK const & K, VectorK & fun) const
    {
      VectorN x_node;
      MatrixK K_mat{K.reshaped(N, S)};
      MatrixK fun_mat;
      for (Integer i{0}; i < S; ++i) {
        x_node = x + K_mat * this->m_tableau.A.row(i).transpose();
        if (!this->m_reverse) {
          fun_mat.col(i) = this->m_system->F(x_node, K_mat.col(i)/h, t + h * this->m_tableau.c(i));
        } else {
          fun_mat.col(i) = this->m_system->F_reverse(x_node, K_mat.col(i)/h, t + h * this->m_tableau.c(i));
        }
      }
      fun = fun_mat.reshaped(N*S, 1);
    }

    /**
    * Compute the Jacobian of the system of equations
    *
    * \f[
    * \begin{cases}
    *   \mathbf{F}_1 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{1j}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_1}{h}, t + h c_1\right) \\
    *   \mathbf{F}_2 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{2j}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_2}{h}, t + h c_2\right) \\[-0.5em]
    *   \vdots \\[-0.5em]
    *   \mathbf{F}_s \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{sj}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_s}{h}, t + h c_s\right)
    * \end{cases} \text{.}
    * \f]
    *
    * with respect to the \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ variables. The \f$ ij \f$-th
    * Jacobian matrix entry is given by
    *
    * \f[
    *   \frac{\partial\mathbf{F}_i}{\partial\tilde{\mathbf{K}}_j}\left(\mathbf{x} + \displaystyle
    *   \sum_{j=1}^s a_{ij}\tilde{\mathbf{K}}_j, \displaystyle\frac{\tilde{\mathbf{K}}_i}{h}, t +
    *   h c_i \right) = a_{ij} \mathbf{JF}_x + \displaystyle\frac{1}{h}\begin{cases}
    *   \mathbf{JF}_{x^{\prime}} & i = j \\
    *   0 & i \neq j
    * \end{cases} \text{,}
    *
    * \f]
    * for \f$ i = 1, 2, \ldots, s \f$ and \f$ j = 1, 2, \ldots, s \f$.
     * \note If \f$ h \to 0 \f$, then the first guess to solve the nonlinear system is given by \f$
    * \tilde{\mathbf{K}} = \mathbf{0} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] h Advancing step \f$ h \f$.
    * \param[in] K Variables \f$ h \mathbf{K} \f$ of the system to be solved.
    * \param[out] jac The Jacobian of system to be solved.
    */
    void irk_jacobian(VectorN const & x, Real const t, Real const h, VectorK const & K, MatrixJ & jac) const
    {
      using Eigen::seqN;

      // Reset the Jacobian matrix
      jac.setZero();

      // Loop through each equation of the system
      MatrixK K_mat{K.reshaped(N, S)};
      Real t_node;
      VectorN x_node, x_dot_node;
      MatrixN JF_x, JF_x_dot;
      auto idx = seqN(0, N), jdx = seqN(0, N);
      for (Integer i{0}; i < S; ++i) {
        t_node = t + h * this->m_tableau.c(i);
        x_node = x + K_mat * this->m_tableau.A.row(i).transpose();

        // Compute the Jacobians with respect to x and x_dot
        x_dot_node = K_mat.col(i) / h;
        if (!this->m_reverse) {
          JF_x     = this->m_system->JF_x(x_node, x_dot_node, t_node);
          JF_x_dot = this->m_system->JF_x_dot(x_node, x_dot_node, t_node);
        } else {
          JF_x     = this->m_system->JF_x_reverse(x_node, x_dot_node, t_node);
          JF_x_dot = this->m_system->JF_x_dot_reverse(x_node, x_dot_node, t_node);
        }

        // Combine the Jacobians with respect to x and x_dot to obtain the Jacobian with respect to K
        idx = seqN(i*N, N);
        for (Integer j{0}; j < S; ++j) {
          jdx = seqN(j*N, N);
          if (i == j) {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x + JF_x_dot / h;
          } else {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x;
          }
        }
      }
    }

    /**
    * Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    * + h_k \f$ using an implicit Runge-Kutta method, given an implicit system of the form \f$
    * \mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0} \f$. If the method is embedded,
    * the step for the next advancing step \f$ h_{k+1}^\star \f$ is also computed according to the
    * error control method.
    * \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    * \param[in] t_old Independent variable (or time) \f$ t_k \f$ at the \f$ k \f$-th step.
    * \param[in] h_old Advancing step \f$ h_k \f$ at the \f$ k \f$-th step.
    * \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    * \param[out] h_new The suggested step \f$ h_{k+1}^\star \f$ for the next advancing step.
    * \param[out] K The \f$ \mathbf{K} \f$ variables of the Runge-Kutta method.
    * \return True if the step is successfully computed, false otherwise.
    */
    bool irk_step(VectorN const & x_old, Real const t_old, Real const h_old, VectorN & x_new,
      Real & h_new, MatrixK & K) const
    {
      VectorK K_vec;
      VectorK K_ini(VectorK::Zero());

      // Check if the solver converged
      if (!this->m_newtonK.solve(
          [this, &x_old, t_old, h_old](VectorK const & K_fun, VectorK & fun)
            {this->irk_function(x_old, t_old, h_old, K_fun, fun);},
          [this, &x_old, t_old, h_old](VectorK const & K_jac, MatrixJ & jac)
            {this->irk_jacobian(x_old, t_old, h_old, K_jac, jac);},
          K_ini, K_vec))
        {return false;}

      // Reshape the K vector to a matrix
      K = K_vec.reshaped(N, S);
      if (!K.allFinite()) {return false;}

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next step
      if (this->m_adaptive && this->m_tableau.is_embedded) {
        VectorN x_emb(x_old + K * this->m_tableau.b_e);
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    /*\
     |   ____ ___ ____  _  __
     |  |  _ \_ _|  _ \| |/ /
     |  | | | | || |_) | ' /
     |  | |_| | ||  _ <| . \
     |  |____/___|_| \_\_|\_\
     |
    \*/

    /**
    * Compute the residual of system to be solved, which is given by the values of the system
    *
    * \f[
    * \begin{array}{l}
    *   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^n a_{nj}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    * \end{array} \text{.}
    * \f]
    *
    * where \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$.
    * \note If \f$ h \to 0 \f$, then the first guess to solve the nonlinear system is given by \f$
    * \tilde{\mathbf{K}} = \mathbf{0} \f$.
    * \param[in] n Equation index \f$ n \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] h Advancing step \f$ h \f$.
    * \param[in] K Variables \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ of the system to be solved.
    * \param[out] fun The residual of system to be solved.
    */
    void dirk_function(Integer n, VectorN const & x, Real const t, Real const h, MatrixK const & K,
      VectorN & fun) const
    {
      using Eigen::all;
      using Eigen::seqN;
      VectorN x_node(x + K(all, seqN(0, n+1)) * this->m_tableau.A(n, seqN(0, n+1)).transpose());
      if (!this->m_reverse) {
        fun = this->m_system->F(x_node, K.col(n)/h, t + h * this->m_tableau.c(n));
      } else {
        fun = this->m_system->F_reverse(x_node, K.col(n)/h, t + h * this->m_tableau.c(n));
      }
    }

    /**
    * Compute the Jacobian of the system of equations
    *
    * \f[
    * \begin{array}{l}
    *   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^n a_{nj}\tilde{\mathbf{K}}_j,
    *   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    * \end{array} \text{.}
    * \f]
    *
    * with respect to the \f$ \tilde{\mathbf{K}}_n = h \mathbf{K}_n \f$ variables. The Jacobian
    * matrix of the \f$ n \f$-th equation with respect to the \f$ \tilde{\mathbf{K}}_n \f$ variables
    * is given by
    *
    * \f[
    *   \frac{\partial\mathbf{F}_n}{\partial\tilde{\mathbf{K}}_n}\left(\mathbf{x} + \displaystyle
    *   \sum_{j=1}^n a_{nj}\tilde{\mathbf{K}}_j, \displaystyle\frac{\tilde{\mathbf{K}}_n}{h},
    *   t + h c_n \right) = a_{nj} \mathbf{JF}_x + \displaystyle\frac{1}{h}\begin{cases}
    *   \mathbf{JF}_{x^{\prime}} & n = j \\
    *   0 & n \neq j
    * \end{cases} \text{,}
    * \f]
    *
    * for \f$ j = 1, 2, \ldots, s \f$.
    * \note If \f$ h \to 0 \f$, then the first guess to solve the nonlinear system is given by \f$
    * \tilde{\mathbf{K}} = \mathbf{0} \f$.
    * \param[in] n Equation index \f$ n \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] h Advancing step \f$ h \f$.
    * \param[in] K Variables \f$ h \mathbf{K} \f$ of the system to be solved.
    * \param[out] jac The Jacobian of system to be solved.
    */
    void dirk_jacobian(Integer n, VectorN const & x, Real const t, Real const h, MatrixK const & K,
      MatrixN & jac) const
    {
      using Eigen::all;
      using Eigen::seqN;
      Real t_node{t + h * this->m_tableau.c(n)};
      VectorN x_node(x + K(all, seqN(0, n+1)) * this->m_tableau.A(n, seqN(0, n+1)).transpose());
      VectorN x_dot_node(K.col(n)/h);
      if (!this->m_reverse) {
        jac = this->m_tableau.A(n,n) * this->m_system->JF_x(x_node, x_dot_node, t_node) +
          this->m_system->JF_x_dot(x_node, x_dot_node, t_node) / h;
      } else {
        jac = this->m_tableau.A(n,n) * this->m_system->JF_x_reverse(x_node, x_dot_node, t_node) +
          this->m_system->JF_x_dot_reverse(x_node, x_dot_node, t_node) / h;
      }
    }

    /**
    * Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    * + h_k \f$ using an diagonally implicit Runge-Kutta method, given an implicit system of the
    * form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0} \f$. If the method is
    * embedded, the step for the next advancing step \f$ h_{k+1}^\star \f$ is also computed
    * according to the error control method.
    * \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    * \param[in] t_old Independent variable (or time) \f$ t_k \f$ at the \f$ k \f$-th step.
    * \param[in] h_old Advancing step \f$ h_k \f$ at the \f$ k \f$-th step.
    * \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    * \param[out] h_new The suggested step \f$ h_{k+1}^\star \f$ for the next advancing step.
    * \param[out] K The \f$ \mathbf{K} \f$ variables of the Runge-Kutta method.
    * \return True if the step is successfully computed, false otherwise.
    */
    bool dirk_step(VectorN const & x_old, Real const t_old, Real const h_old, VectorN & x_new,
      Real & h_new, MatrixK & K) const
    {
      VectorN K_sol;
      VectorN K_ini(VectorN::Zero());

      // Check if the solver converged at each step
      for (Integer n{0}; n < S; ++n) {
        if (this->m_newtonX.solve(
            [this, n, &K, &x_old, t_old, h_old](VectorN const & K_fun, VectorN & fun)
              {K.col(n) = K_fun; this->dirk_function(n, x_old, t_old, h_old, K, fun);},
            [this, n, &K, &x_old, t_old, h_old](VectorN const & K_jac, MatrixN & jac)
              {K.col(n) = K_jac; this->dirk_jacobian(n, x_old, t_old, h_old, K, jac);},
            K_ini, K_sol)) {
            K.col(n) = K_sol;
          } else {
            return false;
          }
      }

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next step
      if (this->m_adaptive && this->m_tableau.is_embedded) {
        VectorN x_emb(x_old + K * this->m_tableau.b_e);
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    /**
    * Compute a step using a generic integration method for a system of the form \f$ \mathbf{F}(
    * \mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0} \f$. The step is automatically selected
    * based on the system properties and the integration method properties.
    * \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    * \param[in] t_old Independent variable (or time) \f$ t_k \f$ at the \f$ k \f$-th step.
    * \param[in] h_old Advancing step \f$ h_k \f$ at the \f$ k \f$-th step.
    * \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    * \param[out] h_new The suggested step \f$ h_{k+1}^\star \f$ for the next advancing step.
    * \param[out] K The \f$ \mathbf{K} \f$ variables of the Runge-Kutta method.
    * \return True if the step is successfully computed, false otherwise.
    */
    bool step(VectorN const & x_old, Real const t_old, Real const h_old, VectorN & x_new, Real & h_new,
      MatrixK & K) const
    {
      #define CMD "Sandals::RungeKutta::step(...): "

      SANDALS_ASSERT(this->m_system->in_domain(x_old, t_old), CMD "in " << this->m_tableau.name <<
        " solver, at t = " << t_old << ", x = " << x_old.transpose() << ", system out of domain.");

      if (this->is_erk() && this->m_system->is_explicit()) {
        return this->erk_explicit_step(x_old, t_old, h_old, x_new, h_new, K);
      } else if (this->is_erk() && this->m_system->is_implicit()) {
        return this->erk_implicit_step(x_old, t_old, h_old, x_new, h_new, K);
      } else if (this->is_dirk()) {
        return this->dirk_step(x_old, t_old, h_old, x_new, h_new, K);
      } else {
        return this->irk_step(x_old, t_old, h_old, x_new, h_new, K);
      }

      #undef CMD
    }

    /**
    * Advance using a generic integration method for a system of the form \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^\prime, t) = \mathbf{0} \f$. The step is automatically selected based on the
    * system properties and the integration method properties. In the advvancing step, the system
    * solution is also projected on the manifold \f$ \mathbf{h}(\mathbf{x}, t) \f$. Substepping
    * may be also used to if the integration step fails with the current step size.
    * \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    * \param[in] t_old Independent variable (or time) \f$ t_k \f$ at the \f$ k \f$-th step.
    * \param[in] h_old Advancing step \f$ h_k \f$ at the \f$ k \f$-th step.
    * \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    * \param[out] h_new The suggested step \f$ h_{k+1}^\star \f$ for the next advancing step.
    * \return True if the step is successfully computed, false otherwise.
    */
    bool advance(VectorN const & x_old, Real const t_old, Real const h_old, VectorN & x_new, Real & h_new) const
    {
      #define CMD "Sandals::RungeKutta::advance(...): "

      // Check step size
      SANDALS_ASSERT(h_old > Real(0.0), CMD "in " << this->m_tableau.name << " solver, h = "<<
        h_old << ", expected > 0.");

      // If the integration step failed, try again with substepping
      MatrixK K;
      if (!this->step(x_old, t_old, h_old, x_new, h_new, K))
      {
        VectorN x_tmp(x_old);
        Real t_tmp{t_old}, h_tmp{h_old / Real(2.0)};

        // Substepping logic
        Integer max_k{this->m_max_substeps * this->m_max_substeps}, k{2};
        Real h_new_tmp;
        while (k > 0) {
          // Calculate the next step with substepping logic
          if (this->step(x_tmp, t_tmp, h_tmp, x_new, h_new_tmp, K)) {

            // Accept the step
            h_tmp = h_new_tmp;

            // If substepping is enabled, double the step size
            if (k > 0 && k < max_k) {
              k -= 1;
              // If the substepping index is even, double the step size
              if (k % 2 == 0) {
                h_tmp = Real(2.0) * h_tmp;
                if (this->m_verbose) {
                  SANDALS_WARNING(CMD "in " << this->m_tableau.name << " solver, at t = " << t_tmp <<
                    ", integration succedded disable one substepping layer.");
                }
              }
            }

            // Check the infinity norm of the solution
            SANDALS_ASSERT(std::isfinite(x_tmp.maxCoeff()), CMD "in " << this->m_tableau.name <<
              " solver, at t = " << t_tmp << ", ||x||_inf = inf, computation interrupted.");

          } else {

            // If the substepping index is too high, abort the integration
            k += 2;
            SANDALS_ASSERT(k < max_k, CMD "in " << this->m_tableau.name << " solver, at t = " <<
              t_tmp << ", integration failed with h = " << h_tmp << ", aborting.");
            return false;

            // Otherwise, try again with a smaller step
            if (this->m_verbose) {SANDALS_WARNING(CMD "in " << this->m_tableau.name << " solver, " <<
              "at t = " << t_tmp << ", integration failed, adding substepping layer.");}
            h_tmp /= Real(2.0);
            continue;
          }

          // Store independent variable (or time)
          t_tmp += h_tmp;
        }

        // Store output states substepping solutions
        x_new = x_tmp;
        h_new = h_tmp;
      }

      // Project intermediate solution on the invariants
      if (this->m_projection) {
        VectorN x_projected;
        if (this->project(x_new, t_old + h_new, x_projected)) {
          x_new = x_projected;
        } else {
          return false;
        }
      }
      return true;

      #undef CMD
    }

    /**
    * Solve the system and calculate the approximate solution over the independent variable  (or
    * time) mesh\f$ \mathbf{t} = \left[ t_1, t_2, \ldots, t_n \right]^\top \f$. The step size is fixed
    * and given by \f$ h = t_{k+1} - t_k \f$.
    * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t} \f$.
    * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    * \param[out] sol The solution of the system over the mesh of independent variable.
    * \return True if the system is successfully solved, false otherwise.
    * \warning Do not use the solution for internal backtracking, as the step callback may directly
    * modify the solution.
    */
    bool solve(VectorX const & t_mesh, VectorN const & ics, Solution<Real, N, M> & sol) const
    {
      using Eigen::last;

      // Instantiate output
      sol.resize(t_mesh.size());

      // Store initial conditions
      sol.t(0)     = t_mesh(0);
      sol.x.col(0) = ics;
      sol.h.col(0) = this->m_system->h(ics, t_mesh(0));

      // Callback on initial conditions
      if (this->m_step_callback) {this->m_step_callback(0, ics, t_mesh(0));}

      // Update the current step
      Integer step{0};
      VectorN x_old_step(ics), x_new_step(ics);
      Real t_step{t_mesh(0)}, h_step{t_mesh(1) - t_mesh(0)}, h_tmp_step{h_step}, h_new_step;
      bool mesh_point_bool, saturation_bool;

      while (true) {
        // Integrate system
        if (!this->advance(x_old_step, t_step, h_step, x_new_step, h_new_step)) {return false;}

        // Update the current step
        t_step += h_step;

        // Saturate the suggested timestep
        mesh_point_bool = std::abs(t_step - t_mesh(step+1)) < SQRT_EPSILON;
        saturation_bool = t_step + h_new_step > t_mesh(step+1) + SQRT_EPSILON;
        if (this->m_adaptive && this->m_tableau.is_embedded && !mesh_point_bool && saturation_bool) {
          h_tmp_step = h_new_step; // Used to store the previous step and keep the integration pace
          h_step     = t_mesh(step+1) - t_step;
        } else {
          h_step = h_new_step;
        }

        // Store solution if the step is a mesh point
        if (!this->m_adaptive || mesh_point_bool) {

          // Update temporaries
          step += 1;
          h_step = h_tmp_step;

          // Update outputs
          sol.t(step)     = t_step;
          sol.x.col(step) = x_new_step;
          sol.h.col(step) = this->m_system->h(x_new_step, t_step);

          // Callback after the step is completed
          if (this->m_step_callback) {this->m_step_callback(step, x_new_step, t_step);}

          // Check if the current step is the last one
          if (std::abs(t_step - t_mesh(last)) < SQRT_EPSILON) {break;}

          // Update the previous step
          x_old_step = x_new_step;
        }
      }
      return true;
    }

    /**
    * Solve the system and calculate the approximate solution over the suggested independent
    * variable mesh \f$ \mathbf{t} = \left[ t_1, t_2, \ldots, t_n \right]^\top \f$, the step size
    * is automatically computed based on the error control method.
    * \param[in] t_mesh Independent variable (or time) mesh \f$ \mathbf{t} \f$.
    * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    * \param[out] sol The solution of the system over the mesh of independent variable.
    * \return True if the system is successfully solved, false otherwise.
    * \warning Do not use the solution for internal backtracking, as the step callback may directly
    * modify the solution.
    */
    bool adaptive_solve(VectorX const & t_mesh, VectorN const & ics, Solution<Real, N, M> & sol) const
    {
      using Eigen::all;
      using Eigen::last;

      #define CMD "Sandals::RungeKutta::adaptive_solve(...): "

      // Check if the adaptive method is enabled and the method is embedded
      if (!this->is_embedded()) {
        SANDALS_WARNING(CMD "the method is not embedded, using solve(...) method.");
        return this->solve(t_mesh, ics, sol);
      } else if (!this->m_adaptive) {
        SANDALS_WARNING(CMD "adaptive method is disabled, using solve(...) method.");
        return this->solve(t_mesh, ics, sol);
      }

      // Instantiate output
      Real t_step{t_mesh(0)}, h_step{t_mesh(1) - t_mesh(0)}, h_new_step, scale{100.0};
      Real h_min{std::max(this->m_min_step, h_step/scale)}, h_max{scale*h_step};
      if (this->m_tableau.is_embedded) {
        Integer safety_length{static_cast<Integer>(std::ceil(std::abs(t_mesh(last) - t_mesh(0))/(2.0*h_min)))};
        sol.resize(safety_length);
      } else {
        sol.resize(t_mesh.size());
      }

      // Store initial conditions
      sol.t(0)     = t_mesh(0);
      sol.x.col(0) = ics;
      sol.h.col(0) = this->m_system->h(ics, t_mesh(0));

      // Callback on initial conditions
      if (this->m_step_callback) {this->m_step_callback(0, ics, t_mesh(0));}

      // Instantiate temporary variables
      Integer step{0};
      VectorN x_old_step(ics), x_new_step(ics);

      while (true) {
        // Integrate system
        this->advance(x_old_step, t_step, h_step, x_new_step, h_new_step);

        // Update the current step
        t_step += h_step;

        // Saturate the suggested timestep
        if (this->m_adaptive && this->m_tableau.is_embedded) {
          h_step = std::max(std::min(h_new_step, h_max), h_min);
        }

        SANDALS_ASSERT(step < sol.size(), CMD "safety length exceeded.");

        // Update temporaries
        step += 1;

        // Update outputs
        sol.t(step)     = t_step;
        sol.x.col(step) = x_new_step;
        sol.h.col(step) = this->m_system->h(x_new_step, t_step);

        // Callback after the step is completed
        if (this->m_step_callback) {this->m_step_callback(step, x_new_step, t_step);}

        // Check if the current step is the last one
        if (std::abs(t_step - t_mesh(last)) < SQRT_EPSILON) {break;}
        else if (t_step + h_step > t_mesh(last)) {h_step = t_mesh(last) - t_step;}

        // Update the previous step
        x_old_step = x_new_step;
      }

      // Resize the output
      sol.conservative_resize(step);

      return true;

      #undef CMD
    }

    /**
    * Project the system solution \f$ \mathbf{x} \f$ on the invariants \f$ \mathbf{h} (\mathbf{x},
    * t) = \mathbf{0} \f$.
    * \param[in] x The initial guess for the states \f$\widetilde{\mathbf{x}} \f$.
    * \param[in] t The independent variable (or time) \f$ t \f$ at which the states are evaluated.
    * \param[out] x_projected The projected states \f$ \mathbf{x} \f$ closest to the invariants
    * manifold \f$ \mathbf{h} (\mathbf{x}, t) = \mathbf{0} \f$.
    * \return True if the solution is successfully projected, false otherwise.
    */
    bool project(VectorN const & x, Real const t, VectorN & x_projected) const
    {
      #define CMD "Sandals::RungeKutta::project(...): "

      // Check if there are any constraints
      x_projected = x;
      if (M > Integer(0)) {

        VectorM h;
        MatrixM Jh_x;
        VectorP b, x_step;
        MatrixP A;
        A.setZero();
        A.template block<N, N>(0, 0) = MatrixN::Identity();
        for (Integer k{0}; k < this->m_max_projection_iterations; ++k) {

          /* Standard projection method
               [A]         {x}    =        {b}
           / I  Jh_x^T \ /   dx   \   / x_t - x_k \
           |           | |        | = |           |
           \ Jh_x    0 / \ lambda /   \    -h     /
          */

          // Evaluate the invariants vector and its Jacobian
          h    = this->m_system->h(x_projected, t);
          Jh_x = this->m_system->Jh_x(x_projected, t);

          // Check if the solution is found
          if (h.norm() < this->m_projection_tolerance) {return true;}

          // Build the linear system
          A.template block<N, M>(0, N) = Jh_x.transpose();
          A.template block<M, N>(N, 0) = Jh_x;
          b.template head<N>() = x - x_projected;
          b.template tail<M>() = -h;

          // Compute the solution of the linear system
          this->m_lu.compute(A);
          SANDALS_ASSERT(this->m_lu.rank() == N+M, CMD "singular Jacobian detected.");
          x_step = this->m_lu.solve(b);

          // Check if the step is too small
          if (x_step.norm() < this->m_projection_tolerance * this->m_projection_tolerance) {return false;}

          // Update the solution
          x_projected.noalias() += x_step(Eigen::seqN(0, N));
        }
        if (this->m_verbose) {SANDALS_WARNING(CMD "maximum number of iterations reached.");}
        return false;
      } else {
        return true;
      }

      #undef CMD
    }

    /**
    * Project the system solution \f$ \mathbf{x} \f$ on the invariants \f$ \mathbf{h} (\mathbf{x},
    * t) = \mathbf{0} \f$.
    * \param[in] x The initial guess for the states \f$\widetilde{\mathbf{x}} \f$.
    * \param[in] t The independent variable (or time) \f$ t \f$ at which the states are evaluated.
    * \param[in] projected_equations The indices of the states to be projected.
    * \param[in] projected_invariants The indices of the invariants to be projected.
    * \param[out] x_projected The projected states \f$ \mathbf{x} \f$ closest to the invariants
    * manifold \f$ \mathbf{h} (\mathbf{x}, t) = \mathbf{0} \f$.
    * \return True if the solution is successfully projected, false otherwise.
    */
    bool project_ics(VectorN const & x, Real const t, std::vector<Integer> const & projected_equations,
      std::vector<Integer> const & projected_invariants, VectorN & x_projected) const
    {
      #define CMD "Sandals::RungeKutta::project_ics(...): "

      Integer X{static_cast<Integer>(projected_equations.size())};
      Integer H{static_cast<Integer>(projected_invariants.size())};

      // Check if there are any constraints
      x_projected = x;
      if (H > Integer(0)) {

        VectorM h;
        MatrixM Jh_x;
        VectorX b(X+H), x_step(X+H);
        MatrixX A(X+H, X+H);
        A.setZero();
        A.block(0, 0, X, X) = MatrixX::Identity(X+H, X+H);
        Eigen::FullPivLU<MatrixX> lu;
        for (Integer k{0}; k < this->m_max_projection_iterations; ++k) {

          /* Standard projection method
               [A]         {x}    =        {b}
           / I  Jh_x^T \ /   dx   \   / x_t - x_k \
           |           | |        | = |           |
           \ Jh_x    0 / \ lambda /   \    -h     /
          */

          // Evaluate the invariants vector and its Jacobian
          h    = this->m_system->h(x_projected, t);
          Jh_x = this->m_system->Jh_x(x_projected, t);

          // Select only the projected invariants
          h    = h(projected_invariants);
          Jh_x = Jh_x(projected_invariants, projected_equations);

          // Check if the solution is found
          if (h.norm() < this->m_projection_tolerance) {return true;}

          // Build the linear system
          A.block(0, X, X, H) = Jh_x.transpose();
          A.block(X, 0, H, X) = Jh_x;
          b.head(X) = x(projected_equations) - x_projected(projected_equations);
          b.tail(H) = -h;

          // Compute the solution of the linear system
          lu.compute(A);
          SANDALS_ASSERT(lu.rank() == X+H, CMD "singular Jacobian detected.");
          x_step = this->m_lu.solve(b);

          // Check if the step is too small
          if (x_step.norm() < this->m_projection_tolerance * this->m_projection_tolerance) {return false;}

          // Update the solution
          x_projected(projected_equations).noalias() += x_step;
        }
        if (this->m_verbose) {SANDALS_WARNING(CMD "maximum number of iterations reached.");}
        return false;
      } else {
        return true;
      }

      #undef CMD
    }

    /**
    * Estimate the order of the Runge-Kutta method.
    * \param[in] t_mesh The vector of time meshes with same initial and final time with *fixed* step.
    * \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    * \param[in] sol The *analytical* solution function.
    * \return The estimated order of the method.
    */
    Real estimate_order(std::vector<VectorX> const & t_mesh, VectorN const & ics, std::function<MatrixX(VectorX)> & sol) const
    {
      using Eigen::last;

      #define CMD "Sandals::RungeKutta::estimate_order(...): "

      SANDALS_ASSERT(t_mesh.size() > Integer(1), CMD "expected at least two time meshes.");

      for (Integer i{0}; i < static_cast<Integer>(t_mesh.size()); ++i) {

        // Check if the time mesh has the same initial and final time
        SANDALS_ASSERT((t_mesh[0](0) - t_mesh[i](0)) < SQRT_EPSILON,
          CMD "expected the same initial time.");
        SANDALS_ASSERT((t_mesh[0](last) - t_mesh[i](last)) < SQRT_EPSILON,
          CMD "expected the same final time.");

        // Check if the mesh step is fixed all over the time mesh
        for (Integer j{1}; j < static_cast<Integer>(t_mesh[i].size()); ++j) {
          SANDALS_ASSERT((t_mesh[i](j) - t_mesh[i](j-1)) - (t_mesh[i](1) - t_mesh[i](0)) < SQRT_EPSILON,
            CMD "expected a fixed step.");
        }
      }

      // Solve the system for each time scale
      Solution<Real, N, M> sol_num;
      MatrixX sol_ana;
      VectorX h_vec(t_mesh.size()), e_vec(t_mesh.size());
      for (Integer i{0}; i < static_cast<Integer>(t_mesh.size()); ++i) {
        SANDALS_ASSERT(this->solve(t_mesh[i], ics, sol_num), CMD "failed to solve the system for " <<
          "the" << i << "-th time mesh.");
        sol_ana = sol(sol_num.t);
        SANDALS_ASSERT(sol_ana.rows() == sol_num.x.rows(),
          CMD "expected the same number of states in analytical solution.");
        SANDALS_ASSERT(sol_ana.cols() == sol_num.x.cols(),
          CMD "expected the same number of steps in analytical solution.");
        h_vec(i) = std::abs(sol_num.t(1) - sol_num.t(0));
        e_vec(i) = (sol_ana - sol_num.x).array().abs().maxCoeff();
      }

      // Compute the order of the method thorugh least squares
      VectorX A(h_vec.array().log());
      VectorX b(e_vec.array().log());
      return ((A.transpose() * A).ldlt().solve(A.transpose() * b))(0);

      #undef CMD
    }

  }; // class RungeKutta

} // namespace Sandals

#endif // SANDALS_RUNGEKUTTA_HH
