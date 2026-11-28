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

#ifndef SANDALS_OPTIMAL_CONTROL_HH
#define SANDALS_OPTIMAL_CONTROL_HH

#include <Sandals.hh>
#include <Sandals/FiniteDifferences.hh>
#include <Sandals/System/Controlled.hh>

namespace Sandals
{

  /*\
   |    ___        _   _                 _  ____            _             _
   |   / _ \ _ __ | |_(_)_ __ ___   __ _| |/ ___|___  _ __ | |_ _ __ ___ | |
   |  | | | | '_ \| __| | '_ ` _ \ / _` | | |   / _ \| '_ \| __| '__/ _ \| |
   |  | |_| | |_) | |_| | | | | | | (_| | | |__| (_) | | | | |_| | | (_) | |
   |   \___/| .__/ \__|_|_| |_| |_|\__,_|_|\____\___/|_| |_|\__|_|  \___/|_|
   |        |_|
  \*/

  /**
  * \brief Class container for an optimal control problem (OCP) of ODEs/DAEs.
  *
  * Class container for an optimal control problem with fixed initial and final times of the type
  *
  * \f[
  * \begin{array}{ll}
  * \min & \phi(\mathbf{x}(t_0), \mathbf{x}(t_f)) + \displaystyle\int_{t_0}^{t_f} L(\mathbf{x}(t),
  * \mathbf{u}(t), t) \, \mathrm{d}t \\
  * \text{s.t.} & \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, \mathbf{u}, t) = \mathbf{0} \\
  * & \mathbf{b}_{\text{lo}} \leq \mathbf{b}(\mathbf{x}(t_0), \mathbf{x}(t_f)) \leq \mathbf{b}_{\text{up}} \\
  * & \mathbf{g}_{\text{lo}} \leq \mathbf{g}(\mathbf{x}, \mathbf{u}, t) \leq \mathbf{g}_{\text{up}} \\
  * & \mathbf{x}_{\text{lo}} \leq \mathbf{x} \leq \mathbf{x}_{\text{up}}, \quad
  * \mathbf{u}_{\text{lo}} \leq \mathbf{u} \leq \mathbf{u}_{\text{up}} \text{.}
  * \end{array}
  * \f]
  *
  * Only the functions values are mandatory when the corresponding terms are present. All first and
  * second derivatives default to central finite differences, override them to provide the exact ones.
  * Every term defaults to "absent" (zero cost, no boundary conditions, no path constraints).
  *
  * The invariants of the system, if any, are enforced by the transcription at every mesh node after
  * the first one (see Sandals::Collocation); the boundary conditions must make the initial states
  * consistent with them.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam U The number of controls.
  * \tparam M The number of invariants of the system.
  */
  template <typename Real, Integer N, Integer U, Integer M = 0>
  class OptimalControlProblem
  {
  public:
    using System   = Controlled<Real, N, U, M>; /**< Controlled system type. */
    using SystemPtr = typename System::Pointer; /**< Controlled system shared pointer type. */
    using VectorF  = Eigen::Vector<Real, N>; /**< States vector type. */
    using VectorU  = Eigen::Vector<Real, U>; /**< Controls vector type. */
    using VectorB  = Eigen::Vector<Real, 2*N>; /**< Boundary states vector type \f$ [\mathbf{x}_0; \mathbf{x}_f] \f$. */
    using MatrixB  = Eigen::Matrix<Real, 2*N, 2*N>; /**< Boundary states Hessian type. */
    using VectorP  = Eigen::Vector<Real, N+U>; /**< Path variables vector type \f$ [\mathbf{x}; \mathbf{u}] \f$. */
    using MatrixP  = Eigen::Matrix<Real, N+U, N+U>; /**< Path variables Hessian type. */
    using VectorX  = Eigen::Vector<Real, Eigen::Dynamic>; /**< Dynamic vector type. */
    using MatrixX  = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; /**< Dynamic matrix type. */

    static constexpr Real INFTY{std::numeric_limits<Real>::infinity()}; /**< Infinity. */

  private:
    std::string m_name;    /**< Name of the optimal control problem. */
    SystemPtr   m_system;  /**< Controlled ODE/DAE system. */
    Real        m_t_ini;   /**< Initial time. */
    Real        m_t_end;   /**< Final time. */
    VectorF     m_x_lower{VectorF::Constant(-INFTY)}; /**< States lower bounds. */
    VectorF     m_x_upper{VectorF::Constant(INFTY)};  /**< States upper bounds. */
    VectorU     m_u_lower{VectorU::Constant(-INFTY)}; /**< Controls lower bounds. */
    VectorU     m_u_upper{VectorU::Constant(INFTY)};  /**< Controls upper bounds. */
    VectorF     m_x_nominal{VectorF::Ones()};         /**< Typical magnitudes of the states. */
    VectorU     m_u_nominal{VectorU::Ones()};         /**< Typical magnitudes of the controls. */

    static VectorF head(VectorB const & z) {return z.template head<N>();}
    static VectorF tail(VectorB const & z) {return z.template tail<N>();}
    static VectorF state(VectorP const & z) {return z.template head<N>();}
    static VectorU control(VectorP const & z) {return z.template tail<U>();}

  public:
    /**
    * Class constructor for the optimal control problem.
    * \param[in] t_name The name of the problem.
    * \param[in] t_system The controlled ODE/DAE system.
    * \param[in] t_t_ini The initial time \f$ t_0 \f$.
    * \param[in] t_t_end The final time \f$ t_f \f$.
    */
    OptimalControlProblem(std::string t_name, SystemPtr t_system, Real const t_t_ini, Real const t_t_end)
      : m_name(t_name), m_system(t_system), m_t_ini(t_t_ini), m_t_end(t_t_end)
    {
      SANDALS_ASSERT(t_t_end > t_t_ini, "Sandals::OptimalControlProblem: final time must be greater than initial time.");
    }

    /**
    * Class destructor.
    */
    virtual ~OptimalControlProblem() {}

    /** \return The name of the problem. */
    std::string const & name() const {return this->m_name;}

    /** \return The controlled ODE/DAE system. */
    SystemPtr const & system() const {return this->m_system;}

    /** \return The initial time \f$ t_0 \f$. */
    Real t_ini() const {return this->m_t_ini;}

    /** \return The final time \f$ t_f \f$. */
    Real t_end() const {return this->m_t_end;}

    /** \return The states lower bounds. */
    VectorF const & x_lower() const {return this->m_x_lower;}

    /** \return The states upper bounds. */
    VectorF const & x_upper() const {return this->m_x_upper;}

    /** \return The controls lower bounds. */
    VectorU const & u_lower() const {return this->m_u_lower;}

    /** \return The controls upper bounds. */
    VectorU const & u_upper() const {return this->m_u_upper;}

    /**
    * Set the states bounds (use infinity for unbounded components).
    * \param[in] t_lower The states lower bounds.
    * \param[in] t_upper The states upper bounds.
    */
    void x_bounds(VectorF const & t_lower, VectorF const & t_upper) {this->m_x_lower = t_lower; this->m_x_upper = t_upper;}

    /**
    * Set the controls bounds (use infinity for unbounded components).
    * \param[in] t_lower The controls lower bounds.
    * \param[in] t_upper The controls upper bounds.
    */
    void u_bounds(VectorU const & t_lower, VectorU const & t_upper) {this->m_u_lower = t_lower; this->m_u_upper = t_upper;}

    /**
    * Set the typical magnitudes of the states and of the controls (default one). The transcriptions
    * scale the variables and the constraints by them, so that the solver works on quantities of order
    * one, and its tolerances are relative to these magnitudes.
    * \param[in] t_x_nominal The states magnitudes (positive).
    * \param[in] t_u_nominal The controls magnitudes (positive).
    */
    void nominal(VectorF const & t_x_nominal, VectorU const & t_u_nominal)
    {
      SANDALS_ASSERT((t_x_nominal.array() > 0.0).all() && (t_u_nominal.array() > 0.0).all(),
        "Sandals::OptimalControlProblem::nominal(...): the magnitudes must be positive.");
      this->m_x_nominal = t_x_nominal;
      this->m_u_nominal = t_u_nominal;
    }

    /** \return The typical magnitudes of the states. */
    VectorF const & x_nominal() const {return this->m_x_nominal;}

    /** \return The typical magnitudes of the controls. */
    VectorU const & u_nominal() const {return this->m_u_nominal;}

    /*\
     |  Mayer term
    \*/

    /**
    * Evaluate the Mayer term \f$ \phi(\mathbf{x}_0, \mathbf{x}_f) \f$ (default zero).
    * \param[in] x_ini Initial states \f$ \mathbf{x}_0 \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_f \f$.
    * \return The Mayer term value.
    */
    virtual Real mayer(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const {return Real(0.0);}

    /**
    * Evaluate the gradient of the Mayer term with respect to \f$ [\mathbf{x}_0; \mathbf{x}_f] \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_0 \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_f \f$.
    * \return The gradient (size \f$ 2N \f$).
    */
    virtual VectorB mayer_gradient(VectorF const & x_ini, VectorF const & x_end) const
    {
      VectorB z; z << x_ini, x_end;
      return finite_differences_gradient<Real, 2*N>([this] (VectorB const & w)
        {return this->mayer(head(w), tail(w));}, z);
    }

    /**
    * Evaluate the Hessian of the Mayer term with respect to \f$ [\mathbf{x}_0; \mathbf{x}_f] \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_0 \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_f \f$.
    * \return The (symmetric) Hessian (size \f$ 2N \times 2N \f$).
    */
    virtual MatrixB mayer_hessian(VectorF const & x_ini, VectorF const & x_end) const
    {
      VectorB z; z << x_ini, x_end;
      return finite_differences_hessian<Real, 2*N>([this] (VectorB const & w)
        {return this->mayer_gradient(head(w), tail(w));}, z);
    }

    /*\
     |  Lagrange term
    \*/

    /**
    * Evaluate the Lagrange term \f$ L(\mathbf{x}, \mathbf{u}, t) \f$ (default zero).
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Lagrange term value.
    */
    virtual Real lagrange(VectorF const & /*x*/, VectorU const & /*u*/, Real const /*t*/) const {return Real(0.0);}

    /**
    * Evaluate the gradient of the Lagrange term with respect to \f$ [\mathbf{x}; \mathbf{u}] \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The gradient (size \f$ N+U \f$).
    */
    virtual VectorP lagrange_gradient(VectorF const & x, VectorU const & u, Real const t) const
    {
      VectorP z; z << x, u;
      return finite_differences_gradient<Real, N+U>([this, t] (VectorP const & w)
        {return this->lagrange(state(w), control(w), t);}, z);
    }

    /**
    * Evaluate the Hessian of the Lagrange term with respect to \f$ [\mathbf{x}; \mathbf{u}] \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The (symmetric) Hessian (size \f$ (N+U) \times (N+U) \f$).
    */
    virtual MatrixP lagrange_hessian(VectorF const & x, VectorU const & u, Real const t) const
    {
      VectorP z; z << x, u;
      return finite_differences_hessian<Real, N+U>([this, t] (VectorP const & w)
        {return this->lagrange_gradient(state(w), control(w), t);}, z);
    }

    /*\
     |  Boundary conditions
    \*/

    /**
    * Get the number of boundary conditions (default zero).
    * \return The number of boundary conditions.
    */
    virtual Integer boundary_conditions_number() const {return 0;}

    /**
    * Evaluate the boundary conditions \f$ \mathbf{b}(\mathbf{x}_0, \mathbf{x}_f) \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_0 \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_f \f$.
    * \return The boundary conditions residual.
    */
    virtual VectorX boundary_conditions(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const
    {return VectorX::Zero(0);}

    /**
    * Get the boundary conditions lower bounds (default zero, i.e., equality conditions).
    * \return The lower bounds (size \f$ n_b \f$).
    */
    virtual VectorX boundary_conditions_lower() const
    {return VectorX::Zero(this->boundary_conditions_number());}

    /**
    * Get the boundary conditions upper bounds (default zero, i.e., equality conditions).
    * \return The upper bounds (size \f$ n_b \f$).
    */
    virtual VectorX boundary_conditions_upper() const
    {return VectorX::Zero(this->boundary_conditions_number());}

    /**
    * Evaluate the Jacobian of the boundary conditions with respect to \f$ [\mathbf{x}_0; \mathbf{x}_f] \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_0 \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_f \f$.
    * \return The Jacobian (size \f$ n_b \times 2N \f$).
    */
    virtual MatrixX boundary_conditions_jacobian(VectorF const & x_ini, VectorF const & x_end) const
    {
      VectorB z; z << x_ini, x_end;
      return finite_differences_jacobian<Real, Eigen::Dynamic, 2*N>([this] (VectorB const & w)
        {return this->boundary_conditions(head(w), tail(w));}, z, this->boundary_conditions_number());
    }

    /**
    * Evaluate the Hessian of the weighted boundary conditions \f$ \boldsymbol{\lambda}^\top\mathbf{b} \f$
    * with respect to \f$ [\mathbf{x}_0; \mathbf{x}_f] \f$.
    * \param[in] x_ini Initial states \f$ \mathbf{x}_0 \f$.
    * \param[in] x_end Final states \f$ \mathbf{x}_f \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$ (size \f$ n_b \f$).
    * \return The (symmetric) Hessian (size \f$ 2N \times 2N \f$).
    */
    virtual MatrixB boundary_conditions_hessian(VectorF const & x_ini, VectorF const & x_end,
      VectorX const & lambda) const
    {
      VectorB z; z << x_ini, x_end;
      return finite_differences_hessian<Real, 2*N>([this, &lambda] (VectorB const & w)
        {return VectorB(this->boundary_conditions_jacobian(head(w), tail(w)).transpose()*lambda);}, z);
    }

    /*\
     |  Path constraints
    \*/

    /**
    * Get the number of path constraints (default zero).
    * \return The number of path constraints.
    */
    virtual Integer path_constraints_number() const {return 0;}

    /**
    * Evaluate the path constraints \f$ \mathbf{g}(\mathbf{x}, \mathbf{u}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The path constraints value.
    */
    virtual VectorX path_constraints(VectorF const & /*x*/, VectorU const & /*u*/, Real const /*t*/) const
    {return VectorX::Zero(0);}

    /**
    * Get the path constraints lower bounds (default \f$ -\infty \f$).
    * \return The lower bounds (size \f$ n_g \f$).
    */
    virtual VectorX path_constraints_lower() const
    {return VectorX::Constant(this->path_constraints_number(), -INFTY);}

    /**
    * Get the path constraints upper bounds (default zero, i.e., \f$ \mathbf{g} \leq \mathbf{0} \f$).
    * \return The upper bounds (size \f$ n_g \f$).
    */
    virtual VectorX path_constraints_upper() const
    {return VectorX::Zero(this->path_constraints_number());}

    /**
    * Evaluate the Jacobian of the path constraints with respect to \f$ [\mathbf{x}; \mathbf{u}] \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian (size \f$ n_g \times (N+U) \f$).
    */
    virtual MatrixX path_constraints_jacobian(VectorF const & x, VectorU const & u, Real const t) const
    {
      VectorP z; z << x, u;
      return finite_differences_jacobian<Real, Eigen::Dynamic, N+U>([this, t] (VectorP const & w)
        {return this->path_constraints(state(w), control(w), t);}, z, this->path_constraints_number());
    }

    /**
    * Evaluate the Hessian of the weighted path constraints \f$ \boldsymbol{\lambda}^\top\mathbf{g} \f$
    * with respect to \f$ [\mathbf{x}; \mathbf{u}] \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$ (size \f$ n_g \f$).
    * \return The (symmetric) Hessian (size \f$ (N+U) \times (N+U) \f$).
    */
    virtual MatrixP path_constraints_hessian(VectorF const & x, VectorU const & u, Real const t,
      VectorX const & lambda) const
    {
      VectorP z; z << x, u;
      return finite_differences_hessian<Real, N+U>([this, t, &lambda] (VectorP const & w)
        {return VectorP(this->path_constraints_jacobian(state(w), control(w), t).transpose()*lambda);}, z);
    }

    /*\
     |  Structural sparsity patterns (default dense)
    \*/

    using MatrixXb = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>; /**< Boolean pattern type. */

    /**
    * \return The structural pattern of the Jacobian of the DAE residual in its arguments
    * \f$ (\mathbf{x}, \mathbf{x}^{\prime}, \mathbf{u}) \f$ (\f$ N \times (2N+U) \f$, true where an entry
    * may be nonzero). The transcriptions only store and compute the entries in the patterns.
    */
    virtual MatrixXb dae_jacobian_pattern() const {return MatrixXb::Constant(N, 2*N+U, true);}

    /** \return The structural pattern of the Hessian of any weighted sum of the DAE residuals in \f$ (\mathbf{x}, \mathbf{x}^{\prime}, \mathbf{u}) \f$. */
    virtual MatrixXb dae_hessian_pattern() const {return MatrixXb::Constant(2*N+U, 2*N+U, true);}

    /** \return The structural pattern of the Hessian of the Lagrange term in \f$ (\mathbf{x}, \mathbf{u}) \f$. */
    virtual MatrixXb lagrange_hessian_pattern() const {return MatrixXb::Constant(N+U, N+U, true);}

    /** \return The structural pattern of the Jacobian of the path constraints in \f$ (\mathbf{x}, \mathbf{u}) \f$. */
    virtual MatrixXb path_jacobian_pattern() const {return MatrixXb::Constant(this->path_constraints_number(), N+U, true);}

    /** \return The structural pattern of the Hessian of any weighted sum of the path constraints in \f$ (\mathbf{x}, \mathbf{u}) \f$. */
    virtual MatrixXb path_hessian_pattern() const {return MatrixXb::Constant(N+U, N+U, true);}

    /** \return The structural pattern of the Jacobian of the invariants in \f$ \mathbf{x} \f$ (\f$ M \times N \f$). */
    virtual MatrixXb invariants_jacobian_pattern() const {return MatrixXb::Constant(M, N, true);}

    /** \return The structural pattern of the second (and contracted third) derivatives of the invariants in \f$ \mathbf{x} \f$. */
    virtual MatrixXb invariants_hessian_pattern() const {return MatrixXb::Constant(N, N, true);}

  }; // class OptimalControlProblem

} // namespace Sandals

#endif // SANDALS_OPTIMAL_CONTROL_HH
