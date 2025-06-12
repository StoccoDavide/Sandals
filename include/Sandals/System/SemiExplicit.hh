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

#ifndef SANDALS_SEMIEXPLICIT_SYSTEM_HH
#define SANDALS_SEMIEXPLICIT_SYSTEM_HH

#include <Sandals.hh>
#include <Sandals/System/Explicit.hh>

namespace Sandals {

  /*\
   |   ____                 _ _____            _ _      _ _
   |  / ___|  ___ _ __ ___ (_) ____|_  ___ __ | (_) ___(_) |_
   |  \___ \ / _ \ '_ ` _ \| |  _| \ \/ / '_ \| | |/ __| | __|
   |   ___) |  __/ | | | | | | |___ >  <| |_) | | | (__| | |_
   |  |____/ \___|_| |_| |_|_|_____/_/\_\ .__/|_|_|\___|_|\__|
   |                                    |_|
  \*/

  /**
  * \brief Class container for the system of semi-explicit ODEs/DAEs.
  *
  * Class container for the system of semi-explicit ordinary differential equations (ODEs) or
  * differential algebraic equations (DAEs)of the type \f$ \mathbf{A}(\mathbf{x}, t)
  * \mathbf{x}^{\prime} = \mathbf{b}(\mathbf{x}, t) \f$, with invariants manifold \f$ \mathbf{h}(
  * \mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam Real const The scalar number type.
  * \tparam N The dimension of the semi-explicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SemiExplicit : public Explicit<Real, N, M>
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */

    using Pointer = std::shared_ptr<SemiExplicit<Real, N, M>>; /**< Shared pointer to a semi-explicit ODE/DAE system. */
    using VectorF = typename Explicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixJF = typename Explicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using MatrixA = typename Explicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using TensorTA = typename std::vector<MatrixJF>; /**< Templetized matrix type. */
    using VectorB = typename Explicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixJB = typename Explicit<Real, N, M>::MatrixJF; /**< Templetized vector type. */
    using Type = typename Explicit<Real, N, M>::Type; /**< System type enumeration. */

  private:
    mutable Eigen::FullPivLU<MatrixA> m_lu; /**< LU decomposition. */

  public:
    /**
    * Class constructor for the semi-explicit ODE/DAE system.
    */
    SemiExplicit() : Explicit<Real, N, M>(Type::SEMIEXPLICIT, "(missing name)") {}

    /**
    * Class constructor for the semi-explicit ODE/DAE system.
    * \param[in] t_name The name of the semi-explicit ODE/DAE system.
    */
    SemiExplicit(std::string t_name) : Explicit<Real, N, M>(Type::SEMIEXPLICIT, t_name) {}

    /**
    * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$
    *
    * \f[
    * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    * \mathbf{A}(\mathbf{x}, t)\mathbf{x}^{\prime} - \mathbf{b}(\mathbf{x}, t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    VectorF F(VectorF const & x, VectorF const & x_dot, Real const t) const override
    {
      return this->A(x, t)*x_dot - this->b(x, t);
    }

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{
    * \prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    * -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    * -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = -\mathbf{A}(\mathbf{x}, t)^{-1}(\mathbf{Jb}_{
    * \mathbf{x}}(\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{x}^\prime) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x(VectorF const & x, VectorF const & x_dot, Real const t) const override
    {
      return -this->Jf_x(x, x_dot, t);
    }

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    * \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    * \mathbf{A}(\mathbf{x}, t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x_dot(VectorF const & x, VectorF const &/*x_dot*/, Real const t) const override
    {
      return this->A(x, t);
    }

    /**
    * Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$ as
    *
    * \f[
    * \mathbf{f}(\mathbf{x}, t) = \mathbf{A}(\mathbf{x}, t)^{-1} \mathbf{b}(\mathbf{x}, t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    */
    VectorF f(VectorF const & x, Real const t) const override
    {
      this->m_lu.compute(this->A(x, t));
      SANDALS_ASSERT(this->m_lu.rank() == N,
        "Sandals:SemiExplicit::f(...): singular mass matrix A detected.");
      return this->m_lu.solve(this->b(x, t));
    }

    /**
    * Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^\prime, t) = \displaystyle\frac{\partial
    * \mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{A}(\mathbf{x}, t)^{-1}(\mathbf{Jb}_{
    * \mathbf{x}}(\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{x}^\prime) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^\prime, t) \f$.
    */
    MatrixJF Jf_x(VectorF const & x, VectorF const & x_dot, Real const t) const
    {
      TensorTA TA_x{this->TA_x(x, t)};
      MatrixJF tAp{MatrixJF::Zero()};
      for (Integer i{0}; i < N; ++i) {tAp.col(i) = TA_x[i] * x_dot;} // TODO: Check if TA_x.size() = N
      this->m_lu.compute(this->A(x, t));
      SANDALS_ASSERT(this->m_lu.rank() == N,
        "Sandals:SemiExplicit::Jf_x(...): singular mass matrix A detected.");
      return this->m_lu.solve(this->Jb_x(x, t) - tAp);
    }

    /**
    * Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = \displaystyle\frac{\partial\mathbf{f}(
    * \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{A}(\mathbf{x}, t)^{-1}( \mathbf{Jb}_{\mathbf{x}}
    * (\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{f}(\mathbf{x}, t)) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    MatrixJF Jf_x(VectorF const & x, Real const t) const override {return this->Jf_x(x, this->f(x, t), t);}

    /**
    * Evaluate the semi-explicit ODE/DAE system mass matrix \f$ \mathbf{A}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{A}(\mathbf{x}, t) \f$.
    */
    virtual MatrixA A(VectorF const & x, Real const t) const = 0;

    /**
    * Evaluate the tensor of the semi-explicit ODE/DAE system mass matrix \f$ \mathbf{A}(\mathbf{x},
    * t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) =
    * \displaystyle\frac{\partial\mathbf{A}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    virtual TensorTA TA_x(VectorF const & x, Real const t) const = 0;

    /**
    * Evaluate the semi-explicit ODE/DAE system right-hand-side function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
    */
    virtual VectorB b(VectorF const & x, Real const t) const = 0;

    /**
    * Evaluate the Jacobian of the semi-explicit ODE/DAE system right-hand-side function \f$ \mathbf{b}
    * (\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) =
    * \displaystyle\frac{\partial\mathbf{b}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    virtual MatrixJB Jb_x(VectorF const & x, Real const t) const = 0;

  }; // class SemiExplicit

  /*\
   |   ____                 _ _____            _ _      _ _ __        __
   |  / ___|  ___ _ __ ___ (_) ____|_  ___ __ | (_) ___(_) |\ \      / / __ __ _ _ __  _ __   ___ _ __
   |  \___ \ / _ \ '_ ` _ \| |  _| \ \/ / '_ \| | |/ __| | __\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
   |   ___) |  __/ | | | | | | |___ >  <| |_) | | | (__| | |_ \ V  V /| | | (_| | |_) | |_) |  __/ |
   |  |____/ \___|_| |_| |_|_|_____/_/\_\ .__/|_|_|\___|_|\__| \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
   |                                    |_|                                     |_|   |_|
  \*/

  /**
  * \brief Class container for the system of semi-explicit ODEs/DAEs wrapper.
  *
  * Class container for the system of semi-explicit ordinary differential equations (ODEs) or
  * differential algebraic equations (DAEs)of the type \f$ \mathbf{A}(\mathbf{x}, t)
  * \mathbf{x}^{\prime} = \mathbf{b}(\mathbf{x}, t) \f$, with invariants manifold \f$ \mathbf{h}(
  * \mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam Real const The scalar number type.
  * \tparam N The dimension of the semi-explicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SemiExplicitWrapper : public SemiExplicit<Real, N, M>
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */

    using Pointer = std::shared_ptr<SemiExplicitWrapper<Real, N, M>>; /**< Shared pointer to a semi-explicit ODE/DAE system. */
    using VectorF = typename SemiExplicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixA = typename SemiExplicit<Real, N, M>::MatrixA; /**< Templetized matrix type. */
    using TensorTA = typename SemiExplicit<Real, N, M>::TensorTA; /**< Templetized matrix type. */
    using VectorB = typename SemiExplicit<Real, N, M>::VectorB; /**< Templetized vector type. */
    using MatrixJB = typename SemiExplicit<Real, N, M>::MatrixJB; /**< Templetized vector type. */
    using VectorH = typename SemiExplicit<Real, N, M>::VectorH; /**< Templetized invariants vector type. */
    using MatrixJH = typename SemiExplicit<Real, N, M>::MatrixJH; /**< Templetized invariants Jacobian type. */
    using FunctionA = std::function<MatrixA(VectorF const &, Real const)>; /**< Function type for the mass matrix. */
    using FunctionTA = std::function<TensorTA(VectorF const &, Real const)>; /**< Function type for the mass matrix. */
    using FunctionB = std::function<VectorB(VectorF const &, Real const)>; /**< Function type for the right-hand-side. */
    using FunctionJB = std::function<MatrixJB(VectorF const &, Real const)>; /**< Function type for the right-hand-side Jacobian. */
    using FunctionH = std::function<VectorH(VectorF const &, Real const)>; /**< Invariants function type. */
    using FunctionJH = std::function<MatrixJH(VectorF const &, Real const)>; /**< Jacobian of the invariants function type. */
    using FunctionID = std::function<bool(VectorF const &, Real const)>; /**< In-domain function type. */

    inline const static FunctionH DefaultH = [](VectorF const &, Real const) {return VectorH::Zero();}; /**< Default mass matrix function. */
    inline const static FunctionJH DefaultJH = [](VectorF const &, Real const) {return MatrixJH::Zero();}; /**< Default system matrix function. */
    inline const static FunctionID DefaultID = [](VectorF const &, Real const) {return true;}; /**< Default in-domain function. */

  private:
    FunctionA m_A{nullptr}; /**< Function for the mass matrix. */
    FunctionTA m_TA_x{nullptr}; /**< Function for the mass matrix. */
    FunctionB m_b{nullptr}; /**< Function for the right-hand-side. */
    FunctionJB m_Jb_x{nullptr}; /**< Function for the right-hand-side Jacobian. */
    FunctionH m_h{nullptr}; /**< Invariants function. */
    FunctionJH m_Jh_x{nullptr}; /**< Jacobian of the invariants function. */
    FunctionID m_in_domain{nullptr}; /**< In-domain function. */

  public:
    /**
    * Class constructor for the semi-explicit ODE/DAE system wrapper.
    * \param[in] t_A The function for the mass matrix.
    * \param[in] t_TA_x The function for the mass matrix.
    * \param[in] t_b The function for the right-hand-side.
    * \param[in] t_Jb_x The function for the right-hand-side Jacobian.
    * \param[in] t_h The invariants function.
    * \param[in] t_Jh_x The Jacobian of the invariants function.
    * \param[in] t_in_domain The in-domain function.
    */
    SemiExplicitWrapper(FunctionA t_A, FunctionTA t_TA_x, FunctionB t_b, FunctionJB t_Jb_x,
      FunctionH t_h = DefaultH, FunctionJH t_Jh_x = DefaultJH, FunctionID t_in_domain = DefaultID)
      : SemiExplicit<Real, N, M>(), m_A(t_A), m_TA_x(t_TA_x), m_b(t_b), m_Jb_x(t_Jb_x), m_h(t_h),
      m_Jh_x(t_Jh_x), m_in_domain(t_in_domain)
    {}

    /**
    * Class constructor for the semi-explicit ODE/DAE system wrapper.
    * \param[in] t_name The name of the semi-explicit ODE/DAE system.
    * \param[in] t_A The function for the mass matrix.
    * \param[in] t_TA_x The function for the mass matrix.
    * \param[in] t_b The function for the right-hand-side.
    * \param[in] t_Jb_x The function for the right-hand-side Jacobian.
    * \param[in] t_h The invariants function.
    * \param[in] t_Jh_x The Jacobian of the invariants function.
    * \param[in] t_in_domain The in-domain function.
    */
    SemiExplicitWrapper(std::string t_name, FunctionA t_A, FunctionTA t_TA_x, FunctionB t_b, FunctionJB t_Jb_x,
      FunctionH t_h = DefaultH, FunctionJH t_Jh_x = DefaultJH, FunctionID t_in_domain = DefaultID)
      : SemiExplicit<Real, N, M>(t_name), m_A(t_A), m_TA_x(t_TA_x), m_b(t_b), m_Jb_x(t_Jb_x), m_h(t_h),
      m_Jh_x(t_Jh_x), m_in_domain(t_in_domain)
    {}

    /**
    * Class destructor for the semi-explicit ODE/DAE system wrapper.
    */
    ~SemiExplicitWrapper() {}

    /**
    * Get the function for the mass matrix.
    * \return The function for the mass matrix.
    */
    FunctionA & A() {return this->m_A;}

    /**
    * Get the function for the mass matrix tensor with respect to the states.
    * \return The function for the mass matrix tensor with respect to the states.
    */
    FunctionTA & TA_x() {return this->m_TA_x;}

    /**
    * Get the function for the right-hand-side.
    * \return The function for the right-hand-side.
    */
    FunctionB & b() {return this->m_b;}

    /**
    * Get the function for the right-hand-side Jacobian.
    * \return The function for the right-hand-side Jacobian.
    */
    FunctionJB & Jb_x() {return this->m_Jb_x;}

    /**
    * Get the system's invariants.
    * \return The system's invariants.
    */
    FunctionH & h() {return this->m_h;}

    /**
    * Get the Jacobian of the system's invariants with respect to the states.
    * \return The Jacobian of the system's invariants with respect to the states.
    */
    FunctionJH & Jh_x() {return this->m_Jh_x;}

    /**
    * Get the in-domain function.
    * \return The in-domain function.
    */
    FunctionID & in_domain() {return this->m_in_domain;}

    /**
    * Evaluate th e semi-explicit ODE/DAE system mass matrix \f$ \mathbf{A}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system mass matrix \f$ \mathbf{A}(\mathbf{x}, t) \f$.
    */
    MatrixA A(VectorF const & x, Real const t) const override
    {
      return this->m_A(x, t);
    }

    /**
    * Evaluate the tensor of the semi-explicit ODE/DAE system mass matrix \f$ \mathbf{A}(\mathbf{x},
    * t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) =
    * \displaystyle\frac{\partial\mathbf{A}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    TensorTA TA_x(VectorF const & x, Real const t) const override
    {
      return this->m_TA_x(x, t);
    }

    /**
    * Evaluate the semi-explicit ODE/DAE system right-hand-side function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
    */
    VectorB b(VectorF const & x, Real const t) const override
    {
      return this->m_b(x, t);
    }

    /**
    * Evaluate the Jacobian of the semi-explicit ODE/DAE system right-hand-side function \f$ \mathbf{b}
    * (\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) =
    * \displaystyle\frac{\partial\mathbf{b}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    MatrixJB Jb_x(VectorF const & x, Real const t) const override
    {
      return this->m_Jb_x(x, t);
    }


    /**
    * Evaluate the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    */
    VectorH h(VectorF const & x, Real const t) const override
    {
      return this->m_h(x, t);
    }

    /**
    * Evaluate the Jacobian of the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jh}_{\mathbf{x}}(\mathbf{x} t) = \displaystyle\frac{\partial\mathbf{h}(\mathbf{x},
    * t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jh}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    MatrixJH Jh_x(VectorF const & x, Real const t) const override
    {
      return this->m_Jh_x(x, t);
    }

    /**
    * Return true if the values \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the
    * domain of the ODE/DAE system.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return True if \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the domain of the ODE/DAE system.
    */
    bool in_domain(VectorF const & x, Real const t) const override
    {
      return this->m_in_domain(x, t);
    }

  }; // class SemiExplicitWrapper

} // namespace Sandals

#endif // SANDALS_SEMIEXPLICIT_SYSTEM_HH
