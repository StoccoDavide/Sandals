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

#ifndef SANDALS_LINEAR_SYSTEM_HH
#define SANDALS_LINEAR_SYSTEM_HH

#include <Sandals.hh>
#include <Sandals/System/Explicit.hh>

namespace Sandals {

  /*\
   |   _     _
   |  | |   (_)_ __   ___  __ _ _ __
   |  | |   | | '_ \ / _ \/ _` | '__|
   |  | |___| | | | |  __/ (_| | |
   |  |_____|_|_| |_|\___|\__,_|_|
   |
  \*/

  /**
  * \brief Class container for the system of linear ODEs/DAEs.
  *
  * Class container for the system of linear ordinary differential equations (ODEs) or differential
  * algebraic equations (DAEs) of the type \f$ \mathbf{E}(t)\mathbf{x}^{\prime} = \mathbf{A}(t)
  * \mathbf{x} + \mathbf{b}(t) \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the semi-explicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class Linear : public Explicit<Real, N, M>
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */
    const Real SQRT_EPSILON{std::sqrt(EPSILON)}; /**< Square root of machine epsilon epsilon static constant value. */ \

    using Pointer = std::shared_ptr<SemiExplicit<Real, N, M>>; /**< Shared pointer to a semi-explicit ODE/DAE system. */
    using VectorF = typename Explicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixJF = typename Explicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using MatrixE = typename Explicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using MatrixA = typename Explicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using VectorB = typename Explicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using Type = typename Explicit<Real, N, M>::Type; /**< System type enumeration. */

  private:
    Eigen::FullPivLU<MatrixE> m_lu; /**< LU decomposition. */

  public:
    /**
    * Class constructor for the linear ODE/DAE system.
    */
    Linear() : Explicit<Real, N, M>(Type::LINEAR, "(missing name)") {}

    /**
    * Class constructor for the linear ODE/DAE system.
    * \param[in] t_name The name of the linear ODE/DAE system.
    */
    Linear(std::string t_name) : Explicit<Real, N, M>(Type::LINEAR, t_name) {}

    /**
    * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    *
    * \f[
    * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{E}(t) \mathbf{x}^{\prime} -
    * \mathbf{A}(t) \mathbf{x} - \mathbf{b}(t)
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    VectorF F(VectorF const &x, VectorF const &x_dot, Real t) const override
    {
      return this->E(t)*x_dot - this->A(t)*x - this->b(t);
    }

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    * t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    * -\mathbf{A}(t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real t) const override {return -this->A(t);}

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    * \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    * \mathbf{E}(t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x_dot(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real t) const override {return this->E(t);}

    /**
    * Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$ as
    *
    * \f[
    * \mathbf{f}(\mathbf{x}, t) = \mathbf{E}(t)^{-1}(\mathbf{A}(t)\mathbf{x} + \mathbf{b}(t))
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    */
    VectorF f(VectorF const &x, Real t) const override
    {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Linear::f(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t)*x + this->b(t));
    }

    /**
    * Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = \displaystyle\frac{\partial\mathbf{f}(
    * \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{E}(t)^{-1} \mathbf{A}(t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    MatrixJF Jf_x(VectorF const &/*x*/, Real t) const override
    {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Linear::Jf_x(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t));
    }

    /**
    * Evaluate the linear ODE/DAE system mass matrix \f$ \mathbf{E}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{E}(t) \f$.
    */
    virtual MatrixE E(Real t) const = 0;

    /**
    * Evaluate the linear ODE/DAE system matrix \f$ \mathbf{A}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{A}(t) \f$.
    */
    virtual MatrixA A(Real t) const = 0;

    /**
    * Evaluate the linear ODE/DAE system vector \f$ \mathbf{b}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{b}(t) \f$.
    */
    virtual VectorB b(Real t) const = 0;

  }; // class Linear

  /*\
   |   _     _                     __        __
   |  | |   (_)_ __   ___  __ _ _ _\ \      / / __ __ _ _ __  _ __   ___ _ __
   |  | |   | | '_ \ / _ \/ _` | '__\ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
   |  | |___| | | | |  __/ (_| | |   \ V  V /| | | (_| | |_) | |_) |  __/ |
   |  |_____|_|_| |_|\___|\__,_|_|    \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
   |                                                   |_|   |_|
  \*/

  /**
  * \brief Class container for the system of linear ODEs/DAEs wrapper.
  *
  * Class container for the system of linear ordinary differential equations (ODEs) or differential
  * algebraic equations (DAEs) of the type \f$ \mathbf{E}(t)\mathbf{x}^{\prime} = \mathbf{A}(t)
  * \mathbf{x} + \mathbf{b}(t) \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the semi-explicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class LinearWrapper : public Linear<Real, N, M>
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */

    using Pointer = std::shared_ptr<LinearWrapper<Real, N, M>>; /**< Shared pointer to a linear ODE/DAE system. */
    using VectorF = typename Linear<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixJF = typename Linear<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using MatrixE = typename Linear<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using MatrixA = typename Linear<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using VectorB = typename Linear<Real, N, M>::VectorF; /**< Templetized vector type. */
    using VectorH = typename Explicit<Real, N, M>::VectorH; /**< Templetized vector type. */
    using MatrixJH = typename Explicit<Real, N, M>::MatrixJH; /**< Templetized matrix type. */
    using FunctionE = std::function<MatrixE(Real)>; /**< Mass matrix function type. */
    using FunctionA = std::function<MatrixA(Real)>; /**< System matrix function type. */
    using FunctionB = std::function<VectorB(Real)>; /**< System vector function type. */
    using FunctionH = std::function<VectorH(VectorF const &, Real)>; /**< Invariants function type. */
    using FunctionJH = std::function<MatrixJH(VectorF const &, Real)>; /**< Jacobian of the invariants function type. */
    using FunctionID = std::function<bool(VectorF const &, Real)>; /**< In-domain function type. */

    inline const static FunctionH DefaultH = [](VectorF const &, Real) {return VectorH::Zero();}; /**< Default mass matrix function. */
    inline const static FunctionJH DefaultJH = [](VectorF const &, Real) {return MatrixJH::Zero();}; /**< Default system matrix function. */
    inline const static FunctionID DefaultID = [](VectorF const &, Real) {return true;}; /**< Default in-domain function. */

  private:
    FunctionE m_E{nullptr}; /**< Mass matrix function. */
    FunctionA m_A{nullptr}; /**< System matrix function. */
    FunctionB m_b{nullptr}; /**< System vector function. */
    FunctionH m_h{nullptr}; /**< Invariants function. */
    FunctionJH m_Jh_x{nullptr}; /**< Jacobian of the invariants function. */
    FunctionID m_in_domain{nullptr}; /**< In-domain function. */

  public:
    /**
    * Class constructor for the linear ODE/DAE system wrapper.
    * \param[in] t_E The mass matrix function.
    * \param[in] t_A The system matrix function.
    * \param[in] t_b The system vector function.
    * \param[in] t_h The invariants function.
    * \param[in] t_Jh_x The Jacobian of the invariants function.
    * \param[in] t_in_domain The in-domain function.
    */
    LinearWrapper(FunctionE t_E, FunctionA t_A, FunctionB t_b, FunctionH t_h = DefaultH, FunctionJH
      t_Jh_x = DefaultJH, FunctionID t_in_domain = DefaultID) : Linear<Real, N, M>("(missing name)"),
      m_E(t_E), m_A(t_A), m_b(t_b), m_h(t_h), m_Jh_x(t_Jh_x), m_in_domain(t_in_domain)
    {}

    /**
    * Class constructor for the linear ODE/DAE system wrapper.
    * \param[in] t_name The name of the linear ODE/DAE system.
    * \param[in] t_E The mass matrix function.
    * \param[in] t_A The system matrix function.
    * \param[in] t_b The system vector function.
    * \param[in] t_h The invariants function.
    * \param[in] t_Jh_x The Jacobian of the invariants function.
    * \param[in] t_in_domain The in-domain function.
    */
    LinearWrapper(std::string t_name, FunctionE t_E, FunctionA t_A, FunctionB t_b, FunctionH t_h = DefaultH,
      FunctionJH t_Jh_x = DefaultJH, FunctionID t_in_domain = DefaultID) : Linear<Real, N, M>(t_name),
      m_E(t_E), m_A(t_A), m_b(t_b), m_h(t_h), m_Jh_x(t_Jh_x), m_in_domain(t_in_domain)
    {}

    /**
    * Class destructor for the linear ODE/DAE system wrapper.
    */
    ~LinearWrapper() {}

    /**
    * Get the mass matrix function.
    * \return The mass matrix function.
    */
    FunctionE & E() {return this->m_E;}

    /**
    * Get the system matrix function.
    * \return The system matrix function.
    */
    FunctionA & A() {return this->m_A;}

    /**
    * Get the system vector function.
    * \return The system vector function.
    */
    FunctionB & b() {return this->m_b;}

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
    * Evaluate the linear ODE/DAE system mass matrix \f$ \mathbf{E}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{E}(t) \f$.
    */
    MatrixE E(Real t) const override
    {
      return this->m_E(t);
    }

    /**
    * Evaluate the linear ODE/DAE system matrix \f$ \mathbf{A}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{A}(t) \f$.
    */
    MatrixA A(Real t) const override
    {
      return this->m_A(t);
    }

    /**
    * Evaluate the linear ODE/DAE system vector \f$ \mathbf{b}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{b}(t) \f$.
    */
    VectorB b(Real t) const override
    {
      return this->m_b(t);
    }

    /**
    * Evaluate the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    */
    VectorH h(VectorF const &x, Real t) const override
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
    MatrixJH Jh_x(VectorF const &x, Real t) const override
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
    bool in_domain(VectorF const &x, Real t) const override
    {
      return this->m_in_domain(x, t);
    }

  }; // class LinearWrapper

} // namespace Sandals

#endif // SANDALS_LINEAR_SYSTEM_HH
