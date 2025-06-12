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

#ifndef SANDALS_IMPLICIT_SYSTEM_HH
#define SANDALS_IMPLICIT_SYSTEM_HH

#include <Sandals.hh>

namespace Sandals
{

  /*\
   |   ___                 _ _      _ _
   |  |_ _|_ __ ___  _ __ | (_) ___(_) |_
   |   | || '_ ` _ \| '_ \| | |/ __| | __|
   |   | || | | | | | |_) | | | (__| | |_
   |  |___|_| |_| |_| .__/|_|_|\___|_|\__|
   |                |_|
  \*/

  /**
  * \brief Class container for the system of implicit ODEs/DAEs.
  *
  * Class container for the system of implicit ordinary differential equations (ODEs) or differential
  * algebraic equations (DAEs) of the type \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
  * \mathbf{0} \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the implicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class Implicit
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */

    using Type = enum class Type : Integer {IMPLICIT = 0, EXPLICIT = 1, SEMIEXPLICIT = 1}; /**< System type enumeration. */
    using Pointer = std::shared_ptr<Implicit<Real, N, M>>; /**< Shared pointer to an implicit ODE system. */
    using VectorF = Eigen::Vector<Real, N>; /**< Templetized vector type. */
    using MatrixJF = Eigen::Matrix<Real, N, N>; /**< Templetized matrix type. */
    using VectorH = Eigen::Vector<Real, M>; /**< Templetized vector type. */
    using MatrixJH = Eigen::Matrix<Real, M, N>; /**< Templetized matrix type. */

  private:
    Type        m_type{Type::IMPLICIT}; /**< ODE/DAE system type. */
    std::string m_name;                 /**< Name of the ODE/DAE system. */

  protected:
    /**
    * Class constructor for an implicit ODE/DAE system.
    * \param[in] t_type The type of the implicit ODE/DAE system.
    * \param[in] t_name The name of the implicit ODE/DAE system.
    */
    Implicit(Type t_type, std::string t_name) : m_type(t_type), m_name(t_name) {}

  public:
    /**
    * Class constructor for an implicit ODE/DAE system.
    */
    Implicit() : m_type(Type::IMPLICIT), m_name("(missing name)") {}

    /**
    * Class constructor for an implicit ODE/DAE system.
    * \param[in] t_name The name of the implicit ODE/DAE system.
    */
    Implicit(std::string t_name) : m_type(Type::IMPLICIT), m_name(t_name) {}

    /**
    * Class destructor for an implicit ODE/DAE system.
    */
    virtual ~Implicit() {}

    /**
    * Get the enumeration type of the ODE/DAE system.
    * \return The enumeration type of the ODE/DAE system.
    */
    Type type() const {return this->m_type;}

    /**
    * Check if the ODE/DAE system is implicit.
    * \return True if the ODE/DAE system is implicit, false otherwise.
    */
    bool is_implicit() const {return this->m_type == Type::IMPLICIT;}

    /**
    * Check if the ODE/DAE system is explicit.
    * \return True if the ODE/DAE system is explicit, false otherwise.
    */
    bool is_explicit() const {return this->m_type == Type::EXPLICIT;}

    /**
    * Check if the ODE/DAE system is semi-explicit.
    * \return True if the ODE/DAE system is semi-explicit, false otherwise.
    */
    bool is_semiexplicit() const {return this->m_type == Type::SEMIEXPLICIT;}

    /**
    * Get the ODE/DAE system name reference.
    * \return The ODE/DAE system name reference.
    */
    std::string & name() {return this->m_name;}

    /**
    * Get the ODE/DAE system name const reference.
    * \return The ODE/DAE system name const reference.
    */
    std::string const & name() const {return this->m_name;}

    /**
    * Get the number of equations of the ODE/DAE system.
    * \return The number of equations of the ODE/DAE system.
    */
    Integer equations_number() const {return N;}

    /**
    * Get the number of invariants of the ODE/DAE system.
    * \return The number of invariants of the ODE/DAE system.
    */
    Integer invariants_number() const {return M;}

    /**
    * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    virtual VectorF F(VectorF const & x, VectorF const & x_dot, Real const t) const = 0;

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle\frac{
    * \partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    virtual MatrixJF JF_x(VectorF const & x, VectorF const & x_dot, Real const t) const = 0;

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}}
    * \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    virtual MatrixJF JF_x_dot(VectorF const & x, VectorF const & x_dot, Real const t) const = 0;

    /**
    * Evaluate the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    */
    virtual VectorH h(VectorF const & x, Real const t) const = 0;

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
    virtual MatrixJH Jh_x(VectorF const & x, Real const t) const = 0;

    /**
    * Return true if the values \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$ is in the
    * domain of the ODE/DAE system.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return True if \f$ \mathbf{F}(\mathbf{x}, t) \f$ is in the domain of the ODE/DAE system.
    */
    virtual bool in_domain(VectorF const & x, Real const t) const = 0;

    /**
    * Time reversal of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{
    * \prime}, t) = -\mathbf{F}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, -t) \f$.
    */
    VectorF F_reverse(VectorF const & x, VectorF const & x_dot, Real const t) const
    {
      return -this->F(x, -x_dot, -t);
    }

    /**
    * Time reversal of the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} = -\mathbf{JF}_{\mathbf{x}}
    * (\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, -\mathbf{x}^{\prime},
    * -t) \f$.
    */
    MatrixJF JF_x_reverse(VectorF const & x, VectorF const & x_dot, Real const t) const
    {
      return -this->JF_x(x, -x_dot, -t);
    }

    /**
    * Time reversal of the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} =
    * -\mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{
    * \prime}, -t) \f$.
    */
    MatrixJF JF_x_dot_reverse(VectorF const & x, VectorF const & x_dot, Real const t) const
    {
      return this->JF_x_dot(x, -x_dot, -t);
    }

  }; // class Implicit

  /*\
   |   ___                 _ _      _ _ __        __
   |  |_ _|_ __ ___  _ __ | (_) ___(_) |\ \      / / __ __ _ _ __  _ __   ___ _ __
   |   | || '_ ` _ \| '_ \| | |/ __| | _ \ \ /\ / / '__/ _` | '_ \| '_ \ / _ \ '__|
   |   | || | | | | | |_) | | | (__| | |_ \ V  V /| | | (_| | |_) | |_) |  __/ |
   |  |___|_| |_| |_| .__/|_|_|\___|_|\__| \_/\_/ |_|  \__,_| .__/| .__/ \___|_|
   |                |_|                                       |_|   |_|
  \*/

  /**
  * \brief Class container for the system of implicit ODEs/DAEs wrapper.
  *
  * Class container for the system of implicit ordinary differential equations (ODEs) or differential
  * algebraic equations (DAEs) of the type \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
  * \mathbf{0} \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the implicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class ImplicitWrapper : public Implicit<Real, N, M>
  {
  public:
    SANDALS_BASIC_CONSTANTS(Real) /**< Basic constants. */

    using Pointer = std::shared_ptr<ImplicitWrapper<Real, N, M>>; /**< Shared pointer to an implicit ODE system. */
    using VectorF = typename Implicit<Real, N, M>::VectorF; /**< Templetized vector type. */
    using MatrixJF = typename Implicit<Real, N, M>::MatrixJF; /**< Templetized matrix type. */
    using VectorH = typename Implicit<Real, N, M>::VectorH; /**< Templetized vector type. */
    using MatrixJH = typename Implicit<Real, N, M>::MatrixJH; /**< Templetized matrix type. */
    using FunctionF = std::function<VectorF(VectorF const &, VectorF const &, Real const)>; /**< Implicit ODE system function type. */
    using FunctionJF = std::function<MatrixJF(VectorF const &, VectorF const &, Real const)>; /**< Jacobian of the ODE system function function type. */
    using FunctionH = std::function<VectorH(VectorF const &, Real const)>; /**< Invariants function type. */
    using FunctionJH = std::function<MatrixJH(VectorF const &, Real const)>; /**< Jacobian of the invariants function type. */
    using FunctionID = std::function<bool(VectorF const &, Real const)>; /**< In-domain function type. */

    inline static const FunctionH DefaultH = [](VectorF const &, Real const) {return VectorH::Zero();}; /**< Default mass matrix function. */
    inline static const FunctionJH DefaultJH = [](VectorF const &, Real const) {return MatrixJH::Zero();}; /**< Default system matrix function. */
    inline static const FunctionID DefaultID = [](VectorF const &, Real const) {return true;}; /**< Default in-domain function. */

  private:
    FunctionF m_F{nullptr}; /**< Implicit ODE system function. */
    FunctionJF m_JF_x{nullptr}; /**< Jacobian of the implicit ODE system function with respect to the states. */
    FunctionJF m_JF_x_dot{nullptr}; /**< Jacobian of the implicit ODE system function with respect to the states derivative. */
    FunctionH m_h{nullptr}; /**< System's invariants. */
    FunctionJH m_Jh_x{nullptr}; /**< Jacobian of the system's invariants with respect to the states. */
    FunctionID m_in_domain{nullptr}; /**< In-domain function. */

  public:
    /**
    * Class constructor for the implicit ODE/DAE system wrapper.
    * \param[in] t_F The implicit ODE system function.
    * \param[in] t_JF_x The Jacobian of the implicit ODE system function with respect to the states.
    * \param[in] t_JF_x_dot The Jacobian of the implicit ODE system function with respect to the states derivative.
    * \param[in] t_h The system's invariants.
    * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to the states.
    * \param[in] t_in_domain The in-domain function.
    */
    ImplicitWrapper(FunctionF t_F, FunctionJF t_JF_x, FunctionJF t_JF_x_dot, FunctionH t_h = DefaultH,
      FunctionJH t_Jh_x = DefaultJH, FunctionID t_in_domain = DefaultID) : Implicit<Real, N, M>(), m_F(t_F), m_JF_x(t_JF_x),
      m_JF_x_dot(t_JF_x_dot), m_h(t_h), m_Jh_x(t_Jh_x), m_in_domain(t_in_domain)
    {}

    /**
    * Class constructor for the implicit ODE/DAE system wrapper.
    * \param[in] t_name The name of the implicit ODE/DAE system.
    * \param[in] t_F The implicit ODE system function.
    * \param[in] t_JF_x The Jacobian of the implicit ODE system function with respect to the states.
    * \param[in] t_JF_x_dot The Jacobian of the implicit ODE system function with respect to the states derivative.
    * \param[in] t_h The implicit ODE system invariants.
    * \param[in] t_Jh_x The Jacobian of the implicit ODE system invariants with respect to the states.
    * \param[in] t_in_domain The in-domain function.
    */
    ImplicitWrapper(std::string t_name, FunctionF t_F, FunctionJF t_JF_x, FunctionJF t_JF_x_dot,
      FunctionH t_h = DefaultH, FunctionJH t_Jh_x = DefaultJH, FunctionID t_in_domain = DefaultID)
    : Implicit<Real, N, M>(t_name), m_F(t_F), m_JF_x(t_JF_x), m_JF_x_dot(t_JF_x_dot), m_h(t_h),
      m_Jh_x(t_Jh_x), m_in_domain(t_in_domain)
    {}

    /**
    * Class destructor for the implicit ODE/DAE system wrapper.
    */
    ~ImplicitWrapper() {}

    /**
    * Get the implicit ODE system function.
    * \return The implicit ODE system function.
    */
    FunctionF & F() {return this->m_F;}

    /**
    * Get the Jacobian of the implicit ODE system function with respect to the states.
    * \return The Jacobian of the implicit ODE system function with respect to the states.
    */
    FunctionJF & JF_x() {return this->m_JF_x;}

    /**
    * Get the Jacobian of the implicit ODE system function with respect to the states derivative.
    * \return The Jacobian of the implicit ODE system function with respect to the states derivative.
    */
    FunctionJF & JF_x_dot() {return this->m_JF_x_dot;}

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
    * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    VectorF F(VectorF const & x, VectorF const & x_dot, Real const t) const override
    {
      return this->m_F(x, x_dot, t);
    }

    /**
    * Evaluate the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle\frac{
    * \partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x(VectorF const & x, VectorF const & x_dot, Real const t) const override
    {
      return this->m_JF_x(x, x_dot, t);
    }

    /**
    * Evaluate the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}}
    * \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x_dot(VectorF const & x, VectorF const & x_dot, Real const t) const override
    {
      return this->m_JF_x_dot(x, x_dot, t);
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
    * Return true if the values \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$ is in the
    * domain of the ODE/DAE system.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return True if \f$ \mathbf{F}(\mathbf{x}, t) \f$ is in the domain of the ODE/DAE system.
    */
    bool in_domain(VectorF const & x, Real const t) const override
    {
      return this->m_in_domain(x, t);
    }

  }; // class ImplicitWrapper

} // namespace Sandals

#endif // SANDALS_IMPLICIT_SYSTEM_HH
