/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_IMPLICIT_SYSTEM_HXX
#define SANDALS_IMPLICIT_SYSTEM_HXX

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

  //! \brief Class container for the system of implicit ODEs/DAEs.
  //!
  //! Class container for the system of inplicit ordinary differential equations (ODEs) or differential
  //! algebraic equations (DAEs) of the type \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
  //! \mathbf{0} \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  //! \tparam N The dimension of the implicit ODE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M>
  class Implicit
  {
  public:
    using Type = enum class Type : Size {IMPLICIT=0, EXPLICIT=1, SEMIEXPLICIT=1}; //!< System type enumeration.
    using Pointer  = std::shared_ptr<Implicit<N, M>>; //!< Shared pointer to an implicit ODE system.
    using VectorF  = Eigen::Vector<Real, N>;          //!< Templetized vector type.
    using MatrixJF = Eigen::Matrix<Real, N, N>;       //!< Templetized matrix type.
    using VectorH  = Eigen::Vector<Real, M>;          //!< Templetized vector type.
    using MatrixJH = Eigen::Matrix<Real, M, N>;       //!< Templetized matrix type.

  private:
    Type        m_type{Type::IMPLICIT}; //!< ODE/DAE system type.
    std::string m_name;                 //!< Name of the ODE/DAE system.

  protected:
    //! Class constructor for an implicit ODE/DAE system.
    //! \param[in] t_type The type of the implicit ODE/DAE system.
    //! \param[in] t_name The name of the implicit ODE/DAE system.
    Implicit(Type t_type, std::string t_name) : m_type(t_type), m_name(t_name) {}

  public:
    //! Class constructor for an implicit ODE/DAE system.
    Implicit() : m_type(Type::IMPLICIT), m_name("(missing name)") {}

    //! Class constructor for an implicit ODE/DAE system.
    //! \param[in] t_name The name of the implicit ODE/DAE system.
    Implicit(std::string t_name) : m_type(Type::IMPLICIT), m_name(t_name) {}

    //! Class destructor for an implicit ODE/DAE system.
    virtual ~Implicit() {}

    //! Get the enumeration type of the ODE/DAE system.
    //! \return The enumeration type of the ODE/DAE system.
    Type type() const {return this->m_type;}

    //! Check if the ODE/DAE system is implicit.
    //! \return True if the ODE/DAE system is implicit, false otherwise.
    bool is_implicit() const {return this->m_type == Type::IMPLICIT;}

    //! Check if the ODE/DAE system is explicit.
    //! \return True if the ODE/DAE system is explicit, false otherwise.
    bool is_explicit() const {return this->m_type == Type::EXPLICIT;}

    //! Check if the ODE/DAE system is semi-explicit.
    //! \return True if the ODE/DAE system is semi-explicit, false otherwise.
    bool is_semiexplicit() const {return this->m_type == Type::SEMIEXPLICIT;}

    //! Get the ODE/DAE system name reference.
    //! \return The ODE/DAE system name reference.
    std::string & name() {return this->m_name;}

    //! Get the ODE/DAE system name const reference.
    //! \return The ODE/DAE system name const reference.
    std::string const & name() const {return this->m_name;}

    //! Get the number of equations of the ODE/DAE system.
    //! \return The number of equations of the ODE/DAE system.
    Size equations_number() const {return N;}

    //! Get the number of invariants of the ODE/DAE system.
    //! \return The number of invariants of the ODE/DAE system.
    Size invariants_number() const {return M;}

    //! Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual VectorF F(VectorF const &x, VectorF const &x_dot, Real t) const = 0;

    //! Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x},
    //! \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle\frac{
    //! \partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual MatrixJF JF_x(VectorF const &x, VectorF const &x_dot, Real t) const = 0;

    //! Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x},
    //! \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    //!
    //! \f[
    //! \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    //! \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}}
    //! \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual MatrixJF JF_x_dot(VectorF const &x, VectorF const &x_dot, Real t) const = 0;

    //! Evaluate the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
    virtual VectorH h(VectorF const &x, Real t) const = 0;

    //! Evaluate the Jacobian of the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$
    //! with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{Jh}_{\mathbf{x}}(\mathbf{x} t) = \displaystyle\frac{\partial\mathbf{h}(\mathbf{x},
    //! t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jh}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual MatrixJH Jh_x(VectorF const &x, Real t) const = 0;

    //! Return true if the values \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$ is in the
    //! domain of the ODE/DAE system.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return True if \f$ \mathbf{F}(\mathbf{x}, t) \f$ is in the domain of the ODE/DAE system.
    virtual bool in_domain(VectorF const &x, Real t) const = 0;

    //! Time reversal of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{
    //! \prime}, t) = -\mathbf{F}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The time-reversed system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, -t) \f$.
    VectorF F_reverse(VectorF const &x, VectorF const &x_dot, Real t) const
    {
      return -this->F(x, -x_dot, -t);
    }

    //! Time reversal of the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    //! \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} = -\mathbf{JF}_{\mathbf{x}}
    //! (\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The time-reversed Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, -\mathbf{x}^{\prime},
    //! -t) \f$.
    MatrixJF JF_x_reverse(VectorF const &x, VectorF const &x_dot, Real t) const
    {
      return -this->JF_x(x, -x_dot, -t);
    }

    //! Time reversal of the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    //! \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} =
    //! -\mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The time-reversed Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{
    //! \prime}, -t) \f$.
    MatrixJF JF_x_dot_reverse(VectorF const &x, VectorF const &x_dot, Real t) const
    {
      return this->JF_x_dot(x, -x_dot, -t);
    }

  }; // class Implicit

} // namespace Sandals

#endif // SANDALS_ImPLICIT_SYSTEM_HXX
