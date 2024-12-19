/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
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

  //! \brief Class container for the system of implicit ODEs/DAEs.
  //!
  //! Class container for the system of inplicit ordinary differential equations (ODEs) or differential
  //! algebraic equations (DAEs) of the type \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
  //! \mathbf{0} \f$.
  template <Size N>
  class Implicit
  {
  public:
    using Pointer = std::shared_ptr<Implicit<N>>; //!< Shared pointer to an implicit ODE system.
    using Vector  = Eigen::Vector<Real, N>;       //!< Templetized vector type.
    using Matrix  = Eigen::Matrix<Real, N, N>;    //!< Templetized matrix type.
    using Type    = enum class Type : Size {
      IMPLICIT = 0, EXPLICIT = 1, SEMI_EXPLICIT = 2
    }; //!< ODE system type enumeration.

  private:
    Type        m_type{Type::IMPLICIT}; //!< ODE system type.
    std::string m_name;                 //! Name of the ODE system.

  protected:
    //! Class constructor for an implicit ODE system.
    //! \param[in] t_type The type of the implicit ODE system.
    //! \param[in] t_name The name of the implicit ODE system.
    Implicit(Type t_type, std::string t_name) : m_type(t_type), m_name(t_name) {}

  public:
    //! Default class constructor for an implicit ODE system.
    //!
    Implicit() : m_type(Type::IMPLICIT), m_name("(missing name)") {}

    //! Class constructor for an implicit ODE system.
    //! \param[in] t_name The name of the implicit ODE system.
    Implicit(std::string t_name) : m_type(Type::IMPLICIT), m_name(t_name) {}

    //! Class destructor for an implicit ODE system.
    virtual ~Implicit() {}

    //! Check if the ODE system is implicit.
    //! \return True if the ODE system is implicit, false otherwise.
    bool is_implicit() const {return this->m_type == Type::IMPLICIT;}

    //! Check if the ODE system is explicit.
    //! \return True if the ODE system is explicit, false otherwise.
    bool is_explicit() const {return this->m_type == Type::EXPLICIT;}

    //! Check if the ODE system is semi-explicit.
    //! \return True if the ODE system is semi-explicit, false otherwise.
    bool is_semiexplicit() const {return this->m_type == Type::SEMI_EXPLICIT;}

    //! Get the ODE system name reference.
    //! \return The ODE system name reference.
    std::string & name() {return this->m_name;}

    //! Get the ODE system name const reference.
    //! \return The ODE system name const reference.
    std::string const & name() const {return this->m_name;}

    //! Get the number of equations of the ODE system.
    //! \return The number of equations of the ODE system.
    Size dimension() const {return N;}

    //! Evaluate the ODE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    //! \param[in] x     States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t     Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual Vector F(Vector const &x, Vector const &x_dot, Real t) const = 0;

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    //! t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle\frac{
    //! \partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} \text{.} \f]
    //!
    //! \param[in] x     States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t     Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual Matrix JF_x(Vector const &x, Vector const &x_dot, Real t) const = 0;

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    //! t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    //! \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}}
    //! \text{.} \f]
    //!
    //! \param[in] x     States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t     Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual Matrix JF_x_dot(Vector const &x, Vector const &x_dot, Real t) const = 0;

    //! Return true if the values \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$ is in the
    //! domain of the ODE system.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return True if \f$ \mathbf{F}(\mathbf{x}, t) \f$ is in the domain of the ODE system.
    virtual bool in_domain(Vector const &x, Real t) const = 0;

  }; // class Implicit

} // namespace Sandals

#endif // SANDALS_ImPLICIT_SYSTEM_HXX
