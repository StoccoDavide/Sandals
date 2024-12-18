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

  //! Class container for the system of implicit Ordinary Differential Equations
  //! (ODEs) of the type \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
  //! \mathbf{0} \f$.
  template <unsigned N>
  class Implicit
  {
  public:
    using ptr  = std::shared_ptr<Implicit<N>>; //!< Shared pointer to an implicit ODE system.
    using vec  = Eigen::Vector<real, N>;       //!< Templetized vector type.
    using mat  = Eigen::Matrix<real, N, N>;    //!< Templetized matrix type.
    using type = enum class type : unsigned {
      IMPLICIT = 0, EXPLICIT = 1,
    }; //!< ODE type enumeration.

  private:
    type        m_type{type::IMPLICIT}; //!< ODE system type.
    std::string m_name;                 //! Name of the ODE system.

  protected:
    //! Class constructor for an implicit ODE system.
    //! \param t_type The type of the implicit ODE system.
    //! \param t_name The name of the implicit ODE system.
    Implicit(type t_type, std::string t_name) : m_type(t_type), m_name(t_name) {}

  public:
    //! Default class constructor for an implicit ODE system.
    Implicit(void) : m_type(type::IMPLICIT), m_name("(missing name)") {}

    //! Class constructor for an implicit ODE system.
    //! \param t_name The name of the implicit ODE system.
    Implicit(std::string t_name) : m_type(type::IMPLICIT), m_name(t_name) {}

    //! Class destructor for an implicit ODE system.
    virtual ~Implicit(void) {}

    //! Check if the ODE system is implicit.
    //! \return True if the ODE system is implicit, false otherwise.
    bool is_implicit(void) const {return this->m_type == type::IMPLICIT;}

    //! Check if the ODE system is explicit.
    //! \return True if the ODE system is explicit, false otherwise.
    bool is_explicit(void) const {return this->m_type == type::EXPLICIT;}

    //! Get the ODE system name reference.
    //! \return The ODE system name reference.
    std::string & name(void) {return this->m_name;}

    //! Get the ODE system name const reference.
    //! \return The ODE system name const reference.
    std::string const & name(void) const {return this->m_name;}

    //! Get the number of equations of the ODE system.
    //! \return The number of equations of the ODE system.
    unsigned dimension(void) const {return N;}

    //! Evaluate the ODE system \f$ \mathbf{F} \f$.
    //! \param x     States \f$ \mathbf{x} \f$.
    //! \param x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param t     Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{F} \f$.
    virtual vec F(vec const &x, vec const &x_dot, real t) const = 0;

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F} \f$
    //! with respect to the states \f$ \mathbf{x} \f$$:
    //!
    //! \f[
    //!   \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //!   \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param x     States \f$ \mathbf{x} \f$.
    //! \param x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param t     Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual mat JF_x(vec const &x, vec const &x_dot, real t) const = 0;

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F} \f$
    //! with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$$:
    //!
    //! \f[
    //!   \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //!   \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} \text{.}
    //! \f]
    //!
    //! \param x     States \f$ \mathbf{x} \f$.
    //! \param x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param t     Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual mat JF_x_dot(vec const &x, vec const &x_dot, real t) const = 0;

    //! Return true if \f$ (\mathbf{x}, t) \f$ is in the domain of the ODE system.
    //! \param x States \f$ \mathbf{x} \f$.
    //! \param t Independent variable \f$ t \f$.
    //! \return True if \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the domain of the ODE system.
    virtual bool in_domain(vec const &x, real t) const = 0;

  }; // class Implicit

} // namespace Sandals

#endif // SANDALS_ImPLICIT_SYSTEM_HXX
