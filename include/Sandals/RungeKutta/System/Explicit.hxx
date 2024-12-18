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

#ifndef SANDALS_EXPLICIT_SYSTEM_HXX
#define SANDALS_EXPLICIT_SYSTEM_HXX

namespace Sandals {

  //! Class container for the system of explicit Ordinary Differential Equations
  //! (ODEs) of the type \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, t) \f$.
  template <unsigned N>
  class Explicit : public Implicit<N>
  {
  public:
    using ptr  = std::shared_ptr<Explicit<N>>; //!< Shared pointer to an explicit ODE system.
    using vec  = typename Implicit<N>::vec;    //!< Templetized vector type.
    using mat  = typename Implicit<N>::mat;    //!< Templetized matrix type.
    using type = typename Implicit<N>::type;   //!< ODE type enumeration.

    //! Default class constructor for an implicit ODE system.
    Explicit(void) : Implicit<N>(type::EXPLICIT, "(missing name)") {}

    //! Class constructor for an explicit ODE system.
    //! \param t_name The name of the explicit ODE system.
    Explicit(std::string t_name) : Implicit<N>(type::EXPLICIT, t_name) {}

    //! Evaluate the explicit ODE system \f$ \mathbf{f} \f$.
    //! \param x States \f$ \mathbf{x} \f$.
    //! \param t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{f} \f$.
    virtual vec f(vec const &x, real t) const = 0;

    //! Evaluate the Jacobian of the explicit ODE system function \f$ \mathbf{F} \f$
    //! with respect to the states \f$ \mathbf{x} \f$$:
    //!
    //! \f[
    //!   \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
    //!   \frac{\partial\mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param x States \f$ \mathbf{x} \f$.
    //! \param t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual mat Jf_x(vec const &x, real t) const = 0;

        //! Evaluate the ODE system \f$ \mathbf{F} \f$.
    //! \param x     States \f$ \mathbf{x} \f$.
    //! \param x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param t     Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{F} \f$.
    virtual vec F(vec const &x, vec const &x_dot, real t) const override
    {
      return x_dot - this->f(x, t);
    }

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
    mat JF_x(vec const &x, vec const &/*x_dot*/, real t) const override
    {
      return -this->Jf_x(x, t);
    }

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
    mat JF_x_dot(vec const &/*x*/, vec const &/*x_dot*/, real /*t*/) const override
    {
      return mat::Identity();
    }

  }; // class Explicit

} // namespace Sandals

#endif // SANDALS_EXPLICIT_SYSTEM_HXX
