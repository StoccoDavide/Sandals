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

  //! \brief Class container for the system of explicit ODEs.
  //!
  //! Class container for the system of explicit ordinary differential equations (ODEs) of the type
  //! \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, t) \f$.
  template <Size N>
  class Explicit : public Implicit<N>
  {
  public:
    using Pointer = std::shared_ptr<Explicit<N>>; //!< Shared pointer to an explicit ODE system.
    using Vector  = typename Implicit<N>::Vector; //!< Templetized vector type.
    using Matrix  = typename Implicit<N>::Matrix; //!< Templetized matrix type.
    using Type    = typename Implicit<N>::Type;   //!< ODE system type enumeration.

    //! Default class constructor for the explicit ODE system.
    //!
    Explicit() : Implicit<N>(Type::EXPLICIT, "(missing name)") {}

    //! Class constructor for the explicit ODE system.
    //! \param[in] t_name The name of the explicit ODE system.
    Explicit(std::string t_name) : Implicit<N>(Type::EXPLICIT, t_name) {}

    //! Evaluate the ODE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    //!
    //! \f[ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{x}^{\prime} - \mathbf{f}(\mathbf{x}, t) \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    Vector F(Vector const &x, Vector const &x_dot, Real t) const override
    {
      return x_dot - this->f(x, t);
    }

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    //! t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //! \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    Matrix JF_x(Vector const &x, Vector const &/*x_dot*/, Real t) const override
    {
      return -this->Jf_x(x, t);
    }

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    //! \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    //! \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    //! \mathbf{I} \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    Matrix JF_x_dot(Vector const &/*x*/, Vector const &/*x_dot*/, Real /*t*/) const override
    {
      return Matrix::Identity();
    }

    //! Evaluate the explicit ODE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    virtual Vector f(Vector const &x, Real t) const = 0;

    //! Evaluate the Jacobian of the explicit ODE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    //! with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
    //! \displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    virtual Matrix Jf_x(Vector const &x, Real t) const = 0;

  }; // class Explicit

} // namespace Sandals

#endif // SANDALS_EXPLICIT_SYSTEM_HXX
