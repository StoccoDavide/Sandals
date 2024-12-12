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
    using ptr = std::shared_ptr<Explicit>; //!< Shared pointer to an explicit ODE system.
    using matN = Implicit<N>; //!< Matrix type.
    using vecN = Implicit<N>; //!< Column vector type.

    //! Class constructor for an explicit ODE system.
    //! \param t_name The name of the explicit ODE system.
    Explicit(std::string t_name) : Implicit<N>(t_name) {}

    //! Evaluate the explicit ODE system \f$ \mathbf{f} \f$.
    //! \param x States \f$ \mathbf{x} \f$.
    //! \param t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{f} \f$.
    virtual vecN f(vecN const &x, real t) const;

    //! Evaluate the Jacobian of the explicit ODE system function \f$ \mathbf{F} \f$
    //! with respect to the states \f$ \mathbf{x} \f$$:
    //!
    //! \f[
    //!   \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
    //!   \dfrac{\partial\mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param x States \f$ \mathbf{x} \f$.
    //! \param t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    virtual matN Jf_x(vecN const &x, real t) const;

  }; // class Explicit

} // namespace Sandals

#endif // SANDALS_EXPLICIT_SYSTEM_HXX
