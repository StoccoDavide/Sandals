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

#ifndef SANDALS_IMPLICITEULER_HXX
#define SANDALS_IMPLICITEULER_HXX

namespace Sandals
{
  //! Butcher tableau for the Implicit Euler method:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   1 & 1 \\
  //!   \hline
  //!     & 1
  //! \end{array} \text{.}
  //! \f]
  class ImplicitEulerTableau : public Tableau<1>
  {
    public:
    using Tableau<1>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<1>::Vector; //!< Templetized vector type.
    using Tableau<1>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the Implicit Euler method.
    ImplicitEulerTableau() {
      this->name  = "ImplicitEuler";
      this->type  = Type::IRK;
      this->order = 1;
      this->A << 1.0;
      this->b << 1.0;
      this->c << 1.0;
    }
  }; // class ImplicitEulerTableau

  //! Class container for the Implicit Euler method.
  template <Size N>
  class ImplicitEuler : public RungeKutta<1, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a Implicit Euler solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    ImplicitEuler() : RungeKutta<1, N>(ImplicitEulerTableau()) {}

    //! Class constructor for a Implicit Euler solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    //! \param[in] t_system The system reference.
    ImplicitEuler(System t_system) : RungeKutta<1, N>(ImplicitEulerTableau(), t_system) {}

  }; // class ImplicitEuler

} // namespace Sandals

#endif // SANDALS_IMPLICITEULER_HXX
