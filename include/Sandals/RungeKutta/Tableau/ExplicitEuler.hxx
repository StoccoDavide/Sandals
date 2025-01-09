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

#ifndef SANDALS_EXPLICITEULER_HXX
#define SANDALS_EXPLICITEULER_HXX

namespace Sandals
{
  //! Butcher tableau for the Explicit Euler method:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   0 & 0 \\
  //!   \hline
  //!     & 1
  //! \end{array} \text{.}
  //! \f]
  class ExplicitEulerTableau : public Tableau<1>
  {
    public:
    using Tableau<1>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<1>::Vector; //!< Templetized vector type.
    using Tableau<1>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the Explicit Euler method.
    ExplicitEulerTableau() {
      this->name  = "ExplicitEuler";
      this->type  = Type::ERK;
      this->order = 1;
      this->A << 0.0;
      this->b << 1.0;
      this->c << 0.0;
    }
  }; // class ExplicitEulerTableau

  //! Class container for the Explicit Euler method.
  template <Size N>
  class ExplicitEuler : public RungeKutta<1, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a Explicit Euler solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    ExplicitEuler() : RungeKutta<1, N>(ExplicitEulerTableau()) {}

    //! Class constructor for a Explicit Euler solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    //! \param[in] t_system The system reference.
    ExplicitEuler(System t_system) : RungeKutta<1, N>(ExplicitEulerTableau(), t_system) {}

  }; // class ExplicitEuler

} // namespace Sandals

#endif // SANDALS_EXPLICITEULER_HXX
