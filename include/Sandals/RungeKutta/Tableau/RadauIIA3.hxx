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

#ifndef SANDALS_RADAUIIA3_HXX
#define SANDALS_RADAUIIA3_HXX

namespace Sandals
{
  //! Butcher tableau for the RadauIIA order 3 method:
  //!
  //! \f[
  //! \begin{array}{c|cc}
  //!   1/3 & 5/12 & -1/12 \\
  //!     1 &  3/4 &   1/4 \\
  //!   \hline
  //!       &  3/4 &   1/4
  //! \end{array} \text{.}
  //! \f]
  class RadauIIA3Tableau : public Tableau<2>
  {
    public:
    using Tableau<2>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<2>::Vector; //!< Templetized vector type.
    using Tableau<2>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the RadauIIA order 3 method.
    RadauIIA3Tableau() {
      this->name  = "RadauIIA3";
      this->type  = Type::IRK;
      this->order = 3;
      this->A << 5.0/12.0, -1.0/12.0,
                 3.0/4.0,  1.0/4.0;
      this->b << 3.0/4.0, 1.0/4.0;
      this->c << 1.0/3.0, 1.0;
    }
  }; // class RadauIIA3Tableau

  //! Class container for the RadauIIA order 3 method.
  template <Size N>
  class RadauIIA3 : public RungeKutta<2, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a RadauIIA order 3 solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    RadauIIA3() : RungeKutta<2, N>(RadauIIA3Tableau()) {}

    //! Class constructor for a RadauIIA order 3 solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    //! \param[in] t_system The system reference.
    RadauIIA3(System t_system) : RungeKutta<2, N>(RadauIIA3Tableau(), t_system) {}

  }; // class RadauIIA3

} // namespace Sandals

#endif // SANDALS_RADAUIIA3_HXX
