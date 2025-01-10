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

#ifndef SANDALS_MTE22_HXX
#define SANDALS_MTE22_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the Minimal Truncation Error order 2 method.
  //!
  //! Butcher tableau for the Minimal Truncation Error order 2 method:
  //!
  //! \f[
  //! \begin{array}{c|cc}
  //!     0 &   0 &   0 \\
  //!   2/3 & 2/3 &   0 \\
  //!   \hline
  //!       & 1/4 & 3/4
  //! \end{array} \text{.}
  //! \f]
  class MTE22Tableau : public Tableau<2>
  {
    public:
    using Tableau<2>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<2>::Vector; //!< Templetized vector type.
    using Tableau<2>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the Minimal Truncation Error order 2 method.
    MTE22Tableau() {
      this->name  = "MTE22";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0,     0.0,
                 2.0/3.0, 0.0;
      this->b << 1.0/4.0, 3.0/4.0;
      this->c << 0.0,     2.0/3.0;
    }
  }; // class MTE22Tableau

  //! Class container for the Minimal Truncation Error order 2 method.
  template <Size N>
  class MTE22 : public RungeKutta<2, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a Minimal Truncation Error order 2 solver given a Tableau reference.
    //!
    MTE22() : RungeKutta<2, N>(MTE22Tableau()) {}

    //! Class constructor for a Minimal Truncation Error order 2 solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    MTE22(System t_system) : RungeKutta<2, N>(MTE22Tableau(), t_system) {}

  }; // class MTE22

} // namespace Sandals

#endif // SANDALS_MTE22_HXX
