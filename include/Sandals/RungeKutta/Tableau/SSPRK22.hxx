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

#ifndef SANDALS_SSPRK22_HXX
#define SANDALS_SSPRK22_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the 2-stage strong-stability preserving Runge-Kutta order 2 method.
  //!
  //! Butcher tableau for the 2-stage strong-stability preserving Runge-Kutta order 2 method:
  //!
  //! \f[
  //! \begin{array}{c|cc}
  //!   0 &   0 &   0 \\
  //!   1 &   1 &   0 \\
  //!   \hline
  //!     & 1/2 & 1/2
  //! \end{array} \text{.}
  //! \f]
  class SSPRK22Tableau : public Tableau<2>
  {
    public:
    using Tableau<2>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<2>::Vector; //!< Templetized vector type.
    using Tableau<2>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the 2-stage strong-stability preserving Runge-Kutta order 2 method.
    SSPRK22Tableau() {
      this->name  = "SSPRK22";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0, 0.0,
                 1.0, 0.0;
      this->b << 1.0/2.0, 1.0/2.0;
      this->c << 0.0,     1.0;
    }
  }; // class SSPRK22Tableau

  //! Class container for the 2-stage strong-stability preserving Runge-Kutta order 2 method.
  template <Size N>
  class SSPRK22 : public RungeKutta<2, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a 2-stage strong-stability preserving Runge-Kutta order 2 solver given
    //! a Tableau reference.
    SSPRK22() : RungeKutta<2, N>(SSPRK22Tableau()) {}

    //! Class constructor for a 2-stage strong-stability preserving Runge-Kutta order 2 solver given
    //! a Tableau reference.
    //! \param[in] t_system The system reference.
    SSPRK22(System t_system) : RungeKutta<2, N>(SSPRK22Tableau(), t_system) {}

  }; // class SSPRK22

} // namespace Sandals

#endif // SANDALS_SSPRK22_HXX
