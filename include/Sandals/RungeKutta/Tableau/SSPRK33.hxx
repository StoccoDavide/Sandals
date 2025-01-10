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

#ifndef SANDALS_SSPRK33_HXX
#define SANDALS_SSPRK33_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the 3-stage strong-stability preserving Runge-Kutta order 3 method.
  //!
  //! Butcher tableau for the 3-stage strong-stability preserving Runge-Kutta order 3 method:
  //!
  //! \f[
  //! \begin{array}{c|ccc}
  //!     0 &   0 &   0 &   0 \\
  //!     1 &   1 &   0 &   0 \\
  //!   1/2 & 1/4 & 1/4 &   0 \\
  //!   \hline
  //!       & 1/6 & 1/6 & 2/3
  //! \end{array} \text{.}
  //! \f]
  class SSPRK33Tableau : public Tableau<3>
  {
    public:
    using Tableau<3>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<3>::Vector; //!< Templetized vector type.
    using Tableau<3>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the 3-stage strong-stability preserving Runge-Kutta order 3 method.
    SSPRK33Tableau() {
      this->name  = "SSPRK33";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,     0.0,
                 1.0,     0.0,     0.0,
                 1.0/4.0, 1.0/4.0, 0.0;
      this->b << 1.0/6.0, 1.0/6.0, 2.0/3.0;
      this->c << 0.0,     1.0,     1.0/2.0;
    }
  }; // class SSPRK33Tableau

  //! Class container for the 3-stage strong-stability preserving Runge-Kutta order 3 method.
  template <Size N>
  class SSPRK33 : public RungeKutta<3, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a 3-stage strong-stability preserving Runge-Kutta order 3 solver given
    //! a Tableau reference.
    SSPRK33() : RungeKutta<3, N>(SSPRK33Tableau()) {}

    //! Class constructor for a 3-stage strong-stability preserving Runge-Kutta order 3 solver given
    //! a Tableau reference.
    //! \param[in] t_system The system reference.
    SSPRK33(System t_system) : RungeKutta<3, N>(SSPRK33Tableau(), t_system) {}

  }; // class SSPRK33

} // namespace Sandals

#endif // SANDALS_SSPRK33_HXX
