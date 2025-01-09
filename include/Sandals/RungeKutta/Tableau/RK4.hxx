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

#ifndef SANDALS_RK4_HXX
#define SANDALS_RK4_HXX

namespace Sandals
{
  //! Butcher tableau for the Runge-Kutta 4 method:
  //!
  //! \f[
  //! \begin{array}{c|cccc}
  //!     0 &   0 &   0 &   0 &   0 \\
  //!   1/2 & 1/2 &   0 &   0 &   0 \\
  //!   1/2 &   0 & 1/2 &   0 &   0 \\
  //!     1 &   0 &   0 &   1 &   0 \\
  //!   \hline
  //!       & 1/6 & 1/3 & 1/3 & 1/6
  //! \end{array} \text{.}
  //! \f]
  class RK4Tableau : public Tableau<4>
  {
    public:
    using Tableau<4>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<4>::Vector; //!< Templetized vector type.
    using Tableau<4>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the RK4 method.
    RK4Tableau() {
      this->name  = "RK4";
      this->type  = Type::ERK;
      this->order = 4;
      this->A << 0.0,     0.0,     0.0,     0.0,
                 1.0/2.0, 0.0,     0.0,     0.0,
                 0.0,     1.0/2.0, 0.0,     0.0,
                 0.0,     0.0,     1.0,     0.0;
      this->b << 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0;
      this->c << 0.0,     1.0/2.0, 1.0/2.0, 1.0;
    }
  }; // class RK4Tableau

  //! Class container for the Runge-Kutta 4 method.
  template <Size N>
  class RK4 : public RungeKutta<4, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a Runge-Kutta 4 solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    RK4() : RungeKutta<4, N>(RK4Tableau()) {}

    //! Class constructor for a Runge-Kutta 4 solver given a Tableau reference.
    //! \param[in] t_tableau The Tableau reference.
    //! \param[in] t_system The system reference.
    RK4(System t_system) : RungeKutta<4, N>(RK4Tableau(), t_system) {}

  }; // class RK4

} // namespace Sandals

#endif // SANDALS_RK4_HXX
