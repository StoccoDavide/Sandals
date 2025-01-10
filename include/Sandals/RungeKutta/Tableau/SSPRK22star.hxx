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

#ifndef SANDALS_SSPRK22STAR_HXX
#define SANDALS_SSPRK22STAR_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the optimal 2-stage strong-stability preserving Runge-Kutta order 2
  //! method.
  //!
  //! Butcher tableau for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method:
  //!
  //! \f[
  //! \begin{array}{c|cc}
  //!                   0 &                 0 & 0 \\
  //!   0.822875655532364 & 0.822875655532364 & 0 \\
  //!   \hline
  //!   & 0.392374781489287 & 0.607625218510713
  //! \end{array} \text{.}
  //! \f]
  class SSPRK22starTableau : public Tableau<2>
  {
    public:
    using Tableau<2>::Type;   //!< Runge-Kutta type enumeration.
    using Tableau<2>::Vector; //!< Templetized vector type.
    using Tableau<2>::Matrix; //!< Templetized matrix type.

    //! Class constructor for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method.
    SSPRK22starTableau() {
      this->name  = "SSPRK22star";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0, 0.0,
                 0.822875655532364, 0.0;
      this->b << 0.392374781489287, 0.607625218510713;
      this->c << 0.0, 0.822875655532364;
    }
  }; // class SSPRK22starTableau

  //! Class container for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method.
  template <Size N>
  class SSPRK22star : public RungeKutta<2, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for an optimal  2-stage strong-stability preserving Runge-Kutta order 2
    //! solver given a Tableau reference.
    SSPRK22star() : RungeKutta<2, N>(SSPRK22starTableau()) {}

    //! Class constructor for an optimal  2-stage strong-stability preserving Runge-Kutta order 2
    //! solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    SSPRK22star(System t_system) : RungeKutta<2, N>(SSPRK22starTableau(), t_system) {}

  }; // class SSPRK22star

} // namespace Sandals

#endif // SANDALS_SSPRK22STAR_HXX
