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

  //! \brief Class container for Runge-Kutta 4 method.
  //! Class container for the Runge-Kutta 4 method. The Butcher tableau is defined
  //! as:
  //!
  //! \f[
  //! \begin{array}{c|ccccc}
  //!     0 &   0 &   0 &   0 &   0 \\
  //!   1/2 & 1/2 &   0 &   0 &   0 \\
  //!   1/2 &   0 & 1/2 &   0 &   0 \\
  //!     1 &   0 &   0 &   1 &   0 \\
  //!   \hline
  //!       & 1/6 & 1/3 & 1/3 & 1/6
  //! \end{array} \text{.}
  //! \f]
  template <unsigned N>
  class RK4 : public RungeKutta<4, N>
  {
  public:
    //! Class constructor for the RK4 method.
    RK4(void)
    {
      Tableau<4> tableau = {
        "RK4", 4, false,
        {0,   0,   0, 0,
        1/2, 0,   0, 0,
        0,   1/2, 0, 0,
        0,   0,   1, 0},
        {1/6, 1/3, 1/3, 1/6},
        {NAN_VEC4},
        {0, 1/2, 1/2, 1}
      };
      this->RungeKutta("RK4", tableau);
    }
  }; // class RK4

} // namespace Sandalsù

#endif // SANDALS_RK4_HXX
