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

#ifndef SANDALS_HEUN3_HXX
#define SANDALS_HEUN3_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the Heun's order 3 method.
  //!
  //! Butcher tableau for the Heun's order 3 method:
  //!
  //! \f[
  //! \begin{array}{c|ccc}
  //!     0 &   0 &   0 &   0 \\
  //!   1/3 & 1/3 &   0 &   0 \\
  //!   2/3 &   0 & 2/3 &   0 \\
  //!   \hline
  //!       & 1/4 &   0 & 3/4
  //! \end{array} \text{.}
  //! \f]
  class Heun3Tableau : public Tableau<3>
  {
  public:
    using Tableau<3>::Type;
    using Tableau<3>::Vector;
    using Tableau<3>::Matrix;

    //! Class constructor for the Heun's order 3 method.
    Heun3Tableau() {
      this->name  = "Heun3";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,     0.0,
                 1.0/3.0, 0.0,     0.0,
                 0.0,     2.0/3.0, 0.0;
      this->b << 1.0/4.0, 0.0,     3.0/4.0;
      this->c << 0.0,     1.0/3.0, 2.0/3.0;
    }
  }; // class Heun3Tableau

  //! Class container for the Heun's order 3 method.
  template <Size N, Size M = 0>
  class Heun3 : public RungeKutta<3, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a Heun's order 3 solver given a Tableau reference.
    //!
    Heun3() : RungeKutta<3, N, M>(Heun3Tableau()) {}

    //! Class constructor for a Heun's order 3 solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    Heun3(System t_system) : RungeKutta<3, N, M>(Heun3Tableau(), t_system) {}

  }; // class Heun3

} // namespace Sandals

#endif // SANDALS_HEUN3_HXX
