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

#ifndef SANDALS_SSPRK43_HXX
#define SANDALS_SSPRK43_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
  //!
  //! Butcher tableau for the 4-stage strong-stability preserving Runge-Kutta order 3 method:
  //!
  //! \f[
  //! \begin{array}{c|cccc}
  //!             0 &           0 &           0 &           0 &          0 \\
  //!   \frac{1}{2} & \frac{1}{2} &           0 &           0 &          0 \\
  //!             1 & \frac{1}{2} & \frac{1}{2} &           0 &          0 \\
  //!   \frac{1}{2} & \frac{1}{6} & \frac{1}{6} & \frac{1}{6} &          0 \\
  //!   \hline
  //!               & \frac{1}{6} & \frac{1}{6} & \frac{1}{6} & \frac{1}{2}
  //! \end{array} \text{.}
  //! \f]
  class SSPRK43Tableau : public Tableau<4>
  {
  public:
    using Tableau<4>::Type;
    using Tableau<4>::Vector;
    using Tableau<4>::Matrix;

    //! Class constructor for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
    SSPRK43Tableau() {
      this->name  = "SSPRK43";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,     0.0,     0.0,
                 1.0/2.0, 0.0,     0.0,     0.0,
                 1.0/2.0, 1.0/2.0, 0.0,     0.0,
                 1.0/6.0, 1.0/6.0, 1.0/6.0, 0.0;
      this->b << 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/2.0;
      this->c << 0.0,     1.0/2.0, 1.0,     1.0/2.0;
    }
  }; // class SSPRK43Tableau

  //! Class container for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
  template <Size N, Size M = 0>
  class SSPRK43 : public RungeKutta<4, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a 4-stage strong-stability preserving Runge-Kutta order 3 solver given
    //! a Tableau reference.
    SSPRK43() : RungeKutta<4, N, M>(SSPRK43Tableau()) {}

    //! Class constructor for a 4-stage strong-stability preserving Runge-Kutta order 3 solver given
    //! a Tableau reference.
    //! \param[in] t_system The system reference.
    SSPRK43(System t_system) : RungeKutta<4, N, M>(SSPRK43Tableau(), t_system) {}

  }; // class SSPRK43

} // namespace Sandals

#endif // SANDALS_SSPRK43_HXX
