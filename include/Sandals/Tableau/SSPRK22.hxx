/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_SSPRK22_HXX
#define SANDALS_SSPRK22_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the 2-stage strong-stability preserving Runge-Kutta (or Heun) order 2 method.
  //!
  //! Butcher tableau for the 2-stage strong-stability preserving Runge-Kutta (or Heun) order 2 method:
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
    using Tableau<2>::Type;
    using Tableau<2>::Vector;
    using Tableau<2>::Matrix;

    //! Class constructor for the 2-stage strong-stability preserving Runge-Kutta (or Heun) order 2 method.
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

  //! Class container for the 2-stage strong-stability preserving Runge-Kutta (or Heun) order 2 method.
  //! \tparam N The dimension of the ODE/DAE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M = 0>
  class SSPRK22 : public RungeKutta<2, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a 2-stage strong-stability preserving Runge-Kutta (or Heun) order 2
    //! solver given a Tableau reference.
    SSPRK22() : RungeKutta<2, N, M>(SSPRK22Tableau()) {}

    //! Class constructor for a 2-stage strong-stability preserving Runge-Kutta (or Heun) order 2
    //! solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    SSPRK22(System t_system) : RungeKutta<2, N, M>(SSPRK22Tableau(), t_system) {}

  }; // class SSPRK22

} // namespace Sandals

#endif // SANDALS_SSPRK22_HXX
