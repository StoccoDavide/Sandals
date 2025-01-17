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

#ifndef SANDALS_HEUN2_HXX
#define SANDALS_HEUN2_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method.
  //!
  //! Butcher tableau for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method:
  //!
  //! \f[
  //! \begin{array}{c|cc}
  //!   0 &   0 &   0 \\
  //!   1 &   1 &   0 \\
  //!   \hline
  //!     & 1/2 & 1/2
  //! \end{array} \text{.}
  //! \f]
  class Heun2Tableau : public Tableau<2>
  {
  public:
    using Tableau<2>::Type;
    using Tableau<2>::Vector;
    using Tableau<2>::Matrix;

    //! Class constructor for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method.
    Heun2Tableau() {
      this->name  = "Heun2";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0, 0.0,
                 1.0, 0.0;
      this->b << 1.0/2.0, 1.0/2.0;
      this->c << 0.0,     1.0;
    }
  }; // class Heun2Tableau

  //! Class container for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2
  //! method.
  //! \tparam N The dimension of the ODE/DAE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M = 0>
  class Heun2 : public RungeKutta<2, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2
    //! solver given a Tableau reference.
    Heun2() : RungeKutta<2, N, M>(Heun2Tableau()) {}

    //! Class constructor for a Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2
    //! solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    Heun2(System t_system) : RungeKutta<2, N, M>(Heun2Tableau(), t_system) {}

  }; // class Heun2

} // namespace Sandals

#endif // SANDALS_HEUN2_HXX
