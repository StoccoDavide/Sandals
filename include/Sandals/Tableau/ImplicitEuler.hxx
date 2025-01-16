/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_IMPLICITEULER_HXX
#define SANDALS_IMPLICITEULER_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the Implicit Euler method.
  //!
  //! Butcher tableau for the Implicit Euler method:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   1 & 1 \\
  //!   \hline
  //!     & 1
  //! \end{array} \text{.}
  //! \f]
  class ImplicitEulerTableau : public Tableau<1>
  {
  public:
    using Tableau<1>::Type;
    using Tableau<1>::Vector;
    using Tableau<1>::Matrix;

    //! Class constructor for the Implicit Euler method.
    ImplicitEulerTableau() {
      this->name  = "ImplicitEuler";
      this->type  = Type::IRK;
      this->order = 1;
      this->A << 1.0;
      this->b << 1.0;
      this->c << 1.0;
    }
  }; // class ImplicitEulerTableau

  //! Class container for the Implicit Euler method.
  //! \tparam N The dimension of the ODE/DAE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M = 0>
  class ImplicitEuler : public RungeKutta<1, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a Implicit Euler solver given a Tableau reference.
    ImplicitEuler() : RungeKutta<1, N, M>(ImplicitEulerTableau()) {}

    //! Class constructor for a Implicit Euler solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    ImplicitEuler(System t_system) : RungeKutta<1, N, M>(ImplicitEulerTableau(), t_system) {}

  }; // class ImplicitEuler

} // namespace Sandals

#endif // SANDALS_IMPLICITEULER_HXX
