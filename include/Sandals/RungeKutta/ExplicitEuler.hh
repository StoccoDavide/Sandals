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

#ifndef SANDALS_EXPLICITEULER_HH
#define SANDALS_EXPLICITEULER_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Explicit Euler method.
  *
  * Butcher tableau for the Explicit Euler method:
  *
  * \f[
  * \begin{array}{c|c}
  *   0 & 0 \\
  *   \hline
  *     & 1
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class ExplicitEulerTableau : public Tableau<Real, 1>
  {
  public:
    using typename Tableau<Real, 1>::Type;
    using typename Tableau<Real, 1>::Vector;
    using typename Tableau<Real, 1>::Matrix;

    /**
    * Class constructor for the Explicit Euler method.
    */
    ExplicitEulerTableau<Real>() {
      this->name  = "ExplicitEuler";
      this->type  = Type::ERK;
      this->order = 1;
      this->A << 0.0;
      this->b << 1.0;
      this->c << 0.0;
    }
  }; // class ExplicitEulerTableau

  /**
  * \brief Class container for the Explicit Euler method.
  *
  * Class container for the Explicit Euler method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class ExplicitEuler : public RungeKutta<Real, 1, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Explicit Euler solver given a Tableau reference.
    */
    ExplicitEuler() : RungeKutta<Real, 1, N, M>(ExplicitEulerTableau<Real>()) {}

    /**
    * Class constructor for a Explicit Euler solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    ExplicitEuler(System t_system) : RungeKutta<Real, 1, N, M>(ExplicitEulerTableau<Real>(), t_system) {}

  }; // class ExplicitEuler

} // namespace Sandals

#endif // SANDALS_EXPLICITEULER_HH
