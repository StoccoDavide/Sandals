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

#ifndef SANDALS_HEUN2_HH
#define SANDALS_HEUN2_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method.
  *
  * Butcher tableau for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method:
  *
  * \f[
  * \begin{array}{c|cc}
  *   0 &   0 &   0 \\
  *   1 &   1 &   0 \\
  *   \hline
  *     & 1/2 & 1/2
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class Heun2Tableau : public Tableau<Real, 2>
  {
  public:
    using typename Tableau<Real, 2>::Type;
    using typename Tableau<Real, 2>::Vector;
    using typename Tableau<Real, 2>::Matrix;

    /**
    * Class constructor for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method.
    */
    Heun2Tableau<Real>() {
      this->name  = "Heun2";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0, 0.0,
                 1.0, 0.0;
      this->b << 1.0/2.0, 1.0/2.0;
      this->c << 0.0,     1.0;
    }
  }; // class Heun2Tableau

  /**
  * \brief Class container for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method.
  *
  * Class container for the Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class Heun2 : public RungeKutta<Real, 2, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2
    * solver given a Tableau reference.
    */
    Heun2() : RungeKutta<Real, 2, N, M>(Heun2Tableau<Real>()) {}

    /**
    * Class constructor for a Heun's (or 2-stage strong-stability preserving Runge-Kutta) order 2
    * solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    Heun2(System t_system) : RungeKutta<Real, 2, N, M>(Heun2Tableau<Real>(), t_system) {}

  }; // class Heun2

} // namespace Sandals

#endif // SANDALS_HEUN2_HH
