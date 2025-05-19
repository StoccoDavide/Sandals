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

#ifndef SANDALS_SSPRK22STAR_HH
#define SANDALS_SSPRK22STAR_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the optimal 2-stage strong-stability preserving Runge-Kutta order 2
  * method.
  *
  * Butcher tableau for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method:
  *
  * \f[
  * \begin{array}{c|cc}
  *                   0 &                 0 & 0 \\
  *   0.822875655532364 & 0.822875655532364 & 0 \\
  *   \hline
  *   & 0.392374781489287 & 0.607625218510713
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class SSPRK22starTableau : public Tableau<Real, 2>
  {
  public:
    using typename Tableau<Real, 2>::Type;
    using typename Tableau<Real, 2>::Vector;
    using typename Tableau<Real, 2>::Matrix;

    /**
    * Class constructor for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method.
    */
    SSPRK22starTableau<Real>() {
      this->name  = "SSPRK22star";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0, 0.0,
                 0.822875655532364, 0.0;
      this->b << 0.392374781489287, 0.607625218510713;
      this->c << 0.0, 0.822875655532364;
    }
  }; // class SSPRK22starTableau

  /**
  * \brief Class container for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method.
  *
  * Class container for the optimal 2-stage strong-stability preserving Runge-Kutta order 2 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SSPRK22star : public RungeKutta<Real, 2, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for an optimal  2-stage strong-stability preserving Runge-Kutta order 2
    * solver given a Tableau reference.
    */
    SSPRK22star() : RungeKutta<Real, 2, N, M>(SSPRK22starTableau<Real>()) {}

    /**
    * Class constructor for an optimal  2-stage strong-stability preserving Runge-Kutta order 2
    * solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    SSPRK22star(System t_system) : RungeKutta<Real, 2, N, M>(SSPRK22starTableau<Real>(), t_system) {}

  }; // class SSPRK22star

} // namespace Sandals

#endif // SANDALS_SSPRK22STAR_HH
