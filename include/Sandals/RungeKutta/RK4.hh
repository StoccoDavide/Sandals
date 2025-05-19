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

#ifndef SANDALS_RK4_HH
#define SANDALS_RK4_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Runge-Kutta order 4 method.
  *
  * Butcher tableau for the Runge-Kutta order 4 method:
  *
  * \f[
  * \begin{array}{c|cccc}
  *     0 &   0 &   0 &   0 &   0 \\
  *   1/2 & 1/2 &   0 &   0 &   0 \\
  *   1/2 &   0 & 1/2 &   0 &   0 \\
  *     1 &   0 &   0 &   1 &   0 \\
  *   \hline
  *       & 1/6 & 1/3 & 1/3 & 1/6
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class RK4Tableau : public Tableau<Real, 4>
  {
  public:
    using typename Tableau<Real, 4>::Type;
    using typename Tableau<Real, 4>::Vector;
    using typename Tableau<Real, 4>::Matrix;

    /**
    * Class constructor for the RK4 method.
    */
    RK4Tableau<Real>() {
      this->name  = "RK4";
      this->type  = Type::ERK;
      this->order = 4;
      this->A << 0.0,     0.0,     0.0,     0.0,
                 1.0/2.0, 0.0,     0.0,     0.0,
                 0.0,     1.0/2.0, 0.0,     0.0,
                 0.0,     0.0,     1.0,     0.0;
      this->b << 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0;
      this->c << 0.0,     1.0/2.0, 1.0/2.0, 1.0;
    }
  }; // class RK4Tableau

  /**
  * \brief Class container for the Runge-Kutta order 4 method.
  *
  * Class container for the Runge-Kutta order 4 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class RK4 : public RungeKutta<Real, 4, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Runge-Kutta order 4 solver given a Tableau reference.
    */
    RK4() : RungeKutta<Real, 4, N, M>(RK4Tableau<Real>()) {}

    /**
    * Class constructor for a Runge-Kutta order 4 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    RK4(System t_system) : RungeKutta<Real, 4, N, M>(RK4Tableau<Real>(), t_system) {}

  }; // class RK4

} // namespace Sandals

#endif // SANDALS_RK4_HH
