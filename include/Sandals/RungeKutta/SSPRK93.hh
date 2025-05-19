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

#ifndef SANDALS_SSPRK93_HH
#define SANDALS_SSPRK93_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the 9-stage strong-stability preserving Runge-Kutta order 3 method.
  *
  * Butcher tableau for the 9-stage strong-stability preserving Runge-Kutta order 3 method:
  *
  * \f[
  * \begin{array}{c|ccccccccc}
  *             0 &           0 &            0 &            0 &            0 &            0 &            0 &           0 &           0 &           0 \\
  *   \frac{1}{6} & \frac{1}{6} &            0 &            0 &            0 &            0 &            0 &           0 &           0 &           0 \\
  *   \frac{1}{3} & \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 &            0 &           0 &           0 &           0 \\
  *   \frac{1}{2} & \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 &           0 &           0 &           0 \\
  *   \frac{2}{3} & \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 &           0 &           0 &           0 \\
  *   \frac{5}{6} & \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 &           0 &           0 &           0 \\
  *   \frac{1}{2} & \frac{1}{6} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} &           0 &           0 &           0 \\
  *   \frac{2}{3} & \frac{1}{6} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{6} &           0 &           0 \\
  *   \frac{5}{6} & \frac{1}{6} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{6} & \frac{1}{6} &           0 \\
  *   \hline
  *               & \frac{1}{6} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{6} & \frac{1}{6} & \frac{1}{6}
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class SSPRK93Tableau : public Tableau<Real, 9>
  {
  public:
    using typename Tableau<Real, 9>::Type;
    using typename Tableau<Real, 9>::Vector;
    using typename Tableau<Real, 9>::Matrix;

    /**
    * Class constructor for the 9-stage strong-stability preserving Runge-Kutta order 3 method.
    */
    SSPRK93Tableau<Real>() {
      this->name  = "SSPRK93";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,     0.0,
                 1.0/6.0, 0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,     0.0,
                 1.0/6.0, 1.0/6.0,  0.0,      0.0,      0.0,      0.0,      0.0,     0.0,     0.0,
                 1.0/6.0, 1.0/6.0,  1.0/6.0,  0.0,      0.0,      0.0,      0.0,     0.0,     0.0,
                 1.0/6.0, 1.0/6.0,  1.0/6.0,  1.0/6.0,  0.0,      0.0,      0.0,     0.0,     0.0,
                 1.0/6.0, 1.0/6.0,  1.0/6.0,  1.0/6.0,  1.0/6.0,  0.0,      0.0,     0.0,     0.0,
                 1.0/6.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 0.0,     0.0,     0.0,
                 1.0/6.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 0.0,     0.0,
                 1.0/6.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 1.0/6.0, 0.0;
      this->b << 1.0/6.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 1.0/6.0, 1.0/6.0;
      this->c << 0.0,     1.0/6.0,  1.0/3.0,  1.0/2.0,  2.0/3.0,  5.0/6.0,  1.0/2.0, 2.0/3.0, 5.0/6.0;
        }
      }; // class SSPRK93Tableau

  /**
  * \brief Class container for the 9-stage strong-stability preserving Runge-Kutta order 3 method.
  *
  * Class container for the 9-stage strong-stability preserving Runge-Kutta order 3 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SSPRK93 : public RungeKutta<Real, 9, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a 9-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    */
    SSPRK93() : RungeKutta<Real, 9, N, M>(SSPRK93Tableau<Real>()) {}

    /**
    * Class constructor for a 9-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    * \param[in] t_system The system reference.
    */
    SSPRK93(System t_system) : RungeKutta<Real, 9, N, M>(SSPRK93Tableau<Real>(), t_system) {}

  }; // class SSPRK93

} // namespace Sandals

#endif // SANDALS_SSPRK93_HH
