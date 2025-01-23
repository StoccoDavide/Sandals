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

#ifndef SANDALS_SSPRK93_HXX
#define SANDALS_SSPRK93_HXX

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
  */
  class SSPRK93Tableau : public Tableau<9>
  {
  public:
    using Tableau<9>::Type;
    using Tableau<9>::Vector;
    using Tableau<9>::Matrix;

    /**
    * Class constructor for the 9-stage strong-stability preserving Runge-Kutta order 3 method.
    */
    SSPRK93Tableau() {
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
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <Integer N, Integer M = 0>
  class SSPRK93 : public RungeKutta<9, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a 9-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    */
    SSPRK93() : RungeKutta<9, N, M>(SSPRK93Tableau()) {}

    /**
    * Class constructor for a 9-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    * \param[in] t_system The system reference.
    */
    SSPRK93(System t_system) : RungeKutta<9, N, M>(SSPRK93Tableau(), t_system) {}

  }; // class SSPRK93

} // namespace Sandals

#endif // SANDALS_SSPRK93_HXX
