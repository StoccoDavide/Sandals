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

#ifndef SANDALS_SSPRK104_HH
#define SANDALS_SSPRK104_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the 10-stage strong-stability preserving Runge-Kutta order 4 method.
  *
  * Butcher tableau for the 10-stage strong-stability preserving Runge-Kutta order 4 method:
  *
  * \f[
  * \begin{array}{c|cccccccccc}
  *             0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 \\
  *   \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 \\
  *   \frac{1}{3} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 &            0 &            0 &            0 &            0 &            0 \\
  *   \frac{1}{2} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 &            0 &            0 &            0 &            0 \\
  *   \frac{2}{3} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 &            0 &            0 &            0 \\
  *   \frac{1}{3} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} &            0 &            0 &            0 &            0 &            0 \\
  *   \frac{1}{2} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} &  \frac{1}{6} &            0 &            0 &            0 &            0 \\
  *   \frac{2}{3} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 &            0 \\
  *   \frac{5}{6} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 &            0 \\
  *             1 & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} & \frac{1}{15} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &  \frac{1}{6} &            0 \\
  *   \hline
  *               & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10} & \frac{1}{10}
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class SSPRK104Tableau : public Tableau<Real, 10>
  {
  public:
    using typename Tableau<Real, 10>::Type;
    using typename Tableau<Real, 10>::Vector;
    using typename Tableau<Real, 10>::Matrix;

    /**
    * Class constructor for the 10-stage strong-stability preserving Runge-Kutta order 4 method.
    */
    SSPRK104Tableau<Real>() {
      this->name  = "SSPRK104";
      this->type  = Type::ERK;
      this->order = 4;
      this->A << 0.0,      0.0,      0.0,      0.0,      0.0,      0.0,     0.0,     0.0,     0.0,     0.0,
                 1.0/6.0,  0.0,      0.0,      0.0,      0.0,      0.0,     0.0,     0.0,     0.0,     0.0,
                 1.0/6.0,  1.0/6.0,  0.0,      0.0,      0.0,      0.0,     0.0,     0.0,     0.0,     0.0,
                 1.0/6.0,  1.0/6.0,  1.0/6.0,  0.0,      0.0,      0.0,     0.0,     0.0,     0.0,     0.0,
                 1.0/6.0,  1.0/6.0,  1.0/6.0,  1.0/6.0,  0.0,      0.0,     0.0,     0.0,     0.0,     0.0,
                 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 0.0,     0.0,     0.0,     0.0,     0.0,
                 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 0.0,     0.0,     0.0,     0.0,
                 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 1.0/6.0, 0.0,     0.0,     0.0,
                 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 0.0,     0.0,
                 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/15.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 0.0;
      this->b << 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0, 1.0/10.0;
      this->c << 0.0,      1.0/6.0,  1.0/3.0,  1.0/2.0,  2.0/3.0,  1.0/3.0, 1.0/2.0, 2.0/3.0, 5.0/6.0, 1.0;
    }
  }; // class SSPRK104Tableau

  /**
  * \brief Class container for the 10-stage strong-stability preserving Runge-Kutta order 4 method.
  *
  * Class container for the 10-stage strong-stability preserving Runge-Kutta order 4 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SSPRK104 : public RungeKutta<Real, 10, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a 10-stage strong-stability preserving Runge-Kutta order 4 solver given
    * a Tableau reference.
    */
    SSPRK104() : RungeKutta<Real, 10, N, M>(SSPRK104Tableau<Real>()) {}

    /**
    * Class constructor for a 10-stage strong-stability preserving Runge-Kutta order 4 solver given
    * a Tableau reference.
    * \param[in] t_system The system reference.
    */
    SSPRK104(System t_system) : RungeKutta<Real, 10, N, M>(SSPRK104Tableau<Real>(), t_system) {}

  }; // class SSPRK104

} // namespace Sandals

#endif // SANDALS_SSPRK104_HH
