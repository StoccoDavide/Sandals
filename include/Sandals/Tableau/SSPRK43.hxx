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

#ifndef SANDALS_SSPRK43_HXX
#define SANDALS_SSPRK43_HXX

namespace Sandals
{
  /**
  * \brief Butcher tableau for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
  *
  * Butcher tableau for the 4-stage strong-stability preserving Runge-Kutta order 3 method:
  *
  * \f[
  * \begin{array}{c|cccc}
  *             0 &           0 &           0 &           0 &          0 \\
  *   \frac{1}{2} & \frac{1}{2} &           0 &           0 &          0 \\
  *             1 & \frac{1}{2} & \frac{1}{2} &           0 &          0 \\
  *   \frac{1}{2} & \frac{1}{6} & \frac{1}{6} & \frac{1}{6} &          0 \\
  *   \hline
  *               & \frac{1}{6} & \frac{1}{6} & \frac{1}{6} & \frac{1}{2}
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class SSPRK43Tableau : public Tableau<Real, 4>
  {
  public:
    using typename Tableau<Real, 4>::Type;
    using typename Tableau<Real, 4>::Vector;
    using typename Tableau<Real, 4>::Matrix;

    /**
    * Class constructor for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
    */
    SSPRK43Tableau<Real>() {
      this->name  = "SSPRK43";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,     0.0,     0.0,
                 1.0/2.0, 0.0,     0.0,     0.0,
                 1.0/2.0, 1.0/2.0, 0.0,     0.0,
                 1.0/6.0, 1.0/6.0, 1.0/6.0, 0.0;
      this->b << 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/2.0;
      this->c << 0.0,     1.0/2.0, 1.0,     1.0/2.0;
    }
  }; // class SSPRK43Tableau

  /**
  * \brief Class container for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
  *
  * Class container for the 4-stage strong-stability preserving Runge-Kutta order 3 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SSPRK43 : public RungeKutta<Real, 4, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a 4-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    */
    SSPRK43() : RungeKutta<Real, 4, N, M>(SSPRK43Tableau<Real>()) {}

    /**
    * Class constructor for a 4-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    * \param[in] t_system The system reference.
    */
    SSPRK43(System t_system) : RungeKutta<Real, 4, N, M>(SSPRK43Tableau<Real>(), t_system) {}

  }; // class SSPRK43

} // namespace Sandals

#endif // SANDALS_SSPRK43_HXX
