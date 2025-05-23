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

#ifndef SANDALS_CHEBYSHEV51_HH
#define SANDALS_CHEBYSHEV51_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the 5-stage Runge-Kutta-Chebyshev order 1 method.
  *
  * Butcher tableau for the 5-stage Runge-Kutta-Chebyshev order 1 method:
  *
  * \f[
  * \begin{array}{c|ccccc}
  * & & & & & \\
  * \frac{1}{25} & \frac{1}{25} & & & & \\
  * \frac{4}{25} & \frac{2}{25} & \frac{2}{25} & & & \\
  * \frac{9}{25} & \frac{3}{25} & \frac{4}{25} & \frac{2}{25} & & \\
  * \frac{16}{25} & \frac{4}{25} & \frac{6}{25} & \frac{4}{25} & \frac{2}{25} & \\
  * \hline & \frac{1}{5} & \frac{8}{25} & \frac{6}{25} & \frac{4}{25} & \frac{2}{25}
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class Chebyshev51Tableau : public Tableau<Real, 5>
  {
  public:
    using typename Tableau<Real, 5>::Type;
    using typename Tableau<Real, 5>::Vector;
    using typename Tableau<Real, 5>::Matrix;

    /**
    * Class constructor for the 5-stage Runge-Kutta-Chebyshev order 1 method.
    */
    Chebyshev51Tableau<Real>() {
      this->name  = "Chebyshev51";
      this->type  = Type::ERK;
      this->order = 1;
      this->A << 0.0,      0.0,      0.0,      0.0,      0.0,
                 1.0/25.0, 0.0,      0.0,      0.0,      0.0,
                 2.0/25.0, 2.0/25.0, 0.0,      0.0,      0.0,
                 3.0/25.0, 4.0/25.0, 2.0/25.0, 0.0,      0.0,
                 4.0/25.0, 6.0/25.0, 4.0/25.0, 2.0/25.0, 0.0;
      this->b << 1.0/5.0, 8.0/25.0, 6.0/25.0, 4.0/25.0, 2.0/25.0;
      this->c << 0.0,     1.0/25.0, 4.0/25.0, 9.0/25.0, 16.0/25.0;
    }
  }; // class Chebyshev51Tableau

  /**
  * \brief Class container for the 5-stage Runge-Kutta-Chebyshev order 1 method.
  *
  * Class container for the 5-stage Runge-Kutta-Chebyshev order 1 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class Chebyshev51 : public RungeKutta<Real, 5, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a 5-stage Runge-Kutta-Chebyshev order 1 solver given a Tableau reference.
    */
    Chebyshev51() : RungeKutta<Real, 5, N, M>(Chebyshev51Tableau<Real>()) {}

    /**
    * Class constructor for a 5-stage Runge-Kutta-Chebyshev order 1 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    Chebyshev51(System t_system) : RungeKutta<Real, 5, N, M>(Chebyshev51Tableau<Real>(), t_system) {}

  }; // class Chebyshev51

} // namespace Sandals

#endif // SANDALS_CHEBYSHEV51_HH
