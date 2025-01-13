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

#ifndef SANDALS_CHEBYSHEV51_HXX
#define SANDALS_CHEBYSHEV51_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the 5-stage Runge-Kutta-Chebyshev order 1 method.
  //!
  //! Butcher tableau for the 5-stage Runge-Kutta-Chebyshev order 1 method:
  //!
  //! \f[
  //! \begin{array}{c|ccccc}
  //! & & & & & \\
  //! \frac{1}{25} & \frac{1}{25} & & & & \\
  //! \frac{4}{25} & \frac{2}{25} & \frac{2}{25} & & & \\
  //! \frac{9}{25} & \frac{3}{25} & \frac{4}{25} & \frac{2}{25} & & \\
  //! \frac{16}{25} & \frac{4}{25} & \frac{6}{25} & \frac{4}{25} & \frac{2}{25} & \\
  //! \hline & \frac{1}{5} & \frac{8}{25} & \frac{6}{25} & \frac{4}{25} & \frac{2}{25}
  //! \end{array} \text{.}
  //! \f]
  class Chebyshev51Tableau : public Tableau<5>
  {
  public:
    using Tableau<5>::Type;
    using Tableau<5>::Vector;
    using Tableau<5>::Matrix;

    //! Class constructor for the 5-stage Runge-Kutta-Chebyshev order 1 method.
    Chebyshev51Tableau() {
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

  //! Class container for the 5-stage Runge-Kutta-Chebyshev order 1 method.
  template <Size N, Size M = 0>
  class Chebyshev51 : public RungeKutta<5, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a 5-stage Runge-Kutta-Chebyshev order 1 solver given a Tableau reference.
    Chebyshev51() : RungeKutta<5, N, M>(Chebyshev51Tableau()) {}

    //! Class constructor for a 5-stage Runge-Kutta-Chebyshev order 1 solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    Chebyshev51(System t_system) : RungeKutta<5, N, M>(Chebyshev51Tableau(), t_system) {}

  }; // class Chebyshev51

} // namespace Sandals

#endif // SANDALS_CHEBYSHEV51_HXX
