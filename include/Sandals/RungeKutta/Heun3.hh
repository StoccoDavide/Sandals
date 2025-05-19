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

#ifndef SANDALS_HEUN3_HH
#define SANDALS_HEUN3_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Heun's order 3 method.
  *
  * Butcher tableau for the Heun's order 3 method:
  *
  * \f[
  * \begin{array}{c|ccc}
  *     0 &   0 &   0 &   0 \\
  *   1/3 & 1/3 &   0 &   0 \\
  *   2/3 &   0 & 2/3 &   0 \\
  *   \hline
  *       & 1/4 &   0 & 3/4
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class Heun3Tableau : public Tableau<Real, 3>
  {
  public:
    using typename Tableau<Real, 3>::Type;
    using typename Tableau<Real, 3>::Vector;
    using typename Tableau<Real, 3>::Matrix;

    /**
    * Class constructor for the Heun's order 3 method.
    */
    Heun3Tableau<Real>() {
      this->name  = "Heun3";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,     0.0,
                 1.0/3.0, 0.0,     0.0,
                 0.0,     2.0/3.0, 0.0;
      this->b << 1.0/4.0, 0.0,     3.0/4.0;
      this->c << 0.0,     1.0/3.0, 2.0/3.0;
    }
  }; // class Heun3Tableau

  /**
  * \brief Class container for the Heun's order 3 method.
  *
  * Class container for the Heun's order 3 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class Heun3 : public RungeKutta<Real, 3, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Heun's order 3 solver given a Tableau reference.
    */
    Heun3() : RungeKutta<Real, 3, N, M>(Heun3Tableau<Real>()) {}

    /**
    * Class constructor for a Heun's order 3 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    Heun3(System t_system) : RungeKutta<Real, 3, N, M>(Heun3Tableau<Real>(), t_system) {}

  }; // class Heun3

} // namespace Sandals

#endif // SANDALS_HEUN3_HH
