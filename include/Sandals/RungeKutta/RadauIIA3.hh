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

#ifndef SANDALS_RADAUIIA3_HH
#define SANDALS_RADAUIIA3_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Radau IIA order 3 method.
  *
  * Butcher tableau for the Radau IIA order 3 method:
  *
  * \f[
  * \begin{array}{c|cc}
  *   1/3 & 5/12 & -1/12 \\
  *     1 &  3/4 &   1/4 \\
  *   \hline
  *       &  3/4 &   1/4
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class RadauIIA3Tableau : public Tableau<Real, 2>
  {
  public:
    using typename Tableau<Real, 2>::Type;
    using typename Tableau<Real, 2>::Vector;
    using typename Tableau<Real, 2>::Matrix;

    /**
    * Class constructor for the Radau IIA order 3 method.
    */
    RadauIIA3Tableau<Real>() {
      this->name  = "RadauIIA3";
      this->type  = Type::IRK;
      this->order = 3;
      this->A << 5.0/12.0, -1.0/12.0,
                 3.0/4.0,  1.0/4.0;
      this->b << 3.0/4.0, 1.0/4.0;
      this->c << 1.0/3.0, 1.0;
    }
  }; // class RadauIIA3Tableau

  /**
  * \brief Class container for the Radau IIA order 3 method.
  *
  * Class container for the Radau IIA order 3 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class RadauIIA3 : public RungeKutta<Real, 2, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Radau IIA order 3 solver given a Tableau reference.
    */
    RadauIIA3() : RungeKutta<Real, 2, N, M>(RadauIIA3Tableau<Real>()) {}

    /**
    * Class constructor for a Radau IIA order 3 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    RadauIIA3(System t_system) : RungeKutta<Real, 2, N, M>(RadauIIA3Tableau<Real>(), t_system) {}

  }; // class RadauIIA3

} // namespace Sandals

#endif // SANDALS_RADAUIIA3_HH
