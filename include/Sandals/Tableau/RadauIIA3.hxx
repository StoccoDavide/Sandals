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

#ifndef SANDALS_RADAUIIA3_HXX
#define SANDALS_RADAUIIA3_HXX

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
  */
  class RadauIIA3Tableau : public Tableau<2>
  {
  public:
    using Tableau<2>::Type;
    using Tableau<2>::Vector;
    using Tableau<2>::Matrix;

    /**
    * Class constructor for the Radau IIA order 3 method.
    */
    RadauIIA3Tableau() {
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
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <Size N, Size M = 0>
  class RadauIIA3 : public RungeKutta<2, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Radau IIA order 3 solver given a Tableau reference.
    */
    RadauIIA3() : RungeKutta<2, N, M>(RadauIIA3Tableau()) {}

    /**
    * Class constructor for a Radau IIA order 3 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    RadauIIA3(System t_system) : RungeKutta<2, N, M>(RadauIIA3Tableau(), t_system) {}

  }; // class RadauIIA3

} // namespace Sandals

#endif // SANDALS_RADAUIIA3_HXX
