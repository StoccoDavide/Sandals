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

#ifndef SANDALS_RALSTON3_HXX
#define SANDALS_RALSTON3_HXX

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Ralston's order 3 method.
  *
  * Butcher tableau for the Ralston's order 3 method:
  *
  * \f[
  * \begin{array}{c|ccc}
  *     0 &   0 &   0 &   0 \\
  *   1/2 & 1/2 &   0 &   0 \\
  *   3/4 &   0 & 3/4 &   0 \\
  *   \hline
  *       & 2/9 & 1/3 & 4/9
  * \end{array} \text{.}
  * \f]
  */
  class Ralston3Tableau : public Tableau<3>
  {
  public:
    using Tableau<3>::Type;
    using Tableau<3>::Vector;
    using Tableau<3>::Matrix;

    /**
    * Class constructor for the Ralston's order 3 method.
    */
    Ralston3Tableau() {
      this->name  = "Ralston3";
      this->type  = Type::ERK;
      this->order = 3;
      this->A << 0.0,     0.0,     0.0,
                 1.0/2.0, 0.0,     0.0,
                 0.0,     3.0/4.0, 0.0;
      this->b << 2.0/9.0, 1.0/3.0, 4.0/9.0;
      this->c << 0.0,     1.0/2.0, 3.0/4.0;
    }
  }; // class Ralston3Tableau

  /**
  * \brief Class container for the Ralston's order 3 method.
  *
  * Class container for the Ralston's order 3 method.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <Integer N, Integer M = 0>
  class Ralston3 : public RungeKutta<3, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Ralston's order 3 solver given a Tableau reference.
    */
    Ralston3() : RungeKutta<3, N, M>(Ralston3Tableau()) {}

    /**
    * Class constructor for a Ralston's order 3 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    Ralston3(System t_system) : RungeKutta<3, N, M>(Ralston3Tableau(), t_system) {}

  }; // class Ralston3

} // namespace Sandals

#endif // SANDALS_RALSTON3_HXX
