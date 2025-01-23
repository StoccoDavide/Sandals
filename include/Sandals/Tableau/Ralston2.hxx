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

#ifndef SANDALS_RALSTON2_HXX
#define SANDALS_RALSTON2_HXX

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Ralston's (or minimal truncation error) order 2 method.
  *
  * Butcher tableau for the Ralston's (or minimal truncation error) order 2 method:
  *
  * \f[
  * \begin{array}{c|cc}
  *     0 &   0 &   0 \\
  *   2/3 & 2/3 &   0 \\
  *   \hline
  *       & 1/4 & 3/4
  * \end{array} \text{.}
  * \f]
  */
  class Ralston2Tableau : public Tableau<2>
  {
  public:
    using Tableau<2>::Type;
    using Tableau<2>::Vector;
    using Tableau<2>::Matrix;

    /**
    * Class constructor for the Ralston's (or minimal truncation error) order 2 method.
    */
    Ralston2Tableau() {
      this->name  = "Ralston2";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0,     0.0,
                 2.0/3.0, 0.0;
      this->b << 1.0/4.0, 3.0/4.0;
      this->c << 0.0,     2.0/3.0;
    }
  }; // class Ralston2Tableau

  /**
  * \brief Class container for the Ralston's (or minimal truncation error) order 2 method.
  *
  * Class container for the Ralston's (or minimal truncation error) order 2 method.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <Integer N, Integer M = 0>
  class Ralston2 : public RungeKutta<2, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Ralston's (or minimal truncation error) order 2 solver given a
    * Tableau reference.
    */
    Ralston2() : RungeKutta<2, N, M>(Ralston2Tableau()) {}

    /**
    * Class constructor for a Ralston's (or minimal truncation error) order 2 solver given a
    * Tableau reference.
    * \param[in] t_system The system reference.
    */
    Ralston2(System t_system) : RungeKutta<2, N, M>(Ralston2Tableau(), t_system) {}

  }; // class Ralston2

} // namespace Sandals

#endif // SANDALS_RALSTON2_HXX
