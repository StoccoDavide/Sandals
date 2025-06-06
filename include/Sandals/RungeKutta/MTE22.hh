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

#ifndef SANDALS_MTE22_HH
#define SANDALS_MTE22_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the minimal truncation error (or Ralston's) order 2 method.
  *
  * Butcher tableau for the minimal truncation error (or Ralston's) order 2 method:
  *
  * \f[
  * \begin{array}{c|cc}
  *     0 &   0 &   0 \\
  *   2/3 & 2/3 &   0 \\
  *   \hline
  *       & 1/4 & 3/4
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class MTE22Tableau : public Tableau<Real, 2>
  {
  public:
    using typename Tableau<Real, 2>::Type;
    using typename Tableau<Real, 2>::Vector;
    using typename Tableau<Real, 2>::Matrix;

    /**
    * Class constructor for the minimal truncation error (or Ralston's) order 2 method.
    */
    MTE22Tableau<Real>() {
      this->name  = "MTE22";
      this->type  = Type::ERK;
      this->order = 2;
      this->A << 0.0,     0.0,
                 2.0/3.0, 0.0;
      this->b << 1.0/4.0, 3.0/4.0;
      this->c << 0.0,     2.0/3.0;
    }
  }; // class MTE22Tableau

  /**
  * \brief Class container for the minimal truncation error (or Ralston's) order 2 method.
  *
  * Class container for the minimal truncation error (or Ralston's) order 2 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class MTE22 : public RungeKutta<Real, 2, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a minimal truncation error order 2 solver given a Tableau reference.
    */
    MTE22() : RungeKutta<Real, 2, N, M>(MTE22Tableau<Real>()) {}

    /**
    * Class constructor for a minimal truncation error order 2 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    MTE22(System t_system) : RungeKutta<Real, 2, N, M>(MTE22Tableau<Real>(), t_system) {}

  }; // class MTE22

} // namespace Sandals

#endif // SANDALS_MTE22_HH
