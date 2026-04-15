/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_IMPLICITEULER_HH
#define SANDALS_IMPLICITEULER_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals {
  /**
   * \brief Butcher tableau for the Implicit Euler method.
   *
   * Butcher tableau for the Implicit Euler method:
   *
   * \f[
   * \begin{array}{c|c}
   *   1 & 1 \\
   *   \hline
   *     & 1
   * \end{array} \text{.}
   * \f]
   * \tparam Real The scalar number type.
   */
  template <typename Real>
  class ImplicitEulerTableau : public Tableau<Real, 1> {
   public:
    using typename Tableau<Real, 1>::Type;
    using typename Tableau<Real, 1>::VectorS;
    using typename Tableau<Real, 1>::MatrixS;

    /**
     * Class constructor for the Implicit Euler method.
     */
    ImplicitEulerTableau<Real>() {
      this->name  = "ImplicitEuler";
      this->type  = Type::IRK;
      this->order = 1;
      this->A << 1.0;
      this->b << 1.0;
      this->c << 1.0;
    }

  };  // class ImplicitEulerTableau

  /**
   * \brief Class container for the Implicit Euler method.
   *
   * Class container for the Implicit Euler method.
   * \tparam Real The scalar number type.
   * \tparam N The dimension of the ODE/DAE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class ImplicitEuler : public RungeKutta<Real, 1, N, M> {
   public:
    using SystemPtr = typename Implicit<Real, N, M>::Pointer;

    /**
     * Class constructor for a Implicit Euler solver given a Tableau reference.
     */
    ImplicitEuler() : RungeKutta<Real, 1, N, M>(ImplicitEulerTableau<Real>()) {}

    /**
     * Class constructor for a Implicit Euler solver given a Tableau reference.
     * \param[in] t_system The system reference.
     */
    ImplicitEuler(SystemPtr t_system)
        : RungeKutta<Real, 1, N, M>(ImplicitEulerTableau<Real>(),
                                    std::move(t_system)) {}

  };  // class ImplicitEuler

}  // namespace Sandals

#endif  // SANDALS_IMPLICITEULER_HH
