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

#ifndef SANDALS_GAUSSLEGENDRE2_HH
#define SANDALS_GAUSSLEGENDRE2_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /*
  Symplectic Runge-Kutta Methods of High Order Based on W-Transformation, Kaifeng Xia, Yuhao Cong
  and Geng Sun. Journal of Applied Analysis ans Computation, Volume 7, Number 3, 2017(8), 1185-1199
  http://www.jaac-online.com/data/article/jaac/preview/pdf/20170325.pdf
  */

  /**
  * \brief Butcher Tableau for the Gauss-Legendre order 2 method.
  *
  * Butcher tableau for the Gauss-Legendre order 2 method:
  *
  * \f[
  * \begin{array}{c|c}
  *   1/2 & 1/2 \\
  *   \hline
  *       & 1
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class GaussLegendre2Tableau : public Tableau<Real, 1>
  {
  public:
    using typename Tableau<Real, 1>::Type;
    using typename Tableau<Real, 1>::Vector;
    using typename Tableau<Real, 1>::Matrix;

    /**
    * Class constructor for the Gauss-Legendre order 2 method.
    */
    GaussLegendre2Tableau<Real>() {
      this->name  = "GaussLegendre2";
      this->type  = Type::IRK;
      this->order = 2;
      this->A << 1.0/2.0;
      this->b << 1.0;
      this->c << 1.0/2.0;
    }
  }; // class GaussLegendre2Tableau

  /**
  * \brief Class container for the Gauss-Legendre order 2 method.
  *
  * Class container for the Gauss-Legendre order 2 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class GaussLegendre2 : public RungeKutta<Real, 1, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Gauss-Legendre order 2 solver given a Tableau reference.
    */
    GaussLegendre2() : RungeKutta<Real, 1, N, M>(GaussLegendre2Tableau<Real>()) {}

    /**
    * Class constructor for a Gauss-Legendre order 2 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    GaussLegendre2(System t_system) : RungeKutta<Real, 1, N, M>(GaussLegendre2Tableau<Real>(), t_system) {}

  }; // class GaussLegendre2

} // namespace Sandals

#endif // SANDALS_GAUSSLEGENDRE2_HH
