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

#ifndef SANDALS_GAUSSLEGENDRE2_HXX
#define SANDALS_GAUSSLEGENDRE2_HXX

namespace Sandals
{
  // Symplectic Runge-Kutta Methods of High Order Based on W-Transformation, Kaifeng Xia, Yuhao Cong
  // and Geng Sun. Journal of Applied Analysis ans Computation, Volume 7, Number 3, 2017(8), 1185-1199
  // http://www.jaac-online.com/data/article/jaac/preview/pdf/20170325.pdf

  //! \brief Butcher Tableau for the Gauss-Legendre order 2 method.
  //!
  //! Butcher tableau for the Gauss-Legendre order 2 method:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   1/2 & 1/2 \\
  //!   \hline
  //!       & 1
  //! \end{array} \text{.}
  //! \f]
  class GaussLegendre2Tableau : public Tableau<1>
  {
  public:
    using Tableau<1>::Type;
    using Tableau<1>::Vector;
    using Tableau<1>::Matrix;

    //! Class constructor for the Gauss-Legendre order 2 method.
    //!
    GaussLegendre2Tableau() {
      this->name  = "GaussLegendre2";
      this->type  = Type::IRK;
      this->order = 2;
      this->A << 1.0/2.0;
      this->b << 1.0;
      this->c << 1.0/2.0;
    }
  }; // class GaussLegendre2Tableau

  //! Class container for the Gauss-Legendre order 2 method.
  //! \tparam N The dimension of the ODE/DAE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M = 0>
  class GaussLegendre2 : public RungeKutta<1, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a Gauss-Legendre order 2 solver given a Tableau reference.
    GaussLegendre2() : RungeKutta<1, N, M>(GaussLegendre2Tableau()) {}

    //! Class constructor for a Gauss-Legendre order 2 solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    GaussLegendre2(System t_system) : RungeKutta<1, N, M>(GaussLegendre2Tableau(), t_system) {}

  }; // class GaussLegendre2

} // namespace Sandals

#endif // SANDALS_GAUSSLEGENDRE2_HXX
