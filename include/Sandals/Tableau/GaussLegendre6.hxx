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

#ifndef SANDALS_GAUSSLEGENDRE6_HXX
#define SANDALS_GAUSSLEGENDRE6_HXX

namespace Sandals
{
  // Symplectic Runge-Kutta Methods of High Order Based on W-Transformation, Kaifeng Xia, Yuhao Cong
  // and Geng Sun. Journal of Applied Analysis ans Computation, Volume 7, Number 3, 2017(8), 1185-1199
  // http://www.jaac-online.com/data/article/jaac/preview/pdf/20170325.pdf

  //! \brief Butcher tableau for the Gauss-Legendre order 6 method.
  //!
  //! Butcher tableau for the Gauss-Legendre order 6 method:
  //!
  //! \f[
  //! \begin{array}{c|ccc}
  //!   1/2-t_1 &     w & z-t_2 & w-t_4 \\
  //!       1/2 & w+t_3 &     z & w-t_3 \\
  //!   1/2+t_1 & w+t_4 & z+t_2 &     w \\
  //!   \hline
  //!         &  5/18 &   4/9 &  5/18
  //! \end{array} \text{,}
  //! \f]
  //!
  //! with \f$ t_1 = \sqrt{15}/10 \f$, \f$ t_2 = \sqrt{15}/15 \f$, \f$ t_3 = \sqrt{15}/24 \f$, \f$
  //! t_4 = \sqrt{15}/30 \f$, \f$ w = 5/36 \f$, and \f$ z = 2/9 \f$.
  class GaussLegendre6Tableau : public Tableau<3>
  {
  public:
    using Tableau<3>::Type;
    using Tableau<3>::Vector;
    using Tableau<3>::Matrix;

    //! Class constructor for the Gauss-Legendre order 6 method.
    GaussLegendre6Tableau() {
      this->name  = "GaussLegendre6";
      this->type  = Type::IRK;
      this->order = 6;
      Real t_1{std::sqrt(15.0)/10.0};
      Real t_2{std::sqrt(15.0)/15.0};
      Real t_3{std::sqrt(15.0)/24.0};
      Real t_4{std::sqrt(15.0)/30.0};
      Real w{5.0/36.0};
      Real z{2.0/9.0};
      this->A << w,           z-t_2,   w-t_4,
                 w+t_3,       z,       w-t_3,
                 w+t_4,       z+t_2,   w;
      this->b << 5.0/18.0,    4.0/9.0, 5.0/18.0;
      this->c << 1.0/2.0-t_1, 1.0/2.0, 1.0/2.0+t_1;
    }
  }; // class GaussLegendre6Tableau

  //! Class container for the Gauss-Legendre order 6 method.
  //! \tparam N The dimension of the ODE/DAE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M = 0>
  class GaussLegendre6 : public RungeKutta<3, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    GaussLegendre6() : RungeKutta<3, N, M>(GaussLegendre6Tableau()) {}

    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    GaussLegendre6(System t_system) : RungeKutta<3, N, M>(GaussLegendre6Tableau(), t_system) {}

  }; // class GaussLegendre6

} // namespace Sandals

#endif // SANDALS_GAUSSLEGENDRE6_HXX
