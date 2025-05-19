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

#ifndef SANDALS_LOBATTOIIIA2_HXX
#define SANDALS_LOBATTOIIIA2_HXX

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Lobatto IIIA order 2 method.
  *
  * Butcher tableau for the Lobatto IIIA order 2 method:
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
  class LobattoIIIA2Tableau : public Tableau<Real, 2>
  {
  public:
    using typename Tableau<Real, 2>::Type;
    using typename Tableau<Real, 2>::Vector;
    using typename Tableau<Real, 2>::Matrix;

    /**
    * Class constructor for the Lobatto IIIA order 2 method.
    */
    LobattoIIIA2Tableau<Real>() {
      this->name  = "LobattoIIIA2";
      this->type  = Type::DIRK;
      this->order = 2;
      this->A << 0.0,      0.0,
                 1.0/2.0,  1.0/2.0;
      this->b << 1.0/2.0, 1.0/2.0;
      this->c << 0.0,     1.0;
    }
  }; // class LobattoIIIA2Tableau

  /**
  * \brief Class container for the Lobatto IIIA order 2 method.
  *
  * Class container for the Lobatto IIIA order 2 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class LobattoIIIA2 : public RungeKutta<Real, 2, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Lobatto IIIA order 2 solver given a Tableau reference.
    */
    LobattoIIIA2() : RungeKutta<Real, 2, N, M>(LobattoIIIA2Tableau<Real>()) {}

    /**
    * Class constructor for a Lobatto IIIA order 2 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    LobattoIIIA2(System t_system) : RungeKutta<Real, 2, N, M>(LobattoIIIA2Tableau<Real>(), t_system) {}

  }; // class LobattoIIIA2

} // namespace Sandals

#endif // SANDALS_LOBATTOIIIA2_HXX
