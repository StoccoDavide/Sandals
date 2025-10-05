/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SANDALS_LOBATTOIIIA4_HH
#define SANDALS_LOBATTOIIIA4_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Lobatto IIIA order 4 method.
  *
  * Butcher tableau for the Lobatto IIIA order 4 method:
  *
  * \f[
  * \begin{array}{c|ccc}
  *     0 &    0 &   0 &     0 \\
  *   1/2 & 5/24 & 1/3 & -1/24 \\
  *     1 &  1/6 & 2/3 &   1/6 \\
  *   \hline
  *       &  1/6 & 2/3 &   1/6 \\
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class LobattoIIIA4Tableau : public Tableau<Real, 3>
  {
  public:
    using typename Tableau<Real, 3>::Type;
    using typename Tableau<Real, 3>::Vector;
    using typename Tableau<Real, 3>::Matrix;

    /**
    * Class constructor for the Lobatto IIIA order 4 method.
    */
    LobattoIIIA4Tableau<Real>() {
      this->name  = "LobattoIIIA4";
      this->type  = Type::DIRK;
      this->order = 4;
      this->A << 0.0,      0.0,      0.0,
                 5.0/24.0, 1.0/3.0, -1.0/24.0,
                 1.0/6.0,  2.0/3.0,  1.0/6.0;
      this->b << 1.0/6.0, 2.0/3.0, 1.0/6.0;
      this->c << 0.0,     1.0/2.0,  1.0;
    }

  }; // class LobattoIIIA4Tableau

  /**
  * \brief Class container for the Lobatto IIIA order 4 method.
  *
  * Class container for the Lobatto IIIA order 4 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class LobattoIIIA4 : public RungeKutta<Real, 3, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Lobatto IIIA order 4 solver given a Tableau reference.
    */
    LobattoIIIA4() : RungeKutta<Real, 3, N, M>(LobattoIIIA4Tableau<Real>()) {}

    /**
    * Class constructor for a Lobatto IIIA order 4 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    LobattoIIIA4(System t_system) : RungeKutta<Real, 3, N, M>(LobattoIIIA4Tableau<Real>(), std::move(t_system)) {}

  }; // class LobattoIIIA4

} // namespace Sandals

#endif // SANDALS_LOBATTOIIIA4_HH
