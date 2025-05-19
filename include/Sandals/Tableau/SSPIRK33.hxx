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

#ifndef SANDALS_SSPIRK33_HXX
#define SANDALS_SSPIRK33_HXX

namespace Sandals
{
  /**
  * \brief Butcher tableau for the 3-stage diagonally-implicit strong-stability preserving Runge-Kuttax order 3 method.
  *
  * Butcher tableau for the 3-stage diagonally-implicit strong-stability preserving Runge-Kutta
  * order 3 method:
  *
  * \f[
  * \begin{array}{c|ccc}
  *   \frac{1}{2}-\frac{\sqrt{2}}{4} & \frac{1}{2}-\frac{\sqrt{2}}{4} &                              0 &                              0 \\
  *                      \frac{1}{2} &             \frac{\sqrt{2}}{4} & \frac{1}{2}-\frac{\sqrt{2}}{4} &                              0 \\
  *   \frac{1}{2}+\frac{\sqrt{2}}{4} &             \frac{\sqrt{2}}{4} &             \frac{\sqrt{2}}{4} & \frac{1}{2}-\frac{\sqrt{2}}{4} \\
  *   \hline
  *                                  &                    \frac{1}{3} &        \frac{1}{3} & \frac{1}{3}
  *   \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class SSPIRK33Tableau : public Tableau<Real, 3>
  {
  public:
    using typename Tableau<Real, 3>::Type;
    using typename Tableau<Real, 3>::Vector;
    using typename Tableau<Real, 3>::Matrix;

    /**
    * Class constructor for the 3-stage diagonally-implicit strong-stability preserving Runge-Kutta
    * order 3 method.
    */
    SSPIRK33Tableau<Real>() {
      this->name  = "SSPIRK33";
      this->type  = Type::DIRK;
      this->order = 3;
      Real t_1{1.0/2.0};
      Real t_2{std::sqrt(2.0)/4.0};
      Real t_3{1.0/3.0};
      Real t_4{t_1-t_2};
      Real t_5{t_1+t_2};
      this->A << t_4, 0.0, 0.0,
                 t_2, t_4, 0.0,
                 t_2, t_2, t_4;
      this->b << t_3, t_3, t_3;
      this->c << t_4, t_1, t_5;
    }
  }; // class SSPIRK33Tableau

  /**
  * \brief Class container for the 3-stage diagonally-implicit strong-stability preserving Runge-Kutta order 3 method.
  *
  * Class container for the 3-stage diagonally-implicit strong-stability preserving Runge-Kutta order 3 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class SSPIRK33 : public RungeKutta<Real, 3, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a 3-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    */
    SSPIRK33() : RungeKutta<Real, 3, N, M>(SSPIRK33Tableau<Real>()) {}

    /**
    * Class constructor for a 3-stage strong-stability preserving Runge-Kutta order 3 solver given
    * a Tableau reference.
    * \param[in] t_system The system reference.
    */
    SSPIRK33(System t_system) : RungeKutta<Real, 3, N, M>(SSPIRK33Tableau<Real>(), t_system) {}

  }; // class SSPIRK33

} // namespace Sandals

#endif // SANDALS_SSPIRK33_HXX
