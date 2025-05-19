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

#ifndef SANDALS_RALSTON4_HH
#define SANDALS_RALSTON4_HH

#include <Sandals.hh>
#include <Sandals/RungeKutta.hh>

namespace Sandals
{
  /**
  * \brief Butcher tableau for the Ralston's order 4 method.
  *
  * Butcher tableau for the Ralston's order 4 method:
  *
  * \f[
  * \begin{array}{c|ccccc}
  *   0           & 0   & 0 & 0 & 0 \\
  *   \frac{2}{5} & \frac{4}{10} & 0 & 0 & 0 \\
  *   \frac{7}{8}-\frac{3}{16}\sqrt{5} &
  *   \frac{357}{256}\sqrt{5}-\frac{2889}{1024} &
  *   \frac{3785}{1024}-\frac{405}{256}\sqrt{5} &
  *    0 &
  *    0 \\
  *   1 &
  *   \frac{1047}{3020}\sqrt{5}-\frac{673}{1208} &
  *   -\frac{975}{2552}-\frac{1523}{1276}\sqrt{5} &
  *   \frac{93408}{48169}+\frac{203968}{240845}\sqrt{5} &
  *   0 \\
  *   \hline
  *     & \frac{263}{1812}+\frac{2}{151}\sqrt{5} &
  *     \frac{125}{3828}-\frac{250}{957}\sqrt{5} &
  *     \frac{3426304}{5924787}+\frac{553984}{1974929}\sqrt{5} &
  *     \frac{10}{41}-\frac{4}{123}\sqrt{5}
  * \end{array} \text{.}
  * \f]
  * \tparam Real The scalar number type.
  */
  template <typename Real>
  class Ralston4Tableau : public Tableau<Real, 4>
  {
  public:
    using typename Tableau<Real, 4>::Type;
    using typename Tableau<Real, 4>::Vector;
    using typename Tableau<Real, 4>::Matrix;

    /**
    * Class constructor for the Ralston's order 4 method.
    */
    Ralston4Tableau<Real>() {
      this->name  = "Ralston4";
      this->type  = Type::ERK;
      this->order = 4;
      Real s_5{std::sqrt(5.0)};
      Real a_21{4.0/10.0};
      Real a_31{(357.0/256.0) * s_5 - 2889.0/1024.0};
      Real a_32{3785.0/1024.0 - (405.0/256.0) * s_5};
      Real a_41{(1047.0/3020.0) * s_5 - 673.0/1208.0};
      Real a_42{-975.0/2552.0 - (1523.0/1276.0) * s_5};
      Real a_43{93408.0/48169.0 + (203968.0/240845.0) * s_5};
      Real b_1{263.0/1812.0 + (2.0/151.0) * s_5};
      Real b_2{125.0/3828.0 - (250.0/957.0) * s_5};
      Real b_3{3426304.0/5924787.0 + (553984.0/1974929.0) * s_5};
      Real b_4{10.0/41.0 - (4.0/123.0) * s_5};
      Real c_2{2.0/5.0};
      Real c_3{7.0/8.0 - (3.0/16.0) * s_5};
      this->A << 0.0,  0.0,  0.0,  0.0,
                 a_21, 0.0,  0.0,  0.0,
                 a_31, a_32, 0.0,  0.0,
                 a_41, a_42, a_43, 0.0;
      this->b << b_1, b_2, b_3, b_4;
      this->c << 0.0, c_2, c_3, 1.0;
    }
  }; // class Ralston4Tableau

  /**
  * \brief Class container for the Ralston's order 4 method.
  *
  * Class container for the Ralston's order 4 method.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <typename Real, Integer N, Integer M = 0>
  class Ralston4 : public RungeKutta<Real, 4, N, M>
  {
  public:
    using System = typename Implicit<Real, N, M>::Pointer; /**< System type. */

    /**
    * Class constructor for a Ralston's order 4 solver given a Tableau reference.
    */
    Ralston4() : RungeKutta<Real, 4, N, M>(Ralston4Tableau<Real>()) {}

    /**
    * Class constructor for a Ralston's order 4 solver given a Tableau reference.
    * \param[in] t_system The system reference.
    */
    Ralston4(System t_system) : RungeKutta<Real, 4, N, M>(Ralston4Tableau<Real>(), t_system) {}

  }; // class Ralston4

} // namespace Sandals

#endif // SANDALS_RALSTON4_HH
