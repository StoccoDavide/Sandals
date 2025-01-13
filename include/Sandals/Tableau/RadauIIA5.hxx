/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_RADAUIIA5_HXX
#define SANDALS_RADAUIIA5_HXX

namespace Sandals
{
  //! \brief Butcher tableau for the Radau IIA order 3 method.
  //!
  //! Butcher tableau for the Radau IIA order 5 method:
  //!
  //! \f[
  //! \begin{array}{c|ccc}
  //!   \frac{2}{5} - \frac{\sqrt{6}}{10} &
  //!   \frac{11}{45} - \frac{7 \sqrt{6}}{360} &
  //!   \frac{37}{225} - \frac{169 \sqrt{6}}{1800}  &
  //!   -\frac{2}{225} + \frac{\sqrt{6}}{75}  \\
  //!   \frac{2}{5} + \frac{\sqrt{6}}{10} &
  //!   \frac{37}{225} + \frac{169 \sqrt{6}}{1800} &
  //!   \frac{11}{45} + \frac{7 \sqrt{6}}{360} &
  //!   -\frac{2}{225} - \frac{\sqrt{6}}{75} \\
  //!   1 &
  //!   \frac{4}{9} - \frac{\sqrt{6}}{36} &
  //!   \frac{4}{9} + \frac{\sqrt{6}}{36} &
  //!   \frac{1}{9} \\
  //!   \hline
  //!     & \frac{4}{9} - \frac{\sqrt{6}}{36} &
  //!       \frac{4}{9} + \frac{\sqrt{6}}{36} &
  //!       \frac{1}{9} \\
  //! \end{array} \text{.}
  //! \f]
  class RadauIIA5Tableau : public Tableau<3>
  {
  public:
    using Tableau<3>::Type;
    using Tableau<3>::Vector;
    using Tableau<3>::Matrix;

    //! Class constructor for the Radau IIA order 5 method.
    RadauIIA5Tableau() {
      this->name  = "RadauIIA5";
      this->type  = Type::IRK;
      this->order = 5;
      Real s6{std::sqrt(6.0)};
      this->A << 11.0/45.0-7.0*s6/360.0,   37.0/225.0-169.0*s6/1800.0, -2.0/225.0+s6/75.0,
                 37.0/225+169.0*s6/1800.0, 11.0/45.0+7.0*s6/360.0,     -2.0/225.0-s6/75.0,
                 4.0/9.0-s6/36.0,          4.0/9.0+s6/36.0,            1.0/9.0;
      this->b << 4.0/9.0-s6/36.0, 4.0/9.0+s6/36, 1.0/9.0;
      this->c << 2.0/5.0-s6/10.0, 2.0/5.0+s6/10, 1.0;
    }
  }; // class RadauIIA5Tableau

  //! Class container for the Radau IIA order 5 method.
  template <Size N, Size M = 0>
  class RadauIIA5 : public RungeKutta<3, N, M>
  {
  public:
    using System = typename Implicit<N, M>::Pointer; //!< System type.

    //! Class constructor for a Radau IIA order 5 solver given a Tableau reference.
    RadauIIA5() : RungeKutta<3, N, M>(RadauIIA5Tableau()) {}

    //! Class constructor for a Radau IIA order 5 solver given a Tableau reference.
    //! \param[in] t_system The system reference.
    RadauIIA5(System t_system) : RungeKutta<3, N, M>(RadauIIA5Tableau(), t_system) {}

  }; // class RadauIIA5

} // namespace Sandals

#endif // SANDALS_RADAUIIA5_HXX
