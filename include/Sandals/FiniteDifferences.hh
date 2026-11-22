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

#ifndef SANDALS_FINITE_DIFFERENCES_HH
#define SANDALS_FINITE_DIFFERENCES_HH

#include <Sandals.hh>

namespace Sandals
{

  /**
  * Compute the Jacobian of a vector function \f$ \mathbf{g}(\mathbf{z}) \f$ by central finite
  * differences, with the step scaled on the magnitude of each component of \f$ \mathbf{z} \f$.
  * \tparam Real The scalar number type.
  * \tparam R The number of rows of the Jacobian (i.e., the size of \f$ \mathbf{g} \f$), can be
  * Eigen::Dynamic.
  * \tparam C The number of columns of the Jacobian (i.e., the size of \f$ \mathbf{z} \f$), can be
  * Eigen::Dynamic.
  * \param[in] g The vector function \f$ \mathbf{g}(\mathbf{z}) \f$.
  * \param[in] z The evaluation point \f$ \mathbf{z} \f$.
  * \param[in] rows The number of rows, needed only if \f$ R \f$ is Eigen::Dynamic.
  * \return The Jacobian matrix, whose \f$ j \f$-th column is \f$ \partial\mathbf{g}/\partial z_j \f$.
  */
  template <typename Real, Integer R, Integer C, typename Function>
  Eigen::Matrix<Real, R, C> finite_differences_jacobian(Function && g, Eigen::Matrix<Real, C, 1> const & z,
    Integer const rows = R)
  {
    static Real const eps{std::cbrt(std::numeric_limits<Real>::epsilon())};
    Eigen::Matrix<Real, R, C> jac(rows, z.size());
    Eigen::Matrix<Real, C, 1> z_p(z), z_m(z);
    for (Integer j{0}; j < static_cast<Integer>(z.size()); ++j) {
      Real const delta{eps*std::max(Real(1.0), std::abs(z(j)))};
      z_p(j) = z(j) + delta;
      z_m(j) = z(j) - delta;
      jac.col(j) = (g(z_p) - g(z_m)) / (z_p(j) - z_m(j));
      z_p(j) = z(j);
      z_m(j) = z(j);
    }
    return jac;
  }

  /**
  * Compute the gradient of a scalar function \f$ g(\mathbf{z}) \f$ by central finite differences.
  * \tparam Real The scalar number type.
  * \tparam C The size of \f$ \mathbf{z} \f$, can be Eigen::Dynamic.
  * \param[in] g The scalar function \f$ g(\mathbf{z}) \f$.
  * \param[in] z The evaluation point \f$ \mathbf{z} \f$.
  * \return The gradient vector \f$ \nabla g(\mathbf{z}) \f$.
  */
  template <typename Real, Integer C, typename Function>
  Eigen::Matrix<Real, C, 1> finite_differences_gradient(Function && g, Eigen::Matrix<Real, C, 1> const & z)
  {
    return finite_differences_jacobian<Real, 1, C>([&g] (Eigen::Matrix<Real, C, 1> const & w)
      {return Eigen::Matrix<Real, 1, 1>(g(w));}, z).transpose();
  }

  /**
  * Compute the symmetric Hessian of a scalar function from its gradient \f$ \nabla g(\mathbf{z}) \f$
  * by central finite differences.
  * \tparam Real The scalar number type.
  * \tparam C The size of \f$ \mathbf{z} \f$, can be Eigen::Dynamic.
  * \param[in] grad The gradient function \f$ \nabla g(\mathbf{z}) \f$.
  * \param[in] z The evaluation point \f$ \mathbf{z} \f$.
  * \return The (symmetrized) Hessian matrix \f$ \nabla^2 g(\mathbf{z}) \f$.
  */
  template <typename Real, Integer C, typename Function>
  Eigen::Matrix<Real, C, C> finite_differences_hessian(Function && grad, Eigen::Matrix<Real, C, 1> const & z)
  {
    Eigen::Matrix<Real, C, C> hess(finite_differences_jacobian<Real, C, C>(grad, z, z.size()));
    return Real(0.5)*(hess + hess.transpose());
  }

} // namespace Sandals

#endif // SANDALS_FINITE_DIFFERENCES_HH
