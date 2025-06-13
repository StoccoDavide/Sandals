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

#ifndef SANDALS_SOLUTION_HH
#define SANDALS_SOLUTION_HH

#include <Sandals.hh>

namespace Sandals {

  /*\
   |   ____        _       _   _
   |  / ___|  ___ | |_   _| |_(_) ___  _ __
   |  \___ \ / _ \| | | | | __| |/ _ \| '_ \
   |   ___) | (_) | | |_| | |_| | (_) | | | |
   |  |____/ \___/|_|\__,_|\__|_|\___/|_| |_|
   |
  \*/

  /**
  * \brief Class container for the numerical solution of a system of ODEs/DAEs.
  *
  * Class container for the numerical solution of a system of Ordinary Differential Equations (ODEs)
  * or Differential Algebraic Equations (DAEs). The solution is stored in a independent variable
  * (or time) vector, states matrix and invariant matrix, where each row of matrix represents the
  * value of a state/invariant at a given mesh node, *i.e.*,
  *
  * \f[
  * \begin{array}{l}
  *   \mathbf{t} &=& \left[ t_1, t_2, \ldots, t_k \right]^\top \text{,} \\
  *   \mathbf{X} &=& \begin{bmatrix} \mathbf{x}_1 \\ \mathbf{x}_2 \\ \vdots \\ \mathbf{x}_k \end{bmatrix}
  *   =  \begin{bmatrix} x_{11} & x_{12} & \cdots & x_{1k} \\
  *                      x_{21} & x_{22} & \cdots & x_{2k} \\
  *                      \vdots & \vdots & \ddots & \vdots \\
  *                      x_{n1} & x_{n2} & \cdots & x_{nk} \end{bmatrix} \text{,} \\
  *   \mathbf{H} &=& \begin{bmatrix} \mathbf{h}_1 \\ \mathbf{h}_2 \\ \vdots \\ \mathbf{h}_k \end{bmatrix}
  *   =  \begin{bmatrix} h_{11} & h_{12} & \cdots & h_{1k} \\
  *                      h_{21} & h_{22} & \cdots & h_{2k} \\
  *                      \vdots & \vdots & \ddots & \vdots \\
  *                      h_{m1} & h_{m2} & \cdots & h_{mk} \end{bmatrix} \text{.}
  * \end{array}
  * \f]
  *
  * Such a solution can be retrieved as standard library vectors or Eigen library vectors.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam M The dimension of the invariants manifold.
  * \note The solution is stored in a row-major format. This allows to extract contiguous
  * vectors from the states and invariants matrices, which is useful for I/O operations.
  */
  template <typename Real, Integer N, Integer M = 0>
  struct Solution
  {
    using Vector  = Eigen::Vector<Real, Eigen::Dynamic>; /**< Templetized vector type. */
    using MatrixN = Eigen::Matrix<Real, N, Eigen::Dynamic, Eigen::RowMajor>; /**< Templetized matrix type. */
    using MatrixM = Eigen::Matrix<Real, M, Eigen::Dynamic, Eigen::RowMajor>; /**< Templetized matrix type. */

    Vector  t; /**< Independent variable (or time) solution vector \f$ \mathbf{t} \f$. */
    MatrixN x; /**< States solution matrix \f$ \mathbf{X} \f$. */
    MatrixM h; /**< Invariants solution matrix \f$ \mathbf{H} \f$. */

    /**
    * Class constructor for the Solution class.
    */
    Solution() : t(0), x(Vector::Zero(N, 0)), h(Vector::Zero(M, 0)) {}

    /**
    * Class constructor for the Solution class given the independent variable (or time) and states sizes.
    * \param[in] size The size of the independent variable (or time) vector and states matrix.
    */
    Solution(Integer const size)
      : t(size), x(Vector::Zero(N, size)), h(Vector::Zero(M, size)) {}

    /**
    * Resize the independent variable (or time) vector and states matrix.
    * \param[in] size The size of the independent variable (or time) vector and states matrix.
    */
    void resize(Integer const size) {
      this->t.resize(size);
      this->x.resize(Eigen::NoChange, size);
      this->h.resize(Eigen::NoChange, size);
    }

    /**
    * Resize conserving the data of the independent variable (or time) vector and states matrix.
    * \param[in] size The size of the independent variable (or time) vector and states matrix.
    */
    void conservative_resize(Integer const size) {
      this->t.conservativeResize(size);
      this->x.conservativeResize(Eigen::NoChange, size);
      this->h.conservativeResize(Eigen::NoChange, size);
    }

    /**
    * Clear the independent variable (or time) vector and states matrix.
    */
    void clear()
    {
      this->t.resize(0);
      this->x.resize(Eigen::NoChange, 0);
      this->h.resize(Eigen::NoChange, 0);
    }

    /**
    * Check if the independent variable (or time) vector and states matrix are empty.
    * \return True if the independent variable (or time) vector and states matrix are empty, false otherwise.
    */
    bool is_empty() const
    {
      return this->t.size() == 0 && this->x.cols()== 0 && this->h.cols() == 0;
    }

    /**
    * Get the size of the independent variable (or time) vector and states matrix.
    * \return The size of the independent variable (or time) vector and states matrix.
    */
    Integer size() const {return this->t.size();}

    /**
    * Get the independent variable (or time) vector \f$ \mathbf{t} \f$ as a standard library vector.
    * \return The independent variable (or time) vector \f$ \mathbf{t} \f$ as a standard library vector.
    */
    std::vector<Real> std_t() const
    {
      return std::vector<Real>(this->t.data(), this->t.data() + this->t.size());
    }

    /**
    * Get the independent variable (or time) vector \f$ \mathbf{t} \f$ as Eigen library vector.
    * \return The independent variable (or time) vector \f$ \mathbf{t} \f$ as Eigen library vector.
    */
    Vector eig_t() const {return this->t;}

    /**
    * Get the \f$ i \f$-th state vector \f$ \mathbf{x}_i \f$ as a standard library vector.
    * \param[in] i The index of the state vector.
    * \return The \f$ i \f$-th state vector as a standard library vector.
    */
    std::vector<Real> std_x(Integer i) const
    {
      Vector tmp(this->x.row(i));
      return std::vector<Real>(tmp.data(), tmp.data() + tmp.size());
    }

    /**
    * Get a map of the states vectors \f$ \mathbf{X} \f$ as standard library vectors.
    * \return A map of the states vectors as standard library vectors.
    */
    std::map<Integer, std::vector<Real>> std_x() const
    {
      std::map<Integer, std::vector<Real>> x_map;
      for (Integer i{0}; i < N; ++i) {x_map[i] = this->std_x(i);}
      return x_map;
    }

    /**
    * Get a map of the states vectors \f$ \mathbf{X} \f$ as Eigen library vectors.
    * \return A map of the states vectors as Eigen library vectors.
    */
    std::map<Integer, Vector> eig_x() const
    {
      std::map<Integer, Vector> x_map;
      for (Integer i{0}; i < N; ++i) {x_map[i] = this->x.row(i);}
      return x_map;
    }

    /**
    * Get a map of the states vectors \f$ \mathbf{X} \f$ as standard vectors.
    * \param[in] names The names of the states vectors.
    * \return A map of the states vectors as standard vectors.
    */
    std::map<std::string, std::vector<Real>> std_x(std::vector<std::string> const & names) const
    {
      std::map<std::string, std::vector<Real>> x_map;
      for (Integer i{0}; i < N; ++i) {x_map[names[i]] = this->std_x(i);}
      return x_map;
    }

    /**
    * Get a map of the states vectors \f$ \mathbf{X} \f$ as Eigen library vectors.
    * \param[in] names The names of the states vectors.
    * \return A map of the states vectors as Eigen library vectors.
    */
    std::map<std::string, Vector> eig_x(std::vector<std::string> const & names) const
    {
      std::map<std::string, Vector> x_map;
      for (Integer i{0}; i < N; ++i) {x_map[names[i]] = this->x.row(i);}
      return x_map;
    }

    /**
    * Get the \f$ i \f$-th invariants vector \f$ \mathbf{h}_i \f$ as a standard library vector.
    * \param[in] i The index of the invariant vector.
    * \return The \f$ i \f$-th invariants vector as a standard library vector.
    */
    std::vector<Real> std_h(Integer const i) const
    {
      Vector tmp(this->h.row(i));
      return std::vector<Real>(tmp.data(), tmp.data() + tmp.size());
    }

    /**
    * Get a map of the invariants vectors \f$ \mathbf{H} \f$ as standard library vectors.
    * \return A map of the invariants vectors as standard library vectors.
    */
    std::map<Integer, std::vector<Real>> std_h() const
    {
      std::map<Integer, std::vector<Real>> h_map;
      for (Integer i{0}; i < M; ++i) {h_map[i] = this->std_h(i);}
      return h_map;
    }

    /**
    * Get a map of the invariants vectors \f$ \mathbf{H} \f$ as Eigen library vectors.
    * \return A map of the invariants vectors as Eigen library vectors.
    */
    std::map<Integer, Vector> eig_h() const
    {
      std::map<Integer, Vector> h_map;
      for (Integer i{0}; i < M; ++i) {h_map[i] = this->h.row(i);}
      return h_map;
    }

    /**
    * Get a map of the invariants vectors \f$ \mathbf{H} \f$ as standard vectors.
    * \param[in] names The names of the invariants vectors.
    * \return A map of the invariants vectors as standard vectors.
    */
    std::map<std::string, std::vector<Real>> std_h(std::vector<std::string> const & names) const
    {
      std::map<std::string, std::vector<Real>> h_map;
      for (Integer i{0}; i < M; ++i) {h_map[names[i]] = this->std_h(i);}
      return h_map;
    }

    /**
    * Get a map of the invariants vectors \f$ \mathbf{H} \f$ as Eigen library vectors.
    * \param[in] names The names of the invariants vectors.
    * \return A map of the invariants vectors as Eigen library vectors.
    */
    std::map<std::string, Vector> eig_h(std::vector<std::string> const & names) const
    {
      std::map<std::string, Vector> h_map;
      for (Integer i{0}; i < M; ++i) {h_map[names[i]] = this->h.row(i);}
      return h_map;
    }
  };

} // namespace Sandals

#endif // SANDALS_SOLUTION_HH
