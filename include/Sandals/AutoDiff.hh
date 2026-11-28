/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_AUTODIFF_HH
#define SANDALS_AUTODIFF_HH

// This header is optional and requires CppAD (https://github.com/coin-or/CppAD).
#include <cppad/cppad.hpp>

#include <Sandals.hh>

// Eigen support for CppAD::AD scalars (the same definitions as <cppad/example/cppad_eigen.hpp>,
// which can be included instead, before Eigen, if its plugin is needed).
#ifndef CPPAD_EXAMPLE_CPPAD_EIGEN_HPP
namespace Eigen
{
  template <class Base>
  struct NumTraits<CppAD::AD<Base>>
  {
    using Real       = CppAD::AD<Base>;
    using NonInteger = CppAD::AD<Base>;
    using Literal    = CppAD::AD<Base>;
    using Nested     = CppAD::AD<Base>;
    enum {
      IsComplex = 0, IsInteger = 0, IsSigned = 1, RequireInitialization = 1,
      ReadCost = 1, AddCost = 2, MulCost = 2
    };
    static CppAD::AD<Base> epsilon() {return CppAD::numeric_limits<CppAD::AD<Base>>::epsilon();}
    static CppAD::AD<Base> dummy_precision() {return 100.0*CppAD::numeric_limits<CppAD::AD<Base>>::epsilon();}
    static CppAD::AD<Base> lowest() {return CppAD::numeric_limits<CppAD::AD<Base>>::min();}
    static CppAD::AD<Base> highest() {return CppAD::numeric_limits<CppAD::AD<Base>>::max();}
    static int digits10() {return CppAD::numeric_limits<CppAD::AD<Base>>::digits10;}
    static int max_digits10() {return CppAD::numeric_limits<CppAD::AD<Base>>::max_digits10;}
    static CppAD::AD<Base> quiet_NaN() {return CppAD::numeric_limits<CppAD::AD<Base>>::quiet_NaN();}
    static CppAD::AD<Base> infinity() {return CppAD::numeric_limits<CppAD::AD<Base>>::infinity();}
  };

  template <typename Base, typename BinaryOp>
  struct ScalarBinaryOpTraits<CppAD::AD<Base>, Base, BinaryOp> {using ReturnType = CppAD::AD<Base>;};

  template <typename Base, typename BinaryOp>
  struct ScalarBinaryOpTraits<Base, CppAD::AD<Base>, BinaryOp> {using ReturnType = CppAD::AD<Base>;};
} // namespace Eigen

namespace CppAD
{
  template <class Base> AD<Base> const & conj(AD<Base> const & x) {return x;}
  template <class Base> AD<Base> const & real(AD<Base> const & x) {return x;}
  template <class Base> AD<Base> imag(AD<Base> const &) {return AD<Base>(0.0);}
  template <class Base> AD<Base> abs2(AD<Base> const & x) {return x*x;}
} // namespace CppAD
#endif

namespace Sandals
{

  /*\
   |   _____
   |  |_   _|_ _ _ __   ___
   |    | |/ _` | '_ \ / _ \
   |    | | (_| | |_) |  __/
   |    |_|\__,_| .__/ \___|
   |            |_|
  \*/

  /**
  * \brief A function \f$ \mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m \f$ recorded by CppAD.
  *
  * The function is recorded once, by calling it with CppAD::AD<double> arguments, and then evaluated
  * with its exact Jacobian and multiplier-weighted Hessian. The recording freezes the operation
  * sequence: data-dependent branches must use CppAD::CondExp (and the constants used by the
  * function are fixed at recording time).
  */
  class Tape
  {
  public:
    using AD      = CppAD::AD<double>;                 /**< First-order AD scalar type. */
    using VectorX = Eigen::Vector<double, Eigen::Dynamic>; /**< Dynamic vector type. */
    using MatrixX = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; /**< Dynamic matrix type. */

  private:
    mutable CppAD::ADFun<double> m_function; /**< Recorded operation sequence. */
    bool m_recorded{false};                  /**< True if the function has been recorded. */

    // Sparse derivatives: the sparsity patterns are computed once (they depend on the operation
    // sequence only, conditional expressions included), then the Jacobian and the Hessian cost one
    // sweep per color of the pattern instead of one per input
    using SizeVector  = std::vector<std::size_t>;
    using Pattern     = CppAD::sparse_rc<SizeVector>;
    using SparseValues = CppAD::sparse_rcv<SizeVector, std::vector<double>>;
    mutable bool                    m_jacobian_ready{false}; /**< The Jacobian pattern is computed. */
    mutable Pattern                 m_jacobian_pattern;      /**< Jacobian sparsity pattern. */
    mutable SparseValues            m_jacobian_values;       /**< Jacobian nonzeros. */
    mutable CppAD::sparse_jac_work  m_jacobian_work;         /**< Jacobian coloring (reused). */
    mutable bool                    m_hessian_ready{false};  /**< The Hessian pattern is computed. */
    mutable Pattern                 m_hessian_pattern;       /**< Hessian sparsity pattern (all the outputs). */
    mutable SparseValues            m_hessian_values;        /**< Hessian nonzeros (lower triangle). */
    mutable CppAD::sparse_hes_work  m_hessian_work;          /**< Hessian coloring (reused). */

    static std::vector<double> to_std(VectorX const & x) {return std::vector<double>(x.data(), x.data() + x.size());}

  public:
    /**
    * Record a function.
    * \param[in] domain The dimension \f$ n \f$ of the domain.
    * \param[in] function The function, mapping a std::vector<AD> of size \f$ n \f$ to a
    * std::vector<AD>.
    * \param[in] point The point at which the function is recorded (any point in its domain).
    */
    template <typename Function>
    void record(Integer const domain, Function && function, double const point = 0.5)
    {
      std::vector<AD> x(domain, AD(point));
      CppAD::Independent(x);
      std::vector<AD> y{function(x)};
      this->m_function.Dependent(x, y);
      this->m_function.optimize();
      this->m_recorded = true;
      this->m_jacobian_ready = this->m_hessian_ready = false;
    }

    /** \return True if the function has been recorded. */
    bool recorded() const {return this->m_recorded;}

    /** \return The dimension of the domain. */
    Integer domain() const {return static_cast<Integer>(this->m_function.Domain());}

    /** \return The dimension of the range. */
    Integer range() const {return static_cast<Integer>(this->m_function.Range());}

    /** \param[in] x The point. \return The value \f$ \mathbf{f}(\mathbf{x}) \f$. */
    VectorX value(VectorX const & x) const
    {
      std::vector<double> const y{this->m_function.Forward(0, to_std(x))};
      return Eigen::Map<VectorX const>(y.data(), y.size());
    }

    /** \param[in] x The point. \return The Jacobian \f$ \partial\mathbf{f}/\partial\mathbf{x} \f$ (\f$ m \times n \f$). */
    /** \return The structural pattern of the Jacobian (\f$ m \times n \f$, true where an entry may be nonzero). */
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> jacobian_pattern() const
    {
      this->prepare_jacobian();
      Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> P(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(this->range(), this->domain(), false));
      for (std::size_t k{0}; k < this->m_jacobian_pattern.nnz(); ++k) {P(this->m_jacobian_pattern.row()[k], this->m_jacobian_pattern.col()[k]) = true;}
      return P;
    }

    /** \return The structural pattern of the Hessian of any weighted sum of the outputs (\f$ n \times n \f$, symmetric). */
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> hessian_pattern() const
    {
      this->prepare_hessian();
      Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> P(Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(this->domain(), this->domain(), false));
      for (std::size_t k{0}; k < this->m_hessian_pattern.nnz(); ++k) {
        P(this->m_hessian_pattern.row()[k], this->m_hessian_pattern.col()[k]) = true;
        P(this->m_hessian_pattern.col()[k], this->m_hessian_pattern.row()[k]) = true;
      }
      return P;
    }

    MatrixX jacobian(VectorX const & x) const
    {
      std::size_t const n{this->m_function.Domain()}, m{this->m_function.Range()};
      this->prepare_jacobian();
      std::vector<double> const xs{to_std(x)};
      this->m_function.sparse_jac_for(1, xs, this->m_jacobian_values, this->m_jacobian_pattern, "cppad", this->m_jacobian_work);
      MatrixX J(MatrixX::Zero(m, n));
      SizeVector const & row{this->m_jacobian_values.row()}, & col{this->m_jacobian_values.col()};
      std::vector<double> const & val{this->m_jacobian_values.val()};
      for (std::size_t k{0}; k < val.size(); ++k) {J(row[k], col[k]) = val[k];}
      return J;
    }

  private:
    void prepare_jacobian() const
    {
      std::size_t const n{this->m_function.Domain()};
      if (!this->m_jacobian_ready) {
        Pattern identity(n, n, n);
        for (std::size_t k{0}; k < n; ++k) {identity.set(k, k, k);}
        this->m_function.for_jac_sparsity(identity, false, false, false, this->m_jacobian_pattern);
        this->m_jacobian_values = SparseValues(this->m_jacobian_pattern);
        this->m_jacobian_work.clear();
        this->m_jacobian_ready = true;
      }
    }

    void prepare_hessian() const
    {
      std::size_t const n{this->m_function.Domain()}, m{this->m_function.Range()};
      if (!this->m_hessian_ready) {
        // Pattern for every output (the multipliers change), lower triangle only for the values. The
        // reverse-mode pattern is used: on tapes with conditional expressions the forward-mode one
        // (for_hes_sparsity) was observed to miss nonzero entries
        std::vector<bool> const select_range(m, true);
        Pattern identity(n, n, n), jacobian, full;
        for (std::size_t k{0}; k < n; ++k) {identity.set(k, k, k);}
        this->m_function.for_jac_sparsity(identity, false, false, true, jacobian); // stored for the reverse sweep
        this->m_function.rev_hes_sparsity(select_range, false, true, full);
        SizeVector const & row{full.row()}, & col{full.col()};
        std::size_t nnz{0};
        for (std::size_t k{0}; k < full.nnz(); ++k) {if (col[k] <= row[k]) {++nnz;}}
        Pattern lower(n, n, nnz);
        for (std::size_t k{0}, e{0}; k < full.nnz(); ++k) {if (col[k] <= row[k]) {lower.set(e++, row[k], col[k]);}}
        this->m_hessian_pattern = full;
        this->m_hessian_values = SparseValues(lower);
        this->m_hessian_work.clear();
        this->m_hessian_ready = true;
      }
    }

  public:

    /**
    * \param[in] x The point.
    * \param[in] w The multipliers (size \f$ m \f$).
    * \return The Hessian \f$ \partial^2(\mathbf{w}^\top\mathbf{f})/\partial\mathbf{x}^2 \f$ (\f$ n \times n \f$).
    */
    MatrixX hessian(VectorX const & x, VectorX const & w) const
    {
      std::size_t const n{this->m_function.Domain()};
      this->prepare_hessian();
      std::vector<double> const xs{to_std(x)}, ws{to_std(w)};
      this->m_function.sparse_hes(xs, ws, this->m_hessian_values, this->m_hessian_pattern, "cppad.symmetric", this->m_hessian_work);
      MatrixX H(MatrixX::Zero(n, n));
      SizeVector const & row{this->m_hessian_values.row()}, & col{this->m_hessian_values.col()};
      std::vector<double> const & val{this->m_hessian_values.val()};
      for (std::size_t k{0}; k < val.size(); ++k) {H(row[k], col[k]) = val[k]; H(col[k], row[k]) = val[k];}
      return H;
    }
  };

  /**
  * Split a std::vector of AD scalars into a fixed-size Eigen vector.
  * \tparam Size The size of the segment.
  * \param[in] v The vector.
  * \param[in] offset The index of the first element of the segment.
  * \return The segment as an Eigen vector.
  */
  template <Integer Size, typename Scalar>
  Eigen::Vector<Scalar, Size> segment(std::vector<Scalar> const & v, Integer const offset)
  {
    Eigen::Vector<Scalar, Size> s;
    for (Integer i{0}; i < Size; ++i) {s(i) = v[offset + i];}
    return s;
  }

  /**
  * Convert an Eigen vector of scalars to a std::vector.
  * \param[in] v The Eigen vector.
  * \return The std::vector.
  */
  template <typename Derived>
  std::vector<typename Derived::Scalar> to_std_vector(Eigen::MatrixBase<Derived> const & v)
  {
    std::vector<typename Derived::Scalar> s(v.size());
    for (Integer i{0}; i < v.size(); ++i) {s[i] = v(i);}
    return s;
  }

} // namespace Sandals

#endif // SANDALS_AUTODIFF_HH
