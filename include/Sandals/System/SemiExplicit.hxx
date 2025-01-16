/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_SEMIEXPLICIT_SYSTEM_HXX
#define SANDALS_SEMIEXPLICIT_SYSTEM_HXX

namespace Sandals {

  /*\
   |   ____                 _ _____            _ _      _ _
   |  / ___|  ___ _ __ ___ (_) ____|_  ___ __ | (_) ___(_) |_
   |  \___ \ / _ \ '_ ` _ \| |  _| \ \/ / '_ \| | |/ __| | __|
   |   ___) |  __/ | | | | | | |___ >  <| |_) | | | (__| | |_
   |  |____/ \___|_| |_| |_|_|_____/_/\_\ .__/|_|_|\___|_|\__|
   |                                    |_|
  \*/

  //! \brief Class container for the system of explicit ODEs/DAEs.
  //!
  //! Class container for the system of semi-explicit ordinary differential equations (ODEs) or
  //! differential algebraic equations (DAEs)of the type \f$ \mathbf{A}(\mathbf{x}, t)
  //! \mathbf{x}^{\prime} = \mathbf{b}(\mathbf{x}, t) \f$, with invariants manifold \f$ \mathbf{h}(
  //! \mathbf{x}, t) = \mathbf{0} \f$.
  //! \tparam N The dimension of the semi-explicit ODE system.
  //! \tparam M The dimension of the invariants manifold.
  template <Size N, Size M = 0>
  class SemiExplicit : public Explicit<N, M>
  {
  public:
    using Pointer  = std::shared_ptr<SemiExplicit<N, M>>; //!< Shared pointer to a semi-explicit ODE/DAE system.
    using VectorF  = typename Explicit<N, M>::VectorF;    //!< Templetized vector type.
    using MatrixJF = typename Explicit<N, M>::MatrixJF;   //!< Templetized matrix type.
    using MatrixA  = typename Explicit<N, M>::MatrixJF;   //!< Templetized matrix type.
    using TensorTA = typename std::vector<MatrixJF>;      //!< Templetized matrix type.
    using VectorB  = typename Explicit<N, M>::VectorF;    //!< Templetized vector type.
    using MatrixJB = typename Explicit<N, M>::MatrixJF;   //!< Templetized vector type.
    using Type     = typename Explicit<N, M>::Type;       //!< System type enumeration.

  private:
    mutable Eigen::FullPivLU<MatrixA> m_lu; //!< LU decomposition.

  public:
    //! Default class constructor for the semi-explicit ODE/DAE system.
    SemiExplicit() : Explicit<N, M>(Type::SEMIEXPLICIT, "(missing name)") {}

    //! Class constructor for the semi-explicit ODE/DAE system.
    //! \param[in] t_name The name of the semi-explicit ODE/DAE system.
    SemiExplicit(std::string t_name) : Explicit<N, M>(Type::SEMIEXPLICIT, t_name) {}

    //! Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$
    //!
    //! \f[
    //! \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //! \mathbf{A}(\mathbf{x}, t)\mathbf{x}^{\prime} - \mathbf{b}(\mathbf{x}, t) \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    VectorF F(VectorF const &x, VectorF const &x_dot, Real t) const override
    {
      return this->A(x, t)*x_dot - this->b(x, t);
    }

    //! Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{
    //! \prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //! \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = -\mathbf{A}(\mathbf{x}, t)^{-1}(\mathbf{Jb}_{
    //! \mathbf{x}}(\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{x}^\prime) \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    MatrixJF JF_x(VectorF const &x, VectorF const &x_dot, Real t) const override
    {
      return -this->Jf_x(x, x_dot, t);
    }

    //! Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    //! \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    //!
    //! \f[
    //! \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    //! \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    //! \mathbf{A}(\mathbf{x}, t) \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    MatrixJF JF_x_dot(VectorF const &x, VectorF const &/*x_dot*/, Real t) const override
    {
      return this->A(x, t);
    }

    //! Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$ as
    //!
    //! \f[
    //! \mathbf{f}(\mathbf{x}, t) = \mathbf{A}(\mathbf{x}, t)^{-1} \mathbf{b}(\mathbf{x}, t) \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    VectorF f(VectorF const &x, Real t) const override
    {
      this->m_lu.compute(this->A(x, t));
      SANDALS_ASSERT(this->m_lu.rank() == N,
        "Sandals:SemiExplicit::f(...): singular mass matrix A detected.");
      return this->m_lu.solve(this->b(x, t));
    }

    //! Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    //! with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^\prime, t) = \displaystyle\frac{\partial
    //! \mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{A}(\mathbf{x}, t)^{-1}(\mathbf{Jb}_{
    //! \mathbf{x}}(\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{x}^\prime) \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^\prime, t) \f$.
    MatrixJF Jf_x(VectorF const &x, VectorF const &x_dot, Real t) const
    {
      TensorTA TA_x{this->TA_x(x, t)};
      MatrixJF tAp{MatrixJF::Zero()};
      for (Size i{0}; i < N; ++i) {tAp.col(i) = TA_x[i] * x_dot;} // TODO: Check if TA_x.size() = N
      this->m_lu.compute(this->A(x, t));
      SANDALS_ASSERT(this->m_lu.rank() == N,
        "Sandals:SemiExplicit::Jf_x(...): singular mass matrix A detected.");
      return this->m_lu.solve(this->Jb_x(x, t) - tAp);
    }

    //! Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    //! with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = \displaystyle\frac{\partial\mathbf{f}(
    //! \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{A}(\mathbf{x}, t)^{-1}( \mathbf{Jb}_{\mathbf{x}}
    //! (\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{f}(\mathbf{x}, t)) \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    MatrixJF Jf_x(VectorF const &x, Real t) const override {return this->Jf_x(x, this->f(x, t), t);}

    //! Evaluate the semi-explicit ODE/DAE system mass matrix \f$ \mathbf{A}(\mathbf{x}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The system function \f$ \mathbf{A}(\mathbf{x}, t) \f$.
    virtual MatrixA A(VectorF const &x, Real t) const = 0;

    //! Evaluate the tensor of the semi-explicit ODE/DAE system mass matrix \f$ \mathbf{A}(\mathbf{x},
    //! t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) =
    //! \displaystyle\frac{\partial\mathbf{A}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    virtual TensorTA TA_x(VectorF const &x, Real t) const = 0;

    //! Evaluate the semi-explicit ODE/DAE system right-hand-side function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The system function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
    virtual VectorB b(VectorF const &x, Real t) const = 0;

    //! Evaluate the Jacobian of the semi-explicit ODE/DAE system right-hand-side function \f$ \mathbf{b}
    //! (\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[
    //! \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) =
    //! \displaystyle\frac{\partial\mathbf{b}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    //! \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable (or time) \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    virtual MatrixJB Jb_x(VectorF const &x, Real t) const = 0;

  }; // class SemiExplicit

} // namespace Sandals

#endif // SANDALS_SEMIEXPLICIT_SYSTEM_HXX
