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

#ifndef SANDALS_LINEAR_SYSTEM_HXX
#define SANDALS_LINEAR_SYSTEM_HXX

namespace Sandals {

  /*\
   |   _     _
   |  | |   (_)_ __   ___  __ _ _ __
   |  | |   | | '_ \ / _ \/ _` | '__|
   |  | |___| | | | |  __/ (_| | |
   |  |_____|_|_| |_|\___|\__,_|_|
   |
  \*/

  /**
  * \brief Class container for the system of linear ODEs/DAEs.
  *
  * Class container for the system of linear ordinary differential equations (ODEs) or differential
  * algebraic equations (DAEs) of the type \f$ \mathbf{E}(t)\mathbf{x}^{\prime} = \mathbf{A}(t)
  * \mathbf{x} + \mathbf{b}(t) \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam N The dimension of the semi-explicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <Integer N, Integer M = 0>
  class Linear : public Explicit<N, M>
  {
  public:
    using Pointer  = std::shared_ptr<SemiExplicit<N, M>>; /**< Shared pointer to a semi-explicit ODE/DAE system. */
    using VectorF  = typename Explicit<N, M>::VectorF;    /**< Templetized vector type. */
    using MatrixJF = typename Explicit<N, M>::MatrixJF;   /**< Templetized matrix type. */
    using MatrixE  = typename Explicit<N, M>::MatrixJF;   /**< Templetized matrix type. */
    using MatrixA  = typename Explicit<N, M>::MatrixJF;   /**< Templetized matrix type. */
    using VectorB  = typename Explicit<N, M>::VectorF;    /**< Templetized vector type. */
    using Type     = typename Explicit<N, M>::Type;       /**< System type enumeration. */

  private:
    Eigen::FullPivLU<MatrixE> m_lu; /**< LU decomposition. */

  public:
    /**
    * Class constructor for the linear ODE/DAE system.
    */
    Linear() : Explicit<N, M>(Type::LINEAR, "(missing name)") {}

    /**
    * Class constructor for the linear ODE/DAE system.
    * \param[in] t_name The name of the linear ODE/DAE system.
    */
    Linear(std::string t_name) : Explicit<N, M>(Type::LINEAR, t_name) {}

    /**
    * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    *
    * \f[
    * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{E}(t) \mathbf{x}^{\prime} -
    * \mathbf{A}(t) \mathbf{x} - \mathbf{b}(t)
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    VectorF F(VectorF const &x, VectorF const &x_dot, Real t) const override
    {
      return this->E(t)*x_dot - this->A(t)*x - this->b(t);
    }

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    * t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    * -\mathbf{A}(t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real t) const override {return -this->A(t);}

    /**
    * Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    * \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    * \mathbf{E}(t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x_dot(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real t) const override {return this->E(t);}

    /**
    * Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$ as
    *
    * \f[
    * \mathbf{f}(\mathbf{x}, t) = \mathbf{E}(t)^{-1}(\mathbf{A}(t)\mathbf{x} + \mathbf{b}(t))
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    */
    VectorF f(VectorF const &x, Real t) const override
    {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Linear::f(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t)*x + this->b(t));
    }

    /**
    * Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = \displaystyle\frac{\partial\mathbf{f}(
    * \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{E}(t)^{-1} \mathbf{A}(t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    MatrixJF Jf_x(VectorF const &/*x*/, Real t) const override
    {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Linear::Jf_x(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t));
    }

    /**
    * Evaluate the linear ODE/DAE system mass matrix \f$ \mathbf{E}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{E}(t) \f$.
    */
    virtual MatrixE E(Real t) const = 0;

    /**
    * Evaluate the linear ODE/DAE system matrix \f$ \mathbf{A}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{A}(t) \f$.
    */
    virtual MatrixA A(Real t) const = 0;

    /**
    * Evaluate the linear ODE/DAE system vector \f$ \mathbf{b}(t) \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{b}(t) \f$.
    */
    virtual VectorB b(Real t) const = 0;

  }; // class Linear

} // namespace Sandals

#endif // SANDALS_LINEAR_SYSTEM_HXX
