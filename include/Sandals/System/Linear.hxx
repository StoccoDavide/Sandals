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

  //! \brief Class container for the system of linear ODEs/DAEs.
  //!
  //! Class container for the system of linear ordinary differential equations (ODEs) or differential
  //! algebraic equations (DAEs) of the type \f$ \mathbf{E}(t)\mathbf{x}^{\prime} = \mathbf{A}(t)
  //! \mathbf{x} + \mathbf{f}(t) \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
  template <Size N, Size M = 0>
  class Linear : public Explicit<N, M>
  {
  public:
    using Pointer = std::shared_ptr<Linear<N>>;   //!< Shared pointer to a linear ODE/DAE system.
    using VectorN = typename Explicit<N, M>::VectorN; //!< Templetized vector type.
    using MatrixN = typename Explicit<N, M>::MatrixN; //!< Templetized matrix type.
    using Tensor  = typename std::vector<MatrixN>; //!< Templetized matrix type.
    using Type    = typename Explicit<N, M>::Type;   //!< System type enumeration.

  private:
    Eigen::FullPivLU<MatrixN> m_lu; //!< LU decomposition.

  public:
    //! Default class constructor for the linear ODE/DAE system.
    //!
    Linear() : Explicit<N, M>(Type::LINEAR, "(missing name)") {}

    //! Class constructor for the linear ODE/DAE system.
    //! \param[in] t_name The name of the linear ODE/DAE system.
    Linear(std::string t_name) : Explicit<N, M>(Type::LINEAR, t_name) {}

    //! Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    //!
    //! \f[ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{E}(t) \mathbf{x}^{\prime} -
    //! \mathbf{A}(t) \mathbf{x} - \mathbf{f}(t) \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    VectorN F(VectorN const &x, VectorN const &x_dot, Real t) const override
    {
      return this->E(t)*x_dot - this->A(t)*x - this->f(t);
    }

    //! Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    //! t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //! \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\mathbf{A}(t) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    MatrixN JF_x(VectorN const &/*x*/, VectorN const &/*x_dot*/, Real t) const override {return -this->A(t);}

    //! Evaluate the Jacobian of the ODE/DAE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    //! \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    //! \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    //! \mathbf{E}(t) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    MatrixN JF_x_dot(VectorN const &/*x*/, VectorN const &/*x_dot*/, Real t) const override {return this->E(t);}

    //! Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$ as
    //!
    //! \f[ \mathbf{f}(\mathbf{x}, t) = \mathbf{E}(t)^{-1}(\mathbf{A}(t)\mathbf{x} + \mathbf{f}(t)) \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    VectorN f(VectorN const &x, Real t) const override
    {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Linear::f(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t)*x + this->f(t));
    }

    //! Evaluate the Jacobian of the explicit ODE/DAE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    //! with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = \displaystyle\frac{\partial\mathbf{f}(
    //! \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{E}(t)^{-1} \mathbf{A}(t) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    MatrixN Jf_x(VectorN const &/*x*/, Real t) const override
    {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(this->m_lu.rank() == N, "Sandals:Linear::Jf_x(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t));
    }

    //! Evaluate the linear ODE/DAE system mass matrix \f$ \mathbf{E}(t) \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{E}(t) \f$.
    virtual MatrixN E(Real t) const = 0;

    //! Evaluate the linear ODE/DAE system matrix \f$ \mathbf{A}(t) \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{A}(t) \f$.
    virtual MatrixN A(Real t) const = 0;

    //! Evaluate the linear ODE/DAE system vector \f$ \mathbf{f}(t) \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{f}(t) \f$.
    virtual VectorN f(Real t) const = 0;

  }; // class Linear

} // namespace Sandals

#endif // SANDALS_LINEAR_SYSTEM_HXX
