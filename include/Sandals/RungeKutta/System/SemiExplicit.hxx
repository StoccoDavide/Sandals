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

#ifndef SANDALS_SEMI_EXPLICIT_SYSTEM_HXX
#define SANDALS_SEMI_EXPLICIT_SYSTEM_HXX

namespace Sandals {

  //! \brief Class container for the system of explicit ODEs.
  //!
  //! Class container for the system of semi-explicit ordinary differential equations (ODEs) of the
  //! type \f$ \mathbf{M}(\mathbf{x}, t)\mathbf{x}^{\prime} = \mathbf{r}(\mathbf{x}, t) \f$.
  template <Size N>
  class SemiExplicit : public Implicit<N>
  {
  public:
    using Pointer = std::shared_ptr<SemiExplicit<N>>; //!< Shared pointer to a semi-explicit ODE system.
    using Vector  = typename Implicit<N>::Vector;     //!< Templetized vector type.
    using Matrix  = typename Implicit<N>::Matrix;     //!< Templetized matrix type.
    using Tensor  = typename std::vector<Matrix>;     //!< Templetized matrix type.
    using Type    = typename Implicit<N>::Type;       //!< ODE system type enumeration.

  private:
    Eigen::FullPivLU<Matrix> m_lu; //!< LU decomposition.

  public:
    //! Default class constructor for the semi-explicit ODE system.
    //!
    SemiExplicit() : Implicit<N>(Type::SEMI_EXPLICIT, "(missing name)") {}

    //! Class constructor for the semi-explicit ODE system.
    //! \param[in] t_name The name of the semi-explicit ODE system.
    SemiExplicit(std::string t_name) : Implicit<N>(Type::SEMI_EXPLICIT, t_name) {}

    //! Evaluate the ODE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    //!
    //! \f[ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{x}^{\prime} - \mathbf{f}(\mathbf{x}, t) \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    Vector F(Vector const &x, Vector const &x_dot, Real t) const override
    {
      return this->M(x, t)*x_dot - this->r(x, t);
    }

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
    //! t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //! \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    //! -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    Matrix JF_x(Vector const &x, Vector const &/*x_dot*/, Real t) const override
    {
      return -this->Jf_x(x, t);
    }

    //! Evaluate the Jacobian of the ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    //! \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    //!
    //! \f[ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    //! \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    //! \mathbf{M}(\mathbf{x}, t) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    Matrix JF_x_dot(Vector const &x, Vector const &/*x_dot*/, Real t) const override
    {
      return this->M(x, t);
    }

    //! Evaluate the explicit ODE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$ as
    //!
    //! \f[ \mathbf{f}(\mathbf{x}, t) = \mathbf{M}(\mathbf{x}, t)^{-1} \mathbf{r}(\mathbf{x}, t) \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    Vector f(Vector const &x, Real t) const override
    {
      this->m_lu.compute(this->M(x, t));
      return this->m_lu.solve(this->r(x, t));
    }

    //! Evaluate the Jacobian of the explicit ODE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    //! with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = \displaystyle\frac{\partial\mathbf{f}(
    //! \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{M}(\mathbf{x}, t)^{-1}( \mathbf{Jr}_{\mathbf{x}}
    //! (\mathbf{x}, t) - \mathbf{TM}_{\mathbf{x}}(\mathbf{x}, t) \mathbf{f}(\mathbf{x}, t)) =
    //! \mathbf{M}(\mathbf{x}, t)^{-1}(\mathbf{Jr}_{\mathbf{x}}(\mathbf{x}, t) - \mathbf{TM}_{\mathbf{x}}
    //! (\mathbf{x}, t) \mathbf{x}^\prime) \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    Matrix Jf_x(Vector const &x, Vector const &x_dot, Real t) const override
    {
      Tensor TM_x{this->TM_x(x, t)};
      Matrix tmp{Matrix::Zero()};
      for (Size i = 0; i < N; ++i) {tmp.col(i) = TM_x[i] * x_dot;} // TODO: Check if TM_x.size() = N
      this->m_lu.compute(this->M(x, t));
      return this->m_lu.solve(this->Jr_x(x, t) - tmp);
    }

    //! Evaluate the semi-explicit ODE system mass matrix \f$ \mathbf{M}(\mathbf{x}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{M}(\mathbf{x}, t) \f$.
    virtual Vector M(Vector const &x, Real t) const = 0;

    //! Evaluate the tensor of the semi-explicit ODE system mass matrix \f$ \mathbf{M}(\mathbf{x},
    //! t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{TM}_{\mathbf{x}}(\mathbf{x}, t) =
    //! \displaystyle\frac{\partial\mathbf{M}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{TM}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    virtual Tensor TM_x(Vector const &x, Real t) const = 0;

    //! Evaluate the semi-explicit ODE system right-hand-side function \f$ \mathbf{r}(\mathbf{x}, t) \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The system function \f$ \mathbf{r}(\mathbf{x}, t) \f$.
    virtual Vector r(Vector const &x, Real t) const = 0;

    //! Evaluate the Jacobian of the semi-explicit ODE system right-hand-side function \f$ \mathbf{r}
    //! (\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    //!
    //! \f[ \mathbf{Jr}_{\mathbf{x}}(\mathbf{x}, t) =
    //! \displaystyle\frac{\partial\mathbf{r}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.} \f]
    //!
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Independent variable \f$ t \f$.
    //! \return The Jacobian \f$ \mathbf{Jr}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    virtual Matrix Jr_x(Vector const &x, Real t) const = 0;

  }; // class SemiExplicit

} // namespace Sandals

#endif // SANDALS_SEMI_EXPLICIT_SYSTEM_HXX
