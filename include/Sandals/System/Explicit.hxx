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

#ifndef SANDALS_EXPLICIT_SYSTEM_HXX
#define SANDALS_EXPLICIT_SYSTEM_HXX

namespace Sandals {

  /*\
   |   _____            _ _      _ _
   |  | ____|_  ___ __ | (_) ___(_) |_
   |  |  _| \ \/ / '_ \| | |/ __| | __|
   |  | |___ >  <| |_) | | | (__| | |_
   |  |_____/_/\_\ .__/|_|_|\___|_|\__|
   |             |_|
  \*/

  /**
  * \brief Class container for the system of explicit ODEs.
  *
  * Class container for the system of explicit ordinary differential equations (ODEs) of the type
  * \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, t) \f$, with invariants manifold \f$ \mathbf{h}
  * (\mathbf{x}, t) = \mathbf{0} \f$.
  * \tparam N The dimension of the explicit ODE system.
  * \tparam M The dimension of the invariants manifold.
  */
  template <Integer N, Integer M>
  class Explicit : public Implicit<N, M>
  {
  public:
    using Pointer  = std::shared_ptr<Explicit<N, M>>;   /**< Shared pointer to an explicit ODE system. */
    using VectorF  = typename Implicit<N, M>::VectorF;  /**< Templetized vector type. */
    using MatrixJF = typename Implicit<N, M>::MatrixJF; /**< Templetized matrix type. */
    using Type     = typename Implicit<N, M>::Type;     /**< System type enumeration. */

  protected:

    /**
    * Class constructor for an explicit ODE/DAE system.
    * \param[in] t_type The type of the explicit ODE/DAE system.
    * \param[in] t_name The name of the explicit ODE/DAE system.
    */
    Explicit(Type t_type, std::string t_name) : Implicit<N, M>(t_type, t_name) {}

  public:

    /**
    * Class constructor for the explicit ODE system.
    */
    Explicit() : Implicit<N, M>(Type::EXPLICIT, "(missing name)") {}

    /**
    * Class constructor for the explicit ODE system.
    * \param[in] t_name The name of the explicit ODE system.
    */
    Explicit(std::string t_name) : Implicit<N, M>(Type::EXPLICIT, t_name) {}

    /**
    * Evaluate the ODE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$
    *
    * \f[
    * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    * \mathbf{x}^{\prime} - \mathbf{f}(\mathbf{x}, t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    VectorF F(VectorF const &x, VectorF const &x_dot, Real t) const override
    {
      return x_dot - this->f(x, t);
    }

    /**
    * Evaluate the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    * -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}} =
    * -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x(VectorF const &x, VectorF const &/*x_dot*/, Real t) const override
    {
      return -this->Jf_x(x, t);
    }

    /**
    * Evaluate the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
    *
    * \f[
    * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \displaystyle
    * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)}{\partial\mathbf{x}^{\prime}} =
    * \mathbf{I} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$.
    */
    MatrixJF JF_x_dot(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const override
    {
      return MatrixJF::Identity();
    }

    /**
    * Evaluate the explicit ODE system function\f$ \mathbf{f}(\mathbf{x}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
    */
    virtual VectorF f(VectorF const &x, Real t) const = 0;

    /**
    * Evaluate the Jacobian of the explicit ODE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} \f$
    *
    * \f[
    * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
    * \displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} \text{.}
    * \f]
    *
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
    */
    virtual MatrixJF Jf_x(VectorF const &x, Real t) const = 0;

    /**
    * Time reversal of the explicit ODE system function \f$ \mathbf{f}(\mathbf{x}, t) = -\mathbf{f}(
    * \mathbf{x}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed system function \f$ \mathbf{f}(\mathbf{x}, -t) \f$.
    */
    VectorF f_reverse(VectorF const &x, Real t) const
    {
      return -this->f(x, -t);
    }

    /**
    * Time reversal of the Jacobian of the explicit ODE system function \f$ \mathbf{f}(\mathbf{x}, t) \f$
    * with respect to the states \f$ \mathbf{x} = -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, -t) \f$.
    */
    MatrixJF Jf_x_reverse(VectorF const &x, Real t) const
    {
      return -this->Jf_x(x, -t);
    }

    /**
    * Time reversal of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{
    * \prime}, t) = -\mathbf{F}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed system function \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, -t) \f$.
    */
    VectorF F_reverse(VectorF const &x, VectorF const &x_dot, Real t) const
    {
      return -x_dot - this->f(x, -t);
    }

    /**
    * Time reversal of the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} = -\mathbf{JF}_{\mathbf{x}}
    * (\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, -\mathbf{x}^{\prime},
    * -t) \f$.
    */
    MatrixJF JF_x_reverse(VectorF const &x, VectorF const &/*x_dot*/, Real t) const
    {
      return -this->Jf_x(x, -t);
    }

    /**
    * Time reversal of the Jacobian of the implicit ODE system function \f$ \mathbf{F}(\mathbf{x},
    * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} =
    * -\mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The time-reversed Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{
    * \prime}, -t) \f$.
    */
    MatrixJF JF_x_dot_reverse(VectorF const &/*x*/, VectorF const &/*x_dot*/, Real /*t*/) const
    {
      return -MatrixJF::Identity();
    }

  }; // class Explicit

} // namespace Sandals

#endif // SANDALS_EXPLICIT_SYSTEM_HXX
