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

#ifndef SANDALS_CONTROLLED_SYSTEM_HH
#define SANDALS_CONTROLLED_SYSTEM_HH

#include <Sandals.hh>
#include <Sandals/FiniteDifferences.hh>

namespace Sandals
{

  /*\
   |    ____            _             _ _          _
   |   / ___|___  _ __ | |_ _ __ ___ | | | ___  __| |
   |  | |   / _ \| '_ \| __| '__/ _ \| | |/ _ \/ _` |
   |  | |__| (_) | | | | |_| | | (_) | | |  __/ (_| |
   |   \____\___/|_| |_|\__|_|  \___/|_|_|\___|\__,_|
   |
  \*/

  /**
  * \brief Class container for a controlled system of implicit ODEs/DAEs.
  *
  * Class container for a controlled system of implicit ordinary differential equations (ODEs) or
  * differential algebraic equations (DAEs) of the type \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
  * \mathbf{u}, t) = \mathbf{0} \f$, where \f$ \mathbf{u} \f$ are the controls. The second
  * derivatives are given as Hessians of the multiplier-weighted function \f$ \boldsymbol{\lambda}^\top
  * \mathbf{F} \f$, split in blocks with respect to \f$ \mathbf{x} \f$, \f$ \mathbf{x}^{\prime} \f$,
  * and \f$ \mathbf{u} \f$. All Hessians default to central finite differences of the Jacobians,
  * override them to provide the exact ones.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam U The number of controls.
  */
  template <typename Real, Integer N, Integer U>
  class Controlled
  {
    static_assert(N > 0 && U > 0, "Sandals::Controlled: N and U must be positive.");

  public:
    using Pointer = std::shared_ptr<Controlled<Real, N, U>>; /**< Shared pointer to a controlled system. */
    using VectorF = Eigen::Vector<Real, N>; /**< States vector type. */
    using VectorU = Eigen::Vector<Real, U>; /**< Controls vector type. */
    using MatrixJF = Eigen::Matrix<Real, N, N>; /**< States Jacobian (and Hessian) matrix type. */
    using MatrixJU = Eigen::Matrix<Real, N, U>; /**< Controls Jacobian (and mixed Hessian) matrix type. */
    using MatrixUU = Eigen::Matrix<Real, U, U>; /**< Controls Hessian matrix type. */

  private:
    std::string m_name; /**< Name of the controlled system. */

  public:
    /**
    * Class constructor for a controlled ODE/DAE system.
    * \param[in] t_name The name of the controlled ODE/DAE system.
    */
    Controlled(std::string t_name = "(missing name)") : m_name(t_name) {}

    /**
    * Class destructor for a controlled ODE/DAE system.
    */
    virtual ~Controlled() {}

    /**
    * Get the controlled ODE/DAE system name const reference.
    * \return The controlled ODE/DAE system name const reference.
    */
    std::string const & name() const {return this->m_name;}

    /**
    * Get the number of equations of the controlled ODE/DAE system.
    * \return The number of equations.
    */
    static constexpr Integer equations_number() {return N;}

    /**
    * Get the number of controls of the controlled ODE/DAE system.
    * \return The number of controls.
    */
    static constexpr Integer controls_number() {return U;}

    /**
    * Evaluate the controlled ODE/DAE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, \mathbf{u},
    * t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{F} \f$.
    */
    virtual VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Jacobian \f$ \partial\mathbf{F}/\partial\mathbf{x} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}} \f$.
    */
    virtual MatrixJF JF_x(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Jacobian \f$ \partial\mathbf{F}/\partial\mathbf{x}^{\prime} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}} \f$.
    */
    virtual MatrixJF JF_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Jacobian \f$ \partial\mathbf{F}/\partial\mathbf{u} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{u}} \f$.
    */
    virtual MatrixJU JF_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{F})/\partial\mathbf{x}^2 \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The (symmetric) Hessian block.
    */
    virtual MatrixJF HF_x_x(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t,
      VectorF const & lambda) const
    {
      return finite_differences_hessian<Real, N>([&] (VectorF const & z)
        {return VectorF(this->JF_x(z, x_dot, u, t).transpose()*lambda);}, x);
    }

    /**
    * Evaluate the mixed Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{F})/\partial\mathbf{x}\,
    * \partial\mathbf{x}^{\prime} \f$, whose \f$ (i,j) \f$ entry is the derivative with respect to
    * \f$ x_i \f$ and \f$ x^{\prime}_j \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The mixed Hessian block.
    */
    virtual MatrixJF HF_x_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t,
      VectorF const & lambda) const
    {
      return finite_differences_jacobian<Real, N, N>([&] (VectorF const & z)
        {return VectorF(this->JF_x(x, z, u, t).transpose()*lambda);}, x_dot);
    }

    /**
    * Evaluate the Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{F})/\partial{\mathbf{x}^{
    * \prime}}^2 \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The (symmetric) Hessian block.
    */
    virtual MatrixJF HF_x_dot_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t,
      VectorF const & lambda) const
    {
      return finite_differences_hessian<Real, N>([&] (VectorF const & z)
        {return VectorF(this->JF_x_dot(x, z, u, t).transpose()*lambda);}, x_dot);
    }

    /**
    * Evaluate the mixed Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{F})/\partial\mathbf{x}\,
    * \partial\mathbf{u} \f$, whose \f$ (i,j) \f$ entry is the derivative with respect to \f$ x_i \f$
    * and \f$ u_j \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The mixed Hessian block.
    */
    virtual MatrixJU HF_x_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t,
      VectorF const & lambda) const
    {
      return finite_differences_jacobian<Real, N, U>([&] (VectorU const & z)
        {return VectorF(this->JF_x(x, x_dot, z, t).transpose()*lambda);}, u);
    }

    /**
    * Evaluate the mixed Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{F})/\partial
    * \mathbf{x}^{\prime}\,\partial\mathbf{u} \f$, whose \f$ (i,j) \f$ entry is the derivative with
    * respect to \f$ x^{\prime}_i \f$ and \f$ u_j \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The mixed Hessian block.
    */
    virtual MatrixJU HF_x_dot_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t,
      VectorF const & lambda) const
    {
      return finite_differences_jacobian<Real, N, U>([&] (VectorU const & z)
        {return VectorF(this->JF_x_dot(x, x_dot, z, t).transpose()*lambda);}, u);
    }

    /**
    * Evaluate the Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{F})/\partial\mathbf{u}^2 \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The (symmetric) Hessian block.
    */
    virtual MatrixUU HF_u_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t,
      VectorF const & lambda) const
    {
      return finite_differences_hessian<Real, U>([&] (VectorU const & z)
        {return VectorU(this->JF_u(x, x_dot, z, t).transpose()*lambda);}, u);
    }

  }; // class Controlled

  /*\
   |    ____            _             _ _          _ _____            _ _      _ _
   |   / ___|___  _ __ | |_ _ __ ___ | | | ___  __| | ____|_  ___ __ | (_) ___(_) |_
   |  | |   / _ \| '_ \| __| '__/ _ \| | |/ _ \/ _` |  _| \ \/ / '_ \| | |/ __| | __|
   |  | |__| (_) | | | | |_| | | (_) | | |  __/ (_| | |___ >  <| |_) | | | (__| | |_
   |   \____\___/|_| |_|\__|_|  \___/|_|_|\___|\__,_|_____/_/\_\ .__/|_|_|\___|_|\__|
   |                                                            |_|
  \*/

  /**
  * \brief Class container for a controlled system of explicit ODEs.
  *
  * Class container for a controlled system of explicit ordinary differential equations (ODEs) of the
  * type \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, \mathbf{u}, t) \f$, i.e., \f$ \mathbf{F} =
  * \mathbf{x}^{\prime} - \mathbf{f}(\mathbf{x}, \mathbf{u}, t) \f$. The Hessians of the weighted
  * function \f$ \boldsymbol{\lambda}^\top\mathbf{f} \f$ default to central finite differences of
  * the Jacobians, override them to provide the exact ones.
  * \tparam Real The scalar number type.
  * \tparam N The dimension of the ODE system.
  * \tparam U The number of controls.
  */
  template <typename Real, Integer N, Integer U>
  class ControlledExplicit : public Controlled<Real, N, U>
  {
  public:
    using typename Controlled<Real, N, U>::VectorF;
    using typename Controlled<Real, N, U>::VectorU;
    using typename Controlled<Real, N, U>::MatrixJF;
    using typename Controlled<Real, N, U>::MatrixJU;
    using typename Controlled<Real, N, U>::MatrixUU;

    /**
    * Class constructor for a controlled explicit ODE system.
    * \param[in] t_name The name of the controlled explicit ODE system.
    */
    ControlledExplicit(std::string t_name = "(missing name)") : Controlled<Real, N, U>(t_name) {}

    /**
    * Evaluate the explicit ODE system \f$ \mathbf{f}(\mathbf{x}, \mathbf{u}, t) \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The system function \f$ \mathbf{f} \f$.
    */
    virtual VectorF f(VectorF const & x, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Jacobian \f$ \partial\mathbf{f}/\partial\mathbf{x} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}} \f$.
    */
    virtual MatrixJF Jf_x(VectorF const & x, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Jacobian \f$ \partial\mathbf{f}/\partial\mathbf{u} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{u}} \f$.
    */
    virtual MatrixJU Jf_u(VectorF const & x, VectorU const & u, Real const t) const = 0;

    /**
    * Evaluate the Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{f})/\partial\mathbf{x}^2 \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The (symmetric) Hessian block.
    */
    virtual MatrixJF Hf_x_x(VectorF const & x, VectorU const & u, Real const t, VectorF const & lambda) const
    {
      return finite_differences_hessian<Real, N>([&] (VectorF const & z)
        {return VectorF(this->Jf_x(z, u, t).transpose()*lambda);}, x);
    }

    /**
    * Evaluate the mixed Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{f})/\partial\mathbf{x}\,
    * \partial\mathbf{u} \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The mixed Hessian block.
    */
    virtual MatrixJU Hf_x_u(VectorF const & x, VectorU const & u, Real const t, VectorF const & lambda) const
    {
      return finite_differences_jacobian<Real, N, U>([&] (VectorU const & z)
        {return VectorF(this->Jf_x(x, z, t).transpose()*lambda);}, u);
    }

    /**
    * Evaluate the Hessian \f$ \partial^2(\boldsymbol{\lambda}^\top\mathbf{f})/\partial\mathbf{u}^2 \f$.
    * \param[in] x States \f$ \mathbf{x} \f$.
    * \param[in] u Controls \f$ \mathbf{u} \f$.
    * \param[in] t Independent variable (or time) \f$ t \f$.
    * \param[in] lambda Multipliers \f$ \boldsymbol{\lambda} \f$.
    * \return The (symmetric) Hessian block.
    */
    virtual MatrixUU Hf_u_u(VectorF const & x, VectorU const & u, Real const t, VectorF const & lambda) const
    {
      return finite_differences_hessian<Real, U>([&] (VectorU const & z)
        {return VectorU(this->Jf_u(x, z, t).transpose()*lambda);}, u);
    }

    VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real const t) const override
    {return x_dot - this->f(x, u, t);}

    MatrixJF JF_x(VectorF const & x, VectorF const &, VectorU const & u, Real const t) const override
    {return -this->Jf_x(x, u, t);}

    MatrixJF JF_x_dot(VectorF const &, VectorF const &, VectorU const &, Real const) const override
    {return MatrixJF::Identity();}

    MatrixJU JF_u(VectorF const & x, VectorF const &, VectorU const & u, Real const t) const override
    {return -this->Jf_u(x, u, t);}

    MatrixJF HF_x_x(VectorF const & x, VectorF const &, VectorU const & u, Real const t,
      VectorF const & lambda) const override
    {return -this->Hf_x_x(x, u, t, lambda);}

    MatrixJF HF_x_x_dot(VectorF const &, VectorF const &, VectorU const &, Real const,
      VectorF const &) const override
    {return MatrixJF::Zero();}

    MatrixJF HF_x_dot_x_dot(VectorF const &, VectorF const &, VectorU const &, Real const,
      VectorF const &) const override
    {return MatrixJF::Zero();}

    MatrixJU HF_x_u(VectorF const & x, VectorF const &, VectorU const & u, Real const t,
      VectorF const & lambda) const override
    {return -this->Hf_x_u(x, u, t, lambda);}

    MatrixJU HF_x_dot_u(VectorF const &, VectorF const &, VectorU const &, Real const,
      VectorF const &) const override
    {return MatrixJU::Zero();}

    MatrixUU HF_u_u(VectorF const & x, VectorF const &, VectorU const & u, Real const t,
      VectorF const & lambda) const override
    {return -this->Hf_u_u(x, u, t, lambda);}

  }; // class ControlledExplicit

} // namespace Sandals

#endif // SANDALS_CONTROLLED_SYSTEM_HH
