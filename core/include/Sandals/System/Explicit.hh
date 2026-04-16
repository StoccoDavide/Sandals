/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the BSD 2-Clause License.        *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * davide.stocco@unitn.it                         enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_SYSTEM_EXPLICIT_HH
#define SANDALS_SYSTEM_EXPLICIT_HH

#include <Sandals.hh>
#include <Sandals/System/Implicit.hh>

namespace Sandals {

  /**
   * \brief Class container for the system of explicit ODEs.
   *
   * Class container for the system of explicit ordinary differential equations
   * (ODEs) of the type
   * \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, t) \f$, with invariants
   * manifold \f$ \mathbf{h}
   * (\mathbf{x}, t) = \mathbf{0} \f$.
   * \tparam Real The scalar number type.
   * \tparam N The dimension of the explicit ODE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class Explicit : public Implicit<Real, N, M> {
   public:
    using Pointer = std::unique_ptr<Explicit<Real, N, M>>;
    using typename Implicit<Real, N, M>::VectorF;
    using typename Implicit<Real, N, M>::MatrixJF;
    using typename Implicit<Real, N, M>::Type;

   protected:
    /**
     * Class constructor for an explicit ODE/DAE system.
     * \param[in] t_type The type of the explicit ODE/DAE system.
     * \param[in] t_name The name of the explicit ODE/DAE system.
     */
    Explicit(Type t_type, std::string t_name)
        : Implicit<Real, N, M>(t_type, t_name) {}

   public:
    /**
     * Class constructor for the explicit ODE system.
     */
    Explicit() : Implicit<Real, N, M>(Type::EXPLICIT, "(undefined name)") {}

    /**
     * Class constructor for the explicit ODE system.
     * \param[in] t_name The name of the explicit ODE system.
     */
    Explicit(std::string t_name)
        : Implicit<Real, N, M>(Type::EXPLICIT, t_name) {}

    /**
     * Evaluate the ODE system \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
     * t) \f$
     *
     * \f[
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \mathbf{x}^{\prime} - \mathbf{f}(\mathbf{x}, t) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     */
    VectorF F(const VectorF &x,
              const VectorF &x_dot,
              const Real t) const override {
      return x_dot - this->f(x, t);
    }

    /**
     * Evaluate the Jacobian of the implicit ODE system function \f$
     * \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
     *
     * \f[
     * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}} =
     * -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}} =
     * -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     */
    MatrixJF JF_x(const VectorF &x,
                  const VectorF & /*x_dot*/,
                  const Real t) const override {
      return -this->Jf_x(x, t);
    }

    /**
     * Evaluate the Jacobian of the implicit ODE system function \f$
     * \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$
     * \mathbf{x}^{\prime} \f$
     *
     * \f[
     * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \displaystyle
     * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}^{\prime}} =
     * \mathbf{I} \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     */
    MatrixJF JF_x_dot(const VectorF & /*x*/,
                      const VectorF & /*x_dot*/,
                      Real /*t*/) const override {
      return MatrixJF::Identity();
    }

    /**
     * Evaluate the explicit ODE system function\f$ \mathbf{f}(\mathbf{x}, t)
     * \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
     */
    virtual VectorF f(const VectorF &x, const Real t) const = 0;

    /**
     * Evaluate the Jacobian of the explicit ODE system function \f$
     * \mathbf{f}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    virtual MatrixJF Jf_x(const VectorF &x, const Real t) const = 0;

    /**
     * Time reversal of the explicit ODE system function \f$
     * \mathbf{f}(\mathbf{x}, t) = -\mathbf{f}(
     * \mathbf{x}, -t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The time-reversed system function \f$ -\mathbf{f}(\mathbf{x}, -t)
     * \f$.
     */
    VectorF f_reverse(const VectorF &x, const Real t) const {
      return -this->f(x, -t);
    }

    /**
     * Time reversal of the Jacobian of the explicit ODE system function \f$
     * \mathbf{f}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$, i.e., \f$ -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, -t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The time-reversed Jacobian \f$
     * -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, -t) \f$.
     */
    MatrixJF Jf_x_reverse(const VectorF &x, const Real t) const {
      return -this->Jf_x(x, -t);
    }

    /**
     * Time reversal of the implicit ODE system function \f$
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{
     * \prime}, t) = -\mathbf{F}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The time-reversed system function \f$ -\mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, -t) \f$.
     */
    VectorF F_reverse(const VectorF &x,
                      const VectorF &x_dot,
                      const Real t) const {
      return -x_dot + this->f(x, -t);
    }

    /**
     * Time reversal of the Jacobian of the implicit ODE system function \f$
     * \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$, i.e.,
     * \f$ -\mathbf{JF}_{\mathbf{x}}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The time-reversed Jacobian \f$
     * -\mathbf{JF}_{\mathbf{x}}(\mathbf{x}, -\mathbf{x}^{\prime}, -t) \f$.
     */
    MatrixJF JF_x_reverse(const VectorF &x,
                          const VectorF & /*x_dot*/,
                          const Real t) const {
      return -this->Jf_x(x, -t);
    }

    /**
     * Time reversal of the Jacobian of the implicit ODE system function \f$
     * \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$ with respect to the states derivative \f$
     * \mathbf{x}^{\prime} =
     * -\mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{\prime}, -t)
     * \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The time-reversed Jacobian \f$
     * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, -\mathbf{x}^{
     * \prime}, -t) \f$.
     */
    MatrixJF JF_x_dot_reverse(const VectorF & /*x*/,
                              const VectorF & /*x_dot*/,
                              Real /*t*/) const {
      return -MatrixJF::Identity();
    }

  };  // class Explicit

  /**
   * \brief Class container for the system of explicit ODEs/DAEs wrapper.
   *
   * Class container for the system of explicit ordinary differential equations
   * (ODEs) or differential algebraic equations (DAEs) of the type \f$
   * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
   * \mathbf{0} \f$, with invariants manifold \f$ \mathbf{h}(\mathbf{x}, t) =
   * \mathbf{0} \f$.
   * \tparam Real The scalar number type.
   * \tparam N The dimension of the explicit ODE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class ExplicitWrapper : public Explicit<Real, N, M> {
   public:
    using typename Explicit<Real, N, M>::VectorF;
    using typename Explicit<Real, N, M>::MatrixJF;
    using typename Explicit<Real, N, M>::VectorH;
    using typename Explicit<Real, N, M>::MatrixJH;
    using typename Explicit<Real, N, M>::TensorTH;

    using Pointer    = std::unique_ptr<ExplicitWrapper<Real, N, M>>;
    using FunctionF  = std::function<VectorF(const VectorF &, const Real)>;
    using FunctionJF = std::function<MatrixJF(const VectorF &, const Real)>;
    using FunctionH  = std::function<VectorH(const VectorF &, const Real)>;
    using FunctionJH = std::function<MatrixJH(const VectorF &, const Real)>;
    using FunctionTH = std::function<TensorTH(const VectorF &, const Real)>;
    using FunctionID = std::function<bool(const VectorF &, const Real)>;

    inline const static FunctionH DefaultH = [](const VectorF &, const Real) {
      return VectorH::Zero();
    }; /**< Default system's invariants. */
    inline const static FunctionJH DefaultJH = [](const VectorF &, const Real) {
      return MatrixJH::Zero();
    }; /**< Default system's invariants Jacobian. */
    inline const static FunctionTH DefaultTH = [](const VectorF &, const Real) {
      TensorTH Th_x;
      for (MatrixJH &m : Th_x) {
        m.setZero();
      }
      return Th_x;
    }; /**< Default system's invariants tensor. */
    inline const static FunctionID DefaultID = [](const VectorF &, const Real) {
      return true;
    }; /**< Default in-domain. */

   private:
    FunctionF m_f{nullptr};     /**< Explicit ODE system function. */
    FunctionJF m_Jf_x{nullptr}; /**< Jacobian of the explicit ODE system
                                   function. */
    FunctionH m_h{nullptr};     /**< System invariants. */
    FunctionJH m_Jh_x{nullptr}; /**< Jacobian of the system's invariants. */
    FunctionTH m_Th_x{nullptr}; /**< Tensor of the system's invariants. */
    FunctionID m_in_domain{nullptr}; /**< In-domain function. */

   public:
    /**
     * Class constructor for the explicit ODE/DAE system wrapper.
     * \param[in] t_f The explicit ODE system function.
     * \param[in] t_Jf_x The Jacobian of the explicit ODE system function with
     * respect to the states.
     * \param[in] t_h The system's invariants.
     * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to
     * the states.
     * * \param[in] t_Th_x The tensor of the system's invariants with respect to
     * the states.
     * \param[in] t_in_domain The in-domain function.
     */
    ExplicitWrapper(FunctionF t_f,
                    FunctionJF t_Jf_x,
                    FunctionH t_h          = DefaultH,
                    FunctionJH t_Jh_x      = DefaultJH,
                    FunctionTH t_Th_x      = DefaultTH,
                    FunctionID t_in_domain = DefaultID)
        : Explicit<Real, N, M>(),
          m_f(t_f),
          m_Jf_x(t_Jf_x),
          m_h(t_h),
          m_Jh_x(t_Jh_x),
          m_Th_x(t_Th_x),
          m_in_domain(t_in_domain) {}

    /**
     * Class constructor for the explicit ODE system wrapper.
     * \param[in] t_name The name of the explicit ODE system.
     * \param[in] t_f The explicit ODE system function.
     * \param[in] t_Jf_x The Jacobian of the explicit ODE system function with
     * respect to the states.
     * \param[in] t_h The system's invariants.
     * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to
     * the states.
     * \param[in] t_Th_x The tensor of the system's invariants with respect to
     * the states.
     * \param[in] t_in_domain The in-domain function.
     */
    ExplicitWrapper(std::string t_name,
                    FunctionF t_f,
                    FunctionJF t_Jf_x,
                    FunctionH t_h          = DefaultH,
                    FunctionJH t_Jh_x      = DefaultJH,
                    FunctionTH t_Th_x      = DefaultTH,
                    FunctionID t_in_domain = DefaultID)
        : Explicit<Real, N, M>(t_name),
          m_f(t_f),
          m_Jf_x(t_Jf_x),
          m_h(t_h),
          m_Jh_x(t_Jh_x),
          m_Th_x(t_Th_x),
          m_in_domain(t_in_domain) {}

    /**
     * Class destructor for the explicit ODE/DAE system wrapper.
     */
    ~ExplicitWrapper() {}

    /**
     * Get the explicit ODE system function.
     * \return The explicit ODE system function.
     */
    FunctionF &f() {
      return this->m_f;
    }

    /**
     * Get the Jacobian of the explicit ODE system function with respect to the
     * states.
     * \return The Jacobian of the explicit ODE system function with respect to
     * the states.
     */
    FunctionJF &Jf_x() {
      return this->m_Jf_x;
    }

    /**
     * Get the system's invariants.
     * \return The system's invariants.
     */
    FunctionH &h() {
      return this->m_h;
    }

    /**
     * Get the Jacobian of the system's invariants with respect to the states.
     * \return The Jacobian of the system's invariants with respect to the
     * states.
     */
    FunctionJH &Jh_x() {
      return this->m_Jh_x;
    }

    /**
     * Get the tensor of the system's invariants with respect to the states.
     * \return The tensor of the system's invariants with respect to the states.
     */
    FunctionTH &Th_x() {
      return this->m_Th_x;
    }

    /**
     * Get the in-domain function.
     * \return The in-domain function.
     */
    FunctionID &in_domain() {
      return this->m_in_domain;
    }

    /**
     * Evaluate the ODE/DAE system \f$ \mathbf{f}(\mathbf{x}, t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
     */
    VectorF f(const VectorF &x, const Real t) const override {
      return this->m_f(x, t);
    }

    /**
     * Evaluate the Jacobian of the explicit ODE system function \f$
     * \mathbf{f}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    MatrixJF Jf_x(const VectorF &x, const Real t) const override {
      return this->m_Jf_x(x, t);
    }

    /**
     * Evaluate the ODE/DAE system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system invariants \f$ \mathbf{h}(\mathbf{x}, t) \f$.
     */
    VectorH h(const VectorF &x, const Real t) const override {
      return this->m_h(x, t);
    }

    /**
     * Evaluate the Jacobian of the ODE/DAE system invariants \f$
     * \mathbf{h}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Jh}_{\mathbf{x}}(\mathbf{x} t) =
     * \displaystyle\frac{\partial\mathbf{h}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jh}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    MatrixJH Jh_x(const VectorF &x, const Real t) const override {
      return this->m_Jh_x(x, t);
    }

    /**
     * Evaluate the tensor of the Jacobian of the ODE/DAE system invariants \f$
     * \mathbf{h}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Th}_{\mathbf{x}\mathbf{x}}(\mathbf{x} t) =
     * \displaystyle\frac{\partial^2\mathbf{h}(\mathbf{x},
     * t)}{\partial\mathbf{x}\partial\mathbf{x}} \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The tensor of the Jacobian \f$
     * \mathbf{Th}_{\mathbf{x}\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    TensorTH Th_x(const VectorF &x, const Real t) const override {
      return this->m_Th_x(x, t);
    }

    /**
     * Return true if the values \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the
     * domain of the ODE/DAE system.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return True if \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the domain of the
     * ODE/DAE system.
     */
    bool in_domain(const VectorF &x, const Real t) const override {
      return this->m_in_domain(x, t);
    }

  };  // class ExplicitWrapper

}  // namespace Sandals

#endif  // SANDALS_SYSTEM_EXPLICIT_HH
