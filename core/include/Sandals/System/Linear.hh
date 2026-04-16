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

#ifndef SANDALS_SYSTEM_LINEAR_HH
#define SANDALS_SYSTEM_LINEAR_HH

#include <Sandals.hh>
#include <Sandals/System/Explicit.hh>

namespace Sandals {

  /**
   * \brief Class container for the system of linear ODEs/DAEs.
   *
   * Class container for the system of linear ordinary differential equations
   * (ODEs) or differential algebraic equations (DAEs) of the type \f$
   * \mathbf{E}(t)\mathbf{x}^{\prime} = \mathbf{A}(t)
   * \mathbf{x} + \mathbf{b}(t) \f$, with invariants manifold \f$
   * \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
   * \tparam Real The scalar number type.
   * \tparam N The dimension of the linear ODE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class Linear : public Explicit<Real, N, M> {
   public:
    using Pointer = std::unique_ptr<Linear<Real, N, M>>;
    using typename Explicit<Real, N, M>::VectorF;
    using typename Explicit<Real, N, M>::MatrixJF;
    using MatrixE = typename Explicit<Real, N, M>::MatrixJF;
    using MatrixA = typename Explicit<Real, N, M>::MatrixJF;
    using VectorB = typename Explicit<Real, N, M>::VectorF;
    using typename Explicit<Real, N, M>::Type;

   private:
    Eigen::FullPivLU<MatrixE> m_lu; /**< LU decomposition. */

   public:
    /**
     * Class constructor for the linear ODE/DAE system.
     */
    Linear() : Explicit<Real, N, M>(Type::LINEAR, "(undefined name)") {}

    /**
     * Class constructor for the linear ODE/DAE system.
     * \param[in] t_name The name of the linear ODE/DAE system.
     */
    Linear(std::string t_name) : Explicit<Real, N, M>(Type::LINEAR, t_name) {}

    /**
     * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     *
     * \f[
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{E}(t)
     * \mathbf{x}^{\prime} -
     * \mathbf{A}(t) \mathbf{x} - \mathbf{b}(t)
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
      return this->E(t) * x_dot - this->A(t) * x - this->b(t);
    }

    /**
     * Evaluate the Jacobian of the ODE/DAE system function \f$
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) \f$ with respect to the
     * states \f$ \mathbf{x} \f$
     *
     * \f[
     * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}} =
     * -\mathbf{A}(t) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     */
    MatrixJF JF_x(const VectorF & /*x*/,
                  const VectorF & /*x_dot*/,
                  const Real t) const override {
      return -this->A(t);
    }

    /**
     * Evaluate the Jacobian of the ODE/DAE system function \f$
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
     * \f$ with respect to the states derivative \f$ \mathbf{x}^{\prime} \f$
     *
     * \f[
     * \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \displaystyle
     * \frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}^{\prime}} =
     * \mathbf{E}(t) \text{.}
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
                      const Real t) const override {
      return this->E(t);
    }

    /**
     * Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x},
     * t) \f$ as
     *
     * \f[
     * \mathbf{f}(\mathbf{x}, t) = \mathbf{E}(t)^{-1}(\mathbf{A}(t)\mathbf{x} +
     * \mathbf{b}(t))
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
     */
    VectorF f(const VectorF &x, const Real t) const override {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(
          this->m_lu.rank() == N,
          "Sandals:Linear::f(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t) * x + this->b(t));
    }

    /**
     * Evaluate the Jacobian of the explicit ODE/DAE system function \f$
     * \mathbf{f}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{f}(
     * \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{E}(t)^{-1} \mathbf{A}(t)
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    MatrixJF Jf_x(const VectorF & /*x*/, const Real t) const override {
      this->m_lu.compute(this->E(t));
      SANDALS_ASSERT(
          this->m_lu.rank() == N,
          "Sandals:Linear::Jf_x(...): singular mass matrix E(t) detected.");
      return this->m_lu.solve(this->A(t));
    }

    /**
     * Evaluate the linear ODE/DAE system mass matrix \f$ \mathbf{E}(t) \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system matrix \f$ \mathbf{E}(t) \f$.
     */
    virtual MatrixE E(const Real t) const = 0;

    /**
     * Evaluate the linear ODE/DAE system matrix \f$ \mathbf{A}(t) \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system matrix \f$ \mathbf{A}(t) \f$.
     */
    virtual MatrixA A(const Real t) const = 0;

    /**
     * Evaluate the linear ODE/DAE system vector \f$ \mathbf{b}(t) \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system vector \f$ \mathbf{b}(t) \f$.
     */
    virtual VectorB b(const Real t) const = 0;

  };  // class Linear

  /**
   * \brief Class container for the system of linear ODEs/DAEs wrapper.
   *
   * Class container for the system of linear ordinary differential equations
   * (ODEs) or differential algebraic equations (DAEs) of the type \f$
   * \mathbf{E}(t)\mathbf{x}^{\prime} = \mathbf{A}(t)
   * \mathbf{x} + \mathbf{b}(t) \f$, with invariants manifold \f$
   * \mathbf{h}(\mathbf{x}, t) = \mathbf{0} \f$.
   * \tparam Real The scalar number type.
   * \tparam N The dimension of the linear ODE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class LinearWrapper : public Linear<Real, N, M> {
   public:
    using Pointer = std::unique_ptr<LinearWrapper<Real, N, M>>;
    using typename Linear<Real, N, M>::VectorF;
    using typename Linear<Real, N, M>::MatrixJF;
    using MatrixE = typename Linear<Real, N, M>::MatrixJF;
    using MatrixA = typename Linear<Real, N, M>::MatrixJF;
    using VectorB = typename Linear<Real, N, M>::VectorF;
    using typename Explicit<Real, N, M>::VectorH;
    using typename Explicit<Real, N, M>::MatrixJH;
    using typename Explicit<Real, N, M>::TensorTH;

    using FunctionE  = std::function<MatrixE(const Real)>;
    using FunctionA  = std::function<MatrixA(const Real)>;
    using FunctionB  = std::function<VectorB(const Real)>;
    using FunctionH  = std::function<VectorH(const VectorF &, const Real)>;
    using FunctionJH = std::function<MatrixJH(const VectorF &, const Real)>;
    using FunctionID = std::function<bool(const VectorF &, const Real)>;
    using FunctionTH = std::function<TensorTH(const VectorF &, const Real)>;

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
    FunctionE m_E{nullptr};     /**< System's mass matrix. */
    FunctionA m_A{nullptr};     /**< System's right-hand-side matrix. */
    FunctionB m_b{nullptr};     /**< System's right-hand-side vector. */
    FunctionH m_h{nullptr};     /**< System's invariants. */
    FunctionJH m_Jh_x{nullptr}; /**< Jacobian of the system's invariants with
                                   respect to the states. */
    FunctionTH m_Th_x{nullptr}; /**< Tensor of the system's invariants with
                                   respect to the states. */
    FunctionID m_in_domain{nullptr}; /**< In-domain function. */

   public:
    /**
     * Class constructor for the linear ODE/DAE system wrapper.
     * \param[in] t_E The mass matrix.
     * \param[in] t_A The system matrix.
     * \param[in] t_b The system vector.
     * \param[in] t_h The system's invariants.
     * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to
     * the states.
     * \param[in] t_Th_x The tensor of the system's invariants with respect to
     * the states.
     * \param[in] t_in_domain The in-domain.
     */
    LinearWrapper(FunctionE t_E,
                  FunctionA t_A,
                  FunctionB t_b,
                  FunctionH t_h          = DefaultH,
                  FunctionJH t_Jh_x      = DefaultJH,
                  FunctionTH t_Th_x      = DefaultTH,
                  FunctionID t_in_domain = DefaultID)
        : Linear<Real, N, M>("(undefined name)"),
          m_E(t_E),
          m_A(t_A),
          m_b(t_b),
          m_h(t_h),
          m_Jh_x(t_Jh_x),
          m_Th_x(t_Th_x),
          m_in_domain(t_in_domain) {}

    /**
     * Class constructor for the linear ODE/DAE system wrapper.
     * \param[in] t_name The name of the linear ODE/DAE system.
     * \param[in] t_E The mass matrix.
     * \param[in] t_A The system matrix.
     * \param[in] t_b The system vector.
     * \param[in] t_h The system's invariants.
     * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to
     * the states.
     * \param[in] t_Th_x The tensor of the system's invariants with respect to
     * the states.
     * \param[in] t_in_domain The in-domain.
     */
    LinearWrapper(std::string t_name,
                  FunctionE t_E,
                  FunctionA t_A,
                  FunctionB t_b,
                  FunctionH t_h          = DefaultH,
                  FunctionJH t_Jh_x      = DefaultJH,
                  FunctionTH t_Th_x      = DefaultTH,
                  FunctionID t_in_domain = DefaultID)
        : Linear<Real, N, M>(t_name),
          m_E(t_E),
          m_A(t_A),
          m_b(t_b),
          m_h(t_h),
          m_Jh_x(t_Jh_x),
          m_Th_x(t_Th_x),
          m_in_domain(t_in_domain) {}

    /**
     * Class destructor for the linear ODE/DAE system wrapper.
     */
    ~LinearWrapper() {}

    /**
     * Get the mass matrix.
     * \return The mass matrix.
     */
    FunctionE &E() {
      return this->m_E;
    }

    /**
     * Get the system matrix.
     * \return The system matrix.
     */
    FunctionA &A() {
      return this->m_A;
    }

    /**
     * Get the system vector.
     * \return The system vector.
     */
    FunctionB &b() {
      return this->m_b;
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
     * Evaluate the linear ODE/DAE system mass matrix \f$ \mathbf{E}(t) \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system matrix \f$ \mathbf{E}(t) \f$.
     */
    MatrixE E(const Real t) const override {
      return this->m_E(t);
    }

    /**
     * Evaluate the linear ODE/DAE system matrix \f$ \mathbf{A}(t) \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system matrix \f$ \mathbf{A}(t) \f$.
     */
    MatrixA A(const Real t) const override {
      return this->m_A(t);
    }

    /**
     * Evaluate the linear ODE/DAE system vector \f$ \mathbf{b}(t) \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system vector \f$ \mathbf{b}(t) \f$.
     */
    VectorB b(const Real t) const override {
      return this->m_b(t);
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

  };  // class LinearWrapper

}  // namespace Sandals

#endif  // SANDALS_SYSTEM_LINEAR_HH
