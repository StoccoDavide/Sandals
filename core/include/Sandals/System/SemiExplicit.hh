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

#ifndef SANDALS_SYSTEM_SEMIEXPLICIT_HH
#define SANDALS_SYSTEM_SEMIEXPLICIT_HH

#include <Sandals.hh>
#include <Sandals/System/Explicit.hh>

namespace Sandals {

  /**
   * \brief Class container for the system of semi-explicit ODEs/DAEs.
   *
   * Class container for the system of semi-explicit ordinary differential
   * equations (ODEs) or differential algebraic equations (DAEs)of the type \f$
   * \mathbf{A}(\mathbf{x}, t)
   * \mathbf{x}^{\prime} = \mathbf{b}(\mathbf{x}, t) \f$, with invariants
   * manifold \f$ \mathbf{h}(
   * \mathbf{x}, t) = \mathbf{0} \f$.
   * \tparam Real const The scalar number type.
   * \tparam N The dimension of the semi-explicit ODE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class SemiExplicit : public Explicit<Real, N, M> {
   public:
    using typename Explicit<Real, N, M>::VectorF;
    using typename Explicit<Real, N, M>::MatrixJF;
    using MatrixA  = typename Explicit<Real, N, M>::MatrixJF;
    using TensorTA = typename std::array<MatrixJF, N>;
    using VectorB  = typename Explicit<Real, N, M>::VectorF;
    using MatrixJB = typename Explicit<Real, N, M>::MatrixJF;
    using typename Explicit<Real, N, M>::Type;

   private:
    mutable Eigen::FullPivLU<MatrixA> m_lu; /**< LU decomposition. */

   public:
    /**
     * Class constructor for the semi-explicit ODE/DAE system.
     */
    SemiExplicit()
        : Explicit<Real, N, M>(Type::SEMIEXPLICIT, "(undefined name)") {}

    /**
     * Class constructor for the semi-explicit ODE/DAE system.
     * \param[in] t_name The name of the semi-explicit ODE/DAE system.
     */
    SemiExplicit(std::string t_name)
        : Explicit<Real, N, M>(Type::SEMIEXPLICIT, t_name) {}

    /**
     * Evaluate the ODE/DAE system \f$ \mathbf{F}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$
     *
     * \f[
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \mathbf{A}(\mathbf{x}, t)\mathbf{x}^{\prime} - \mathbf{b}(\mathbf{x}, t)
     * \text{.}
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
      return this->A(x, t) * x_dot - this->b(x, t);
    }

    /**
     * Evaluate the Jacobian of the ODE/DAE system function \f$
     * \mathbf{F}(\mathbf{x}, \mathbf{x}^{
     * \prime}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
     *
     * \f[
     * \mathbf{JF}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
     * \displaystyle\frac{\partial\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}} =
     * -\displaystyle\frac{\partial\mathbf{f}(\mathbf{x}, \mathbf{x}^{\prime},
     * t)}{\partial\mathbf{x}} =
     * -\mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) = -\mathbf{A}(\mathbf{x},
     * t)^{-1}(\mathbf{Jb}_{
     * \mathbf{x}}(\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t)
     * \mathbf{x}^\prime) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     */
    MatrixJF JF_x(const VectorF &x,
                  const VectorF &x_dot,
                  const Real t) const override {
      return -this->Jf_x(x, x_dot, t);
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
     * \mathbf{A}(\mathbf{x}, t) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{JF}_{\mathbf{x}^{\prime}}(\mathbf{x},
     * \mathbf{x}^{\prime}, t) \f$.
     */
    MatrixJF JF_x_dot(const VectorF &x,
                      const VectorF & /*x_dot*/,
                      const Real t) const override {
      return this->A(x, t);
    }

    /**
     * Evaluate the explicit ODE/DAE system function\f$ \mathbf{f}(\mathbf{x},
     * t) \f$ as
     *
     * \f[
     * \mathbf{f}(\mathbf{x}, t) = \mathbf{A}(\mathbf{x}, t)^{-1}
     * \mathbf{b}(\mathbf{x}, t) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{f}(\mathbf{x}, t) \f$.
     */
    VectorF f(const VectorF &x, const Real t) const override {
      this->m_lu.compute(this->A(x, t));
      SANDALS_ASSERT(this->m_lu.rank() == N,
                     "Sandals:SemiExplicit::f(...): singular mass matrix A.");
      return this->m_lu.solve(this->b(x, t));
    }

    /**
     * Evaluate the Jacobian of the explicit ODE/DAE system function \f$
     * \mathbf{f}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, \mathbf{x}^\prime, t) =
     * \displaystyle\frac{\partial
     * \mathbf{f}(\mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{A}(\mathbf{x},
     * t)^{-1}(\mathbf{Jb}_{
     * \mathbf{x}}(\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t)
     * \mathbf{x}^\prime) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] x_dot States derivative \f$ \mathbf{x}^{\prime} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x},
     * \mathbf{x}^\prime, t) \f$.
     */
    MatrixJF Jf_x(const VectorF &x, const VectorF &x_dot, const Real t) const {
      TensorTA TA_x{this->TA_x(x, t)};
      MatrixJF tAp{MatrixJF::Zero()};
      for (Integer i{0}; i < N; ++i) {
        tAp.col(i) = TA_x[i] * x_dot;
      }  // TODO: Check if TA_x.size() = N
      this->m_lu.compute(this->A(x, t));
      SANDALS_ASSERT(
          this->m_lu.rank() == N,
          "Sandals:SemiExplicit::Jf_x(...): singular mass matrix A.");
      return this->m_lu.solve(this->Jb_x(x, t) - tAp);
    }

    /**
     * Evaluate the Jacobian of the explicit ODE/DAE system function \f$
     * \mathbf{f}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{f}(
     * \mathbf{x}, t)}{\partial\mathbf{x}} = \mathbf{A}(\mathbf{x}, t)^{-1}(
     * \mathbf{Jb}_{\mathbf{x}}
     * (\mathbf{x}, t) - \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t)
     * \mathbf{f}(\mathbf{x}, t)) \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jf}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    MatrixJF Jf_x(const VectorF &x, const Real t) const override {
      return this->Jf_x(x, this->f(x, t), t);
    }

    /**
     * Evaluate the semi-explicit ODE/DAE system mass matrix \f$
     * \mathbf{A}(\mathbf{x}, t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{A}(\mathbf{x}, t) \f$.
     */
    virtual MatrixA A(const VectorF &x, const Real t) const = 0;

    /**
     * Evaluate the tensor of the semi-explicit ODE/DAE system mass matrix \f$
     * \mathbf{A}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{A}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    virtual TensorTA TA_x(const VectorF &x, const Real t) const = 0;

    /**
     * Evaluate the semi-explicit ODE/DAE system right-hand-side function \f$
     * \mathbf{b}(\mathbf{x}, t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
     */
    virtual VectorB b(const VectorF &x, const Real t) const = 0;

    /**
     * Evaluate the Jacobian of the semi-explicit ODE/DAE system right-hand-side
     * function \f$ \mathbf{b}
     * (\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
     *
     * \f[
     * \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{b}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    virtual MatrixJB Jb_x(const VectorF &x, const Real t) const = 0;

  };  // class SemiExplicit

  /**
   * \brief Class container for the system of semi-explicit ODEs/DAEs wrapper.
   *
   * Class container for the system of semi-explicit ordinary differential
   * equations (ODEs) or differential algebraic equations (DAEs)of the type \f$
   * \mathbf{A}(\mathbf{x}, t)
   * \mathbf{x}^{\prime} = \mathbf{b}(\mathbf{x}, t) \f$, with invariants
   * manifold \f$ \mathbf{h}(
   * \mathbf{x}, t) = \mathbf{0} \f$.
   * \tparam Real const The scalar number type.
   * \tparam N The dimension of the semi-explicit ODE system.
   * \tparam M The dimension of the invariants manifold.
   */
  template <typename Real, Integer N, Integer M = 0>
  class SemiExplicitWrapper : public SemiExplicit<Real, N, M> {
   public:
    using typename SemiExplicit<Real, N, M>::VectorF;
    using typename SemiExplicit<Real, N, M>::MatrixA;
    using typename SemiExplicit<Real, N, M>::TensorTA;
    using typename SemiExplicit<Real, N, M>::VectorB;
    using typename SemiExplicit<Real, N, M>::MatrixJB;
    using typename SemiExplicit<Real, N, M>::VectorH;
    using typename SemiExplicit<Real, N, M>::MatrixJH;
    using typename SemiExplicit<Real, N, M>::TensorTH;

    using Pointer    = std::unique_ptr<SemiExplicitWrapper<Real, N, M>>;
    using FunctionA  = std::function<MatrixA(const VectorF &, const Real)>;
    using FunctionTA = std::function<TensorTA(const VectorF &, const Real)>;
    using FunctionB  = std::function<VectorB(const VectorF &, const Real)>;
    using FunctionJB = std::function<MatrixJB(const VectorF &, const Real)>;
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
    FunctionA m_A{nullptr};     /**< System's mass matrix. */
    FunctionTA m_TA_x{nullptr}; /**< Jacobian of the system's mass matrix with
                                   respect to the states. */
    FunctionB m_b{nullptr};     /**< System's Right-hand-side. */
    FunctionJB m_Jb_x{nullptr}; /**< Jacobian of the system's right-hand-side
                                   with respect to the states. */
    FunctionH m_h{nullptr};     /**< System's invariants. */
    FunctionJH m_Jh_x{nullptr}; /**< Jacobian of the system's invariants with
                                   respect to the states. */
    FunctionTH m_Th_x{nullptr}; /**< Tensor of the system's invariants with
                                   respect to the states. */
    FunctionID m_in_domain{nullptr}; /**< In-domain function. */

   public:
    /**
     * Class constructor for the semi-explicit ODE/DAE system wrapper.
     * \param[in] t_A The systems's mass matrix.
     * \param[in] t_TA_x The Jacobian of the systems's mass matrix with respect
     * to the states.
     * \param[in] t_b The systems's right-hand-side.
     * \param[in] t_Jb_x The Jacobian of the right-hand-side with respect to the
     * states.
     * \param[in] t_h The system's invariants.
     * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to
     * the states.
     * \param[in] t_Th_x The tensor of the system's invariants with respect to
     * the states.
     * \param[in] t_in_domain The in-domain.
     */
    SemiExplicitWrapper(FunctionA t_A,
                        FunctionTA t_TA_x,
                        FunctionB t_b,
                        FunctionJB t_Jb_x,
                        FunctionH t_h          = DefaultH,
                        FunctionJH t_Jh_x      = DefaultJH,
                        FunctionTH t_Th_x      = DefaultTH,
                        FunctionID t_in_domain = DefaultID)
        : SemiExplicit<Real, N, M>(),
          m_A(t_A),
          m_TA_x(t_TA_x),
          m_b(t_b),
          m_Jb_x(t_Jb_x),
          m_h(t_h),
          m_Jh_x(t_Jh_x),
          m_Th_x(t_Th_x),
          m_in_domain(t_in_domain) {}

    /**
     * Class constructor for the semi-explicit ODE/DAE system wrapper.
     * \param[in] t_name The name of the semi-explicit ODE/DAE system.
     * \param[in] t_A The systems's mass matrix.
     * \param[in] t_TA_x The Jacobian of the systems's mass matrix with respect
     * to the states.
     * \param[in] t_b The systems's right-hand-side.
     * \param[in] t_Jb_x The Jacobian of the right-hand-side with respect to the
     * states.
     * \param[in] t_h The system's invariants.
     * \param[in] t_Jh_x The Jacobian of the system's invariants with respect to
     * the states.
     * \param[in] t_Th_x The tensor of the system's invariants with respect to
     * the states.
     * \param[in] t_in_domain The in-domain.
     */
    SemiExplicitWrapper(std::string t_name,
                        FunctionA t_A,
                        FunctionTA t_TA_x,
                        FunctionB t_b,
                        FunctionJB t_Jb_x,
                        FunctionH t_h          = DefaultH,
                        FunctionJH t_Jh_x      = DefaultJH,
                        FunctionTH t_Th_x      = DefaultTH,
                        FunctionID t_in_domain = DefaultID)
        : SemiExplicit<Real, N, M>(t_name),
          m_A(t_A),
          m_TA_x(t_TA_x),
          m_b(t_b),
          m_Jb_x(t_Jb_x),
          m_h(t_h),
          m_Jh_x(t_Jh_x),
          m_Th_x(t_Th_x),
          m_in_domain(t_in_domain) {}

    /**
     * Class destructor for the semi-explicit ODE/DAE system wrapper.
     */
    ~SemiExplicitWrapper() {}

    /**
     * Get the the mass matrix.
     * \return The the mass matrix.
     */
    FunctionA &A() {
      return this->m_A;
    }

    /**
     * Get the Jacobian of the mass matrix tensor.
     * \return The Jacobian of the mass matrix tensor.
     */
    FunctionTA &TA_x() {
      return this->m_TA_x;
    }

    /**
     * Get the right-hand-side.
     * \return The right-hand-side.
     */
    FunctionB &b() {
      return this->m_b;
    }

    /**
     * Get the Jacobian of the right-hand-side.
     * \return The Jacobian of the right-hand-side.
     */
    FunctionJB &Jb_x() {
      return this->m_Jb_x;
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
     * Evaluate the semi-explicit ODE/DAE system mass matrix \f$
     * \mathbf{A}(\mathbf{x}, t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system mass matrix \f$ \mathbf{A}(\mathbf{x}, t) \f$.
     */
    MatrixA A(const VectorF &x, const Real t) const override {
      return this->m_A(x, t);
    }

    /**
     * Evaluate the tensor of the semi-explicit ODE/DAE system mass matrix \f$
     * \mathbf{A}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{A}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{TA}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    TensorTA TA_x(const VectorF &x, const Real t) const override {
      return this->m_TA_x(x, t);
    }

    /**
     * Evaluate the semi-explicit ODE/DAE system right-hand-side function \f$
     * \mathbf{b}(\mathbf{x}, t) \f$.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The system function \f$ \mathbf{b}(\mathbf{x}, t) \f$.
     */
    VectorB b(const VectorF &x, const Real t) const override {
      return this->m_b(x, t);
    }

    /**
     * Evaluate the Jacobian of the semi-explicit ODE/DAE system right-hand-side
     * function \f$ \mathbf{b}
     * (\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x} \f$
     *
     * \f[
     * \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{b}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Jb}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    MatrixJB Jb_x(const VectorF &x, const Real t) const override {
      return this->m_Jb_x(x, t);
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
     * Evaluate the tensor of the ODE/DAE system invariants \f$
     * \mathbf{h}(\mathbf{x}, t) \f$ with respect to the states \f$ \mathbf{x}
     * \f$
     *
     * \f[
     * \mathbf{Th}_{\mathbf{x}}(\mathbf{x}, t) =
     * \displaystyle\frac{\partial\mathbf{h}(\mathbf{x}, t)}{\partial\mathbf{x}}
     * \text{.}
     * \f]
     *
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return The Jacobian \f$ \mathbf{Th}_{\mathbf{x}}(\mathbf{x}, t) \f$.
     */
    TensorTH Th_x(const VectorF &x, const Real t) const override {
      return this->m_Th_x(x, t);
    }

    /**
     * Return true if the values \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the
     * domain of the ODE/DAE system.
     * \param[in] x States \f$ \mathbf{x} \f$.
     * \param[in] t Independent variable (or time) \f$ t \f$.
     * \return True if \f$ \mathbf{f}(\mathbf{x}, t) \f$ is in the domain of
     * the ODE/DAE system.
     */
    bool in_domain(const VectorF &x, const Real t) const override {
      return this->m_in_domain(x, t);
    }

  };  // class SemiExplicitWrapper

}  // namespace Sandals

#endif  // SANDALS_SYSTEM_SEMIEXPLICIT_HH
