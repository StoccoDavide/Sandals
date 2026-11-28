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

#ifndef SANDALS_OPTIMAL_CONTROL_AD_HH
#define SANDALS_OPTIMAL_CONTROL_AD_HH

#include <Sandals/AutoDiff.hh>
#include <Sandals/System/Controlled.hh>
#include <Sandals/OptimalControl.hh>

namespace Sandals
{

  /*\
   |    ___   ____ ____       _    ____
   |   / _ \ / ___|  _ \     / \  |  _ \
   |  | | | | |   | |_) |   / _ \ | | | |
   |  | |_| | |___|  __/   / ___ \| |_| |
   |   \___/ \____|_|     /_/   \_\____/
   |
  \*/

  /**
  * \brief Optimal control problem of a DAE with invariants, with automatic derivatives.
  *
  * The derived class writes the model once, as templates on the scalar type T:
  *
  * \code
  * template <typename T> Eigen::Vector<T, N> dae(x, x_dot, u, t) const;      // F(x, x', u, t) (required)
  * template <typename T> Eigen::Vector<T, M> invariants(x, t) const;         // h(x, t)
  * template <typename T> T running_cost(x, u, t) const;                      // Lagrange term L
  * template <typename T> T terminal_cost(x_ini, x_end) const;                // Mayer term phi
  * template <typename T> Eigen::Vector<T, NB> boundary(x_ini, x_end) const;  // boundary conditions b
  * template <typename T> Eigen::Vector<T, NG> path(x, u, t) const;           // path constraints g
  * \endcode
  *
  * with arguments of type Eigen::Vector<T, ...> and T for the time. The optional terms default to
  * absent. Each function is recorded once by CppAD (at the first use), and all the first and second
  * derivatives needed by the solvers (and the third derivatives of the invariants) are computed
  * exactly. Math functions must be called unqualified (e.g., `sin(x)`, not `std::sin(x)`), so that the
  * CppAD overloads are found, and data-dependent branches must use CppAD::CondExp.
  *
  * The class implements the Sandals::OptimalControlProblem and Sandals::Controlled interfaces, so it
  * can be solved by any transcription (multiple shooting, collocation).
  * \tparam Derived The derived class (CRTP).
  * \tparam N The number of states.
  * \tparam U The number of controls.
  * \tparam M The number of invariants.
  * \tparam NB The number of boundary conditions.
  * \tparam NG The number of path constraints.
  */
  template <typename Derived, Integer N, Integer U, Integer M = 0, Integer NB = 0, Integer NG = 0>
  class OptimalControlProblemAD : public OptimalControlProblem<double, N, U, M>
  {
  public:
    using Base    = OptimalControlProblem<double, N, U, M>; /**< Optimal control problem interface. */
    using AD      = CppAD::AD<double>;                       /**< AD scalar type. */
    using VectorX = Eigen::Vector<double, Eigen::Dynamic>;   /**< Dynamic vector type. */
    using MatrixX = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; /**< Dynamic matrix type. */
    using typename Base::VectorF;
    using typename Base::VectorU;
    using typename Base::VectorB;
    using typename Base::MatrixB;
    using typename Base::VectorP;
    using typename Base::MatrixP;

    template <typename T> using StateVector    = Eigen::Vector<T, N>;  /**< States of scalar type T. */
    template <typename T> using ControlVector  = Eigen::Vector<T, U>;  /**< Controls of scalar type T. */
    template <typename T> using InvariantVector = Eigen::Vector<T, M>; /**< Invariants of scalar type T. */
    template <typename T> using BoundaryVector = Eigen::Vector<T, NB>; /**< Boundary conditions of scalar type T. */
    template <typename T> using PathVector     = Eigen::Vector<T, NG>; /**< Path constraints of scalar type T. */

    /*\
     |  Default (absent) optional terms, hidden by the derived class when defined
    \*/

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const &, T const &) const {return InvariantVector<T>::Zero();}

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const &, T const &) const {return T(0.0);}

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const &) const {return T(0.0);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const &, StateVector<T> const &) const {return BoundaryVector<T>::Zero();}

    template <typename T>
    PathVector<T> path(StateVector<T> const &, ControlVector<T> const &, T const &) const {return PathVector<T>::Zero();}

  private:
    /** The controlled DAE system, evaluated through the recorded model. */
    class System;

    /** The recorded model functions. */
    struct Tapes
    {
      Tape dae;            // [x, x', u, t]  -> F
      Tape invariants;     // [x, t]         -> h
      Tape curvature;      // [x, t, mu, w]  -> mu^T Jh_x(x, t) w (for the third derivatives of h)
      Tape running_cost;   // [x, u, t]      -> L
      Tape terminal_cost;  // [x_ini, x_end] -> phi
      Tape boundary;       // [x_ini, x_end] -> b
      Tape path;           // [x, u, t]      -> g
    };

    mutable std::unique_ptr<Tapes> m_tapes; /**< Recorded model (at the first use). */

  public:
    /**
    * Class constructor.
    * \param[in] t_name The name of the problem.
    * \param[in] t_t_ini The initial time.
    * \param[in] t_t_end The final time.
    */
    OptimalControlProblemAD(std::string t_name, double const t_t_ini, double const t_t_end)
      : Base(t_name, std::make_shared<System>(*this), t_t_ini, t_t_end)
    {}

    /** Discard the recorded model, e.g., after changing a parameter of the derived class. */
    void rerecord() {this->m_tapes.reset();}

    /*\
     |  Structural sparsity patterns, from the recorded functions
    \*/

    using typename Base::MatrixXb;

    MatrixXb dae_jacobian_pattern() const override {return this->tapes().dae.jacobian_pattern().leftCols(2*N + U);}

    MatrixXb dae_hessian_pattern() const override {return this->tapes().dae.hessian_pattern().topLeftCorner(2*N + U, 2*N + U);}

    MatrixXb lagrange_hessian_pattern() const override {return this->tapes().running_cost.hessian_pattern().topLeftCorner(N + U, N + U);}

    MatrixXb path_jacobian_pattern() const override
    {
      if constexpr (NG == 0) {return MatrixXb(0, N + U);}
      else {return this->tapes().path.jacobian_pattern().leftCols(N + U);}
    }

    MatrixXb path_hessian_pattern() const override
    {
      if constexpr (NG == 0) {return MatrixXb::Constant(N + U, N + U, false);}
      else {return this->tapes().path.hessian_pattern().topLeftCorner(N + U, N + U);}
    }

    MatrixXb invariants_jacobian_pattern() const override
    {
      if constexpr (M == 0) {return MatrixXb(0, N);}
      else {return this->tapes().invariants.jacobian_pattern().leftCols(N);}
    }

    MatrixXb invariants_hessian_pattern() const override
    {
      if constexpr (M == 0) {return MatrixXb::Constant(N, N, false);}
      else {
        // second derivatives of h, and the contraction of its third derivatives (the curvature tape)
        MatrixXb P{this->tapes().invariants.hessian_pattern().topLeftCorner(N, N)};
        MatrixXb const T{this->tapes().curvature.hessian_pattern().topLeftCorner(N, N)};
        return P.array() || T.array();
      }
    }

    /*\
     |  Optimal control problem interface
    \*/

    double mayer(VectorF const & x_ini, VectorF const & x_end) const override
    {return this->tapes().terminal_cost.value(join(x_ini, x_end))(0);}

    VectorB mayer_gradient(VectorF const & x_ini, VectorF const & x_end) const override
    {return this->tapes().terminal_cost.jacobian(join(x_ini, x_end)).transpose();}

    MatrixB mayer_hessian(VectorF const & x_ini, VectorF const & x_end) const override
    {return this->tapes().terminal_cost.hessian(join(x_ini, x_end), VectorX::Ones(1));}

    double lagrange(VectorF const & x, VectorU const & u, double const t) const override
    {return this->tapes().running_cost.value(join(x, u, t))(0);}

    VectorP lagrange_gradient(VectorF const & x, VectorU const & u, double const t) const override
    {return this->tapes().running_cost.jacobian(join(x, u, t)).transpose().template topRows<N+U>();}

    MatrixP lagrange_hessian(VectorF const & x, VectorU const & u, double const t) const override
    {return this->tapes().running_cost.hessian(join(x, u, t), VectorX::Ones(1)).template topLeftCorner<N+U, N+U>();}

    Integer boundary_conditions_number() const override {return NB;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & x_end) const override
    {
      if constexpr (NB == 0) {return VectorX(0);}
      else {return this->tapes().boundary.value(join(x_ini, x_end));}
    }

    MatrixX boundary_conditions_jacobian(VectorF const & x_ini, VectorF const & x_end) const override
    {return this->tapes().boundary.jacobian(join(x_ini, x_end));}

    MatrixB boundary_conditions_hessian(VectorF const & x_ini, VectorF const & x_end, VectorX const & lambda) const override
    {return this->tapes().boundary.hessian(join(x_ini, x_end), lambda);}

    Integer path_constraints_number() const override {return NG;}

    VectorX path_constraints(VectorF const & x, VectorU const & u, double const t) const override
    {
      if constexpr (NG == 0) {return VectorX(0);}
      else {return this->tapes().path.value(join(x, u, t));}
    }

    MatrixX path_constraints_jacobian(VectorF const & x, VectorU const & u, double const t) const override
    {return this->tapes().path.jacobian(join(x, u, t)).leftCols(N+U);}

    MatrixP path_constraints_hessian(VectorF const & x, VectorU const & u, double const t, VectorX const & lambda) const override
    {return this->tapes().path.hessian(join(x, u, t), lambda).template topLeftCorner<N+U, N+U>();}

  private:
    Derived const & derived() const {return static_cast<Derived const &>(*this);}

    /** Stack vectors and scalars into a single vector. */
    template <typename... Args>
    static VectorX join(Args const & ... args)
    {
      VectorX v(Integer((size(args) + ...)));
      Integer offset{0};
      (put(v, offset, args), ...);
      return v;
    }
    template <typename V> static Integer size(V const & v) {return static_cast<Integer>(v.size());}
    static Integer size(double) {return 1;}
    template <typename V> static void put(VectorX & v, Integer & offset, V const & a) {v.segment(offset, a.size()) = a; offset += a.size();}
    static void put(VectorX & v, Integer & offset, double a) {v(offset++) = a;}

    /** Record the model functions at the first use. */
    Tapes const & tapes() const
    {
      if (!this->m_tapes) {
        this->m_tapes = std::make_unique<Tapes>();
        this->record(*this->m_tapes);
      }
      return *this->m_tapes;
    }

    /** Record every model function of the derived class. */
    void record(Tapes & t) const
    {
      Derived const & model{this->derived()};

      t.dae.record(2*N + U + 1, [&model] (std::vector<AD> const & v) {
        return to_std_vector(model.template dae<AD>(segment<N>(v, 0), segment<N>(v, N), segment<U>(v, 2*N), v[2*N + U]));
      });

      t.running_cost.record(N + U + 1, [&model] (std::vector<AD> const & v) {
        return std::vector<AD>{model.template running_cost<AD>(segment<N>(v, 0), segment<U>(v, N), v[N + U])};
      });

      t.terminal_cost.record(2*N, [&model] (std::vector<AD> const & v) {
        return std::vector<AD>{model.template terminal_cost<AD>(segment<N>(v, 0), segment<N>(v, N))};
      });

      if constexpr (NB > 0) {
        t.boundary.record(2*N, [&model] (std::vector<AD> const & v) {
          return to_std_vector(model.template boundary<AD>(segment<N>(v, 0), segment<N>(v, N)));
        });
      }

      if constexpr (NG > 0) {
        t.path.record(N + U + 1, [&model] (std::vector<AD> const & v) {
          return to_std_vector(model.template path<AD>(segment<N>(v, 0), segment<U>(v, N), v[N + U]));
        });
      }

      if constexpr (M > 0) {
        t.invariants.record(N + 1, [&model] (std::vector<AD> const & v) {
          return to_std_vector(model.template invariants<AD>(segment<N>(v, 0), v[N]));
        });
        this->record_curvature(t.curvature);
      }
    }

    /**
    * Record the function \f$ \psi(\mathbf{x}, t, \boldsymbol{\mu}, \mathbf{w}) = \boldsymbol{\mu}^\top
    * \mathbf{Jh}_{\mathbf{x}}(\mathbf{x}, t)\,\mathbf{w} \f$, whose Hessian in \f$ \mathbf{x} \f$ gives
    * the third derivatives of the invariants. The Jacobian of the invariants is itself computed by
    * CppAD, recorded on AD<AD<double>> scalars (nested taping).
    */
    void record_curvature(Tape & curvature) const
    {
      using AD2 = CppAD::AD<AD>;
      Derived const & model{this->derived()};

      // Invariants on AD<double> base, to differentiate them while recording psi
      std::vector<AD2> v2(N + 1, AD2(0.5));
      CppAD::Independent(v2);
      std::vector<AD2> h2{to_std_vector(model.template invariants<AD2>(segment<N>(v2, 0), v2[N]))};
      CppAD::ADFun<AD> invariants(v2, h2);

      curvature.record(N + 1 + M + N, [&invariants] (std::vector<AD> const & v) {
        std::vector<AD> const xt(v.begin(), v.begin() + N + 1);
        std::vector<AD> const Jh{invariants.Jacobian(xt)}; // M x (N+1), row-major
        AD psi(0.0);
        for (Integer j{0}; j < M; ++j) {
          for (Integer a{0}; a < N; ++a) {psi += v[N + 1 + j] * Jh[j*(N + 1) + a] * v[N + 1 + M + a];}
        }
        return std::vector<AD>{psi};
      });
    }
  };

  /*\
   |  The controlled DAE system of an OptimalControlProblemAD
  \*/

  template <typename Derived, Integer N, Integer U, Integer M, Integer NB, Integer NG>
  class OptimalControlProblemAD<Derived, N, U, M, NB, NG>::System : public Controlled<double, N, U, M>
  {
    using Base = Controlled<double, N, U, M>;
    using typename Base::VectorF;
    using typename Base::VectorU;
    using typename Base::VectorH;
    using typename Base::MatrixJF;
    using typename Base::MatrixJU;
    using typename Base::MatrixUU;
    using typename Base::MatrixJH;

    // Layout of the input of the DAE tape: [x, x', u, t]
    static constexpr Integer X{0}, X_DOT{N}, CONTROLS{2*N}, TIME{2*N + U}, INPUTS{2*N + U + 1};

    OptimalControlProblemAD const & m_ocp; /**< The problem owning the recorded model. */

    // The last DAE Jacobian and Hessian, reused by the calls for the different blocks
    mutable VectorX m_jacobian_input, m_hessian_input;
    mutable MatrixX m_jacobian, m_hessian;

  public:
    System(OptimalControlProblemAD const & t_ocp) : Base("ControlledAD"), m_ocp(t_ocp) {}

    VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t) const override
    {return this->dae().value(join(x, x_dot, u, t));}

    MatrixJF JF_x(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t) const override
    {return this->jacobian(x, x_dot, u, t).template middleCols<N>(X);}

    MatrixJF JF_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t) const override
    {return this->jacobian(x, x_dot, u, t).template middleCols<N>(X_DOT);}

    MatrixJU JF_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t) const override
    {return this->jacobian(x, x_dot, u, t).template middleCols<U>(CONTROLS);}

    MatrixJF HF_x_x(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const override
    {return this->hessian(x, x_dot, u, t, l).template block<N, N>(X, X);}

    MatrixJF HF_x_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const override
    {return this->hessian(x, x_dot, u, t, l).template block<N, N>(X, X_DOT);}

    MatrixJF HF_x_dot_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const override
    {return this->hessian(x, x_dot, u, t, l).template block<N, N>(X_DOT, X_DOT);}

    MatrixJU HF_x_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const override
    {return this->hessian(x, x_dot, u, t, l).template block<N, U>(X, CONTROLS);}

    MatrixJU HF_x_dot_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const override
    {return this->hessian(x, x_dot, u, t, l).template block<N, U>(X_DOT, CONTROLS);}

    MatrixUU HF_u_u(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const override
    {return this->hessian(x, x_dot, u, t, l).template block<U, U>(CONTROLS, CONTROLS);}

    VectorH h(VectorF const & x, double const t) const override
    {
      if constexpr (M == 0) {return VectorH::Zero();}
      else {return this->m_ocp.tapes().invariants.value(join(x, t));}
    }

    MatrixJH Jh_x(VectorF const & x, double const t) const override
    {
      if constexpr (M == 0) {return MatrixJH::Zero();}
      else {return this->m_ocp.tapes().invariants.jacobian(join(x, t)).leftCols(N);}
    }

    MatrixJF Hh_x_x(VectorF const & x, double const t, VectorH const & mu) const override
    {
      if constexpr (M == 0) {return MatrixJF::Zero();}
      else {return this->m_ocp.tapes().invariants.hessian(join(x, t), mu).template topLeftCorner<N, N>();}
    }

    MatrixJF Th_x_x(VectorF const & x, double const t, VectorH const & mu, VectorF const & w) const override
    {
      if constexpr (M == 0) {return MatrixJF::Zero();}
      else {return this->m_ocp.tapes().curvature.hessian(join(x, t, mu, w), VectorX::Ones(1)).template topLeftCorner<N, N>();}
    }

  private:
    Tape const & dae() const {return this->m_ocp.tapes().dae;}

    MatrixX const & jacobian(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t) const
    {
      VectorX const input{join(x, x_dot, u, t)};
      if (input.size() != this->m_jacobian_input.size() || input != this->m_jacobian_input) {
        this->m_jacobian = this->dae().jacobian(input);
        this->m_jacobian_input = input;
      }
      return this->m_jacobian;
    }

    MatrixX const & hessian(VectorF const & x, VectorF const & x_dot, VectorU const & u, double const t, VectorF const & l) const
    {
      VectorX const input{join(x, x_dot, u, t, l)};
      if (input.size() != this->m_hessian_input.size() || input != this->m_hessian_input) {
        this->m_hessian = this->dae().hessian(input.head(INPUTS), l);
        this->m_hessian_input = input;
      }
      return this->m_hessian;
    }
  };

} // namespace Sandals

#endif // SANDALS_OPTIMAL_CONTROL_AD_HH
