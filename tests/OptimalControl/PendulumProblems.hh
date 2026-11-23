/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PENDULUM_PROBLEMS_HH
#define TESTS_PENDULUM_PROBLEMS_HH

#include "Sandals.hh"
#include "Sandals/System/Controlled.hh"
#include "Sandals/OptimalControl.hh"

// Swing-up of a pendulum driven by a horizontal force u on the bob, written in Cartesian coordinates
// with different differentiation indices, and in the angle (reference ODE).
//
// The index-3 formulation, with states (x, y, vx, vy, lambda), is
//   x'    = vx
//   y'    = vy
//   m vx' = -2 x lambda + u
//   m vy' = -2 y lambda - m g
//   0     = x^2 + y^2 - L^2                    (position level)
// Differentiating the constraint once or twice lowers the index:
//   index 2:  0 = x vx + y vy                  (velocity level)
//   index 1:  0 = vx^2 + vy^2 + x vx' + y vy'  (acceleration level)
//   index 0:  lambda eliminated, states (x, y, vx, vy)
// The constraints lost by the index reduction are the invariants of the reduced formulations:
//   h0 = x^2 + y^2 - L^2 (index 2, 1, 0)
//   h1 = x vx + y vy     (index 1, 0)
// Each reduced formulation is templated on the number of invariants M, so that it can be solved
// with (M > 0) or without (M = 0) the projection on the invariants.
namespace TestOCP
{
  using namespace Sandals;
  using Real = double;

  /*\
   |  Shared data
  \*/

  // Indices of the Cartesian states
  enum PendulumState : Integer {X = 0, Y = 1, VX = 2, VY = 3, LAMBDA = 4};

  // Physical parameters
  struct PendulumParameters
  {
    Real m{1.0};  // Mass (kg)
    Real L{1.0};  // Length (m)
    Real g{9.81}; // Gravity (m/s^2)

    // Hanging at rest, with the static value of lambda (when the state includes it)
    template <Integer N>
    Eigen::Vector<Real, N> hanging_at_rest() const
    {
      Eigen::Vector<Real, N> x(Eigen::Vector<Real, N>::Zero());
      x(Y) = -L;
      if constexpr (N > LAMBDA) {x(LAMBDA) = m*g/(2.0*L);}
      return x;
    }
  };

  // The first M of the invariants (h0 = x^2 + y^2 - L^2, h1 = x vx + y vy) on N states. They are
  // quadratic in the states, so their third derivatives vanish.
  template <Integer N, Integer M>
  struct PendulumInvariants
  {
    static_assert(M >= 0 && M <= 2, "PendulumInvariants: at most two invariants.");
    using VectorF  = Eigen::Vector<Real, N>;
    using VectorH  = Eigen::Vector<Real, M>;
    using MatrixJF = Eigen::Matrix<Real, N, N>;
    using MatrixJH = Eigen::Matrix<Real, M, N>;

    static VectorH h(VectorF const & x, Real const L)
    {
      VectorH h;
      if constexpr (M > 0) {h(0) = x(X)*x(X) + x(Y)*x(Y) - L*L;}
      if constexpr (M > 1) {h(1) = x(X)*x(VX) + x(Y)*x(VY);}
      return h;
    }

    static MatrixJH Jh_x(VectorF const & x)
    {
      MatrixJH J(MatrixJH::Zero());
      if constexpr (M > 0) {
        J(0, X) = 2.0*x(X);
        J(0, Y) = 2.0*x(Y);
      }
      if constexpr (M > 1) {
        J(1, X)  = x(VX);
        J(1, Y)  = x(VY);
        J(1, VX) = x(X);
        J(1, VY) = x(Y);
      }
      return J;
    }

    static MatrixJF Hh_x_x(VectorH const & mu)
    {
      MatrixJF H(MatrixJF::Zero());
      if constexpr (M > 0) {
        H(X, X) = 2.0*mu(0);
        H(Y, Y) = 2.0*mu(0);
      }
      if constexpr (M > 1) {
        H(X, VX) = H(VX, X) = mu(1);
        H(Y, VY) = H(VY, Y) = mu(1);
      }
      return H;
    }
  };

  // Implicit pendulum formulations share the four dynamics equations (index 3, 2, 1) and differ in
  // the last one, the constraint. The derived classes implement the constraint and its derivatives.
  template <Integer M>
  class PendulumImplicitBase : public Controlled<Real, 5, 1, M>, public PendulumParameters
  {
  public:
    using Base = Controlled<Real, 5, 1, M>;
    using typename Base::VectorF;
    using typename Base::VectorU;
    using typename Base::VectorH;
    using typename Base::MatrixJF;
    using typename Base::MatrixJU;
    using typename Base::MatrixUU;
    using typename Base::MatrixJH;
    using Invariants = PendulumInvariants<5, M>;

    PendulumImplicitBase(std::string const & name) : Base(name) {}

    // Constraint equation (last row of F) and its derivatives
    virtual Real constraint(VectorF const & x, VectorF const & x_dot) const = 0;
    virtual VectorF constraint_x(VectorF const & x, VectorF const & x_dot) const = 0;
    virtual VectorF constraint_x_dot(VectorF const & x, VectorF const & x_dot) const = 0;
    virtual MatrixJF constraint_x_x() const = 0;     // constant for all the formulations
    virtual MatrixJF constraint_x_x_dot() const = 0; // constant for all the formulations

    VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real /*t*/) const override
    {
      VectorF F;
      F(0) = x_dot(X) - x(VX);
      F(1) = x_dot(Y) - x(VY);
      F(2) = m*x_dot(VX) + 2.0*x(X)*x(LAMBDA) - u(0);
      F(3) = m*x_dot(VY) + 2.0*x(Y)*x(LAMBDA) + m*g;
      F(4) = this->constraint(x, x_dot);
      return F;
    }

    MatrixJF JF_x(VectorF const & x, VectorF const & x_dot, VectorU const &, Real) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(0, VX)     = -1.0;
      J(1, VY)     = -1.0;
      J(2, X)      = 2.0*x(LAMBDA);
      J(2, LAMBDA) = 2.0*x(X);
      J(3, Y)      = 2.0*x(LAMBDA);
      J(3, LAMBDA) = 2.0*x(Y);
      J.row(4)     = this->constraint_x(x, x_dot).transpose();
      return J;
    }

    MatrixJF JF_x_dot(VectorF const & x, VectorF const & x_dot, VectorU const &, Real) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(0, X)  = 1.0;
      J(1, Y)  = 1.0;
      J(2, VX) = m;
      J(3, VY) = m;
      J.row(4) = this->constraint_x_dot(x, x_dot).transpose();
      return J;
    }

    MatrixJU JF_u(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {
      MatrixJU J(MatrixJU::Zero());
      J(2, 0) = -1.0;
      return J;
    }

    // Nonlinear terms: 2 x lambda, 2 y lambda (equations 2, 3) and the constraint (equation 4)
    MatrixJF HF_x_x(VectorF const &, VectorF const &, VectorU const &, Real, VectorF const & w) const override
    {
      MatrixJF H(w(4)*this->constraint_x_x());
      H(X, LAMBDA) += 2.0*w(2);
      H(LAMBDA, X) += 2.0*w(2);
      H(Y, LAMBDA) += 2.0*w(3);
      H(LAMBDA, Y) += 2.0*w(3);
      return H;
    }

    MatrixJF HF_x_x_dot(VectorF const &, VectorF const &, VectorU const &, Real, VectorF const & w) const override
    {return w(4)*this->constraint_x_x_dot();}

    MatrixJF HF_x_dot_x_dot(VectorF const &, VectorF const &, VectorU const &, Real, VectorF const &) const override
    {return MatrixJF::Zero();}

    MatrixJU HF_x_u(VectorF const &, VectorF const &, VectorU const &, Real, VectorF const &) const override
    {return MatrixJU::Zero();}

    MatrixJU HF_x_dot_u(VectorF const &, VectorF const &, VectorU const &, Real, VectorF const &) const override
    {return MatrixJU::Zero();}

    MatrixUU HF_u_u(VectorF const &, VectorF const &, VectorU const &, Real, VectorF const &) const override
    {return MatrixUU::Zero();}

    // Invariants
    VectorH h(VectorF const & x, Real) const override {return Invariants::h(x, L);}
    MatrixJH Jh_x(VectorF const & x, Real) const override {return Invariants::Jh_x(x);}
    MatrixJF Hh_x_x(VectorF const &, Real, VectorH const & mu) const override {return Invariants::Hh_x_x(mu);}
    MatrixJF Th_x_x(VectorF const &, Real, VectorH const &, VectorF const &) const override {return MatrixJF::Zero();}
  };

  /*\
   |  Index 3: 0 = x^2 + y^2 - L^2
  \*/

  class PendulumIndex3 : public PendulumImplicitBase<0>
  {
  public:
    PendulumIndex3() : PendulumImplicitBase<0>("PendulumIndex3") {}

    Real constraint(VectorF const & x, VectorF const &) const override
    {return x(X)*x(X) + x(Y)*x(Y) - L*L;}

    VectorF constraint_x(VectorF const & x, VectorF const &) const override
    {
      VectorF d(VectorF::Zero());
      d(X) = 2.0*x(X);
      d(Y) = 2.0*x(Y);
      return d;
    }

    VectorF constraint_x_dot(VectorF const &, VectorF const &) const override {return VectorF::Zero();}

    MatrixJF constraint_x_x() const override
    {
      MatrixJF H(MatrixJF::Zero());
      H(X, X) = 2.0;
      H(Y, Y) = 2.0;
      return H;
    }

    MatrixJF constraint_x_x_dot() const override {return MatrixJF::Zero();}
  };

  /*\
   |  Index 2: 0 = x vx + y vy, invariant h0
  \*/

  template <Integer M>
  class PendulumIndex2 : public PendulumImplicitBase<M>
  {
    static_assert(M == 0 || M == 1, "PendulumIndex2: the only invariant is the position constraint.");

  public:
    using typename PendulumImplicitBase<M>::VectorF;
    using typename PendulumImplicitBase<M>::MatrixJF;

    PendulumIndex2() : PendulumImplicitBase<M>("PendulumIndex2") {}

    Real constraint(VectorF const & x, VectorF const &) const override
    {return x(X)*x(VX) + x(Y)*x(VY);}

    VectorF constraint_x(VectorF const & x, VectorF const &) const override
    {
      VectorF d(VectorF::Zero());
      d(X)  = x(VX);
      d(Y)  = x(VY);
      d(VX) = x(X);
      d(VY) = x(Y);
      return d;
    }

    VectorF constraint_x_dot(VectorF const &, VectorF const &) const override {return VectorF::Zero();}

    MatrixJF constraint_x_x() const override
    {
      MatrixJF H(MatrixJF::Zero());
      H(X, VX) = H(VX, X) = 1.0;
      H(Y, VY) = H(VY, Y) = 1.0;
      return H;
    }

    MatrixJF constraint_x_x_dot() const override {return MatrixJF::Zero();}
  };

  /*\
   |  Index 1: 0 = vx^2 + vy^2 + x vx' + y vy', invariants h0 and h1
  \*/

  template <Integer M>
  class PendulumIndex1 : public PendulumImplicitBase<M>
  {
  public:
    using typename PendulumImplicitBase<M>::VectorF;
    using typename PendulumImplicitBase<M>::MatrixJF;

    PendulumIndex1() : PendulumImplicitBase<M>("PendulumIndex1") {}

    Real constraint(VectorF const & x, VectorF const & x_dot) const override
    {return x(VX)*x(VX) + x(VY)*x(VY) + x(X)*x_dot(VX) + x(Y)*x_dot(VY);}

    VectorF constraint_x(VectorF const & x, VectorF const & x_dot) const override
    {
      VectorF d(VectorF::Zero());
      d(X)  = x_dot(VX);
      d(Y)  = x_dot(VY);
      d(VX) = 2.0*x(VX);
      d(VY) = 2.0*x(VY);
      return d;
    }

    VectorF constraint_x_dot(VectorF const & x, VectorF const &) const override
    {
      VectorF d(VectorF::Zero());
      d(VX) = x(X);
      d(VY) = x(Y);
      return d;
    }

    MatrixJF constraint_x_x() const override
    {
      MatrixJF H(MatrixJF::Zero());
      H(VX, VX) = 2.0;
      H(VY, VY) = 2.0;
      return H;
    }

    // Entry (i, j) is the derivative with respect to x_i and x'_j
    MatrixJF constraint_x_x_dot() const override
    {
      MatrixJF H(MatrixJF::Zero());
      H(X, VX) = 1.0;
      H(Y, VY) = 1.0;
      return H;
    }
  };

  /*\
   |  Index 0: lambda eliminated, invariants h0 and h1
   |    lambda = (m |v|^2 + x u - m g y) / (2 (x^2 + y^2))
   |    vx'    = (-2 x lambda + u) / m
   |    vy'    = (-2 y lambda - m g) / m
   |  Exact Jacobians, finite differences Hessians (defaults).
  \*/

  template <Integer M>
  class PendulumIndex0 : public ControlledExplicit<Real, 4, 1, M>, public PendulumParameters
  {
  public:
    using Base = ControlledExplicit<Real, 4, 1, M>;
    using typename Base::VectorF;
    using typename Base::VectorU;
    using typename Base::VectorH;
    using typename Base::MatrixJF;
    using typename Base::MatrixJU;
    using typename Base::MatrixJH;
    using Invariants = PendulumInvariants<4, M>;

    PendulumIndex0() : Base("PendulumIndex0") {}

    // The multiplier lambda as a function of the states and the control
    Real lambda(VectorF const & x, VectorU const & u) const
    {
      Real const r2{x(X)*x(X) + x(Y)*x(Y)};
      Real const v2{x(VX)*x(VX) + x(VY)*x(VY)};
      return (m*v2 + x(X)*u(0) - m*g*x(Y)) / (2.0*r2);
    }

    // Gradient of lambda with respect to [x, y, vx, vy, u]
    Eigen::Vector<Real, 5> lambda_gradient(VectorF const & x, VectorU const & u) const
    {
      Real const r2{x(X)*x(X) + x(Y)*x(Y)};
      Real const l{this->lambda(x, u)};
      Eigen::Vector<Real, 5> d;
      d(X)  = (u(0) - 4.0*l*x(X)) / (2.0*r2);
      d(Y)  = (-m*g - 4.0*l*x(Y)) / (2.0*r2);
      d(VX) = m*x(VX) / r2;
      d(VY) = m*x(VY) / r2;
      d(4)  = x(X) / (2.0*r2);
      return d;
    }

    VectorF f(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      Real const l{this->lambda(x, u)};
      VectorF f;
      f(X)  = x(VX);
      f(Y)  = x(VY);
      f(VX) = (-2.0*x(X)*l + u(0)) / m;
      f(VY) = (-2.0*x(Y)*l - m*g) / m;
      return f;
    }

    MatrixJF Jf_x(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      Real const l{this->lambda(x, u)};
      Eigen::Vector<Real, 4> const dl{this->lambda_gradient(x, u).template head<4>()};
      MatrixJF J(MatrixJF::Zero());
      J(X, VX)  = 1.0;
      J(Y, VY)  = 1.0;
      J.row(VX) = -2.0*x(X)*dl.transpose() / m;
      J.row(VY) = -2.0*x(Y)*dl.transpose() / m;
      J(VX, X) -= 2.0*l / m;
      J(VY, Y) -= 2.0*l / m;
      return J;
    }

    MatrixJU Jf_u(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      Real const dl_du{this->lambda_gradient(x, u)(4)};
      MatrixJU J(MatrixJU::Zero());
      J(VX, 0) = (1.0 - 2.0*x(X)*dl_du) / m;
      J(VY, 0) = -2.0*x(Y)*dl_du / m;
      return J;
    }

    // Invariants
    VectorH h(VectorF const & x, Real) const override {return Invariants::h(x, L);}
    MatrixJH Jh_x(VectorF const & x, Real) const override {return Invariants::Jh_x(x);}
    MatrixJF Hh_x_x(VectorF const &, Real, VectorH const & mu) const override {return Invariants::Hh_x_x(mu);}
    MatrixJF Th_x_x(VectorF const &, Real, VectorH const &, VectorF const &) const override {return MatrixJF::Zero();}
  };

  /*\
   |  Swing-up problem, common to all the Cartesian formulations
  \*/

  // From hanging at rest to upright at rest in time T:
  //   min  int_0^T u^2/2 dt,  s.t. |u| <= u_max
  // Boundary conditions:
  //  - all the initial states are fixed, since the first mesh node is neither a collocation point
  //    nor projected on the invariants; when lambda is a state, its node value only shifts the stage
  //    variables (lambda' never appears), so it is a gauge, fixed to its static value;
  //  - at the final node the constraint (or the invariants) fix the radial components, so only the
  //    tangential ones x(T) = 0 and vx(T) = 0 are imposed as equalities (imposing more would violate
  //    the LICQ). These also hold hanging at rest, so the inequality y(T) >= 0 selects the upright
  //    position (it is inactive at the solution, so it does not affect the LICQ).
  template <typename System>
  class PendulumSwingUpOCP : public OptimalControlProblem<Real, System::equations_number(), 1, System::invariants_number()>
  {
  public:
    static constexpr Integer N{System::equations_number()};
    static constexpr Integer M{System::invariants_number()};
    using Base = OptimalControlProblem<Real, N, 1, M>;
    using typename Base::VectorF;
    using typename Base::VectorU;
    using typename Base::VectorP;
    using typename Base::MatrixP;
    using typename Base::MatrixB;
    using typename Base::VectorX;
    using typename Base::MatrixX;

  private:
    VectorF m_x_ini; // Hanging at rest

  public:
    PendulumSwingUpOCP(Real const T, Real const u_max = std::numeric_limits<Real>::infinity())
      : Base(std::string("SwingUp") + System().name(), std::make_shared<System>(), 0.0, T)
    {
      this->u_bounds(VectorU(-u_max), VectorU(u_max));
      this->m_x_ini = PendulumParameters().template hanging_at_rest<N>();
    }

    Real lagrange(VectorF const &, VectorU const & u, Real) const override
    {
      return 0.5*u(0)*u(0);
    }

    VectorP lagrange_gradient(VectorF const &, VectorU const & u, Real) const override
    {
      VectorP gradient(VectorP::Zero());
      gradient(N) = u(0);
      return gradient;
    }

    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {
      MatrixP hessian(MatrixP::Zero());
      hessian(N, N) = 1.0;
      return hessian;
    }

    Integer boundary_conditions_number() const override {return N + 3;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & x_end) const override
    {
      VectorX b(N + 3);
      b.head(N) = x_ini - this->m_x_ini;
      b(N)      = x_end(X);
      b(N + 1)  = x_end(VX);
      b(N + 2)  = x_end(Y);
      return b;
    }

    // Equalities, except y(T) >= 0
    VectorX boundary_conditions_upper() const override
    {
      VectorX upper(VectorX::Zero(N + 3));
      upper(N + 2) = std::numeric_limits<Real>::infinity();
      return upper;
    }

    MatrixX boundary_conditions_jacobian(VectorF const &, VectorF const &) const override
    {
      MatrixX jacobian(MatrixX::Zero(N + 3, 2*N)); // columns: [x_ini, x_end]
      jacobian.topLeftCorner(N, N).setIdentity();
      jacobian(N, N + X)      = 1.0;
      jacobian(N + 1, N + VX) = 1.0;
      jacobian(N + 2, N + Y)  = 1.0;
      return jacobian;
    }

    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override
    {
      return MatrixB::Zero();
    }
  };

  using PendulumIndex3OCP = PendulumSwingUpOCP<PendulumIndex3>;

  /*\
   |  Reference: pendulum in the angle theta from the downward vertical (x = L sin(theta),
   |  y = -L cos(theta)). The horizontal force has tangential component u cos(theta):
   |    theta' = omega
   |    omega' = (-g sin(theta) + u cos(theta) / m) / L
  \*/

  class PendulumAngle : public ControlledExplicit<Real, 2, 1>, public PendulumParameters
  {
  public:
    PendulumAngle() : ControlledExplicit<Real, 2, 1>("PendulumAngle") {}

    VectorF f(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      Real const theta{x(0)}, omega{x(1)};
      VectorF f;
      f(0) = omega;
      f(1) = (-g*std::sin(theta) + u(0)*std::cos(theta)/m)/L;
      return f;
    }

    MatrixJF Jf_x(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      Real const theta{x(0)};
      MatrixJF J(MatrixJF::Zero());
      J(0, 1) = 1.0;
      J(1, 0) = (-g*std::cos(theta) - u(0)*std::sin(theta)/m)/L;
      return J;
    }

    MatrixJU Jf_u(VectorF const & x, VectorU const & /*u*/, Real /*t*/) const override
    {
      Real const theta{x(0)};
      MatrixJU J;
      J(0, 0) = 0.0;
      J(1, 0) = std::cos(theta)/(m*L);
      return J;
    }

    // Only omega' is nonlinear, in theta and u
    MatrixJF Hf_x_x(VectorF const & x, VectorU const & u, Real /*t*/, VectorF const & lambda) const override
    {
      Real const theta{x(0)};
      MatrixJF H(MatrixJF::Zero());
      H(0, 0) = lambda(1)*(g*std::sin(theta) - u(0)*std::cos(theta)/m)/L;
      return H;
    }

    MatrixJU Hf_x_u(VectorF const & x, VectorU const & /*u*/, Real /*t*/, VectorF const & lambda) const override
    {
      Real const theta{x(0)};
      MatrixJU H(MatrixJU::Zero());
      H(0, 0) = -lambda(1)*std::sin(theta)/(m*L);
      return H;
    }

    MatrixUU Hf_u_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixUU::Zero();}
  };

  // Swing-up in the angle: theta(0) = 0, omega(0) = 0, theta(T) = pi, omega(T) = 0.
  class PendulumAngleOCP : public OptimalControlProblem<Real, 2, 1>
  {
  public:
    PendulumAngleOCP(Real const T, Real const u_max = std::numeric_limits<Real>::infinity())
      : OptimalControlProblem<Real, 2, 1>("PendulumAngleOCP", std::make_shared<PendulumAngle>(), 0.0, T)
    {
      this->u_bounds(VectorU(-u_max), VectorU(u_max));
    }

    Real lagrange(VectorF const &, VectorU const & u, Real) const override
    {
      return 0.5*u(0)*u(0);
    }

    VectorP lagrange_gradient(VectorF const &, VectorU const & u, Real) const override
    {
      return VectorP(0.0, 0.0, u(0));
    }

    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {
      MatrixP hessian(MatrixP::Zero());
      hessian(2, 2) = 1.0;
      return hessian;
    }

    Integer boundary_conditions_number() const override {return 4;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & x_end) const override
    {
      VectorX b(4);
      b(0) = x_ini(0);
      b(1) = x_ini(1);
      b(2) = x_end(0) - M_PI;
      b(3) = x_end(1);
      return b;
    }

    MatrixX boundary_conditions_jacobian(VectorF const &, VectorF const &) const override
    {
      return MatrixX::Identity(4, 4);
    }

    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override
    {
      return MatrixB::Zero();
    }
  };

} // namespace TestOCP

#endif // TESTS_PENDULUM_PROBLEMS_HH
