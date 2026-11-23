/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_OPTIMAL_CONTROL_PROBLEMS_HH
#define TESTS_OPTIMAL_CONTROL_PROBLEMS_HH

#include "Sandals.hh"
#include "Sandals/System/Controlled.hh"
#include "Sandals/OptimalControl.hh"

// Optimal control test problems. All derivatives are exact, except where stated otherwise (to
// exercise the finite differences defaults). Exact derivatives are required to converge IPOPT to
// tight tolerances (1e-12).
namespace TestOCP
{
  using namespace Sandals;
  using Real = double;

  static constexpr Real INFTY{std::numeric_limits<Real>::infinity()};

  /*\
   |  __     __          ____            ____       _
   |  \ \   / /_ _ _ __ |  _ \  ___ _ __|  _ \ ___ | |
   |   \ \ / / _` | '_ \| | | |/ _ \ '__| |_) / _ \| |
   |    \ V / (_| | | | | |_| |  __/ |  |  __/ (_) | |
   |     \_/ \__,_|_| |_|____/ \___|_|  |_|   \___/|_|
   |
  \*/

  // Van der Pol oscillator with additive control:
  //   x0' = x1
  //   x1' = (1 - x0^2) x1 - x0 + u
  class VanDerPol : public ControlledExplicit<Real, 2, 1>
  {
  public:
    VanDerPol() : ControlledExplicit<Real, 2, 1>("VanDerPol") {}

    VectorF f(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      VectorF f;
      f(0) = x(1);
      f(1) = (1.0 - x(0)*x(0))*x(1) - x(0) + u(0);
      return f;
    }

    MatrixJF Jf_x(VectorF const & x, VectorU const & /*u*/, Real /*t*/) const override
    {
      MatrixJF J;
      J(0, 0) = 0.0;                     J(0, 1) = 1.0;
      J(1, 0) = -2.0*x(0)*x(1) - 1.0;    J(1, 1) = 1.0 - x(0)*x(0);
      return J;
    }

    MatrixJU Jf_u(VectorF const & /*x*/, VectorU const & /*u*/, Real /*t*/) const override
    {
      MatrixJU J;
      J(0, 0) = 0.0;
      J(1, 0) = 1.0;
      return J;
    }

    // Only f1 is nonlinear: d2 f1/dx0^2 = -2 x1, d2 f1/dx0 dx1 = -2 x0.
    MatrixJF Hf_x_x(VectorF const & x, VectorU const & /*u*/, Real /*t*/, VectorF const & lambda) const override
    {
      MatrixJF H;
      H(0, 0) = -2.0*lambda(1)*x(1);    H(0, 1) = -2.0*lambda(1)*x(0);
      H(1, 0) = -2.0*lambda(1)*x(0);    H(1, 1) = 0.0;
      return H;
    }

    MatrixJU Hf_x_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJU::Zero();}
    MatrixUU Hf_u_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixUU::Zero();}
  };

  // min  10 |x(T)|^2 + int_0^T (x0^2 + x1^2 + u^2) dt
  // s.t. x(0) = (1, 0), x1^2 + u^2 <= 4, |u| <= 1.5
  class VanDerPolOCP : public OptimalControlProblem<Real, 2, 1>
  {
  public:
    VanDerPolOCP(Real const T = 5.0)
      : OptimalControlProblem<Real, 2, 1>("VanDerPolOCP", std::make_shared<VanDerPol>(), 0.0, T)
    {
      this->u_bounds(VectorU(-1.5), VectorU(1.5));
    }

    // Mayer term on [x_ini; x_end]
    Real mayer(VectorF const & /*x_ini*/, VectorF const & x_end) const override
    {
      return 10.0*x_end.squaredNorm();
    }

    VectorB mayer_gradient(VectorF const & /*x_ini*/, VectorF const & x_end) const override
    {
      VectorB gradient(VectorB::Zero());
      gradient.tail<2>() = 20.0*x_end;
      return gradient;
    }

    MatrixB mayer_hessian(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
    {
      MatrixB hessian(MatrixB::Zero());
      hessian(2, 2) = 20.0;
      hessian(3, 3) = 20.0;
      return hessian;
    }

    // Lagrange term on [x; u]
    Real lagrange(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return x.squaredNorm() + u.squaredNorm();
    }

    VectorP lagrange_gradient(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      VectorP gradient;
      gradient << 2.0*x, 2.0*u;
      return gradient;
    }

    MatrixP lagrange_hessian(VectorF const & /*x*/, VectorU const & /*u*/, Real /*t*/) const override
    {
      return 2.0*MatrixP::Identity();
    }

    // Boundary conditions x(0) = (1, 0)
    Integer boundary_conditions_number() const override {return 2;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & /*x_end*/) const override
    {
      return x_ini - VectorF(1.0, 0.0);
    }

    MatrixX boundary_conditions_jacobian(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
    {
      MatrixX jacobian(MatrixX::Zero(2, 4));
      jacobian(0, 0) = 1.0;
      jacobian(1, 1) = 1.0;
      return jacobian;
    }

    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override
    {
      return MatrixB::Zero();
    }

    // Path constraint x1^2 + u^2 <= 4
    Integer path_constraints_number() const override {return 1;}

    VectorX path_constraints(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return VectorX::Constant(1, x(1)*x(1) + u(0)*u(0));
    }

    VectorX path_constraints_upper() const override {return VectorX::Constant(1, 4.0);}

    MatrixX path_constraints_jacobian(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      MatrixX jacobian(1, 3);
      jacobian << 0.0, 2.0*x(1), 2.0*u(0);
      return jacobian;
    }

    MatrixP path_constraints_hessian(VectorF const &, VectorU const &, Real, VectorX const & lambda) const override
    {
      MatrixP hessian(MatrixP::Zero());
      hessian(1, 1) = 2.0*lambda(0);
      hessian(2, 2) = 2.0*lambda(0);
      return hessian;
    }
  };

  // Same problem, defined by the functions values only: every derivative of the problem uses the
  // finite differences defaults (the control bounds are left to the caller).
  class VanDerPolOCPValuesOnly : public OptimalControlProblem<Real, 2, 1>
  {
  public:
    VanDerPolOCPValuesOnly(Real const T = 5.0)
      : OptimalControlProblem<Real, 2, 1>("VanDerPolOCPValuesOnly", std::make_shared<VanDerPol>(), 0.0, T)
    {}

    Real mayer(VectorF const & /*x_ini*/, VectorF const & x_end) const override
    {
      return 10.0*x_end.squaredNorm();
    }

    Real lagrange(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return x.squaredNorm() + u.squaredNorm();
    }

    Integer boundary_conditions_number() const override {return 2;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & /*x_end*/) const override
    {
      return x_ini - VectorF(1.0, 0.0);
    }

    Integer path_constraints_number() const override {return 1;}

    VectorX path_constraints(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return VectorX::Constant(1, x(1)*x(1) + u(0)*u(0));
    }

    VectorX path_constraints_upper() const override {return VectorX::Constant(1, 4.0);}
  };

  /*\
   |   ___           _                _   ____    _    _____
   |  |_ _|_ __   __| | _____  __    / | |  _ \  / \  | ____|
   |   | || '_ \ / _` |/ _ \ \/ /____| | | | | |/ _ \ |  _|
   |   | || | | | (_| |  __/>  <_____| | | |_| / ___ \| |___
   |  |___|_| |_|\__,_|\___/_/\_\    |_| |____/_/   \_\_____|
   |
  \*/

  // Index-1 DAE with differential state x and algebraic state z:
  //   x' = -x + z + u
  //   0  = z - x^2
  // The Hessians are left to the finite differences defaults on purpose.
  class SemiLinearDAE : public Controlled<Real, 2, 1>
  {
  public:
    SemiLinearDAE() : Controlled<Real, 2, 1>("SemiLinearDAE") {}

    VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real /*t*/) const override
    {
      Real const state{x(0)}, algebraic{x(1)};
      VectorF F;
      F(0) = x_dot(0) + state - algebraic - u(0);
      F(1) = algebraic - state*state;
      return F;
    }

    MatrixJF JF_x(VectorF const & x, VectorF const & /*x_dot*/, VectorU const & /*u*/, Real /*t*/) const override
    {
      MatrixJF J;
      J(0, 0) = 1.0;          J(0, 1) = -1.0;
      J(1, 0) = -2.0*x(0);    J(1, 1) = 1.0;
      return J;
    }

    MatrixJF JF_x_dot(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(0, 0) = 1.0; // the algebraic state has no derivative
      return J;
    }

    MatrixJU JF_u(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {
      MatrixJU J;
      J(0, 0) = -1.0;
      J(1, 0) = 0.0;
      return J;
    }
  };

  // min  int_0^T (x^2 + u^2) dt
  // s.t. x(0) = x0, z(0) = x(0)^2
  // The initial conditions must be consistent: the algebraic state is fixed by z(0) = x(0)^2,
  // otherwise z(0) is a free direction of the NLP.
  class SemiLinearDAEOCP : public OptimalControlProblem<Real, 2, 1>
  {
    Real m_x0;

  public:
    SemiLinearDAEOCP(Real const x0 = 0.5, Real const T = 2.0)
      : OptimalControlProblem<Real, 2, 1>("SemiLinearDAEOCP", std::make_shared<SemiLinearDAE>(), 0.0, T),
        m_x0(x0)
    {}

    Real lagrange(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return x(0)*x(0) + u(0)*u(0);
    }

    VectorP lagrange_gradient(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return VectorP(2.0*x(0), 0.0, 2.0*u(0));
    }

    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {
      MatrixP hessian(MatrixP::Zero());
      hessian(0, 0) = 2.0;
      hessian(2, 2) = 2.0;
      return hessian;
    }

    Integer boundary_conditions_number() const override {return 2;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & /*x_end*/) const override
    {
      VectorX b(2);
      b(0) = x_ini(0) - this->m_x0;
      b(1) = x_ini(1) - x_ini(0)*x_ini(0);
      return b;
    }

    MatrixX boundary_conditions_jacobian(VectorF const & x_ini, VectorF const & /*x_end*/) const override
    {
      MatrixX jacobian(MatrixX::Zero(2, 4));
      jacobian(0, 0) = 1.0;
      jacobian(1, 0) = -2.0*x_ini(0);
      jacobian(1, 1) = 1.0;
      return jacobian;
    }

    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const & lambda) const override
    {
      MatrixB hessian(MatrixB::Zero());
      hessian(0, 0) = -2.0*lambda(1);
      return hessian;
    }
  };

  // Equivalent ODE, obtained by eliminating z = x^2:
  //   x' = -x + x^2 + u
  class SemiLinearODE : public ControlledExplicit<Real, 1, 1>
  {
  public:
    SemiLinearODE() : ControlledExplicit<Real, 1, 1>("SemiLinearODE") {}

    VectorF f(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return VectorF(-x(0) + x(0)*x(0) + u(0));
    }

    MatrixJF Jf_x(VectorF const & x, VectorU const & /*u*/, Real /*t*/) const override
    {
      return MatrixJF(-1.0 + 2.0*x(0));
    }

    MatrixJU Jf_u(VectorF const &, VectorU const &, Real) const override {return MatrixJU(1.0);}

    MatrixJF Hf_x_x(VectorF const &, VectorU const &, Real, VectorF const & lambda) const override
    {
      return MatrixJF(2.0*lambda(0));
    }

    MatrixJU Hf_x_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJU::Zero();}
    MatrixUU Hf_u_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixUU::Zero();}
  };

  class SemiLinearODEOCP : public OptimalControlProblem<Real, 1, 1>
  {
    Real m_x0;

  public:
    SemiLinearODEOCP(Real const x0 = 0.5, Real const T = 2.0)
      : OptimalControlProblem<Real, 1, 1>("SemiLinearODEOCP", std::make_shared<SemiLinearODE>(), 0.0, T),
        m_x0(x0)
    {}

    Real lagrange(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return x(0)*x(0) + u(0)*u(0);
    }

    VectorP lagrange_gradient(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return VectorP(2.0*x(0), 2.0*u(0));
    }

    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {
      return 2.0*MatrixP::Identity();
    }

    Integer boundary_conditions_number() const override {return 1;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & /*x_end*/) const override
    {
      return VectorX::Constant(1, x_ini(0) - this->m_x0);
    }

    MatrixX boundary_conditions_jacobian(VectorF const &, VectorF const &) const override
    {
      MatrixX jacobian(MatrixX::Zero(1, 2));
      jacobian(0, 0) = 1.0;
      return jacobian;
    }

    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override
    {
      return MatrixB::Zero();
    }
  };

  /*\
   |   ____              _     _        ___       _                       _
   |  |  _ \  ___  _   _| |__ | | ___  |_ _|_ __ | |_ ___  __ _ _ __ __ _| |_ ___  _ __
   |  | | | |/ _ \| | | | '_ \| |/ _ \  | || '_ \| __/ _ \/ _` | '__/ _` | __/ _ \| '__|
   |  | |_| | (_) | |_| | |_) | |  __/  | || | | | ||  __/ (_| | | | (_| | || (_) | |
   |  |____/ \___/ \__,_|_.__/|_|\___| |___|_| |_|\__\___|\__, |_|  \__,_|\__\___/|_|
   |                                                      |___/
  \*/

  // Double integrator:
  //   x0' = x1
  //   x1' = u
  class DoubleIntegrator : public ControlledExplicit<Real, 2, 1>
  {
  public:
    DoubleIntegrator() : ControlledExplicit<Real, 2, 1>("DoubleIntegrator") {}

    VectorF f(VectorF const & x, VectorU const & u, Real /*t*/) const override
    {
      return VectorF(x(1), u(0));
    }

    MatrixJF Jf_x(VectorF const &, VectorU const &, Real) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(0, 1) = 1.0;
      return J;
    }

    MatrixJU Jf_u(VectorF const &, VectorU const &, Real) const override
    {
      MatrixJU J;
      J(0, 0) = 0.0;
      J(1, 0) = 1.0;
      return J;
    }

    MatrixJF Hf_x_x(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJF::Zero();}
    MatrixJU Hf_x_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJU::Zero();}
    MatrixUU Hf_u_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixUU::Zero();}
  };

  // min  int_0^1 u^2/2 dt,  s.t. x(0) = x_ini, x(1) = x_end
  //  - rest-to-rest x(0) = (0, 0), x(1) = (1, 0): optimal u = 6 - 12 t, J* = 6;
  //  - Bryson-Denham x(0) = (0, 1), x(1) = (0, -1) with x0 <= 1/9: J* = 4.
  class DoubleIntegratorOCP : public OptimalControlProblem<Real, 2, 1>
  {
    VectorF m_x_ini;
    VectorF m_x_end;

  public:
    DoubleIntegratorOCP(VectorF const & x_ini, VectorF const & x_end)
      : OptimalControlProblem<Real, 2, 1>("DoubleIntegratorOCP", std::make_shared<DoubleIntegrator>(), 0.0, 1.0),
        m_x_ini(x_ini), m_x_end(x_end)
    {}

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
      b << x_ini - this->m_x_ini, x_end - this->m_x_end;
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

#endif // TESTS_OPTIMAL_CONTROL_PROBLEMS_HH
