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

namespace TestOCP
{
  using namespace Sandals;
  using Real = double;

  /*\
   |  Van der Pol oscillator: x0' = x1, x1' = (1 - x0^2) x1 - x0 + u (exact derivatives)
  \*/

  class VanDerPol : public ControlledExplicit<Real, 2, 1>
  {
  public:
    VanDerPol() : ControlledExplicit<Real, 2, 1>("VanDerPol") {}
    VectorF f(VectorF const & x, VectorU const & u, Real) const override
    {return VectorF(x(1), (1.0 - x(0)*x(0))*x(1) - x(0) + u(0));}
    MatrixJF Jf_x(VectorF const & x, VectorU const &, Real) const override
    {MatrixJF J; J << 0.0, 1.0, -2.0*x(0)*x(1) - 1.0, 1.0 - x(0)*x(0); return J;}
    MatrixJU Jf_u(VectorF const &, VectorU const &, Real) const override
    {return MatrixJU(0.0, 1.0);}
    MatrixJF Hf_x_x(VectorF const & x, VectorU const &, Real, VectorF const & l) const override
    {MatrixJF H; H << -2.0*l(1)*x(1), -2.0*l(1)*x(0), -2.0*l(1)*x(0), 0.0; return H;}
    MatrixJU Hf_x_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJU::Zero();}
    MatrixUU Hf_u_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixUU::Zero();}
  };

  // min 10 |x(T)|^2 + int (x0^2 + x1^2 + u^2) dt, x(0) = (1, 0), x1^2 + u^2 <= 4, |u| <= 1.5
  class VanDerPolOCP : public OptimalControlProblem<Real, 2, 1>
  {
  public:
    VanDerPolOCP(Real const T = 5.0) : OptimalControlProblem<Real, 2, 1>("VanDerPolOCP", std::make_shared<VanDerPol>(), 0.0, T)
    {this->u_bounds(VectorU(-1.5), VectorU(1.5));}

    Real mayer(VectorF const &, VectorF const & x_end) const override {return 10.0*x_end.squaredNorm();}
    VectorB mayer_gradient(VectorF const &, VectorF const & x_end) const override
    {VectorB g(VectorB::Zero()); g.tail<2>() = 20.0*x_end; return g;}
    MatrixB mayer_hessian(VectorF const &, VectorF const &) const override
    {MatrixB H(MatrixB::Zero()); H.bottomRightCorner<2, 2>().diagonal().setConstant(20.0); return H;}

    Real lagrange(VectorF const & x, VectorU const & u, Real) const override {return x.squaredNorm() + u.squaredNorm();}
    VectorP lagrange_gradient(VectorF const & x, VectorU const & u, Real) const override
    {VectorP g; g << 2.0*x, 2.0*u; return g;}
    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override {return 2.0*MatrixP::Identity();}

    Integer boundary_conditions_number() const override {return 2;}
    VectorX boundary_conditions(VectorF const & x_ini, VectorF const &) const override
    {return x_ini - VectorF(1.0, 0.0);}
    MatrixX boundary_conditions_jacobian(VectorF const &, VectorF const &) const override
    {MatrixX J(MatrixX::Zero(2, 4)); J.leftCols(2).setIdentity(); return J;}
    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override
    {return MatrixB::Zero();}

    Integer path_constraints_number() const override {return 1;}
    VectorX path_constraints(VectorF const & x, VectorU const & u, Real) const override
    {return VectorX::Constant(1, x(1)*x(1) + u(0)*u(0));}
    VectorX path_constraints_upper() const override {return VectorX::Constant(1, 4.0);}
    MatrixX path_constraints_jacobian(VectorF const & x, VectorU const & u, Real) const override
    {MatrixX J(1, 3); J << 0.0, 2.0*x(1), 2.0*u(0); return J;}
    MatrixP path_constraints_hessian(VectorF const &, VectorU const &, Real, VectorX const & l) const override
    {MatrixP H(MatrixP::Zero()); H(1, 1) = 2.0*l(0); H(2, 2) = 2.0*l(0); return H;}
  };

  // Same problem with values only: every derivative uses the finite differences defaults.
  class VanDerPolOCPValuesOnly : public OptimalControlProblem<Real, 2, 1>
  {
  public:
    VanDerPolOCPValuesOnly(Real const T = 5.0)
      : OptimalControlProblem<Real, 2, 1>("VanDerPolOCPValuesOnly", std::make_shared<VanDerPol>(), 0.0, T) {}
    Real mayer(VectorF const &, VectorF const & x_end) const override {return 10.0*x_end.squaredNorm();}
    Real lagrange(VectorF const & x, VectorU const & u, Real) const override {return x.squaredNorm() + u.squaredNorm();}
    Integer boundary_conditions_number() const override {return 2;}
    VectorX boundary_conditions(VectorF const & x_ini, VectorF const &) const override {return x_ini - VectorF(1.0, 0.0);}
    Integer path_constraints_number() const override {return 1;}
    VectorX path_constraints(VectorF const & x, VectorU const & u, Real) const override
    {return VectorX::Constant(1, x(1)*x(1) + u(0)*u(0));}
    VectorX path_constraints_upper() const override {return VectorX::Constant(1, 4.0);}
  };

  /*\
   |  Index-1 DAE: x' = -x + z + u, 0 = z - x^2 (states (x, z), finite differences Hessians)
  \*/

  class SemiLinearDAE : public Controlled<Real, 2, 1>
  {
  public:
    SemiLinearDAE() : Controlled<Real, 2, 1>("SemiLinearDAE") {}
    VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real) const override
    {return VectorF(x_dot(0) + x(0) - x(1) - u(0), x(1) - x(0)*x(0));}
    MatrixJF JF_x(VectorF const & x, VectorF const &, VectorU const &, Real) const override
    {MatrixJF J; J << 1.0, -1.0, -2.0*x(0), 1.0; return J;}
    MatrixJF JF_x_dot(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {MatrixJF J; J << 1.0, 0.0, 0.0, 0.0; return J;}
    MatrixJU JF_u(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {return MatrixJU(-1.0, 0.0);}
  };

  // min int (x^2 + u^2) dt, x(0) = x0, on [0, T]. The initial conditions must be consistent: the
  // algebraic state is fixed by z(0) = x(0)^2, otherwise z(0) is a free direction of the NLP.
  class SemiLinearDAEOCP : public OptimalControlProblem<Real, 2, 1>
  {
    Real m_x0;
  public:
    SemiLinearDAEOCP(Real const x0 = 0.5, Real const T = 2.0)
      : OptimalControlProblem<Real, 2, 1>("SemiLinearDAEOCP", std::make_shared<SemiLinearDAE>(), 0.0, T), m_x0(x0) {}
    Real lagrange(VectorF const & x, VectorU const & u, Real) const override {return x(0)*x(0) + u(0)*u(0);}
    VectorP lagrange_gradient(VectorF const & x, VectorU const & u, Real) const override {return VectorP(2.0*x(0), 0.0, 2.0*u(0));}
    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {MatrixP H(MatrixP::Zero()); H(0, 0) = 2.0; H(2, 2) = 2.0; return H;}
    Integer boundary_conditions_number() const override {return 2;}
    VectorX boundary_conditions(VectorF const & x_ini, VectorF const &) const override
    {return VectorX((VectorX(2) << x_ini(0) - this->m_x0, x_ini(1) - x_ini(0)*x_ini(0)).finished());}
  };

  // Equivalent ODE: x' = -x + x^2 + u
  class SemiLinearODE : public ControlledExplicit<Real, 1, 1>
  {
  public:
    SemiLinearODE() : ControlledExplicit<Real, 1, 1>("SemiLinearODE") {}
    VectorF f(VectorF const & x, VectorU const & u, Real) const override {return VectorF(-x(0) + x(0)*x(0) + u(0));}
    MatrixJF Jf_x(VectorF const & x, VectorU const &, Real) const override {return MatrixJF(-1.0 + 2.0*x(0));}
    MatrixJU Jf_u(VectorF const &, VectorU const &, Real) const override {return MatrixJU(1.0);}
  };

  class SemiLinearODEOCP : public OptimalControlProblem<Real, 1, 1>
  {
    Real m_x0;
  public:
    SemiLinearODEOCP(Real const x0 = 0.5, Real const T = 2.0)
      : OptimalControlProblem<Real, 1, 1>("SemiLinearODEOCP", std::make_shared<SemiLinearODE>(), 0.0, T), m_x0(x0) {}
    Real lagrange(VectorF const & x, VectorU const & u, Real) const override {return x(0)*x(0) + u(0)*u(0);}
    Integer boundary_conditions_number() const override {return 1;}
    VectorX boundary_conditions(VectorF const & x_ini, VectorF const &) const override {return VectorX::Constant(1, x_ini(0) - this->m_x0);}
  };

  /*\
   |  Double integrator x0' = x1, x1' = u, min int u^2/2
  \*/

  class DoubleIntegrator : public ControlledExplicit<Real, 2, 1>
  {
  public:
    DoubleIntegrator() : ControlledExplicit<Real, 2, 1>("DoubleIntegrator") {}
    VectorF f(VectorF const & x, VectorU const & u, Real) const override {return VectorF(x(1), u(0));}
    MatrixJF Jf_x(VectorF const &, VectorU const &, Real) const override {MatrixJF J; J << 0.0, 1.0, 0.0, 0.0; return J;}
    MatrixJU Jf_u(VectorF const &, VectorU const &, Real) const override {return MatrixJU(0.0, 1.0);}
    MatrixJF Hf_x_x(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJF::Zero();}
    MatrixJU Hf_x_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixJU::Zero();}
    MatrixUU Hf_u_u(VectorF const &, VectorU const &, Real, VectorF const &) const override {return MatrixUU::Zero();}
  };

  // Rest-to-rest: x(0) = (0, 0), x(1) = (1, 0), optimal u = 6 - 12 t, J* = 6.
  // Bryson-Denham: x(0) = (0, 1), x(1) = (0, -1), x0 <= 1/9, J* = 4.
  class DoubleIntegratorOCP : public OptimalControlProblem<Real, 2, 1>
  {
    VectorF m_x_ini, m_x_end;
  public:
    DoubleIntegratorOCP(VectorF const & x_ini, VectorF const & x_end)
      : OptimalControlProblem<Real, 2, 1>("DoubleIntegratorOCP", std::make_shared<DoubleIntegrator>(), 0.0, 1.0),
        m_x_ini(x_ini), m_x_end(x_end) {}
    Real lagrange(VectorF const &, VectorU const & u, Real) const override {return 0.5*u(0)*u(0);}
    VectorP lagrange_gradient(VectorF const &, VectorU const & u, Real) const override {return VectorP(0.0, 0.0, u(0));}
    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {MatrixP H(MatrixP::Zero()); H(2, 2) = 1.0; return H;}
    Integer boundary_conditions_number() const override {return 4;}
    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & x_end) const override
    {VectorX b(4); b << x_ini - this->m_x_ini, x_end - this->m_x_end; return b;}
    MatrixX boundary_conditions_jacobian(VectorF const &, VectorF const &) const override {return MatrixX::Identity(4, 4);}
    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override {return MatrixB::Zero();}
  };

} // namespace TestOCP

#endif // TESTS_OPTIMAL_CONTROL_PROBLEMS_HH
