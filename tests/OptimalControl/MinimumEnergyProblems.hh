/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_MINIMUM_ENERGY_PROBLEMS_HH
#define TESTS_MINIMUM_ENERGY_PROBLEMS_HH

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/System/Implicit.hh"
#include "Sandals/System/Controlled.hh"
#include "Sandals/OptimalControl.hh"

// Minimum energy problem (Example 1.0.1 in Gerdts, Optimal Control of ODEs and DAEs, 2023):
//   min  1/2 int_0^1 u^2 dt
//   s.t. y1' = u - y4,  y2' = u,  y3' = -y2,  0 = y1 + y3   (index 2, y4 algebraic)
//        y1(0) = y1(1) = 0,  y2(0) = -y2(1) = 1,  y3(0) = 0
// The optimal control is u = -2, with J* = 2, and
//   y1 = t - t^2,  y2 = 1 - 2 t,  y3 = t^2 - t,  y4 = 2 t - 3.
namespace TestOCP
{
  using namespace Sandals;
  using Real = double;

  /*\
   |  Indirect approach: index-0 boundary value problem with invariants
  \*/

  // Necessary conditions of the index-0 formulation (the constraint differentiated twice before
  // applying the minimum principle). States x = (y1, y2, y3, y4, l1, l2, l3, l4), control
  // u = -l1 - l2 + l4, and
  //   y1' - u + y4 = 0,  y2' - u = 0,  y3' + y2 = 0,  y4' - u' + u = 0,
  //   l1' = 0,  l2' - l3 = 0,  l3' = 0,  l4' - l1 = 0,
  // with the invariants lost by the index reduction h = (y1 + y3, u - y2 - y4).
  // Exact solution: u = -2, l = (0, 2, 0, 0).
  class MinimumEnergyIndex0 : public Implicit<Real, 8, 2>
  {
  public:
    enum : Integer {Y1 = 0, Y2 = 1, Y3 = 2, Y4 = 3, L1 = 4, L2 = 5, L3 = 6, L4 = 7};

    MinimumEnergyIndex0() : Implicit<Real, 8, 2>("MinimumEnergyIndex0") {}

    // Control and its gradient with respect to the states
    static Real control(VectorF const & x) {return -x(L1) - x(L2) + x(L4);}
    static VectorF control_gradient()
    {
      VectorF d(VectorF::Zero());
      d(L1) = -1.0;
      d(L2) = -1.0;
      d(L4) = 1.0;
      return d;
    }

    VectorF F(VectorF const & x, VectorF const & x_dot, Real /*t*/) const override
    {
      Real const u{control(x)};
      Real const u_dot{control_gradient().dot(x_dot)};
      VectorF F;
      F(0) = x_dot(Y1) - u + x(Y4);
      F(1) = x_dot(Y2) - u;
      F(2) = x_dot(Y3) + x(Y2);
      F(3) = x_dot(Y4) - u_dot + u;
      F(4) = x_dot(L1);
      F(5) = x_dot(L2) - x(L3);
      F(6) = x_dot(L3);
      F(7) = x_dot(L4) - x(L1);
      return F;
    }

    MatrixJF JF_x(VectorF const &, VectorF const &, Real) const override
    {
      VectorF const du{control_gradient()};
      MatrixJF J(MatrixJF::Zero());
      J.row(0)  = -du.transpose();
      J(0, Y4) += 1.0;
      J.row(1)  = -du.transpose();
      J(2, Y2)  = 1.0;
      J.row(3)  = du.transpose();
      J(5, L3)  = -1.0;
      J(7, L1)  = -1.0;
      return J;
    }

    MatrixJF JF_x_dot(VectorF const &, VectorF const &, Real) const override
    {
      MatrixJF J(MatrixJF::Identity());
      J.row(3) -= control_gradient().transpose(); // y4' - u'
      return J;
    }

    VectorH h(VectorF const & x, Real /*t*/) const override
    {
      return VectorH(x(Y1) + x(Y3), control(x) - x(Y2) - x(Y4));
    }

    MatrixJH Jh_x(VectorF const &, Real) const override
    {
      MatrixJH J(MatrixJH::Zero());
      J(0, Y1)  = 1.0;
      J(0, Y3)  = 1.0;
      J.row(1)  = control_gradient().transpose();
      J(1, Y2) -= 1.0;
      J(1, Y4) -= 1.0;
      return J;
    }

    bool in_domain(VectorF const &, Real) const override {return true;}

    // Exact solution
    static VectorF solution(Real const t)
    {
      VectorF x(VectorF::Zero());
      x(Y1) = t - t*t;
      x(Y2) = 1.0 - 2.0*t;
      x(Y3) = t*t - t;
      x(Y4) = 2.0*t - 3.0;
      x(L2) = 2.0;
      return x;
    }
  };

  // Boundary conditions:
  //   y1(0) = 0, y2(0) = 1, y3(0) = 0, y1(1) = 0, y2(1) = -1, l4(0) = 0, l3(1) = 0, l4(1) = 0
  template <typename Integrator>
  class MinimumEnergyIndex0BVP : public Problem<Real, 8, 2, Integrator>
  {
    using S = MinimumEnergyIndex0;

  public:
    using typename Problem<Real, 8, 2, Integrator>::VectorF;
    using typename Problem<Real, 8, 2, Integrator>::MatrixJF;

    MinimumEnergyIndex0BVP(std::shared_ptr<Integrator> integrator)
      : Problem<Real, 8, 2, Integrator>("MinimumEnergyIndex0BVP", std::make_shared<MinimumEnergyIndex0>(), integrator)
    {}

    VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
    {
      VectorF b;
      b(0) = x_ini(S::Y1);
      b(1) = x_ini(S::Y2) - 1.0;
      b(2) = x_ini(S::Y3);
      b(3) = x_end(S::Y1);
      b(4) = x_end(S::Y2) + 1.0;
      b(5) = x_ini(S::L4);
      b(6) = x_end(S::L3);
      b(7) = x_end(S::L4);
      return b;
    }

    MatrixJF Jb_x_ini(VectorF const &, VectorF const &) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(0, S::Y1) = 1.0;
      J(1, S::Y2) = 1.0;
      J(2, S::Y3) = 1.0;
      J(5, S::L4) = 1.0;
      return J;
    }

    MatrixJF Jb_x_end(VectorF const &, VectorF const &) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(3, S::Y1) = 1.0;
      J(4, S::Y2) = 1.0;
      J(6, S::L3) = 1.0;
      J(7, S::L4) = 1.0;
      return J;
    }
  };

  /*\
   |  Direct approach: index-2 optimal control problem
  \*/

  // States (y1, y2, y3, y4), with y4 algebraic: linear system, so the default Hessians are exact.
  class MinimumEnergyIndex2 : public Controlled<Real, 4, 1>
  {
  public:
    MinimumEnergyIndex2() : Controlled<Real, 4, 1>("MinimumEnergyIndex2") {}

    VectorF F(VectorF const & x, VectorF const & x_dot, VectorU const & u, Real /*t*/) const override
    {
      VectorF F;
      F(0) = x_dot(0) - u(0) + x(3);
      F(1) = x_dot(1) - u(0);
      F(2) = x_dot(2) + x(1);
      F(3) = x(0) + x(2);
      return F;
    }

    MatrixJF JF_x(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J(0, 3) = 1.0;
      J(2, 1) = 1.0;
      J(3, 0) = 1.0;
      J(3, 2) = 1.0;
      return J;
    }

    MatrixJF JF_x_dot(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {
      MatrixJF J(MatrixJF::Zero());
      J.diagonal() << 1.0, 1.0, 1.0, 0.0; // y4 is algebraic
      return J;
    }

    MatrixJU JF_u(VectorF const &, VectorF const &, VectorU const &, Real) const override
    {
      MatrixJU J(MatrixJU::Zero());
      J(0, 0) = -1.0;
      J(1, 0) = -1.0;
      return J;
    }
  };

  // Boundary conditions of the problem, plus y4(0) = -3: the node value of the algebraic y4 only
  // shifts the stage variables (y4' never appears), so it is a gauge, fixed to its exact value.
  class MinimumEnergyIndex2OCP : public OptimalControlProblem<Real, 4, 1>
  {
  public:
    MinimumEnergyIndex2OCP()
      : OptimalControlProblem<Real, 4, 1>("MinimumEnergyIndex2OCP", std::make_shared<MinimumEnergyIndex2>(), 0.0, 1.0)
    {}

    Real lagrange(VectorF const &, VectorU const & u, Real) const override {return 0.5*u(0)*u(0);}

    VectorP lagrange_gradient(VectorF const &, VectorU const & u, Real) const override
    {
      VectorP gradient(VectorP::Zero());
      gradient(4) = u(0);
      return gradient;
    }

    MatrixP lagrange_hessian(VectorF const &, VectorU const &, Real) const override
    {
      MatrixP hessian(MatrixP::Zero());
      hessian(4, 4) = 1.0;
      return hessian;
    }

    Integer boundary_conditions_number() const override {return 6;}

    VectorX boundary_conditions(VectorF const & x_ini, VectorF const & x_end) const override
    {
      VectorX b(6);
      b << x_ini(0), x_ini(1) - 1.0, x_ini(2), x_ini(3) + 3.0, x_end(0), x_end(1) + 1.0;
      return b;
    }

    MatrixX boundary_conditions_jacobian(VectorF const &, VectorF const &) const override
    {
      MatrixX J(MatrixX::Zero(6, 8)); // columns: [x_ini, x_end]
      J(0, 0) = 1.0;
      J(1, 1) = 1.0;
      J(2, 2) = 1.0;
      J(3, 3) = 1.0;
      J(4, 4) = 1.0;
      J(5, 5) = 1.0;
      return J;
    }

    MatrixB boundary_conditions_hessian(VectorF const &, VectorF const &, VectorX const &) const override
    {
      return MatrixB::Zero();
    }
  };

} // namespace TestOCP

#endif // TESTS_MINIMUM_ENERGY_PROBLEMS_HH
