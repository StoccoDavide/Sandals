/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef TESTS_GERDTS_PROCESS_HH
#define TESTS_GERDTS_PROCESS_HH

#include "DAEProblemsAD.hh"

// Process engineering example: the semi-batch chemical reactor (index-1 DAE), from M. Gerdts,
// Optimal Control of ODEs and DAEs, 2nd edition, De Gruyter, 2023 (the example numbers refer to the
// book). Conventions as in DAEProblemsAD.hh: the node values of the algebraic states at t = 0 are
// gauges, fixed by a boundary condition.
namespace TestOCP
{

  /*\
   |  Semi-batch chemical reactor (Example 8.2.8): reactions A + B -> C and B + C -> D, fed with B at
   |  the rate FB and cooled with the power Q, maximizing the product C:
   |    min  -MC(tf) + 1e-2 int_0^20 (FB^2 + Q^2) dt
   |    s.t. MA' = -V r1,  MB' = FB - V (r1 + r2),  MC' = V (r1 - r2),  MD' = V r2,
   |         H' = 20 FB - Q + V (r1 + 75 r2),
   |         0 = H - sum_i Mi (alpha_i (TR - Tref) + beta_i/2 (TR^2 - Tref^2))     (index 1)
   |  with the reactor temperature TR algebraic, the volume V = sum_i Mi/rho_i, the concentrations
   |  Ci = Mi/V, and the rates r1 = A1 exp(-E1/TR) CA CB, r2 = A2 exp(-E2/TR) CB CC.
  \*/

  class ChemicalReactor : public OptimalControlProblemAD<ChemicalReactor, 6, 2, 0, 6>
  {
  public:
    static constexpr Real T_REF{298.0};
    static constexpr Real ALPHA[4]{0.1723, 0.2, 0.16, 0.155};
    static constexpr Real BETA[4]{0.000473, 0.0005, 0.00055, 0.000323};
    static constexpr Real RHO[4]{11250.0, 16000.0, 10400.0, 10000.0};
    static constexpr Real A1{0.8}, A2{0.002}, E1{3000.0}, E2{2400.0};

    ChemicalReactor() : OptimalControlProblemAD("ChemicalReactor", 0.0, 20.0)
    {
      VectorF nominal;
      nominal << 1.0e4, 10.0, 1.0e3, 1.0, 1.0e5, 100.0; // masses, energy, temperature
      this->nominal(nominal, VectorU(10.0, 0.01));
    }

    // The temperature TR from the energy balance, a quadratic equation (positive root)
    static Real temperature(Eigen::Vector<Real, 4> const & mass, Real const H)
    {
      Real a{0.0}, b{0.0};
      for (Integer i{0}; i < 4; ++i) {a += mass(i)*BETA[i]/2.0; b += mass(i)*ALPHA[i];}
      Real const c{-H - a*T_REF*T_REF - b*T_REF};
      return (-b + std::sqrt(b*b - 4.0*a*c))/(2.0*a);
    }

    // Consistent initial states
    static VectorF initial_states()
    {
      VectorF x;
      x << 9000.0, 0.0, 0.0, 0.0, 152509.97, 0.0;
      x(5) = temperature(x.head<4>(), x(4));
      return x;
    }

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      T const TR{x(5)}, FB{u(0)}, Q{u(1)};
      T V{0.0}, energy{0.0};
      for (Integer i{0}; i < 4; ++i) {
        V += x(i)/RHO[i];
        energy += x(i)*(ALPHA[i]*(TR - T_REF) + BETA[i]/2.0*(TR*TR - T_REF*T_REF));
      }
      T const CA{x(0)/V}, CB{x(1)/V}, CC{x(2)/V};
      T const r1{A1*exp(-E1/TR)*CA*CB}, r2{A2*exp(-E2/TR)*CB*CC};
      StateVector<T> F;
      F(0) = x_dot(0) + V*r1;
      F(1) = x_dot(1) - FB + V*(r1 + r2);
      F(2) = x_dot(2) - V*(r1 - r2);
      F(3) = x_dot(3) - V*r2;
      F(4) = x_dot(4) - 20.0*FB + Q - V*(r1 + 75.0*r2);
      F(5) = x(4) - energy;
      return F;
    }

    template <typename T>
    T running_cost(StateVector<T> const &, ControlVector<T> const & u, T const &) const {return 1.0e-2*(u(0)*u(0) + u(1)*u(1));}

    template <typename T>
    T terminal_cost(StateVector<T> const &, StateVector<T> const & x_end) const {return -x_end(2);}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {
      VectorF const start{initial_states()};
      BoundaryVector<T> b;
      for (Integer i{0}; i < 6; ++i) {b(i) = x_ini(i) - start(i);} // TR(0) is a gauge
      return b;
    }
  };

} // namespace TestOCP

#endif // TESTS_GERDTS_PROCESS_HH
