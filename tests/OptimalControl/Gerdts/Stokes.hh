/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_GERDTS_STOKES_HH
#define TESTS_GERDTS_STOKES_HH

#include "DAEProblemsAD.hh"

// Distributed control of the 2D Stokes flow, discretized in space by the method of lines (an index-2
// DAE), from M. Gerdts, Optimal Control of ODEs and DAEs, 2nd edition, De Gruyter, 2023 (Examples
// 1.0.4, 1.1.11 and 3.1.12).
namespace TestOCP
{

  /*\
   |  The instationary Stokes equations v_t = Laplace(v) - grad(p) + u, 0 = div(v) on the unit square,
   |  with v = 0 on the boundary and at t = 0, discretized on the grid of step h = 1/G (the velocities
   |  and pressures at the (G-1)^2 interior points, the pressure set to zero on the last row and column)
   |  by the finite differences of Example 1.1.11:
   |    v' = A v + B p + u,  0 = B^T v                                   (index 2, B^T B non-singular)
   |  with A the 5-point Laplacian of each component and B = -grad_h (forward differences), so that
   |  B^T = div_h (backward differences). The control u is distributed on every velocity component, and
   |    min  1/2 int_0^T |v - vd|^2 dt + alpha/2 int_0^T |u|^2 dt,  T = 2,  alpha = 1e-5
   |  tracks the desired flow vd = (-q(t, x) q_y(t, y), q(t, y) q_x(t, x)), q(t, z) = (1 - z)^2 (1 -
   |  cos(2 pi z t)) (Example 3.1.12), which is divergence free.
   |  Index 1: the pressure from the time derivative of the constraint, 0 = B^T v' (a discrete Poisson
   |  equation for p), with the invariants h = B^T v.
   |  States (v1, v2, p) with the unknowns ordered row by row, controls (u1, u2).
  \*/

  template <Integer G>
  struct StokesGrid
  {
    static constexpr Integer n{G - 1};     // interior points per direction
    static constexpr Integer NV{2*n*n};    // velocity components
    static constexpr Integer NP{n*n};      // pressures
    static constexpr Integer N{NV + NP};   // states
  };

  template <Integer G, Integer Index>
  class Stokes : public OptimalControlProblemAD<Stokes<G, Index>, StokesGrid<G>::N, StokesGrid<G>::NV,
    Index == 2 ? 0 : StokesGrid<G>::NP, StokesGrid<G>::N>
  {
    using Grid = StokesGrid<G>;
    using Base = OptimalControlProblemAD<Stokes<G, Index>, Grid::N, Grid::NV, Index == 2 ? 0 : Grid::NP, Grid::N>;
    using MatrixVV = Eigen::Matrix<Real, Grid::NV, Grid::NV>;
    using MatrixVP = Eigen::Matrix<Real, Grid::NV, Grid::NP>;

    static constexpr Integer n{Grid::n}, NV{Grid::NV}, NP{Grid::NP};
    static constexpr Real h{1.0/G};

    MatrixVV m_A{MatrixVV::Zero()}; // discrete Laplacian
    MatrixVP m_B{MatrixVP::Zero()}; // minus the discrete gradient

    // Index of the interior point (i, j), i, j = 1, ..., n
    static Integer point(Integer const i, Integer const j) {return (i - 1) + (j - 1)*n;}

  public:
    template <typename T> using StateVector     = typename Base::template StateVector<T>;
    template <typename T> using ControlVector   = typename Base::template ControlVector<T>;
    template <typename T> using InvariantVector = typename Base::template InvariantVector<T>;
    template <typename T> using BoundaryVector  = typename Base::template BoundaryVector<T>;

    static constexpr Real ALPHA{1.0e-5};

    Stokes() : Base("Stokes" + std::to_string(G) + "Index" + std::to_string(Index), 0.0, 2.0)
    {
      for (Integer c{0}; c < 2; ++c) {       // velocity component
        for (Integer j{1}; j <= n; ++j) {
          for (Integer i{1}; i <= n; ++i) {
            Integer const k{c*NP + point(i, j)};
            m_A(k, k) = -4.0/(h*h);          // neighbours on the boundary are zero
            if (i > 1) {m_A(k, k - 1) = 1.0/(h*h);}
            if (i < n) {m_A(k, k + 1) = 1.0/(h*h);}
            if (j > 1) {m_A(k, k - n) = 1.0/(h*h);}
            if (j < n) {m_A(k, k + n) = 1.0/(h*h);}
            // -grad p: forward difference in x (c = 0) or y (c = 1), p = 0 beyond the last point
            m_B(k, point(i, j)) = 1.0/h;
            if (c == 0 && i < n) {m_B(k, point(i + 1, j)) = -1.0/h;}
            if (c == 1 && j < n) {m_B(k, point(i, j + 1)) = -1.0/h;}
          }
        }
      }
    }

    // The desired flow at the interior points
    template <typename T>
    static Eigen::Vector<T, NV> desired(T const & t)
    {
      auto q   = [&t] (Real const z) {return (1.0 - z)*(1.0 - z)*(1.0 - cos(2.0*M_PI*z*t));};
      auto q_z = [&t] (Real const z) {return -2.0*(1.0 - z)*(1.0 - cos(2.0*M_PI*z*t)) + (1.0 - z)*(1.0 - z)*2.0*M_PI*t*sin(2.0*M_PI*z*t);};
      Eigen::Vector<T, NV> vd;
      for (Integer j{1}; j <= n; ++j) {
        for (Integer i{1}; i <= n; ++i) {
          Real const x{i*h}, y{j*h};
          vd(point(i, j))      = -q(x)*q_z(y);
          vd(NP + point(i, j)) = q(y)*q_z(x);
        }
      }
      return vd;
    }

    // The squared distance of the desired flow from the discretely divergence-free fields (its part
    // B (B^T B)^-1 B^T vd along the discrete gradients), which no feasible flow can track: half its
    // time integral bounds the optimal cost from below.
    Real unreachable(Real const t) const
    {
      Eigen::Vector<Real, NV> const vd{desired(t)};
      Eigen::Vector<Real, NP> const z{(m_B.transpose()*m_B).ldlt().solve(m_B.transpose()*vd)};
      return (m_B*z).squaredNorm();
    }

    // Divergence of the velocity field (the constraint B^T v)
    template <typename T>
    Eigen::Vector<T, NP> divergence(Eigen::Vector<T, NV> const & v) const {return m_B.transpose().template cast<T>()*v;}

    template <typename T>
    StateVector<T> dae(StateVector<T> const & x, StateVector<T> const & x_dot, ControlVector<T> const & u, T const &) const
    {
      Eigen::Vector<T, NV> const v{x.template head<NV>()}, v_dot{x_dot.template head<NV>()};
      StateVector<T> F;
      F.template head<NV>() = v_dot - m_A.template cast<T>()*v - m_B.template cast<T>()*x.template tail<NP>() - u;
      F.template tail<NP>() = Index == 2 ? divergence(v) : divergence(v_dot);
      return F;
    }

    template <typename T>
    InvariantVector<T> invariants(StateVector<T> const & x, T const &) const
    {
      InvariantVector<T> c;
      if constexpr (Index == 1) {c = divergence<T>(x.template head<NV>());}
      return c;
    }

    template <typename T>
    T running_cost(StateVector<T> const & x, ControlVector<T> const & u, T const & t) const
    {return 0.5*(x.template head<NV>() - desired(t)).squaredNorm() + 0.5*ALPHA*u.squaredNorm();}

    template <typename T>
    BoundaryVector<T> boundary(StateVector<T> const & x_ini, StateVector<T> const &) const
    {return x_ini;} // at rest, and the pressure gauge p(0) = 0
  };

} // namespace TestOCP

#endif // TESTS_GERDTS_STOKES_HH
