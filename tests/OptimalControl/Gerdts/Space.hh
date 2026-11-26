/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_GERDTS_SPACE_HH
#define TESTS_GERDTS_SPACE_HH

#include "DAEProblemsAD.hh"
#include "Sandals/System/Explicit.hh"

// Space example: the Arenstorf orbit of a satellite in the Earth-Moon system, from M. Gerdts, Optimal
// Control of ODEs and DAEs, 2nd edition, De Gruyter, 2023 (Example 4.4.1), as an explicit ODE for the
// Runge-Kutta integrators of the library.
namespace TestOCP
{

  /*\
   |  Arenstorf orbit (Example 4.4.1): the restricted three-body problem in the rotating frame of the
   |  Earth (mass 1 - mu at x = -mu) and the Moon (mass mu at x = 1 - mu),
   |    x'' = x + 2 y' - (1 - mu) (x + mu)/D1 - mu (x - 1 + mu)/D2,
   |    y'' = y - 2 x' - (1 - mu) y/D1 - mu y/D2,
   |    D1 = ((x + mu)^2 + y^2)^(3/2),  D2 = ((x - 1 + mu)^2 + y^2)^(3/2),
   |  with mu = 0.012277471. From (x, y, x', y') = (0.994, 0, 0, -2.001585106379) the orbit is periodic
   |  with period 17.065216560158, after close passages by the Moon: the book uses it to show that the
   |  classic Runge-Kutta method needs N = 10000 fixed steps to close it, against a few thousand
   |  evaluations with a step-size control.
  \*/

  struct ArenstorfModel
  {
    static constexpr Real MU{0.012277471};
    static constexpr Real X0{0.994};
    static constexpr Real VY0{-2.001585106379};
    static constexpr Real PERIOD{17.065216560158};

    using System = Sandals::ExplicitWrapper<Real, 4, 0>;

    // The initial states (x, y, x', y') of the periodic orbit
    static Eigen::Vector4d initial_states() {return Eigen::Vector4d(X0, 0.0, 0.0, VY0);}

    // The explicit ODE x' = f(x) for the Runge-Kutta integrators (the Jacobian by central differences)
    static std::shared_ptr<System> system()
    {
      auto f = [] (System::VectorF const & s, Real) {return dynamics<Real>(s);};
      auto Jf_x = [] (System::VectorF const & s, Real) {
        System::MatrixJF J;
        for (Integer j{0}; j < 4; ++j) {
          System::VectorF s_p(s), s_m(s);
          s_p(j) += 1.0e-7;
          s_m(j) -= 1.0e-7;
          J.col(j) = (dynamics<Real>(s_p) - dynamics<Real>(s_m))/2.0e-7;
        }
        return J;
      };
      return std::make_shared<System>(f, Jf_x);
    }

    // The time derivatives of (x, y, x', y')
    template <typename T>
    static Eigen::Vector<T, 4> dynamics(Eigen::Vector<T, 4> const & s)
    {
      T const x{s(0)}, y{s(1)}, vx{s(2)}, vy{s(3)};
      T const r1{sqrt((x + MU)*(x + MU) + y*y)}, r2{sqrt((x - 1.0 + MU)*(x - 1.0 + MU) + y*y)};
      T const D1{r1*r1*r1}, D2{r2*r2*r2};
      return Eigen::Vector<T, 4>(vx, vy,
        x + 2.0*vy - (1.0 - MU)*(x + MU)/D1 - MU*(x - 1.0 + MU)/D2,
        y - 2.0*vx - (1.0 - MU)*y/D1 - MU*y/D2);
    }
  };

} // namespace TestOCP

#endif // TESTS_GERDTS_SPACE_HH
