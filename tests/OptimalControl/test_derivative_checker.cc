/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// IPOPT's derivative checker (first and second order) on the multiple shooting transcription of every
// test problem: the gradient, the constraints Jacobian, and the Hessian of the Lagrangian, through
// the Runge-Kutta substeps and the node least squares, against IPOPT's own finite differences.

#include <fstream>

#include <gtest/gtest.h>

#include "Sandals/Ipopt.hh"
#include "Sandals/RungeKutta/RadauIIA5.hh"
#include "Sandals/RungeKutta/RK4.hh"

#include "OptimalControlProblems.hh"
#include "PendulumProblems.hh"
#include "DAEProblemsAD.hh"
#include "HardProblemsAD.hh"
#include "GerdtsProblemsAD.hh"

using namespace Sandals;
using namespace TestOCP;

using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;

using MatrixX = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

template <typename NLP>
MatrixX dense_jacobian(NLP const & nlp, VectorX const & z)
{
  MatrixX J(MatrixX::Zero(nlp.constraints_number(), nlp.variables_number()));
  VectorX const v{nlp.jacobian_values(z)};
  for (Integer e{0}; e < nlp.jacobian_nnz(); ++e) {J(nlp.jacobian_rows()[e], nlp.jacobian_cols()[e]) += v(e);}
  return J;
}

// Central-difference step, relative to the magnitude of the variable (as IPOPT's checker does).
Real step(Real const z) {return 1.0e-6*std::max(Real(1.0), std::abs(z));}

// Constraints Jacobian and objective gradient against central differences of the values.
template <typename NLP>
Real jacobian_error(NLP const & nlp, VectorX const & z)
{
  MatrixX J_fd(nlp.constraints_number() + 1, z.size()), J(nlp.constraints_number() + 1, z.size());
  J << nlp.objective_gradient(z).transpose(), dense_jacobian(nlp, z);
  for (Integer j{0}; j < z.size(); ++j) {
    Real const delta{step(z(j))};
    VectorX z_p(z), z_m(z);
    z_p(j) += delta;
    z_m(j) -= delta;
    J_fd.col(j) << (nlp.objective(z_p) - nlp.objective(z_m)) / (2.0*delta), (nlp.constraints(z_p) - nlp.constraints(z_m)) / (2.0*delta);
  }
  return (J - J_fd).norm() / std::max(Real(1.0), J_fd.norm());
}

// Hessian of the Lagrangian against central differences of the exact gradient and Jacobian.
template <typename NLP>
Real lagrangian_hessian_error(NLP const & nlp, VectorX const & z, Real const sigma, VectorX const & lambda)
{
  auto gradient = [&] (VectorX const & w) {return VectorX(sigma*nlp.objective_gradient(w) + dense_jacobian(nlp, w).transpose()*lambda);};
  MatrixX H(MatrixX::Zero(z.size(), z.size()));
  VectorX const v{nlp.hessian_values(z, sigma, lambda)};
  for (Integer e{0}; e < nlp.hessian_nnz(); ++e) {
    Integer const r{nlp.hessian_rows()[e]}, c{nlp.hessian_cols()[e]};
    H(r, c) += v(e);
    if (r != c) {H(c, r) += v(e);}
  }
  MatrixX H_fd(z.size(), z.size());
  for (Integer j{0}; j < z.size(); ++j) {
    Real const delta{step(z(j))};
    VectorX z_p(z), z_m(z);
    z_p(j) += delta;
    z_m(j) -= delta;
    H_fd.col(j) = (gradient(z_p) - gradient(z_m)) / (2.0*delta);
  }
  return (H - H_fd).norm() / std::max(Real(1.0), H_fd.norm());
}

// The NLP in the variables and constraints scaled by the nominal magnitudes of the problem (as IPOPT
// sees it with its user scaling), so that the forward differences of IPOPT's checker act on quantities
// of order one: on a badly scaled model (states up to 1e5) their round-off would hide errors.
template <typename NLP>
class ScaledNLP
{
  NLP & m_nlp;
  VectorX m_z, m_c; // scaling factors of the variables and of the constraints

  VectorX unscale(VectorX const & w) const {return w.cwiseQuotient(m_z);}

public:
  explicit ScaledNLP(NLP & nlp) : m_nlp(nlp) {nlp.scaling(m_z, m_c);}

  Integer variables_number() const {return m_nlp.variables_number();}
  Integer constraints_number() const {return m_nlp.constraints_number();}
  Integer jacobian_nnz() const {return m_nlp.jacobian_nnz();}
  Integer hessian_nnz() const {return m_nlp.hessian_nnz();}
  auto const & jacobian_rows() const {return m_nlp.jacobian_rows();}
  auto const & jacobian_cols() const {return m_nlp.jacobian_cols();}
  auto const & hessian_rows() const {return m_nlp.hessian_rows();}
  auto const & hessian_cols() const {return m_nlp.hessian_cols();}

  void bounds(VectorX & z_lower, VectorX & z_upper, VectorX & c_lower, VectorX & c_upper) const
  {
    m_nlp.bounds(z_lower, z_upper, c_lower, c_upper);
    z_lower = z_lower.cwiseProduct(m_z); z_upper = z_upper.cwiseProduct(m_z);
    c_lower = c_lower.cwiseProduct(m_c); c_upper = c_upper.cwiseProduct(m_c);
  }
  VectorX guess() const {return m_nlp.guess().cwiseProduct(m_z);}
  void guess(VectorX const & w) {m_nlp.guess(unscale(w));}

  Real objective(VectorX const & w) const {return m_nlp.objective(unscale(w));}
  VectorX objective_gradient(VectorX const & w) const {return m_nlp.objective_gradient(unscale(w)).cwiseQuotient(m_z);}
  VectorX constraints(VectorX const & w) const {return m_nlp.constraints(unscale(w)).cwiseProduct(m_c);}
  VectorX jacobian_values(VectorX const & w) const
  {
    VectorX v{m_nlp.jacobian_values(unscale(w))};
    for (Integer e{0}; e < v.size(); ++e) {v(e) *= m_c(m_nlp.jacobian_rows()[e]) / m_z(m_nlp.jacobian_cols()[e]);}
    return v;
  }
  VectorX hessian_values(VectorX const & w, Real const sigma, VectorX const & lambda) const
  {
    VectorX v{m_nlp.hessian_values(unscale(w), sigma, lambda.cwiseProduct(m_c))};
    for (Integer e{0}; e < v.size(); ++e) {v(e) /= m_z(m_nlp.hessian_rows()[e]) * m_z(m_nlp.hessian_cols()[e]);}
    return v;
  }
};

using Radau = RadauIIA5Tableau<Real>;
using ExplicitRK4 = RK4Tableau<Real>;

// Run the second-order derivative checker at the simulated starting point, moved away from zero
// components (where forward differences of quadratic terms are inaccurate). Returns the checker report.
template <typename TableauType, Integer N, Integer U, Integer M, typename Problem>
std::string derivative_check(std::shared_ptr<Problem> problem, Integer const intervals, Integer const substeps,
  Eigen::Vector<Real, N> const & x_ini, Real const u_guess, Real const epsilon, std::string const & order,
  Real const tolerance)
{
  constexpr Integer S{TableauType::Vector::RowsAtCompileTime};
  MultipleShooting<Real, S, N, U, M> nlp(problem, TableauType(), intervals, substeps, epsilon);
  EXPECT_TRUE(nlp.simulate_guess(x_ini, [u_guess] (Real) {return Eigen::Vector<Real, U>::Constant(u_guess).eval();}));
  VectorX z{nlp.guess()};
  for (Integer i{0}; i < z.size(); ++i) {z(i) += 0.01 + 0.005*std::sin(Real(i));}
  nlp.guess(z);

  // Central differences of the exact gradient and Jacobian (truncation O(delta^2))
  std::srand(5);
  VectorX const lambda{VectorX::Random(nlp.constraints_number())};
  EXPECT_LT(jacobian_error(nlp, z), 1.0e-7) << problem->name();
  EXPECT_LT(lagrangian_hessian_error(nlp, z, 0.7, lambda), 1.0e-6) << problem->name();

  std::string const file{"derivative_check_" + problem->name() + ".txt"};
  IpoptOptions opts;
  opts.integers["print_level"] = 0;
  opts.strings["sb"] = "yes";
  opts.integers["max_iter"] = 0;
  opts.strings["derivative_test"] = order;
  opts.numerics["derivative_test_tol"] = tolerance;
  // Forward-difference step balancing truncation (delta f'') and round-off (eps |c| / delta): with
  // constraints of order 10, 1e-7 leaves a round-off floor of ~2e-8 on exactly-zero entries, above
  // the checker's absolute floor (~ tol^2)
  opts.numerics["derivative_test_perturbation"] = 1.0e-6;
  opts.numerics["point_perturbation_radius"] = 0.0;
  opts.strings["output_file"] = file;
  opts.integers["file_print_level"] = 5;
  Real J;
  if (nlp.scaled()) {
    ScaledNLP<decltype(nlp)> scaled(nlp);
    ipopt_solve(scaled, opts, J);
  } else {
    ipopt_solve(nlp, opts, J);
  }
  std::ifstream log(file);
  return std::string((std::istreambuf_iterator<char>(log)), std::istreambuf_iterator<char>());
}

// The checker report must state that it found no errors (otherwise show the flagged entries).
void expect_no_errors(std::string const & report, std::string const & name)
{
  bool const clean{report.find("No errors detected by derivative checker") != std::string::npos};
  std::string flagged;
  std::istringstream lines(report);
  for (std::string line; std::getline(lines, line);) {if (line.rfind("*", 0) == 0) {flagged += line + "\n";}}
  EXPECT_TRUE(clean) << name << ":\n" << flagged.substr(0, 3000);
}

// Check one problem (the report is named after the problem).
template <typename TableauType, Integer N, Integer U, Integer M, typename Problem>
void check(std::shared_ptr<Problem> problem, Integer const intervals, Integer const substeps,
  Eigen::Vector<Real, N> const & x_ini, Real const u_guess, Real const epsilon = 0.0,
  std::string const & order = "second-order", Real const tolerance = 1.0e-3)
{
  expect_no_errors(derivative_check<TableauType, N, U, M>(problem, intervals, substeps, x_ini, u_guess, epsilon, order, tolerance),
    problem->name() + " (" + order + ")");
}

TEST(DerivativeChecker, HandCodedProblems) {
  PendulumParameters const p;
  check<Radau, 2, 1, 0>(std::make_shared<VanDerPolOCP>(1.0), 4, 2, Eigen::Vector2d(1.0, 0.0), 0.2);
  check<Radau, 2, 1, 0>(std::make_shared<SemiLinearDAEOCP>(), 4, 2, Eigen::Vector2d(0.5, 0.25), 0.2);
  check<Radau, 5, 1, 2>(std::make_shared<PendulumSwingUpOCP<PendulumIndex1<2>>>(2.0), 4, 2, p.hanging_at_rest<5>(), 0.5);
  check<Radau, 5, 1, 1>(std::make_shared<PendulumSwingUpOCP<PendulumIndex2<1>>>(2.0), 4, 2, p.hanging_at_rest<5>(), 0.5, 1.0e-3);
}

TEST(DerivativeChecker, DAEProblems) {
  Eigen::Vector4d const me(0.0, 1.0, 0.0, -3.0);
  Eigen::Vector<Real, 5> gp;
  gp << 1.0, 0.0, 0.0, 0.0, 0.0;
  check<Radau, 4, 1, 0>(std::make_shared<MinimumEnergyIndex2>(), 4, 2, me, 0.3);
  check<Radau, 4, 1, 1>(std::make_shared<MinimumEnergyIndex1>(), 4, 2, me, 0.3);
  // Index 3: the stage Jacobians are ill-conditioned (~ h^-3), which amplifies the round-off of the exact
  // gradients beyond what forward differences of them can resolve; the Hessian is checked below with
  // central differences instead
  check<Radau, 5, 1, 0>(std::make_shared<GerdtsPendulum<3>>(), 4, 2, gp, 0.3, 0.0, "first-order");
  check<Radau, 5, 1, 1>(std::make_shared<GerdtsPendulum<2>>(), 4, 2, gp, 0.3);
  check<Radau, 5, 1, 2>(std::make_shared<GerdtsPendulum<1>>(), 4, 2, gp, 0.3);
  check<Radau, 4, 1, 2>(std::make_shared<GerdtsPendulum<0>>(), 4, 2, gp.head<4>().eval(), 0.3);
  Eigen::Vector<Real, 7> as(Eigen::Vector<Real, 7>::Zero());
  as(0) = 1.0;
  check<Radau, 7, 3, 1>(std::make_shared<AttitudeSlew>(), 3, 2, as, 0.2);
  // Kepler, double pendulum and cart-pole (Cartesian): large curvature (1/r^3 gravity over long
  // intervals, stiff rod constraints), so forward differences are accurate to ~5e-3 only; the central
  // differences above confirm their derivatives to 1e-7 (first order) and 1e-6 (second order)
  Real const curved{1.0e-2};
  auto kepler = std::make_shared<KeplerOrbit<2>>(0.5, 1);
  check<ExplicitRK4, 4, 1, 2>(kepler, 4, 3, kepler->periapsis(), 0.1, 0.0, "second-order", curved);
  DoublePendulumParameters const dp;
  Eigen::Vector<Real, 10> dc(Eigen::Vector<Real, 10>::Zero());
  dc(1) = -dp.L1; dc(3) = -dp.L1 - dp.L2; dc(8) = (dp.m1 + dp.m2)*dp.g/(2.0*dp.L1); dc(9) = dp.m2*dp.g/(2.0*dp.L2);
  check<Radau, 10, 1, 4>(std::make_shared<DoublePendulumCartesian>(), 3, 2, dc, 0.3, 0.0, "second-order", curved);
  check<Radau, 4, 1, 0>(std::make_shared<DoublePendulumAngles>(), 3, 2, Eigen::Vector4d::Zero().eval(), 0.3);
  check<ExplicitRK4, 1, 1, 0>(std::make_shared<StiffTracking>(), 4, 10, Eigen::Vector<Real, 1>(0.0), 0.5);
}

TEST(DerivativeChecker, Index3HessianCentralDifferences) {
  Eigen::Vector<Real, 5> gp;
  gp << 1.0, 0.0, 0.0, 0.0, 0.0;
  MultipleShooting<Real, 3, 5, 1, 0> nlp(std::make_shared<GerdtsPendulum<3>>(), Radau(), 4, 2);
  ASSERT_TRUE(nlp.simulate_guess(gp, [] (Real) {return Eigen::Vector<Real, 1>(0.3);}));
  std::srand(5);
  VectorX const lambda{VectorX::Random(nlp.constraints_number())};
  EXPECT_LT(lagrangian_hessian_error(nlp, nlp.guess(), 0.7, lambda), 1.0e-6);
}

TEST(DerivativeChecker, HardProblems) {
  check<Radau, 4, 1, 0>(std::make_shared<GoddardRocket>(), 4, 2, Eigen::Vector4d(1.0, 0.0, 1.0, 0.2), 1.0);
  check<Radau, 4, 1, 0>(std::make_shared<OrbitRaisingPolar>(), 4, 2, Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), 0.1);
  check<Radau, 4, 1, 0>(std::make_shared<OrbitRaisingCartesian>(), 4, 2, Eigen::Vector4d(1.0, 0.0, 0.0, 1.0), 1.6);
  check<Radau, 4, 1, 0>(std::make_shared<CartPoleAngles>(), 4, 2, Eigen::Vector4d(0.0, 0.3, 0.0, 0.0), 1.0);
  Eigen::Vector<Real, 7> cp;
  cp << 0.0, 0.5*std::sin(0.3), -0.5*std::cos(0.3), 0.0, 0.0, 0.0, 0.3*9.81*std::cos(0.3)/(2.0*0.5);
  check<Radau, 7, 1, 2>(std::make_shared<CartPoleCartesian>(), 4, 2, cp, 1.0, 0.0, "second-order", 1.0e-2);
  check<Radau, 1, 1, 0>(std::make_shared<Hypersensitive>(5.0), 4, 2, Eigen::Vector<Real, 1>(1.5), 0.2);
  check<Radau, 4, 1, 1>(std::make_shared<Brachistochrone>(), 4, 2, Eigen::Vector4d(0.0, 0.0, 0.0, 0.6), 0.8);
}

TEST(DerivativeChecker, GerdtsAlgebraicCost) {
  check<Radau, 3, 1, 0>(std::make_shared<AlgebraicCost<2>>(1.0), 4, 2, Eigen::Vector3d(1.0, 0.0, 0.0), 0.3);
  check<Radau, 3, 1, 1>(std::make_shared<AlgebraicCost<1>>(1.0), 4, 2, Eigen::Vector3d(1.0, 0.0, 0.0), 0.3);
}

TEST(DerivativeChecker, GerdtsBangSingularBang) {
  check<Radau, 3, 1, 0>(std::make_shared<BangSingularBang<2>>(), 3, 2, Eigen::Vector3d(1.0, 1.0, 0.0), 0.3);
  check<Radau, 3, 1, 1>(std::make_shared<BangSingularBang<1>>(), 3, 2, Eigen::Vector3d(1.0, 1.0, 0.0), 0.3);
}

TEST(DerivativeChecker, GerdtsPendulumGGL) {
  Eigen::Vector<Real, 6> ggl(Eigen::Vector<Real, 6>::Zero());
  ggl(0) = 1.0;
  check<Radau, 6, 1, 0>(std::make_shared<GerdtsPendulumGGL>(), 4, 2, ggl, 0.3);
}

TEST(DerivativeChecker, GerdtsSkiDriver) {
  SkiDriver<3> const ski;
  check<Radau, 9, 2, 1>(std::make_shared<SkiDriver<2>>(), 3, 1, ski.initial_states(8.5), 1.0);
  check<Radau, 9, 2, 2>(std::make_shared<SkiDriver<1>>(), 3, 1, ski.initial_states(8.5), 1.0);
}

TEST(DerivativeChecker, GerdtsKostyukova) {
  Eigen::Vector<Real, 5> kost;
  kost << 1.0, 2.0/3.0, 1.0, 2.0/3.0, 0.0;
  check<Radau, 5, 2, 0>(std::make_shared<Kostyukova<2>>(), 4, 1, kost, 0.3);
  check<Radau, 5, 2, 1>(std::make_shared<Kostyukova<1>>(), 4, 1, kost, 0.3);
}

// Badly scaled models (states up to 1e5): IPOPT's checker runs on the NLP scaled by their nominal
// magnitudes, with the default tolerance. For the landing, the forward differences of the gradients
// have a round-off floor of ~1e-6 on second derivatives of order 1e-3 (the time-scaled dynamics over
// 700 s): the Hessian is checked by the central differences in derivative_check (1e-6) instead
TEST(DerivativeChecker, GerdtsEmergencyLanding) {
  check<Radau, 7, 2, 0>(std::make_shared<EmergencyLanding>(), 3, 1, EmergencyLanding::initial_states(700.0), 0.1, 0.0,
    "first-order");
  check<Radau, 8, 1, 1>(std::make_shared<EmergencyLandingPath<2>>(), 3, 1, EmergencyLandingPath<2>::initial_states(700.0), 0.1,
    0.0, "first-order");
  check<Radau, 8, 1, 2>(std::make_shared<EmergencyLandingPath<1>>(), 3, 1, EmergencyLandingPath<1>::initial_states(700.0), 0.1,
    0.0, "first-order");
}

TEST(DerivativeChecker, GerdtsDegenerateAndUnsmoothed) {
  check<Radau, 2, 1, 0>(std::make_shared<Degenerate<2>>(), 4, 2, Eigen::Vector2d(0.0, 0.0), 0.3);
  check<Radau, 2, 1, 1>(std::make_shared<Degenerate<1>>(), 4, 2, Eigen::Vector2d(0.0, 0.0), 0.3);
  check<Radau, 4, 2, 0>(std::make_shared<KostyukovaUnsmoothed>(), 4, 1, Eigen::Vector4d(1.0, 2.0/3.0, 1.0, 0.0), 0.3);
}

TEST(DerivativeChecker, GerdtsHalfCar) {
  // The first interval [0, 0.75] contains the pothole, where the road terms are active. As for the
  // landing, the forward differences of the gradients have a round-off floor of ~1e-6 on second
  // derivatives of order 1e-4: the Hessian is checked by the central differences (1e-6)
  check<Radau, 10, 2, 2>(std::make_shared<HalfCar<2>>(), 4, 2, HalfCar<2>::initial_states(), 0.4, 0.0, "first-order");
  check<Radau, 10, 2, 4>(std::make_shared<HalfCar<1>>(), 4, 2, HalfCar<1>::initial_states(), 0.4, 0.0, "first-order");
}

TEST(DerivativeChecker, GerdtsDocking) {
  check<Radau, 10, 1, 2>(std::make_shared<Docking<2>>(), 3, 2, DockingModel::docked_states(), 0.0);
  check<Radau, 10, 1, 4>(std::make_shared<Docking<1>>(), 3, 2, DockingModel::docked_states(), 0.0);
}

TEST(DerivativeChecker, GerdtsStokes) {
  // The 3 x 3 grid (8 velocities, 4 pressures): the same code as on the finer grids, at a fraction of
  // the cost of the dense finite differences
  using Grid = StokesGrid<3>;
  Eigen::Vector<Real, Grid::N> const rest(Eigen::Vector<Real, Grid::N>::Zero());
  check<Radau, Grid::N, Grid::NV, 0>(std::make_shared<Stokes<3, 2>>(), 3, 1, rest, 0.3);
  check<Radau, Grid::N, Grid::NV, Grid::NP>(std::make_shared<Stokes<3, 1>>(), 3, 1, rest, 0.3);
}

TEST(DerivativeChecker, GerdtsChemicalReactor) {
  check<Radau, 6, 2, 0>(std::make_shared<ChemicalReactor>(), 3, 2, ChemicalReactor::initial_states(), 1.0);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
