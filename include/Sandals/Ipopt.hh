/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_IPOPT_HH
#define SANDALS_IPOPT_HH

// This header is optional and requires IPOPT (https://coin-or.github.io/Ipopt) to be available.
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>

#include <Sandals/Collocation.hh>

namespace Sandals
{

  /*\
   |   ___                  _
   |  |_ _|_ __   ___  _ __ | |_
   |   | || '_ \ / _ \| '_ \| __|
   |   | || |_) | (_) | |_) | |_
   |  |___| .__/ \___/| .__/ \__|
   |      |_|         |_|
  \*/

  /**
  * \brief IPOPT interface for the collocation transcription of an optimal control problem.
  *
  * Thin adapter exposing a Sandals::Collocation as an Ipopt::TNLP, with exact first and second
  * derivatives. Infinite bounds are mapped to IPOPT's infinity (\f$ \pm 10^{19} \f$ by default).
  * \tparam S The number of stages of the Runge-Kutta method.
  * \tparam N The dimension of the ODE/DAE system.
  * \tparam U The number of controls.
  * \tparam M The number of invariants.
  */
  template <Integer S, Integer N, Integer U, Integer M = 0>
  class IpoptCollocation : public Ipopt::TNLP
  {
  public:
    using Real = Ipopt::Number; /**< Scalar number type (IPOPT is double precision). */
    using NLP = Collocation<Real, S, N, U, M>; /**< Collocation transcription type. */
    using VectorX = typename NLP::VectorX; /**< Dynamic vector type. */

  private:
    NLP &   m_nlp;                       /**< Collocation transcription. */
    Real    m_infinity{1.0e19};          /**< Value used by IPOPT as infinity. */
    VectorX m_solution;                  /**< Primal solution. */
    VectorX m_multipliers;               /**< Constraints multipliers at the solution. */
    Real    m_objective{0.0};            /**< Objective at the solution. */
    Ipopt::SolverReturn m_status{Ipopt::UNASSIGNED}; /**< Solver status. */

    /** View an IPOPT array as an Eigen vector. */
    static Eigen::Map<VectorX const> map(Ipopt::Index const n, Real const * x) {return Eigen::Map<VectorX const>(x, n);}

    /**
    * Run an evaluation, reporting failure to IPOPT (which then backtracks) if it throws or produces
    * non-finite values.
    */
    template <typename Evaluation>
    static bool evaluate(Evaluation && evaluation)
    {
      try {
        return evaluation();
      } catch (std::exception const &) {
        return false;
      }
    }

  public:
    /**
    * Class constructor.
    * \param[in] t_nlp The collocation transcription (its guess is used as starting point).
    */
    IpoptCollocation(NLP & t_nlp) : m_nlp(t_nlp) {}

    /** \param[in] t_infinity The value used by IPOPT as infinity (option nlp_upper_bound_inf). */
    void infinity(Real const t_infinity) {this->m_infinity = t_infinity;}

    /** \return The primal solution. */
    VectorX const & solution() const {return this->m_solution;}

    /** \return The constraints multipliers at the solution. */
    VectorX const & multipliers() const {return this->m_multipliers;}

    /** \return The objective value at the solution. */
    Real objective() const {return this->m_objective;}

    /** \return The solver status. */
    Ipopt::SolverReturn status() const {return this->m_status;}

    bool get_nlp_info(Ipopt::Index & n, Ipopt::Index & m, Ipopt::Index & nnz_jac_g, Ipopt::Index & nnz_h_lag,
      IndexStyleEnum & index_style) override
    {
      n = this->m_nlp.variables_number();
      m = this->m_nlp.constraints_number();
      nnz_jac_g = this->m_nlp.jacobian_nnz();
      nnz_h_lag = this->m_nlp.hessian_nnz();
      index_style = TNLP::C_STYLE;
      return true;
    }

    bool get_bounds_info(Ipopt::Index n, Real * x_l, Real * x_u, Ipopt::Index m, Real * g_l, Real * g_u) override
    {
      VectorX z_l, z_u, c_l, c_u;
      this->m_nlp.bounds(z_l, z_u, c_l, c_u);
      Real const inf{this->m_infinity};
      auto clip = [inf] (Real const v) {return std::isinf(v) ? (v > 0 ? inf : -inf) : v;};
      for (Ipopt::Index i{0}; i < n; ++i) {x_l[i] = clip(z_l(i)); x_u[i] = clip(z_u(i));}
      for (Ipopt::Index i{0}; i < m; ++i) {g_l[i] = clip(c_l(i)); g_u[i] = clip(c_u(i));}
      return true;
    }

    bool get_starting_point(Ipopt::Index n, bool init_x, Real * x, bool init_z, Real * z_L, Real * z_U,
      Ipopt::Index m, bool init_lambda, Real * lambda) override
    {
      if (init_x) {Eigen::Map<VectorX>(x, n) = this->m_nlp.guess();}
      if (init_z) {std::fill(z_L, z_L + n, 0.0); std::fill(z_U, z_U + n, 0.0);}
      if (init_lambda) {std::fill(lambda, lambda + m, 0.0);}
      return true;
    }

    bool eval_f(Ipopt::Index n, Real const * x, bool /*new_x*/, Real & obj_value) override
    {
      return evaluate([&] {
        obj_value = this->m_nlp.objective(map(n, x));
        return std::isfinite(obj_value);
      });
    }

    bool eval_grad_f(Ipopt::Index n, Real const * x, bool /*new_x*/, Real * grad_f) override
    {
      return evaluate([&] {
        Eigen::Map<VectorX> gradient(grad_f, n);
        gradient = this->m_nlp.objective_gradient(map(n, x));
        return gradient.allFinite();
      });
    }

    bool eval_g(Ipopt::Index n, Real const * x, bool /*new_x*/, Ipopt::Index m, Real * g) override
    {
      return evaluate([&] {
        Eigen::Map<VectorX> constraints(g, m);
        constraints = this->m_nlp.constraints(map(n, x));
        return constraints.allFinite();
      });
    }

    bool eval_jac_g(Ipopt::Index n, Real const * x, bool /*new_x*/, Ipopt::Index /*m*/, Ipopt::Index nele_jac,
      Ipopt::Index * iRow, Ipopt::Index * jCol, Real * values) override
    {
      // First call: sparsity pattern only
      if (values == nullptr) {
        std::copy(this->m_nlp.jacobian_rows().begin(), this->m_nlp.jacobian_rows().end(), iRow);
        std::copy(this->m_nlp.jacobian_cols().begin(), this->m_nlp.jacobian_cols().end(), jCol);
        return true;
      }
      return evaluate([&] {
        Eigen::Map<VectorX> jacobian(values, nele_jac);
        jacobian = this->m_nlp.jacobian_values(map(n, x));
        return jacobian.allFinite();
      });
    }

    bool eval_h(Ipopt::Index n, Real const * x, bool /*new_x*/, Real obj_factor, Ipopt::Index m, Real const * lambda,
      bool /*new_lambda*/, Ipopt::Index nele_hess, Ipopt::Index * iRow, Ipopt::Index * jCol, Real * values) override
    {
      // First call: sparsity pattern only (lower triangle)
      if (values == nullptr) {
        std::copy(this->m_nlp.hessian_rows().begin(), this->m_nlp.hessian_rows().end(), iRow);
        std::copy(this->m_nlp.hessian_cols().begin(), this->m_nlp.hessian_cols().end(), jCol);
        return true;
      }
      return evaluate([&] {
        Eigen::Map<VectorX> hessian(values, nele_hess);
        hessian = this->m_nlp.hessian_values(map(n, x), obj_factor, VectorX(map(m, lambda)));
        return hessian.allFinite();
      });
    }

    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, Real const * x, Real const * /*z_L*/,
      Real const * /*z_U*/, Ipopt::Index m, Real const * /*g*/, Real const * lambda, Real obj_value,
      Ipopt::IpoptData const * /*ip_data*/, Ipopt::IpoptCalculatedQuantities * /*ip_cq*/) override
    {
      this->m_status = status;
      this->m_solution = map(n, x);
      this->m_multipliers = map(m, lambda);
      this->m_objective = obj_value;
    }

  }; // class IpoptCollocation

  /**
  * \brief Options passed to IPOPT, grouped by type.
  */
  struct IpoptOptions
  {
    std::map<std::string, std::string> strings;  /**< String options, e.g., {"linear_solver", "mumps"}. */
    std::map<std::string, double>      numerics; /**< Numeric options, e.g., {"tol", 1e-10}. */
    std::map<std::string, int>         integers; /**< Integer options, e.g., {"print_level", 0}. */
  };

  /**
  * Solve the collocation transcription of an optimal control problem with IPOPT, starting from the
  * current guess of the transcription. The default convergence tolerance is \f$ 10^{-12} \f$ (IPOPT
  * option tol), which can be overridden through the options. On return the guess is updated with the last iterate, so that
  * successive calls (e.g., on refined problems) can be warm-started.
  * \param[in,out] nlp The collocation transcription.
  * \param[in] options IPOPT options.
  * \param[out] objective The objective value at the last iterate.
  * \return The IPOPT application return status.
  */
  template <Integer S, Integer N, Integer U, Integer M>
  Ipopt::ApplicationReturnStatus ipopt_solve(Collocation<double, S, N, U, M> & nlp, IpoptOptions const & options,
    double & objective)
  {
    #define CMD "Sandals::ipopt_solve(...): "
    Ipopt::SmartPtr<IpoptCollocation<S, N, U, M>> tnlp{new IpoptCollocation<S, N, U, M>(nlp)};
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app{IpoptApplicationFactory()};
    app->RethrowNonIpoptException(false);
    // Tight default convergence tolerance, unless overridden by the options
    SANDALS_ASSERT(app->Options()->SetNumericValue("tol", 1.0e-12), CMD "invalid option tol.");
    for (auto const & [name, value] : options.strings) {
      SANDALS_ASSERT(app->Options()->SetStringValue(name, value), CMD "invalid option " << name << ".");
    }
    for (auto const & [name, value] : options.numerics) {
      SANDALS_ASSERT(app->Options()->SetNumericValue(name, value), CMD "invalid option " << name << ".");
    }
    for (auto const & [name, value] : options.integers) {
      SANDALS_ASSERT(app->Options()->SetIntegerValue(name, value), CMD "invalid option " << name << ".");
    }
    #undef CMD
    Ipopt::ApplicationReturnStatus status{app->Initialize()};
    if (status != Ipopt::Solve_Succeeded) {return status;}
    status = app->OptimizeTNLP(Ipopt::SmartPtr<Ipopt::TNLP>(GetRawPtr(tnlp)));
    objective = tnlp->objective();
    if (tnlp->solution().size() == nlp.variables_number()) {nlp.guess(tnlp->solution());}
    return status;
  }

} // namespace Sandals

#endif // SANDALS_IPOPT_HH
