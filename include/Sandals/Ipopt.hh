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
#include <type_traits>
#include <utility>

#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>

#include <Sandals/Collocation.hh>
#include <Sandals/MultipleShooting.hh>

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
  * \brief Progress of one IPOPT iteration.
  */
  struct IpoptIteration
  {
    Integer iteration;            /**< Iteration number. */
    double  objective;            /**< Objective value. */
    double  primal_infeasibility; /**< Constraints violation (unscaled). */
    double  dual_infeasibility;   /**< Optimality error (scaled). */
    double  barrier;              /**< Barrier parameter. */
    bool    restoration;          /**< True during the feasibility restoration phase. */
  };

  /**
  * \brief IPOPT interface for a sparse nonlinear program (NLP).
  *
  * Thin adapter exposing an NLP as an Ipopt::TNLP. The NLP class must provide (in double precision):
  * - variables_number(), constraints_number(), jacobian_nnz(), hessian_nnz();
  * - jacobian_rows(), jacobian_cols(), hessian_rows(), hessian_cols() (lower triangle);
  * - bounds(z_lower, z_upper, c_lower, c_upper), with infinite bounds as infinity;
  * - guess() and guess(z), the starting point and its setter;
  * - objective(z), objective_gradient(z), constraints(z), jacobian_values(z), and
  *   hessian_values(z, sigma, lambda).
  * Sandals::Collocation and Sandals::MultipleShooting are such classes. Infinite bounds are mapped
  * to IPOPT's infinity (\f$ \pm 10^{19} \f$ by default).
  * \tparam NLP The nonlinear program type.
  */
  /** True if the NLP provides scaling factors, scaled() and scaling(z_scale, c_scale). */
  template <typename NLP, typename = void>
  struct HasScaling : std::false_type {};

  template <typename NLP>
  struct HasScaling<NLP, std::void_t<decltype(std::declval<NLP const &>().scaled())>> : std::true_type {};

  template <typename NLP>
  class IpoptAdapter : public Ipopt::TNLP
  {
  public:
    using Real = Ipopt::Number; /**< Scalar number type (IPOPT is double precision). */
    using VectorX = Eigen::Vector<Real, Eigen::Dynamic>; /**< Dynamic vector type. */

  private:
    NLP &   m_nlp;                       /**< Nonlinear program. */
    Real    m_infinity{1.0e19};          /**< Value used by IPOPT as infinity. */
    VectorX m_solution;                  /**< Primal solution. */
    VectorX m_multipliers;               /**< Constraints multipliers at the solution. */
    Real    m_objective{0.0};            /**< Objective at the solution. */
    Ipopt::SolverReturn m_status{Ipopt::UNASSIGNED}; /**< Solver status. */
    std::vector<IpoptIteration> m_history;           /**< Progress of the iterations. */

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
    * \param[in] t_nlp The nonlinear program (its guess is used as starting point).
    */
    IpoptAdapter(NLP & t_nlp) : m_nlp(t_nlp) {}

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

    /** \return The progress of the iterations. */
    std::vector<IpoptIteration> const & history() const {return this->m_history;}

    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Real obj_value, Real inf_pr,
      Real inf_du, Real mu, Real /*d_norm*/, Real /*regularization_size*/, Real /*alpha_du*/, Real /*alpha_pr*/,
      Ipopt::Index /*ls_trials*/, Ipopt::IpoptData const * /*ip_data*/, Ipopt::IpoptCalculatedQuantities * /*ip_cq*/) override
    {
      this->m_history.push_back({iter, obj_value, inf_pr, inf_du, mu, mode == Ipopt::RestorationPhaseMode});
      return true;
    }

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

    bool get_scaling_parameters(Real & obj_scaling, bool & use_x_scaling, Ipopt::Index n, Real * x_scaling,
      bool & use_g_scaling, Ipopt::Index m, Real * g_scaling) override
    {
      obj_scaling = 1.0;
      use_x_scaling = use_g_scaling = false;
      if constexpr (HasScaling<NLP>::value) {
        VectorX z_scale, c_scale;
        this->m_nlp.scaling(z_scale, c_scale);
        Eigen::Map<VectorX>(x_scaling, n) = z_scale;
        Eigen::Map<VectorX>(g_scaling, m) = c_scale;
        use_x_scaling = use_g_scaling = true;
        // The objective as IPOPT's gradient-based scaling would, on the scaled variables: its largest
        // gradient component at the starting point is brought to at most one
        Real const gradient{this->m_nlp.objective_gradient(this->m_nlp.guess()).cwiseQuotient(z_scale).template lpNorm<Eigen::Infinity>()};
        if (std::isfinite(gradient) && gradient > 1.0) {obj_scaling = 1.0/gradient;}
      }
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

  }; // class IpoptAdapter

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
  * Solve a nonlinear program (e.g., the multiple shooting or collocation transcription of an optimal
  * control problem) with IPOPT, starting from the current guess of the nonlinear program. The default convergence tolerance is \f$ 10^{-12} \f$ (IPOPT
  * option tol), which can be overridden through the options. On return the guess is updated with the last iterate, so that
  * successive calls (e.g., on refined problems) can be warm-started.
  * \param[in,out] nlp The nonlinear program.
  * \param[in] options IPOPT options.
  * \param[out] objective The objective value at the last iterate.
  * \param[out] history If not null, the progress of the iterations.
  * \param[out] multipliers If not null, the constraints multipliers at the last iterate (with the
  * sign convention of IPOPT, i.e., of the Lagrangian \f$ f + \boldsymbol{\lambda}^\top \mathbf{c} \f$).
  * \return The IPOPT application return status.
  */
  template <typename NLP>
  Ipopt::ApplicationReturnStatus ipopt_solve(NLP & nlp, IpoptOptions const & options, double & objective,
    std::vector<IpoptIteration> * history = nullptr, Eigen::Vector<double, Eigen::Dynamic> * multipliers = nullptr)
  {
    #define CMD "Sandals::ipopt_solve(...): "
    Ipopt::SmartPtr<IpoptAdapter<NLP>> tnlp{new IpoptAdapter<NLP>(nlp)};
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app{IpoptApplicationFactory()};
    app->RethrowNonIpoptException(false);
    // Tight default convergence tolerance, unless overridden by the options
    SANDALS_ASSERT(app->Options()->SetNumericValue("tol", 1.0e-12), CMD "invalid option tol.");
    // The NLP's own scaling (e.g., nominal magnitudes of the states), unless chosen otherwise
    if constexpr (HasScaling<NLP>::value) {
      if (nlp.scaled() && options.strings.count("nlp_scaling_method") == 0) {
        SANDALS_ASSERT(app->Options()->SetStringValue("nlp_scaling_method", "user-scaling"), CMD "invalid option nlp_scaling_method.");
      }
    }
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
    if (history != nullptr) {*history = tnlp->history();}
    if (multipliers != nullptr) {*multipliers = tnlp->multipliers();}
    if (tnlp->solution().size() == nlp.variables_number()) {nlp.guess(tnlp->solution());}
    return status;
  }

} // namespace Sandals

#endif // SANDALS_IPOPT_HH
