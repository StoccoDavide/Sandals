/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *
 *                                                                           *
 * The Sandals project is distributed under the GNU GPLv3.                   *
 *                                                                           *
 * Davide Stocco                                           Enrico Bertolazzi *
 * University of Trento                                 University of Trento *
 * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#pragma once

#ifndef SANDALS_RUNGEKUTTA_HXX
#define SANDALS_RUNGEKUTTA_HXX

namespace Sandals {

  //! \brief Class container for the Butcher tableau of a Runge-Kutta solver.
  //! Class container for the Butcher tableau of a Runge-Kutta solver. The
  //! Butcher tableau is defined as:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   \mathbf{c} & \mathbf{A} \\ \hline
  //!              & \mathbf{b} \\
  //!              & \\mathbf{b}_{e}
  //! \end{array} \text{,}
  //! \f]
  //!
  //! where \f$ \mathbf{A} \f$ is the Runge-Kutta matrix (lower triangular matrix):
  //!
  //! \f[
  //! \mathbf{A} = \begin{bmatrix}
  //!   a_{11} & a_{12} & \dots  & a_{1s} \\
  //!   a_{21} & a_{22} & \dots  & a_{2s} \\
  //!   \vdots & \vdots & \ddots & \vdots \\
  //!   a_{s1} & a_{s2} & \dots  & a
  //! \end{bmatrix} \text{,}
  //! \f]
  //!
  //! \f$ \mathbf{b} \f$ is the Runge-Kutta weights vector relative to a method
  //! of order \f$ p \f$ (row vector):
  //!
  //! \f[
  //!   \mathbf{b} = \left[ b_1, b_2, \dots, b_s \right] \text{,}
  //! \f]
  //!
  //! \f$ \\mathbf{b}_{e} \f$ is the (optional) embedded Runge-Kutta weights
  //! vector relative to a method of order \f$ \hat{p} \f$
  //! (usually \f$ \hat{p} = p−1 \f$ or \f$ \hat{p} = p+1 \f$) (row vector):
  //!
  //! \f[
  //!   \\mathbf{b}_{e} = \left[ \hat{b}_1, \hat{b}_2, \dots, \hat{b}_s \right] \text{,}
  //! \f]
  //!
  //! and \f$ \mathbf{c} \f$ is the Runge-Kutta nodes vector (column vector):
  //!
  //! \f[
  //!   \mathbf{c} = \left[ c_1, c_2, \dots, c_s \right]^{\top} \text{.}
  //! \f]
  template <unsigned S>
  struct Tableau
  {
    using vec = Eigen::Vector<real, S>;    //!< Vector type.
    using mat = Eigen::Matrix<real, S, S>; //!< Matrix type.
    using type = enum class type : unsigned {
      ERK = 0, IRK = 1, DIRK = 2
    }; //!< Runge-Kutta type enumeration.

    std::string name;        //!< Name of the method.
    type        rk_type;     //!< Runge-Kutta type.
    unsigned    order;       //!< Order of the method.
    bool        is_embedded; //!< Embedded method boolean.
    mat         A;           //!< Matrix \f$ \mathbf{A} \f$.
    vec         b;           //!< Weights vector \f$ \mathbf{b} \f$ (row vector).
    vec         b_e;         //!< Embedded weights vector \f$ \mathbf{b}_{e} \f$ (row vector).
    vec         c;           //!< Nodes vector \f$ \mathbf{c} \f$ (column vector).
  };

  //! \brief Class container for Runge-Kutta solvers of the system of ODEs.
  //! Class container for Runge-Kutta solvers of the system of Ordinary
  //! Differential Equations (ODEs). The user must define the Butcher tableau
  //! of the solver of order \f$ p \f$ and \f$ s \f$ stages, which defined as:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   \mathbf{c} & \mathbf{A} \\ \hline
  //!              & \mathbf{b} \\
  //!              & \\mathbf{b}_{e}
  //! \end{array} \text{,}
  //! \f]
  //!
  //! where \f$ \mathbf{A} \f$ is the Runge-Kutta matrix (lower triangular matrix):
  //!
  //! \f[
  //! \mathbf{A} = \begin{bmatrix}
  //!   a_{11} & a_{12} & \dots  & a_{1s} \\
  //!   a_{21} & a_{22} & \dots  & a_{2s} \\
  //!   \vdots & \vdots & \ddots & \vdots \\
  //!   a_{s1} & a_{s2} & \dots  & a_{ss}
  //! \end{bmatrix} \text{,}
  //! \f]
  //!
  //! \f$ \mathbf{b} \f$ is the Runge-Kutta weights vector relative to a method
  //! of order \f$ p \f$ (row vector):
  //!
  //! \f[
  //!   \mathbf{b} = \left[ b_1, b_2, \dots, b_s \right] \text{,}
  //! \f]
  //!
  //! \f$ \\mathbf{b}_{e} \f$ is the (optional) embedded Runge-Kutta weights
  //! vector relative to a method of order \f$ \hat{p} \f$
  //! (usually \f$ \hat{p} = p−1 \f$ or \f$ \hat{p} = p+1 \f$) (row vector):
  //!
  //! \f[
  //!   \\mathbf{b}_{e} = \left[ \hat{b}_1, \hat{b}_2, \dots, \hat{b}_s \right] \text{,}
  //! \f]
  //!
  //! and \f$ \mathbf{c} \f$ is the Runge-Kutta nodes vector (column vector):
  //!
  //! \f[
  //!   \mathbf{c} = \left[ c_1, c_2, \dots, c_s \right]^{\top} \text{.}
  //! \f]
  //
  template <unsigned S, unsigned N>
  class RungeKutta
  {
    using vecK = Eigen::Vector<real, N*S>;      //!< Templetized vector type.
    using matK = Eigen::Matrix<real, N, S>;     //!< Templetized matrix type.
    using matJ = Eigen::Matrix<real, N*S, N*S>; //!< Templetized matrix type.
    using vecS = typename Tableau<S>::vec;      //!< Templetized vector type.
    using matS = typename Tableau<S>::mat;      //!< Templetized matrix type.
    using vecN = typename Implicit<N>::vec;     //!< Templetized vector type.
    using matN = typename Implicit<N>::mat;     //!< Templetized matrix type.

  public:
    using system   = typename Implicit<N>::ptr;              //!< Shared pointer to an implicit ODE system.
    using type     = typename Tableau<S>::type;              //!< Runge-Kutta type enumeration.
    using time     = Eigen::Vector<real, Eigen::Dynamic>;    //!< Templetized vector type for the time.
    using solution = Eigen::Matrix<real, N, Eigen::Dynamic>; //!< Templetized matrix type for the solution.

  private:
    Newton<N*S> m_newton;                //!< Newton solver for implicit methods.
    Tableau<S>  m_tableau;               //!< Butcher tableau of the Runge-Kutta method.
    system      m_system;                //!< ODE system object pointer.
    real        m_abs_tol{EPSILON_LOW};  //!< Absolute tolerance for adaptive step.
    real        m_rel_tol{EPSILON_LOW};  //!< Relative tolerance for adaptive step.
    real        m_sft_fac{0.9};          //!< Safety factor for adaptive step.
    real        m_min_sft_fac{0.2};      //!< Minimum safety factor for adaptive step.
    real        m_max_sft_fac{1.5};      //!< Maximum safety factor for adaptive step.
    real        m_d_t_min{EPSILON_HIGH}; //!< minimum step for advancing;
    bool        m_adaptive_step{false};  //!< Aadaptive step mode boolean.
    bool        m_verbose{false};        //!< Verbose mode boolean.
    unsigned    m_order;                 //!< Order of the solver.
    bool        m_is_embedded;           //!< Boolean to check if the method is embedded.

  public:
    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    //! \param t_tableau The Tableau reference.
    RungeKutta(Tableau<S> const &t_tableau) : m_tableau(t_tableau) {
      //SANDALS_WARNING(
      //  "RungeKutta::RungeKutta(...): the order of the method " << t_tableau.name << " is " <<
      //  this->tableau_order(this->m_tableau.A, this->m_tableau.b, this->m_tableau.c));
    }

    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    //! \param t_tableau The Tableau reference.
    //! \param t_system The ODE system shared pointer.
    RungeKutta(Tableau<S> const &t_tableau, system t_system)
      : m_tableau(t_tableau), m_system(t_system) {
      //SANDALS_WARNING(
      //  "RungeKutta::RungeKutta(...): the order of the method " << t_tableau.name << " is " <<
      //  this->tableau_order(this->m_tableau.A, this->m_tableau.b, this->m_tableau.c));
    }

    //! Check if the solver is explicit, implicit or diagonal implicit
    //! \return The type of the solver.
    type rk_type(void) const {return this->m_tableau.rk_type;}

    //! Cheack if the solver is explicit.
    //! \return True if the solver is explicit, false otherwise.
    bool is_explicit(void) const {return this->m_tableau.rk_type == type::ERK;}

    //! Cheack if the solver is implicit.
    //! \return True if the solver is implicit, false otherwise.
    bool is_implicit(void) const {return this->m_tableau.rk_type == type::IRK;}

    //! Cheack if the solver is diagonal implicit.
    //! \return True if the solver is diagonal implicit, false otherwise.
    bool is_diagonal_implicit(void) const {return this->m_tableau.rk_type == type::DIRK;}

    //! Get the Tableau reference.
    //! \return The Tableau reference.
    Tableau<S> & tableau(void) {return this->m_tableau;}

    //! Get the Tableau const reference.
    //! \return The Tableau const reference.
    Tableau<S> const & tableau(void) const {return this->m_tableau;}

    //! Get the stages number of the method.
    //! \return The stages number of the method.
    unsigned stages(void) const {return S;}

    //! Get the name of the method.
    //! \return The name of the method.
    std::string name(void) const {return this->m_tableau.name;}

    //! Get the order of the method.
    //! \return The order of the method.
    unsigned order(void) const {return this->m_tableau.order;}

    //! Check if the method is embedded.
    //! \return True if the method is embedded, false otherwise.
    bool is_embedded() const {return this->m_tableau.is_embedded;}

    //! Get the matrix \f$ \mathbf{A} \f$.
    //! \return The matrix \f$ \mathbf{A} \f$.
    matS A(void) const {return this->m_tableau.A;}

    //! Get the weights vector \f$ \mathbf{b} \f$.
    //! \return The weights vector \f$ \mathbf{b} \f$.
    vecS b(void) const {return this->m_tableau.b;}

    //! Get the embedded weights vector \f$ \\mathbf{b}_{e} \f$.
    //! \return The embedded weights vector \f$ \\mathbf{b}_{e} \f$.
    vecS b_embedded(void) const {return this->m_tableau.b_e;}

    //! Get the nodes vector \f$ \mathbf{c} \f$.
    //! \return The nodes vector \f$ \mathbf{c} \f$.
    vecS c(void) const {return c;}

    //! Get the ODE system shared pointer.
    //! \return The ODE system shared pointer.
    system ode(void) {return this->m_system;}

    //! Set the ODE system shared pointer.
    //! \param t_system The ODE system shared pointer.
    void ode(system t_system) {this->m_system = t_system;}

    //! Get the adaptive step absolute tolerance.
    //! \return The adaptive step absolute tolerance.
    real absolute_tolerance(void) {return this->m_abs_tol;}

    //! Get the adaptive step absolute tolerance reference.
    //! \param t_abs_tol The adaptive step absolute tolerance.
    void absolute_tolerance(real t_abs_tol) {this->m_abs_tol = t_abs_tol;}

    //! Get the adaptive step relative tolerance.
    //! \return The adaptive step relative tolerance.
    real relative_tolerance(void) {return this->m_rel_tol;}

    //! Get the adaptive step relative tolerance reference.
    //! \param t_rel_tol The adaptive step relative tolerance.
    void relative_tolerance(real t_rel_tol) {this->m_rel_tol = t_rel_tol;}

    //! Get the safety factor for adaptive step.
    //! \return The safety factor for adaptive step.
    real & safety_factor(void) {return this->m_sft_fac;}

    //! Set safety factor for adaptive step.
    //! \param t_sft_fac The safety factor for adaptive step.
    void safety_factor(real t_sft_fac) {this->m_sft_fac = t_sft_fac;}

    //! Get the minimum safety factor for adaptive step.
    //! \return The minimum safety factor for adaptive step.
    real & min_safety_factor(void) {return this->m_min_sft_fac;}

    //! Set the minimum safety factor for adaptive step.
    //! \param t_min_sft_fac The minimum safety factor for adaptive step.
    void min_safety_factor(real t_min_sft_fac) {this->m_min_sft_fac = t_min_sft_fac;}

    //! Get the maximum safety factor for adaptive step.
    //! \return The maximum safety factor for adaptive step.
    real & max_safety_factor(void) {return this->m_max_sft_fac;}

    //! Set the maximum safety factor for adaptive step.
    //! \param t_max_sft_fac The maximum safety factor for adaptive step.
    void max_safety_factor(real t_max_sft_fac) {this->m_max_sft_fac = t_max_sft_fac;}

    //! Enable verbose mode.
    void enable_verbose_mode(void) {this->m_verbose = true;}

    //! Disable verbose mode.
    void disable_verbose_mode(void) {this->m_verbose = false;}

    //! Enable adaptive step mode.
    void enable_adaptive_step(void) {this->m_adaptive_step = true;}

    //! Disable adaptive step mode.
    void disable_adaptive_step(void) {this->m_adaptive_step = false;}

    //! Compute adaptive time step for the next advancing step according to the
    //! error control method. The error control method used is the local truncation
    //! error method, which is based on the following formula:
    //!
    //! \f[
    //! e = \sqrt{\dfrac{1}{n} \displaystyle\sum_{i=1}{n}\left(\dfrac
    //!   {\mathbf{x} - \hat{\mathbf{x}}}
    //!   {s c_i}
    //! \right)^2}
    //! \f]
    //!
    //! where \f$ \mathbf{x} \f$ is the approximation of the states at computed
    //! with higher order method of \f$ p \f$, and \f$ \hat{\mathbf{x}} \f$ is the
    //! approximation of the states at computed with lower order method of \f$
    //! \hat{p} \f$. To compute the suggested time step for the next advancing step
    //! \f$ \Delta t_{k+1} \f$, The error is compared to \f$ 1 \f$ in order to find
    //! an optimal step size. From the error behaviour \f$ e \approx Ch^{q+1} \f$
    //! and from \f$ 1 \approx Ch_{opt}^{q+1} \f$ (where \f$ q = \min(p,\hat{p}) \f$)
    //! the optimal step size is obtained as:
    //!
    //! \f[
    //! h_{opt} = h \left( \dfrac{1}{e} \right)^{\frac{1}{q+1}}
    //! \f]
    //!
    //! We multiply the previous quation by a safety factor \f$ f \f$, usually
    //! \f$ f = 0.8 \f$, \f$ 0.9 \f$, \f$ (0.25)^{1/(q+1)} \f$, or \f$ (0.38)^{1/(q+1)} \f$,
    //! so that the error will be acceptable the next time with high probability.
    //! Further, \f$ h \f$ is not allowed to increase nor to decrease too fast.
    //! So we put:
    //!
    //! \f[
    //! h_{new} = h \min\left(f_{max}, \max\left(f_{max}, f \left(\dfrac{1}{e}\right)^{\frac{1}{q+1}}
    //! \right) \right)
    //! \f]
    //!
    //! for the new step size. Then, if \f$ e \leq 1 \f$, the computed step is
    //! accepted and the solution is advanced to \f$ \mathbf{x} \f$ and a new step
    //! is tried with \f$ h_{new} \f$ as step size. Else, the step is rejected
    //! and the computations are repeated with the new step size \f$ h_{new} \f$.
    //! Typially, \f$ f \f$ is set in the interval \f$ [0.8, 0.9] \f$,
    //! \f$ f_{max} \f$ is set in the interval \f$ [1.5, 5] \f$, and \f$ f_{min} \f$
    //! is set in the interval \f$ [0.1, 0.2] \f$.
    //!
    //! \param x_h Approximation of the states at \f$ k+1 \f$-th time step \f$
    //!            \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$ with higher order method.
    //! \param x_l Approximation of the states at \f$ k+1 \f$-th time step \f$
    //!            \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$ with lower order method.
    //! \param d_t Actual advancing time step \f$ \Delta t\f$.
    //! \return The suggested time step for the next advancing step \f$ \Delta
    //!         t_{k+1} \f$.
    real estimate_step(vecN const &x_h, vecN const &x_l, real d_t) const {
      return d_t * std::min(this->m_max_sft_fac, std::max(this->m_min_sft_fac,
        this->m_sft_fac * ((x_h - x_l) / (
        this->m_abs_tol + this->m_rel_tol*std::max(x_h.cwiseAbs().maxCoeff(), x_l.cwiseAbs().maxCoeff())
        )).cwiseAbs().maxCoeff()));
    }

    //! Print the Runge-Kutta method information.
    //! \param os The output stream.
    void info(std::ostream &os) const {
      os
        << "Runge-Kutta method:\t" << this->name() << std::endl
        << "\t- order:\t" << this->order() << std::endl
        << "\t- stages:\t" << this->stages() << std::endl
        << "\t- type:\t" << this->rk_type() << std::endl
        << "\t- embedded:\t" << this->is_embedded() << std::endl;
    }

    //! Compute an integration step using the explicit Runge-Kutta method for a
    //! system of the form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t)
    //! = \mathbf{0} \f$. Notice that if the Runge-Kutta method and the ODE
    //! system are both explicit then the \f$ \mathbf{K}_i \f$ variables can be
    //! straightforwardly calculated.
    //!
    //! **Solution Algorithm**
    //!
    //! Consider a Runge-Kutta method, written for a system of the
    //! form \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, t) \f$:
    //!
    //! \f[
    //! \begin{array}{l}
    //! \mathbf{K}_i = \mathbf{f} \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j,
    //!   \, t_k + c_i \Delta t
    //!   \right), \qquad i = 1, 2, \ldots, s \\
    //! \mathbf{x}_{k+1} = \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s b_j
    //! \mathbf{K}_j \, ,
    //! \end{array}
    //! \f]
    //!
    //! Beacuse of the nature of the matrix \f$ \mathbf{A} \f$ (lower triangular)
    //! the \f$ s\f$ stages for a generic explicit Runge-Kutta method take the
    //! form:
    //!
    //! \f[
    //! \mathbf{K}_i = \mathbf{f} \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{i-1} a_{ij}
    //!   \mathbf{K}_j, \, t_k + c_i \Delta t
    //!   \right), \qquad i = 1, 2, \ldots, s.
    //! \f]
    //!
    //! Then the explicit Runge-Kutta method for an implicit system of the form
    //! \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{0} \f$ can be
    //! written as:
    //!
    //! \f[
    //! \begin{array}{l}
    //! \mathbf{F}_i \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{i-1} a_{ij}
    //!     \mathbf{K}_j, \, \mathbf{K}_i, \, t_k + c_i \Delta t
    //! \right) = \mathbf{0}, \qquad i = 1, 2, \ldots, s \\
    //! \mathbf{x}_{k+1} = \mathbf{x}_k + \displaystyle\sum_{j=1}^s b_j \mathbf{K}_j.
    //! \end{array}
    //! \f]
    //!
    //! It is important to notice that the system of \f$ s \f$ equations
    //! \f$ \mathbf{F}_i \f$ is a triangular system (which may be non-linear in
    //! the \f$ \mathbf{K}_i \f$ variables), so it can be solved using forward
    //! substitution and the solution of the system is the vector \f$ \mathbf{K}
    //! \f$. Thus, the final system to be solved is the following:
    //!
    //! \f[
    //! \left\{\begin{array}{l}
    //! \mathbf{F}_1 \left(
    //!   \mathbf{x}_k, \, \mathbf{K}_1, \, t_k + c_1 \Delta t
    //! \right) = \mathbf{0} \\
    //! \mathbf{F}_2 \left(
    //!   \mathbf{x}_k + \Delta t \, a_{21} \mathbf{K}_1, \,
    //!   \mathbf{K}_2, \, t_k + c_2 \Delta t
    //! \right) = \mathbf{0} \\
    //! ~~ \vdots \\
    //! \mathbf{F}_s \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{s-1} a_{sj}
    //!   \mathbf{K}_j, \, \mathbf{K}_s, \, t_k + c_s \Delta t
    //! \right) = \mathbf{0}
    //! \end{array}\right.
    //! \f]
    //!
    //! The \f$ \mathbf{K}_i \f$ variable are computed using the Newton's method.
    //!
    //! **Note**
    //!
    //! Another approach is to directly solve the whole system of equations by
    //! Newton'smethod. In other words, the system of equations is solved
    //! iteratively by computing the Jacobian matrixes of the system and using
    //! them to compute the solution. This approach is used in the implicit
    //! Runge-Kutta method. For this reason, a Butcher tableau relative to an
    //! explicit Runge-Kutta method can also be used in the `ImplicitRungeKutta` class.
    //!
    //! The suggested time step for the next advancing step \f$ \Delta t_{k+1} \f$,
    //! is the same as the input time step \f$ \Delta t \f$ since in the explicit
    //! Runge-Kutta method the time step is not modified through any error control
    //! method.
    //!
    //! \param x_old States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param t_old Time step \f$ t_k \f$.
    //! \param t_step Advancing time step \f$ \Delta t\f$.
    //! \param x_new The approximation of the states at \f$ k+1 \f$-th time step
    //!            \f$ \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$
    //! \param t_est The suggested time step for the next advancing step
    //!                 \f$ \Delta t_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool explicit_step(vecN const &x_old, real t_old, real t_step, vecN &x_new, real &t_est)
    {
      using Eigen::all;
      using Eigen::seqN;

      // Compute the K variables in the case of an explicit method and explicit system
      real t_node;
      vecN x_node;
      matK K;
      for (unsigned i = 0; i < S; ++i)
      {
        // FIXME: optimize this loop
        t_node = t_old + this->m_tableau.c(i) * t_step;
        x_node = x_old + K(all, seqN(0, i)) * this->m_tableau.A(i, seqN(0, i)).transpose();
        K.col(i) = t_step * static_cast<Explicit<N> const *>(this->m_system.get())->f(x_node, t_node);
        std::cout << "K[" << i << "] = " << K.col(i).transpose() << std::endl;
      }
      if (!K.allFinite()) {return false;}

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_is_embedded) {
        // FIXME: optimize this loop
        vecN x_e = x_old + K * this->m_tableau.b_e;
        t_est = this->estimate_step(x_new, x_e, t_step);
      }
      return true;
    }

    //! Compute the residual of system to be solved:
    //!
    //! \f[
    //! \mathbf{F}_i\left(\mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{s}
    //!   a_{ij} \mathbf{K}_j, \, \mathbf{K}_i, \, t_k + c_i \Delta t
    //! \right) = \mathbf{0}.
    //! \f]
    //!
    //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param K_k Variables \f$ \mathbf{K} \f$ of the system to be solved.
    //! \param t_k Time step \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t \f$.
    //! \param fun The residual of system to be solved.
    void implicit_function(vecN const &x_k, vecK const &K_k, real t_k, real d_t, vecK &fun) const
    {
      // Loop through each equation of the system
      real t_i;
      vecN x_i;
      matK K{K_k.reshaped(N, S)}; // FIXME: can i avoid reshaping?
      matK fun_tmp;
      for (unsigned i = 0; i < S; ++i) {
        // FIXME: optimize this loop
        t_i = t_k + this->m_tableau.c(i) * d_t;
        x_i = x_k + K * this->m_tableau.A.row(i).transpose();
        fun_tmp.col(i) = this->m_system->F(x_i, K.col(i) / d_t, t_i);
      }
      fun = fun_tmp.reshaped(N*S, 1);
    }

    //! Compute the Jacobian of the system of equations:
    //!
    //! \f[
    //! \mathbf{F}_i\left(\mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s
    //!   a_{ij} \mathbf{K}_j, \, \mathbf{K}_i, \, t_k + c_i \Delta t
    //! \right) = \mathbf{0}
    //! \f]
    //!
    //! to be solved in the \f$ \mathbf{K} \f$ variable:
    //!
    //! \f[
    //! \dfrac{\partial \mathbf{F}_i}{\partial \mathbf{K}_i} \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s a_{ij} \mathbf{K}_j,
    //!   \, \mathbf{K}_i, \, t_k + c_i \Delta t
    //! \right)
    //! \f]
    //!
    //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param K_k Variables \f$ \mathbf{K} \f$ of the system to be solved.
    //! \param t_k Time step \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t \f$.
    //! \param jac The Jacobian of the system of equations to be solved.
    void implicit_jacobian(vecN const &x_k, vecK const &K_k, real t_k, real d_t, matJ &jac) const
    {
      using Eigen::all;
      using Eigen::seqN;

      // Loop through each equation of the system
      matK K{K_k.reshaped(N, S)};
      vecN x_dot_i;
      matN JF_x, JF_x_dot;
      auto idx = seqN(0, N-1);
      auto jdx = idx;
      real t_i;
      vecN x_i;
      for (unsigned i = 0; i < S; ++i) {
        t_i = t_k + this->m_tableau.c(i) * d_t;
        x_i = x_k + K * this->m_tableau.A.row(i).transpose();

        // Compute the Jacobians with respect to x and x_dot
        x_dot_i  = K.col(i) / d_t;
        JF_x     = this->m_system->JF_x(x_i, x_dot_i, t_i);
        JF_x_dot = this->m_system->JF_x_dot(x_i, x_dot_i, t_i);

        // Derivative of F(x_i, K(:,i)/d_t, t_i)
        jdx = seqN(0, N-1);
        for (unsigned j = 0; j < S; ++j) {
          // Combine the Jacobians with respect to x and x_dot to obtain the
          // Jacobian with respect to K
          if (i == j) {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x + JF_x_dot / d_t;
          } else {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x;
          }
          jdx = seqN(j*N, (j+1)*N-1);
        }
        idx = seqN(i*N, (i+1)*N-1);
      }
    }

    //! Compute an integration step using the implicit Runge-Kutta method for a
    //! system of the form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) =
    //! \mathbf{0} \f$.
    //!
    //! **Solution Algorithm**
    //!
    //! Consider a Runge-Kutta method, written for a system of the form
    //! \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, \mathbf{v}, t) \f$:
    //!
    //! \f[
    //! \begin{array}{l}
    //! \mathbf{K}_i = \mathbf{f} \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{s} a_{ij}
    //!   \mathbf{K}_j,
    //!   \, t_k + c_i \Delta t
    //!   \right), \qquad i = 1, 2, \ldots, s \\
    //! \mathbf{x}_{k+1} = \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s b_j
    //! \mathbf{K}_j \, ,
    //! \end{array}
    //! \f]
    //!
    //! Then the implicit Runge-Kutta method for an implicit system of the form
    //! \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{0} \f$ can be written
    //! as:
    //!
    //! \f[
    //! \begin{array}{l}
    //! \mathbf{F}_i \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s a_{ij}
    //!     \mathbf{K}_j, \, \mathbf{K}_i, \, t_k + c_i \Delta t
    //! \right) = \mathbf{0}, \qquad i = 1, 2, \ldots, s \\
    //! \mathbf{x}_{k+1} = \mathbf{x}_k + \displaystyle\sum_{j=1}^s b_j \mathbf{K}_j.
    //! \end{array}
    //! \f]
    //!
    //! Thus, the final system to be solved is the following:
    //!
    //! \f[
    //! \left\{\begin{array}{l}
    //! \mathbf{F}_1 \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s a_{1j}
    //!   \mathbf{K}_j, \, \mathbf{K}_1, \, t_k + c_1 \Delta t
    //! \right) = \mathbf{0} \\
    //! \mathbf{F}_2 \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s a_{2j}
    //!   \mathbf{K}_j, \, \mathbf{K}_2, \, t_k + c_2 \Delta t
    //! \right) = \mathbf{0} \\
    //! ~~ \vdots \\
    //! \mathbf{F}_s \left(
    //!   \mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s a_{sj}
    //!   \mathbf{K}_j, \, \mathbf{K}_s, \, t_k + c_s \Delta t
    //! \right) = \mathbf{0}
    //! \end{array}\right.
    //! \f]
    //!
    //! The \f$ \mathbf{K} \f$ variables are computed using the Newton's method.
    //!
    //! The suggested time step for the next advancing step
    //! \f$ \Delta t_{k+1} \f$, is the same as the input time step
    //! \f$ \Delta t \f$ since in the implicit Runge-Kutta method the time step
    //! is not modified through any error control method.
    //!
    //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param t_k Time step \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t\f$.
    //! \param x_o The approximation of the states at \f$ k+1 \f$-th time step
    //!            \f$ \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$
    //! \param d_t_star The suggested time step for the next advancing step
    //!                 \f$ \Delta t_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool implicit_step(vecN const &x_k, real t_k, real d_t, vecN &x_o, real &d_t_star)
    {
      std::cout << "implicit_step: ";
      // Check if the solver converged
      vecK K;
      if (!this->m_newton.solve(
        [this, &x_k, t_k, d_t](vecK const &K, vecK &fun) {this->implicit_function(x_k, K, t_k, d_t, fun);},
        [this, &x_k, t_k, d_t](vecK const &K, matJ &jac) {this->implicit_jacobian(x_k, K, t_k, d_t, jac);},
        vecK::Zero(), K)) {return false;}

      // Perform the step and obtain x_k+1
      x_o = x_k + K.reshaped(N, S) * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_tableau.is_embedded) {
        vecN x_e{x_k + K.reshaped(N, S) * this->m_tableau.b_e};
        d_t_star = this->estimate_step(x_o, x_e, d_t);
      }
      std::cout << "DONE" << std::endl;
      return true;
    }

    //! Compute a step using a generic integration method for a system of the
    //! form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, \mathbf{v}, t) = \mathbf{0}
    //! \f$. The step is based on the following formula:
    //!
    //! \f[
    //! \mathbf{x}_{k+1}(t_{k}+\Delta t) = \mathbf{x}_k(t_{k}) +
    //! \mathcal{S}(\mathbf{x}_k(t_k), \mathbf{x}^{\prime}_k(t_k), t_k, \Delta t)
    //! \f]
    //!
    //! where \f$ \mathcal{S} \f$ is the generic advancing step of the solver.
    //!
    //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param t_k Time step \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t\f$.
    //! \param x_o The approximation of the states at \f$ k+1 \f$-th time step
    //!              \f$ \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$
    //! \param d_t_star The suggested time step for the next advancing step
    //!                \f$ \Delta t_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool step(vecN const &x_k, real t_k, real d_t, vecN &x_o, real &d_t_star)
    {
      if (this->is_explicit() && this->m_system->is_explicit()) {
        return this->explicit_step(x_k, t_k, d_t, x_o, d_t_star);
      } else {
        return this->implicit_step(x_k, t_k, d_t, x_o, d_t_star);
      }
    }

    //! Advance using a generic integration method for a system of the form
    //! \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, \mathbf{v}, t) = \mathbf{0} \f$.
    //! The step is based on the following formula:
    //!
    //! \f[
    //! \mathbf{x}_{k+1}(t_{k}+\Delta t) = \mathbf{x}_k(t_{k}) +
    //! \mathcal{S}(\mathbf{x}_k(t_k), \mathbf{x}^{\prime}_k(t_k), t_k, \Delta t)
    //! \f]
    //!
    //! where \f$ \mathcal{S} \f$ is the generic advancing step of the solver.
    //! In the advvancing step, the system solution is also projected on the
    //! manifold \f$ \mathcal{h}(\mathbf{x}, \mathbf{v}, t) \f$. Substepping is
    //! also used to ensure that the solution is accurate.
    //!
    //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param t_k Time step \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t\f$.
    //! \param x_new The approximation of the states at \f$ k+1 \f$-th time step
    //!              \f$ \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$
    //! \param d_t_star The suggested time step for the next advancing step
    //!                 \f$ \Delta t_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool advance(vecN const &x_k, real t_k, real d_t, vecN &x_new, real &d_t_star)
    {
      #define CMD "Indigo.RungeKutta.advance(...): "

      // Check step size
      SANDALS_ASSERT(d_t > real(0.0), CMD "in " << this->m_tableau.name <<
        "solver, d_t = "<< d_t << ", expected > 0.");

      // If the integration step failed, try again with substepping
      if (!this->step(x_k, t_k, d_t, x_new, d_t_star))
      {
        vecN x_tmp   = x_k;
        real t_tmp   = t_k;
        real d_t_tmp = d_t / real(2.0);

        // Substepping logic
        unsigned max_substeps = 10;
        unsigned max_k = max_substeps * max_substeps;
        unsigned k = 2;
        real d_t_star_tmp;
        while (k > 0)
        {
          // Calculate the next time step with substepping logic
          if (this->step(x_tmp, t_tmp, d_t_tmp, x_new, d_t_star_tmp)) {

            // Accept the step
            d_t_tmp = d_t_star_tmp;

            // If substepping is enabled, double the step size
            if (k > 0 && k < max_k) {
              k = k - 1;
              // If the substepping index is even, double the step size
              if (k % 2 == 0) {
                d_t_tmp = real(2.0) * d_t_tmp;
                if (this->m_verbose) {
                  SANDALS_WARNING(CMD "in " << this->m_tableau.name << " solver, at t = " <<
                  t_tmp << ", integration succedded disable one substepping layer.");
                }
              }
            }

            // Check the infinity norm of the solution
            SANDALS_ASSERT(std::isfinite(x_tmp.maxCoeff()), CMD "in " <<
            this->m_tableau.name << " solver, at t = " << t_tmp <<
            ", ||x||_inf = inf, computation interrupted.");

          } else {

            // If the substepping index is too high, abort the integration
            k = k + 2;
            SANDALS_ASSERT(k < max_k, CMD "in " << this->m_tableau.name <<
              " solver, at t = " << t_tmp << ", integration failed with d_t = " <<
              d_t << ", aborting.");
            return false;

            // Otherwise, try again with a smaller step
            if (this->m_verbose) {
              SANDALS_WARNING(CMD "in " << this->m_tableau.name << " solver, at t = " <<
                t_tmp << ", integration failed, adding substepping layer.");
            }
            d_t_tmp = d_t_tmp / real(2.0);
            continue;

          }

          // Store time solution
          t_tmp = t_tmp + d_t_tmp;
        }

        // Store output states substepping solutions
        x_new    = x_tmp;
        d_t_star = d_t_tmp;
        return true;
      } else {
        return true;
      }
      #undef CMD
    }

    //! Solve the system and calculate the approximate solution over the mesh of
    //! time points.
    //!
    //! \param t Time mesh points \f$ \mathbf{t} = \left[ t_0, t_1, \ldots, t_n
    //!          \right]^T \f$.
    //! \param ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    //! \param sol A matrix \f$ \left[(\mathrm{size}(\mathbf{x}) \times \mathrm{size}
    //!            (\mathbf{t})\right] \f$ containing the approximated solution
    //!            over the mesh of time points.
    //! \return True if the system is successfully solved, false otherwise.
    bool solve(vecD const &t, vecN const &ics, solution &sol)
    {
      using Eigen::last;

      // Store first step
      sol.resize(ics.size(), t.size());
      sol.col(0) = ics;
      unsigned s{0};

      // Update the current step
      vecN x_s{ics};
      real t_s{t(0)};
      real d_t_s{t(1) - t(0)};
      real d_t_tmp{d_t_s}, d_t_star;
      bool mesh_point_bool, saturation_bool;

      while (true) {
        // Integrate system
        if (!this->advance(sol.col(s), t_s, d_t_s, x_s, d_t_star)) {return false;}

        // Update the current step
        t_s = t_s + d_t_s;

        // Saturate the suggested timestep
        mesh_point_bool = std::abs(t_s - t(s+1)) < SQRT_EPSILON;
        saturation_bool = t_s + d_t_star > t(s+1) + SQRT_EPSILON;
        if (this->m_adaptive_step && !mesh_point_bool && saturation_bool) {
          d_t_tmp = d_t_star;
          d_t_s   = t(s+1) - t_s;
        } else {
          d_t_s = d_t_star;
        }

        // Store solution if the step is a mesh point
        if (!this->m_adaptive_step || mesh_point_bool) {
          // Update temporaries
          s     += 1;
          d_t_s = d_t_tmp;

          // Update outputs
          sol.col(s) = x_s;

          // Check if the current step is the last one
          if (std::abs(t_s - t(last)) < SQRT_EPSILON) {break;}
        }
      }
      return true;
    }

    //! Solve the system and calculate the approximate solution over the
    //! suggested mesh of time points with adaptive step size control.
    //! \param t Time mesh points \f$ \mathbf{t} = \left[ t_0, t_1, \ldots, t_n
    //!          \right]^T \f$.
    //! \param ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    //! \param sol A matrix \f$ \left[(\mathrm{size}(\mathbf{x}) \times \mathrm{size}
    //!           (\mathbf{t})\right] \f$ containing the approximated solution
    //!           over the mesh of time points.
    //! \return True if the system is successfully solved, false otherwise.
    bool adaptive_solve(vecD t, vecN const &ics, solution &sol) const
    {
      using Eigen::last;

      #define CMD "Indigo.RungeKutta.adaptive_solve(...): "

      // Collect optional arguments
      real d_t{t(1) - t(0)}, d_t_star;
      real scale{100.0};
      real t_min{std::max(this->m_d_t_min, d_t/scale)};
      real t_max{scale*d_t};

      SANDALS_ASSERT(t_max > t_min && t_min > real(0.0), CMD "invalid time bounds detected.");
      d_t = std::max(std::min(d_t, t_max), t_min);

      // Instantiate output
      unsigned safety_length = std::ceil(real(1.5)/this->m_min_sft_fac) * t.size();
      vecN t_o(safety_length);

      // Store first step
      t_o(0)     = t(0);
      sol.col(0) = ics;

      // Instantiate temporary variables
      unsigned s{0}; // Current step

      while (true) {
        // Integrate system
        this->advance(sol.col(s), t_o(s), d_t, sol.col(s+1), d_t_star);

        // Saturate the suggested timestep
        d_t = std::max(std::min(d_t_star, t_max), t_min);

        // Store solution
        t_o(s+1) = t_o(s) + d_t;

        // Check if the current step is the last one
        if (t_o(s+1) + d_t > t(last)) {break;}

        // Update steps counter
        s += 1;
      }

      // Resize the output
      t_o.conservativeResize(s-1);
      sol.conservativeResize(Eigen::NoChange, s-1);

      return true;

      #undef CMD
    }

    //! Check Butcher tableau consistency for an explicit Runge-Kutta method.
    //! \param tbl.A   Matrix \f$ \mathbf{A} \f$.
    //! \param tbl.b   Weights vector \f$ \mathbf{b} \f$.
    //! \param tbl.b_e [optional] Embedded weights vector \f$ \\mathbf{b}_{e} \f$.
    //! \param tbl.c   Nodes vector \f$ \mathbf{c} \f$.
    //! \return True if the Butcher tableau is consistent, false otherwise.
    //[out,order,e_order] = check_tableau( this, tbl )

    //! Check the order of a Runge-Kutta tableau according to the conditions taken from:
    //! *A family of embedded Runge-Kutta formulae*, J. R. Dormand and P. J. Prince,
    //! Journal of Computational and Applied Mathematics, volume 6(1), 1980.
    //! Doi: [10.1016/0771-0509(80)90013-3](https://doi.org/10.1016/0771-0509(80)90013-3)
    //! \param tableau The Runge-Kutta tableau to be checked.
    //! \return The order of the Runge-Kutta tableau.
    /*unsigned tableau_order(matS const &A, vecS const &b, vecS const &c) const
    {
      #define CMD "Indigo.RungeKutta.tableau_order(...): "

      // Temporary variables initialization
      real tol{std::pow(EPSILON, real(2.0/3.0))};
      vecN one{vecN::Ones()};
      vecS Ac{A*c};
      matS bA{(b*A).transpose()};
      vecS err{A*one - c};
      unsigned order{0};
      std::string msg{""};

      // Check consistency
      SANDALS_ASSERT(std::max(std::abs(err)) > tol,
        CMD "consistency precheck failed, A*[1] - c = " << err << " ≠ 0.");

      // Check order 1
      if (std::abs(b.sum() - 1) > tol) {
        SANDALS_WARNING(CMD "order 1 check failed, found sum(b) == " << b.sum() << " ≠ 1.");
        return order;
      }

      order = 1; // Order 1 is the highest order that can be checked

      // Check order 2
      auto bc{b.array() * c.array()};
      if (std::abs(bc.sum() - 1/2) > tol) {
        SANDALS_WARNING(CMD "order 2 check failed, sum(b*c) = " << bc.sum() << " ≠ 1/2.");
        return order;
      }

      order = 2; // Order 2 is the highest order that can be checked

      // Check order 3
      auto bc2{b.array() * (c.pow(2)).array()};
      if (std::abs(bc2.sum() - 1/3) > tol) {
        SANDALS_WARNING(CMD "order 3 check failed, sum(b*c^2) = " << bc2.sum() << " ≠ 1/3.");
        return order;
      }

      auto bAc{b.array() * Ac.array()};
      if (std::abs(bAc.sum() - 1/6) > tol) {
        SANDALS_WARNING(CMD "order 3 check failed, sum(b*d) = " << bAc.sum() << " ≠ 1/6.");
        return order;
      }

      order = 3; // Order 3 is the highest order that can be checked

      // Check order 4
      auto bc3{b.array() * (c.pow(3)).array()};
      if (std::abs(bc3.sum() - 1/4) > tol) {
        SANDALS_WARNING(CMD "order 4 check failed, sum(b*c^3) = " << bc3.sum() << " ≠ 1/4.");
        return order;
      }

      auto cAc{c.array() * Ac.array()};
      vecS bcAc{b.cwiseProduct(cAc)};
      if (std::abs(bcAc.sum() - 1/8) > tol) {
        SANDALS_WARNING(CMD "order 4 check failed, sum(b*c*A*c) = " << bcAc.sum() << " ≠ 1/8.");
        return order;
      }

      auto bAc2{bA.array() * (c.pow(2)).array()};
      if (std::abs(bAc2.sum() - 1/12) > tol) {
        SANDALS_WARNING(CMD "order 4 check failed, sum(b*A*c^2) = " << bAc2.sum() << " ≠ 1/12.");
        return order;
      }

      auto bAAc{bA.array() * Ac.array()};
      if (std::abs(bAAc.sum() - 1/24) > tol) {
        SANDALS_WARNING(CMD "order 4 check failed, sum(b*A*A*c) = " << bAAc.sum() << " ≠ 1/24.");
        return order;
      }

      order = 4; // Order 4 is the highest order that can be checked

      // Check order 5
      auto bc4{b.array() * (c.pow(4)).array()};
      if (std::abs(bc4.sum() - 1/5) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*c^4) = " << bc4.sum() << " ≠ 1/5.");
        return order;
      }

      auto bc2Ac{bc2.array() * Ac.array()};
      if (std::abs(bc2Ac.sum() - 1/10) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*c^2*A*c) = " << bc2Ac.sum() << " ≠ 1/10.");
        return order;
      }

      vecS bAcAc{(b.cwiseProduct(Ac)).cwiseProduct(Ac)};
      if (std::abs(bAcAc.sum() - 1/20) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*A*c*A*c) = " << bAcAc.sum() << " ≠ 1/20.");
        return order;
      }

      matS Ac2{A*(c.pow(2))};
      matS bcAc2{bc.cwiseProduct(Ac2)};
      if (std::abs(bcAc2.sum() - 1/15) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*c*A*c^2) = " << bcAc2.sum() << " ≠ 1/15.");
        return order;
      }

      matS Ac3{A*(c.pow(3))};
      matS bAc3{b.cwiseProduct(Ac3)};
      if (std::abs(bAc3.sum() - 1/20) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*A*c^3) = " << bAc3.sum() << " ≠ 1/20.");
        return order;
      }

      bAcAc = bA.cwiseProduct(c.cwiseProduct(Ac));
      if (std::abs(bAcAc.sum() - 1/40) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*A*c*A*c) = " << bAcAc.sum() << " ≠ 1/40.");
        return order;
      }

      matS bAAc2{bA.cwiseProduct(Ac2)};
      if (std::abs(bAAc2.sum() - 1/60) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*A*c*A*c) = " << bAAc2.sum() << " ≠ 1/60.");
        return order;
      }

      matS AAc{A*Ac};
      matS bAAAc{bA.cwiseProduct(AAc)};
      if (std::abs(bAAAc.sum() - 1/120) > tol) {
        SANDALS_WARNING(CMD "order 5 check failed, sum(b*A*c*A*c) = " << bAAAc.sum() << " ≠ 1/120.");
        return order;
      }

      order = 5; // Order 5 is the highest order that can be checked

      // Check order 6
      vecS bc5{b.cwiseProduct(c.pow(5))};
      if (std::abs(bc5.sum() - 1/6) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(b*c^5) = " << bc5.sum() << " ≠ 1/6.");
        return order;
      }

      vecS bc3Ac{bc3.cwiseProduct(Ac)};
      if (std::abs(bc3Ac.sum() - 1/12) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bc3Ac) = " << bc3Ac.sum() << " ≠ 1/12.");
        return order;
      }

      vecS bcAcAc{bc.cwiseProduct(Ac).pow(2)};
      if (std::abs(bcAcAc.sum() - 1/24) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bcAcAc) = " << bc3Ac.sum() << " ≠ 1/24.");
        return order;
      }

      vecS bc2Ac2{bc2.cwiseProduct(Ac2)};
      if (std::abs(bc2Ac2.sum() - 1/18) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bc2Ac2) = " << bc2Ac2.sum() << " ≠ 1/18.");
        return order;
      }

      vecS bAc2Ac{b.cwiseProduct(Ac2.cwiseProduct(Ac))};
      if (std::abs(bAc2Ac.sum() - 1/36) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAc2Ac) = " << bAc2Ac.sum() << " ≠ 1/36.");
        return order;
      }

      vecS bcAc3{bc.cwiseProduct(Ac3)};
      if (std::abs(bcAc3.sum() - 1/24) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bcAc3) = " << bcAc3.sum() << " ≠ 1/24.");
        return order;
      }

      vecS Ac4{A*c.pow(4)};
      vecS bAc4{b.cwiseProduct(Ac4)};
      if (std::abs(bAc4.sum() - 1/30) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAc4) = " << bAc4.sum() << " ≠ 1/30.");
        return order;
      }

      vecS bc2A{A.transpose()*bc2};
      vecS bc2AAc{bc2A.cwiseProduct(Ac)};
      if (std::abs(bc2AAc.sum() - 1/36) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bc2AAc) = " << bc2AAc.sum() << " ≠ 1/36.");
        return order;
      }

      matS bAcAAc = bAc.cwiseProduct(A)*Ac;
      if (std::abs(bAcAAc.sum() - 1/72) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAcAAc) = " << bAcAAc.sum() << " ≠ 1/72.");
        return order;
      }

      vecS bcA{A.transpose()*bc};
      bcAcAc = bcA.cwiseProduct(cAc);
      if (std::abs(bcAcAc.sum() - 1/48) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bcAcAc) = " << bcAcAc.sum() << " ≠ 1/48.");
        return order;
      }

      bAc2Ac = bA.cwiseProduct(c.pow(2)).cwiseProduct(Ac);
      if (std::abs(bAc2Ac.sum() - 1/60) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAc2Ac) = " << bAc2Ac.sum() << " ≠ 1/60.");
        return order;
      }

      vecS bAAcAc{bA.cwiseProduct(Ac.pow(2))};
      if (std::abs(bAAcAc.sum() - 1/120) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAAcAc) = " << bAAcAc.sum() << " ≠ 1/120.");
        return order;
      }

      vecS bcAAc2{bcA.cwiseProduct(Ac2)};
      if (std::abs(bcAAc2.sum() - 1/72) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bcAAc2) = " << bcAAc2.sum() << " ≠ 1/72.");
        return order;
      }

      vecS bAcAc2{bA.cwiseProduct(c).cwiseProduct(Ac2)};
      if (std::abs(bAcAc2.sum() - 1/90) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAcAc2) = " << bAcAc2.sum() << " ≠ 1/90.");
        return order;
      }

      vecS bAAc3{bA.cwiseProduct(Ac3)};
      if (std::abs(bAAc3.sum() - 1/120) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAAc3) = " << bAAc3.sum() << " ≠ 1/120.");
        return order;
      }

      vecS bcAAAc{bcA.cwiseProduct(A)*Ac};
      if (std::abs(bcAAAc.sum() - 1/144) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bcAAAc) = " << bcAAAc.sum() << " ≠ 1/144.");
        return order;
      }

      bAcAAc = (bA.cwiseProduct(c)).cwiseProduct(A)*Ac;
      if (std::abs(bAcAAc.sum() - 1/180) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAcAAc) = " << bAcAAc.sum() << " ≠ 1/180.");
        return order;
      }

      bAAcAc = bA.cwiseProduct(A)*(cAc);
      if (std::abs(bAAcAc.sum() - 1/240) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAAcAc) = " << bAAcAc.sum() << " ≠ 1/240.");
        return order;
      }

      vecS bAAAc2{bA.cwiseProduct(A)*Ac2};
      if (std::abs(bAAAc2.sum() - 1/360) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAAcAc) = " << bAAAc2.sum() << " ≠ 1/360.");
        return order;
      }

      vecS bAAAAc{bA.cwiseProduct(A)*(A*Ac)};
      if (std::abs(bAAAAc.sum() - 1/720) > tol) {
        SANDALS_WARNING(CMD "order 6 check failed, sum(bAAcAc) = " << bAAAAc.sum() << " ≠ 1/720.");
        return order;
      }

      order = 6; // Order 6 is the highest order that can be checked
      return order;
    }*/

  }; // class RungeKutta

} // namespace Sandals

#endif // SANDALS_RUNGEKUTTA_HXX
