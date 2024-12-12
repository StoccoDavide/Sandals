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
    using vecS = Eigen::Matrix<real, S, 1>; //!< Vector type.
    using matS = Eigen::Matrix<real, S, S>; //!< Matrix type.

    std::string name;        //!< Name of the method.
    unsigned    order;       //!< Order of the method.
    bool        is_embedded; //!< Embedded method boolean.
    matS        A;           //!< Matrix \f$ \mathbf{A} \f$.
    vecS        b;           //!< Weights vector \f$ \mathbf{b} \f$ (row vector).
    vecS        b_e;         //!< Embedded weights vector \f$ \mathbf{b}_{e} \f$ (row vector).
    vecS        c;           //!< Nodes vector \f$ \mathbf{c} \f$ (column vector).
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
    using matK = Eigen::Matrix<real, N, S>;  //!< Matrix type.

  public:
    using vecS = typename Tableau<S>::vecS;   //!< Vector type.
    using matS = typename Tableau<S>::matS;   //!< Matrix type.
    using vecN = typename Implicit<N>::vecN;  //!< Vector type.
    using matN = typename Implicit<N>::matN;  //!< Matrix type.
    using odeN = typename Implicit<N>::ptr;   //!< Matrix type.
    using nlsN = NonlinearSolver<vecN, matN>; //!< Non-linear system solver.
    using RKType = enum class RKType : integer {
      ERK = 0, IRK = 1, DIRK = 2
    }; //!< Runge-Kutta type enumeration.

  private:
    nlsN       m_nlsolver;             //!< Non-linear system solver.
    RKType     m_rk_type;              //!< Runge-Kutta type (explicit, implicit, diagonal implicit).
    Tableau<S> m_tableau;              //!< Butcher tableau of the Runge-Kutta method.
    odeN       m_sys;                  //!< ODE system object pointer.
    real       m_abs_tol{EPSILON_LOW}; //!< Absolute tolerance for adaptive step.
    real       m_rel_tol{EPSILON_LOW}; //!< Relative tolerance for adaptive step.
    real       m_sfy_fac{0.9};         //!< Safety factor for adaptive step.
    real       m_min_sfy_fac{0.2};     //!< Minimum safety factor for adaptive step.
    real       m_max_sfy_fac{1.5};     //!< Maximum safety factor for adaptive step.
    real       m_dt_min{EPSILON_HIGH}; //!< minimum step for advancing;
    bool       m_adaptive_step{false}; //!< Aadaptive step mode boolean.
    bool       m_verbose{false};       //!< Verbose mode boolean.
    unsigned   m_order;                //!< Order of the solver.
    bool       m_is_embedded;          //!< Boolean to check if the method is embedded.

  public:
    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    //! \param t_tableau The Tableau reference.
    RungeKutta(Tableau<S> t_tableau) : m_tableau(t_tableau) {}

    //! Class constructor for a Runge-Kutta solver given a namne and a Tableau.
    //! \param t_name The name of the solver.
    //! \param t_tableau The Tableau reference.
    RungeKutta(std::string t_name, Tableau<S> t_tableau) : m_tableau(t_tableau) {
      this->m_tableau.name = t_name;
    }

    //! Class constructor for a Runge-Kutta solver given all Tableau elements.
    //! \param t_name The name of the solver.
    //! \param t_A The matrix \f$ \mathbf{A} \f$ (lower triangular matrix).
    //! \param t_b The weights vector \f$ \mathbf{b} \f$.
    //! \param t_b_e The embedded weights vector \f$ \\mathbf{b}_{e} \f$.
    //! \param t_c The nodes vector \f$ \mathbf{c} \f$.
    RungeKutta(std::string t_name, matS const &t_A,vecS const &t_b, vecS const &t_b_e,
      vecS const &t_c) : m_tableau({t_name, t_A, t_b, t_b_e, t_c}) {}

    //! Check if the solver is explicit, implicit or diagonal implicit
    //! \return The type of the solver.
    RKType rk_type(void) const {return this->m_rk_type;}

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

    //! Get the ODE system pointer.
    //! \return The ODE system pointer.
    odeN ode(void) {return this->m_sys;}

    //! Set the ODE system pointer.
    //! \param t_ode_ptr The ODE system pointer.
    void ode(odeN t_ode_ptr) {this->m_sys = t_ode_ptr;}

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
    real & safety_factor(void) {return this->m_sfy_fac;}

    //! Set safety factor for adaptive step.
    //! \param t_sfy_fac The safety factor for adaptive step.
    void safety_factor(real t_sfy_fac) {this->m_sfy_fac = t_sfy_fac;}

    //! Get the minimum safety factor for adaptive step.
    //! \return The minimum safety factor for adaptive step.
    real & min_safety_factor(void) {return this->m_min_sfy_fac;}

    //! Set the minimum safety factor for adaptive step.
    //! \param t_min_sfy_fac The minimum safety factor for adaptive step.
    void min_safety_factor(real t_min_sfy_fac) {this->m_min_sfy_fac = t_min_sfy_fac;}

    //! Get the maximum safety factor for adaptive step.
    //! \return The maximum safety factor for adaptive step.
    real & max_safety_factor(void) {return this->m_max_sfy_fac;}

    //! Set the maximum safety factor for adaptive step.
    //! \param t_max_sfy_fac The maximum safety factor for adaptive step.
    void max_safety_factor(real t_max_sfy_fac) {this->m_max_sfy_fac = t_max_sfy_fac;}

    //! Enable verbose mode.
    void enable_verbose_mode(void) {this->m_verbose = true;}

    //! Disable verbose mode.
    void disable_verbose_mode(void) {this->m_verbose = false;}

    //! Enable adaptive step mode.
    void enable_adaptive_step(void) {this->m_adaptive_step = true;}

    //! Disable adaptive step mode.
    void disable_adaptive_step(void) {this->m_adaptive_step = false;}

    //! Error estimation function for adaptive step (override in derived classes
    //! for specific error estimation).
    //! \param x_h The state vector computed with higher order method.
    //! \param x_l The state vector computed with lower order method.
    //! \return The error estimation.
    real estimate_error(vecN const &x_h, vecN const &x_l) const {
      return ((x_h - x_l) / (
          this->m_abs_tol + this->m_rel_tol*std::max(std::abs(x_h), std::abs(x_l))
        )).abs().maxCoeff();
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
    //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
    //! \param t_k Time step \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t\f$.
    //! \param x_o The approximation of the states at \f$ k+1 \f$-th time step
    //!              \f$ \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$
    //! \param d_t_star The suggested time step for the next advancing step
    //!                 \f$ \Delta t_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool explicit_step(vecS const &x_k, real t_k, real d_t, vecS &x_o, real &d_t_star) const
    {
      using Eigen::all;
      using Eigen::seqN;

      // Compute the K variables in the case of an explicit method and explicit system
      real t_s;
      vecS x_s;
      matK K{matK::Zero()};
      for (unsigned i = 0; i < S; ++i)
      {
        t_s = t_k + m_tableau.c(i) * d_t;
        x_s = x_k + K(all, seqN(0, i)) * m_tableau.A(i, seqN(0, i)).transpose();
        K.col(i) = d_t * m_sys->f(x_s, t_s);
      }
      if (!K.allFinite()) {return false;}

      // Perform the step and obtain the next state
      x_o = x_k + K * this->m_tableau.b.transpose();

      // Adapt next time step
      if (this->m_adaptive_step && this->m_is_embedded) {
        vecS x_e = x_k + K * this->m_tableau.b_e.transpose();
        d_t_star = d_t * std::min(this->m_max_sfy_fac, std::max(this->m_min_sfy_fac,
          this->m_safety_factor*this->estimate_error(x_o, x_e)
        ));
      }
      return true;
    }

    //! Compute a step using a generic integration method for a system of the
    //! form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}', \mathbf{v}, t) = \mathbf{0}
    //! \f$. The step is based on the following formula:
    //!
    //! \f[
    //! \mathbf{x}_{k+1}(t_{k}+\Delta t) = \mathbf{x}_k(t_{k}) +
    //! \mathcal{S}(\mathbf{x}_k(t_k), \mathbf{x}'_k(t_k), t_k, \Delta t)
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
    bool step(vecS const &x_k, real t_k, real d_t, vecS &x_o, real &d_t_star) const
    {
      if (this->is_explicit() && this->m_sys->is_explicit()) {
        return this->explicit_step(x_k, t_k, d_t, x_o, d_t_star);
      } else {
        return this->implicit_step(x_k, t_k, d_t, x_o, d_t_star);
      }
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
      //! \param i   Index of the step to be computed.
      //! \param x_i States value at \f$ i \f$-th node.
      //! \param K   Variable \f$ \mathbf{K} \f$ of the system to be solved.
      //! \param t_k Time step \f$ t_k \f$.
      //! \param dt  Advancing time step \f$ \Delta t\f$.
      //!
      //! \return The Jacobian of the system of equations to be solved.
      // out = K_jacobian( this, K_in, x_k, t_k, dt )

      //! Compute the residual of system to be solved:
      //!
      //! \f[
      //!   \mathbf{F}_i\left(\mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^{s}
      //!   a_{ij} \mathbf{K}_j, \, \mathbf{K}_i, \, t_k + c_i \Delta t
      //! \right) = \mathbf{0}.
      //! \f]
      //!
      //! \param x_k States value at \f$ k \f$-th time step \f$ \mathbf{x}(t_k) \f$.
      //! \param K   Variable \f$ \mathbf{K} \f$ of the system to be solved.
      //! \param t_k Time step \f$ t_k \f$.
      //! \param dt  Advancing time step \f$ \Delta t\f$.
      //!
      //! \return The residual of system to be solved.
      // out = K_residual( this, K_in, x_k, t_k, dt )

      //! Solve the \f$ i \f$-th implicit step of the system to find the
      //! \f$ \mathbf{K} \f$ variables:
      //!
      //! \f[
      //! \mathbf{F}_i\left(\mathbf{x}_k + \Delta t \displaystyle\sum_{j=1}^s
      //!   a_{ij} \mathbf{K}_j, \, \mathbf{K}_i, \, t_k + c_i \Delta t
      //! \right) = \mathbf{0}
      //! \f]
      //!
      //! by Newton's method.
      //!
      //! \param x_k States value at \f$ k \f$-th time step \f$
      //!            \mathbf{x}(t_k) \f$.
      //! \param K   Initial guess for the \f$ \mathbf{K} \f$ variables to be
      //!            found.
      //! \param t_k Time step \f$ t_k \f$.
      //! \param dt  Advancing time step \f$ \Delta t\f$.
      //!
      //! \return The \f$ \mathbf{K} \f$ variables of system to be solved and the
      //!         error control flag.
      // [K, ierr] = K_solve( this, x_k, t_k, dt )

      //! The suggested time step for the next advancing step
      //! \f$ \Delta t_{k+1} \f$, is the same as the input time step
      //! \f$ \Delta t \f$ since in the implicit Runge-Kutta method the time step
      //! is not modified through any error control method.
      //!
      //! \param x_k States value at \f$ k \f$-th time step
      //!            \f$ \mathbf{x}(t_k) \f$.
      //! \param t_k Time step \f$ t_k \f$.
      //! \param dt  Advancing time step \f$ \Delta t\f$.
      //!
      //! \return The approximation of the states at \f$ k+1 \f$-th time step \f$
      //!         \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$, the approximation of the
      //!         states derivatives at \f$ k+1 \f$-th time step
      //!         \f$ \mathbf{x}^{\prime}_{k+1} (t_{k}+\Delta t) \f$, the suggested time
      //!         step for the next advancing step \f$ \Delta t_{k+1} \f$, and the
      //!         error control flag.
      //
      // [x_o, d_t_star, ierr] = step( this, x_k, t_k, dt )

      //! Solve the system and calculate the approximate solution over
      //! the mesh of time points.
      //!
      //! \param t   Time mesh points
      //!            \f$ \mathbf{t} = \left[ t_0, t_1, \ldots, t_n \right]^{\top} \f$.
      //! \param x_0 Initial states value \f$ \mathbf{x}(t_0) \f$.
      //!
      //! \return A matrix
      //!         \f$ \left[(\mathrm{size}(\mathbf{x}) \times \mathrm{size}(\mathbf{t})\right] \f$
      //!         containing the approximated solution
      //!         \f$ \mathbf{x}(t) \f$ and
      //!         \f$ \mathbf{x}^\prime(t) \f$
      //!         of the system.
      // [x_o, t] = solve( this, t, x_0 )

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
      //! h_{new} = h \min \left( f_{max}, \max \left( f_{max}, f \left(
      //!   \dfrac{1}{e} \right)^{\frac{1}{q+1}}
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
      //! \param dt  Actual advancing time step \f$ \Delta t\f$.
      //!
      //! \return The suggested time step for the next advancing step \f$ \Delta
      //!         t_{k+1} \f$.
      // out = estimate_step( this, x_h, x_l, dt )




      //! Check Butcher tableau consistency for an explicit Runge-Kutta method.
      //! \param tbl.A   Matrix \f$ \mathbf{A} \f$.
      //! \param tbl.b   Weights vector \f$ \mathbf{b} \f$.
      //! \param tbl.b_e [optional] Embedded weights vector \f$ \\mathbf{b}_{e} \f$.
      //! \param tbl.c   Nodes vector \f$ \mathbf{c} \f$.
      //! \return True if the Butcher tableau is consistent, false otherwise.
      //[out,order,e_order] = check_tableau( this, tbl )

      //[order,msg] = tableau_order( this, A, b, c )
  }; // class RungeKutta

} // namespace Sandals

#endif // SANDALS_RUNGEKUTTA_HXX
