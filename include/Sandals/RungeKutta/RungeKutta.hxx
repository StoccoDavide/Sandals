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
  //!
  //! Class container for the Butcher tableau of a Runge-Kutta solver. The
  //! Butcher tableau is defined as:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   \mathbf{c} & \mathbf{A} \\ \hline
  //!              & \mathbf{b} \\
  //!              & \mathbf{b}_{e}
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
  //! \f$ \mathbf{b}_{e} \f$ is the (optional) embedded Runge-Kutta weights
  //! vector relative to a method of order \f$ \hat{p} \f$
  //! (usually \f$ \hat{p} = p−1 \f$ or \f$ \hat{p} = p+1 \f$) (row vector):
  //!
  //! \f[
  //!   \mathbf{b}_{e} = \left[ \hat{b}_1, \hat{b}_2, \dots, \hat{b}_s \right] \text{,}
  //! \f]
  //!
  //! and \f$ \mathbf{c} \f$ is the Runge-Kutta nodes vector (column vector):
  //!
  //! \f[
  //!   \mathbf{c} = \left[ c_1, c_2, \dots, c_s \right]^{\top} \text{.}
  //! \f]
  template <Size S>
  struct Tableau
  {
    using Vector = Eigen::Vector<Real, S>;    //!< Vector type.
    using Matrix = Eigen::Matrix<Real, S, S>; //!< Matrix type.
    using Type = enum class type : Size {
      ERK = 0, IRK = 1, DIRK = 2
    }; //!< Runge-Kutta type enumeration.

    std::string name;        //!< Name of the method.
    type        rk_type;     //!< Runge-Kutta type.
    Size        order;       //!< Order of the method.
    bool        is_embedded; //!< Embedded method boolean.
    Matrix      A;           //!< Matrix \f$ \mathbf{A} \f$.
    Vector      b;           //!< Weights vector \f$ \mathbf{b} \f$ (row vector).
    Vector      b_e;         //!< Embedded weights vector \f$ \mathbf{b}_{e} \f$ (row vector).
    Vector      c;           //!< Nodes vector \f$ \mathbf{c} \f$ (column vector).
  };


  template <Size N>
  struct Solution
  {
    using Vector = Eigen::Vector<Real, Eigen::Dynamic>;    //!< Templetized matrix type.
    using Matrix = Eigen::Matrix<Real, N, Eigen::Dynamic>; //!< Templetized matrix type.

    Vector time;   //!< Time vector.
    Matrix states; //!< States solution matrix.

    //! Default class constructor for the Solution class.
    Solution() : time(0), states(Vector::Zero(N, 0)) {}

    //! Class constructor for the Solution class given the time and states sizes.
    //! \param size The size of the time and states vectors.
    Solution(Size size)
      : time(size), states(Vector::Zero(N, size)) {}

    //! Resize the time and states vectors.
    //! \param size The size of the time and states vectors.
    void resize(Size size) {
      this->time.resize(size);
      this->states.resize(N, size);
    }

    //! Resize conserving the data of the time and states vectors.
    //! \param size The size of the time and states vectors.
    void conservative_resize(Size size) {
      this->time.conservativeResize(size);
      this->states.conservativeResize(Eigen::NoChange, size);
    }

    //! Clear the time and states vectors.
    void clear() {
      this->time.resize(0);
      this->states.resize(N, 0);
    }

    //! Check if the time and states vectors are empty.
    //! \return True if the time and states vectors are empty, false otherwise.
    bool is_empty() const {return this->time.size() == 0 && this->states.cols()== 0 ;}

    //! Get the size of the time and states vectors.
    //! \return The size of the time and states vectors.
    Size size() const {return this->time.size();}

    //! Get the time vector as a standard vector.
    //! \return The time vector as a standard vector.
    std::vector<Real> std_time() const {
      return std::vector<Real>(this->time.data(), this->time.data() + this->time.size());
    }

    //! Get the states i-th state vector as a standard vector.
    //! \param i The index of the state vector.
    //! \return The states i-th state vector as a standard vector.
    std::vector<Real> std_state(Size i) const {
      Vector tmp(this->states.row(i));
      return std::vector<Real>(tmp.data(), tmp.data() + tmp.size());
    }

    //! Get a map of the states vectors as standard vectors.
    //! \return A map of the states vectors as standard vectors.
    std::map<Size, std::vector<Real>> std_states() const {
      std::map<Size, std::vector<Real>> states_map;
      for (Size i = 0; i < N; ++i) {
        states_map[i] = this->std_states(i);
      }
      return states_map;
    }
    //! Get a map of the states vectors as Eigen vectors.
    //! \return A map of the states vectors as Eigen vectors.
    std::map<Size, Vector> eig_states() const {
      std::map<Size, Vector> states_map;
      for (Size i = 0; i < N; ++i) {
        states_map[i] = this->states.row(i);
      }
      return states_map;
    }

    //! Get a map of the states vectors as standard vectors.
    //! \param names The names of the states vectors.
    //! \return A map of the states vectors as standard vectors.
    std::map<std::string, std::vector<Real>> std_states(std::vector<std::string> names) const {
      std::map<std::string, std::vector<Real>> states_map;
      for (Size i = 0; i < N; ++i) {
        states_map[names[i]] = this->std_state(i);
      }
      return states_map;
    }

    //! Get a map of the states vectors as Eigen vectors.
    //! \param names The names of the states vectors.
    //! \return A map of the states vectors as Eigen vectors.
    std::map<std::string, Vector> eig_states(std::vector<std::string> names) const {
      std::map<std::string, Vector> states_map;
      for (Size i = 0; i < N; ++i) {
        states_map[names[i]] = this->states.row(i);
      }
      return states_map;
    }
  };


  //! \brief Class container for Runge-Kutta solvers of the system of ODEs.
  //!
  //! Class container for Runge-Kutta solvers of the system of Ordinary
  //! Differential Equations (ODEs). The user must define the Butcher tableau
  //! of the solver of order \f$ p \f$ and \f$ s \f$ stages, which defined as:
  //!
  //! \f[
  //! \begin{array}{c|c}
  //!   \mathbf{c} & \mathbf{A} \\ \hline
  //!              & \mathbf{b} \\
  //!              & \mathbf{b}_{e}
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
  //! \f$ \mathbf{b}_{e} \f$ is the (optional) embedded Runge-Kutta weights
  //! vector relative to a method of order \f$ \hat{p} \f$
  //! (usually \f$ \hat{p} = p−1 \f$ or \f$ \hat{p} = p+1 \f$) (row vector):
  //!
  //! \f[
  //!   \mathbf{b}_{e} = \left[ \hat{b}_1, \hat{b}_2, \dots, \hat{b}_s \right] \text{,}
  //! \f]
  //!
  //! and \f$ \mathbf{c} \f$ is the Runge-Kutta nodes vector (column vector):
  //!
  //! \f[
  //!   \mathbf{c} = \left[ c_1, c_2, \dots, c_s \right]^{\top} \text{.}
  //! \f]
  //
  template <Size S, Size N>
  class RungeKutta
  {
    using VectorK = Eigen::Vector<Real, N*S>;      //!< Templetized vector type.
    using MatrixK = Eigen::Matrix<Real, N, S>;     //!< Templetized matrix type.
    using MatrixJ = Eigen::Matrix<Real, N*S, N*S>; //!< Templetized matrix type.
    using VectorS = typename Tableau<S>::Vector;   //!< Templetized vector type.
    using MatrixS = typename Tableau<S>::Matrix;   //!< Templetized matrix type.
    using VectorN = typename Implicit<N>::Vector;  //!< Templetized vector type.
    using MatrixN = typename Implicit<N>::Matrix;  //!< Templetized matrix type.

  public:
    using System   = typename Implicit<N>::Pointer;       //!< Shared pointer to an implicit ODE system.
    using Type     = typename Tableau<S>::Type;           //!< Runge-Kutta type enumeration.
    using Time     = Eigen::Vector<Real, Eigen::Dynamic>; //!< Templetized vector type for the time.
    using Solution = Solution<N>;                         //!< Templetized structure type for solution.

  private:
    Newton<N*S> m_newton;                //!< Newton solver for implicit methods.
    Tableau<S>  m_tableau;               //!< Butcher tableau of the Runge-Kutta method.
    System      m_system;                //!< ODE system object pointer.
    Real        m_abs_tol{EPSILON_LOW};  //!< Absolute tolerance for adaptive step.
    Real        m_rel_tol{EPSILON_LOW};  //!< Relative tolerance for adaptive step.
    Real        m_sft_fac{0.9};          //!< Safety factor for adaptive step.
    Real        m_min_sft_fac{0.2};      //!< Minimum safety factor for adaptive step.
    Real        m_max_sft_fac{1.5};      //!< Maximum safety factor for adaptive step.
    Real        m_d_t_min{EPSILON_HIGH}; //!< minimum step for advancing;
    bool        m_adaptive_step{false};  //!< Aadaptive step mode boolean.
    bool        m_verbose{false};        //!< Verbose mode boolean.
    Size        m_order;                 //!< Order of the solver.
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
    RungeKutta(Tableau<S> const &t_tableau, System t_system)
      : m_tableau(t_tableau), m_system(t_system) {
      //SANDALS_WARNING(
      //  "RungeKutta::RungeKutta(...): the order of the method " << t_tableau.name << " is " <<
      //  this->tableau_order(this->m_tableau.A, this->m_tableau.b, this->m_tableau.c));
    }

    //! Check if the solver is explicit, implicit or diagonal implicit
    //! \return The type of the solver.
    Type rk_type() const {return this->m_tableau.rk_type;}

    //! Cheack if the solver is explicit.
    //! \return True if the solver is explicit, false otherwise.
    bool is_explicit() const {return this->m_tableau.rk_type == Type::ERK;}

    //! Cheack if the solver is implicit.
    //! \return True if the solver is implicit, false otherwise.
    bool is_implicit() const {return this->m_tableau.rk_type == Type::IRK;}

    //! Cheack if the solver is diagonal implicit.
    //! \return True if the solver is diagonal implicit, false otherwise.
    bool is_diagonal_implicit() const {return this->m_tableau.rk_type == Type::DIRK;}

    //! Get the Tableau reference.
    //! \return The Tableau reference.
    Tableau<S> & tableau() {return this->m_tableau;}

    //! Get the Tableau const reference.
    //! \return The Tableau const reference.
    Tableau<S> const & tableau() const {return this->m_tableau;}

    //! Get the stages number of the method.
    //! \return The stages number of the method.
    Size stages() const {return S;}

    //! Get the name of the method.
    //! \return The name of the method.
    std::string name() const {return this->m_tableau.name;}

    //! Get the order of the method.
    //! \return The order of the method.
    Size order() const {return this->m_tableau.order;}

    //! Check if the method is embedded.
    //! \return True if the method is embedded, false otherwise.
    bool is_embedded() const {return this->m_tableau.is_embedded;}

    //! Get the matrix \f$ \mathbf{A} \f$.
    //! \return The matrix \f$ \mathbf{A} \f$.
    MatrixS A() const {return this->m_tableau.A;}

    //! Get the weights vector \f$ \mathbf{b} \f$.
    //! \return The weights vector \f$ \mathbf{b} \f$.
    VectorS b() const {return this->m_tableau.b;}

    //! Get the embedded weights vector \f$ \mathbf{b}_{e} \f$.
    //! \return The embedded weights vector \f$ \mathbf{b}_{e} \f$.
    VectorS b_embedded() const {return this->m_tableau.b_e;}

    //! Get the nodes vector \f$ \mathbf{c} \f$.
    //! \return The nodes vector \f$ \mathbf{c} \f$.
    VectorS c() const {return c;}

    //! Get the ODE system shared pointer.
    //! \return The ODE system shared pointer.
    System ode() {return this->m_system;}

    //! Set the ODE system shared pointer.
    //! \param t_system The ODE system shared pointer.
    void ode(System t_system) {this->m_system = t_system;}

    //! Get the adaptive step absolute tolerance.
    //! \return The adaptive step absolute tolerance.
    Real absolute_tolerance() {return this->m_abs_tol;}

    //! Get the adaptive step absolute tolerance reference.
    //! \param t_abs_tol The adaptive step absolute tolerance.
    void absolute_tolerance(Real t_abs_tol) {this->m_abs_tol = t_abs_tol;}

    //! Get the adaptive step relative tolerance.
    //! \return The adaptive step relative tolerance.
    Real relative_tolerance() {return this->m_rel_tol;}

    //! Get the adaptive step relative tolerance reference.
    //! \param t_rel_tol The adaptive step relative tolerance.
    void relative_tolerance(Real t_rel_tol) {this->m_rel_tol = t_rel_tol;}

    //! Get the safety factor for adaptive step.
    //! \return The safety factor for adaptive step.
    Real & safety_factor() {return this->m_sft_fac;}

    //! Set safety factor for adaptive step.
    //! \param t_sft_fac The safety factor for adaptive step.
    void safety_factor(Real t_sft_fac) {this->m_sft_fac = t_sft_fac;}

    //! Get the minimum safety factor for adaptive step.
    //! \return The minimum safety factor for adaptive step.
    Real & min_safety_factor() {return this->m_min_sft_fac;}

    //! Set the minimum safety factor for adaptive step.
    //! \param t_min_sft_fac The minimum safety factor for adaptive step.
    void min_safety_factor(Real t_min_sft_fac) {this->m_min_sft_fac = t_min_sft_fac;}

    //! Get the maximum safety factor for adaptive step.
    //! \return The maximum safety factor for adaptive step.
    Real & max_safety_factor() {return this->m_max_sft_fac;}

    //! Set the maximum safety factor for adaptive step.
    //! \param t_max_sft_fac The maximum safety factor for adaptive step.
    void max_safety_factor(Real t_max_sft_fac) {this->m_max_sft_fac = t_max_sft_fac;}

    //! Enable verbose mode.
    void enable_verbose_mode() {this->m_verbose = true;}

    //! Disable verbose mode.
    void disable_verbose_mode() {this->m_verbose = false;}

    //! Enable adaptive step mode.
    void enable_adaptive_step() {this->m_adaptive_step = true;}

    //! Disable adaptive step mode.
    void disable_adaptive_step() {this->m_adaptive_step = false;}

    //! Compute adaptive time step for the next advancing step according to the
    //! error control method. The error control method used is the local truncation
    //! error method, which is based on the following formula:
    //!
    //! \f[
    //! e = \sqrt{\frac{1}{n} \displaystyle\sum_{i=1}{n}\left(\frac
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
    //! h_{opt} = h \left( \frac{1}{e} \right)^{\frac{1}{q+1}}
    //! \f]
    //!
    //! We multiply the previous quation by a safety factor \f$ f \f$, usually
    //! \f$ f = 0.8 \f$, \f$ 0.9 \f$, \f$ (0.25)^{1/(q+1)} \f$, or \f$ (0.38)^{1/(q+1)} \f$,
    //! so that the error will be acceptable the next time with high probability.
    //! Further, \f$ h \f$ is not allowed to increase nor to decrease too fast.
    //! So we put:
    //!
    //! \f[
    //! h_{new} = h \min\left(f_{max}, \max\left(f_{max}, f \left(\frac{1}{e}\right)^{\frac{1}{q+1}}
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
    Real estimate_step(VectorN const &x_h, VectorN const &x_l, Real d_t) const
    {
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
    //! \f$ \mathbf{F}_i \f$ is a triangular system (which may be nonlinear in
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
    //! \param t_step Advancing time step \f$ \Delta t\f$.
    //! \param x_out The approximation of the states at \f$ k+1 \f$-th time step
    //!            \f$ \mathbf{x_{k+1}}(t_{k}+\Delta t) \f$
    //! \param t_adap The suggested adaptive time step for the next advancing step
    //!               \f$ \Delta t_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool explicit_step(VectorN const &x_k, Real t_k, Real t_step, VectorN &x_out, Real &t_adap) const
    {
      using Eigen::all;
      using Eigen::seqN;

      // Compute the K variables in the case of an explicit method and explicit system
      Real t_i;
      VectorN x_i;
      MatrixK K;
      for (Size i = 0; i < S; ++i)
      {
        t_i = t_k + this->m_tableau.c(i) * t_step;
        x_i = x_k + K(all, seqN(0, i)) * this->m_tableau.A(i, seqN(0, i)).transpose();
        K.col(i) = t_step * static_cast<Explicit<N> const *>(this->m_system.get())->f(x_i, t_i);
      }
      if (!K.allFinite()) {return false;}

      // Perform the step and obtain the next state
      x_out = x_k + K * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_is_embedded) {
        VectorN x_e = x_k + K * this->m_tableau.b_e;
        t_adap = this->estimate_step(x_out, x_e, t_step);
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
    void implicit_function(VectorN const &x_k, VectorK const &K_k, Real t_k, Real d_t, VectorK &fun) const
    {
      fun.setZero();

      // Loop through each equation of the system
      Real t_i;
      VectorN x_i;
      MatrixK K{K_k.reshaped(N, S)}; // FIXME: can i avoid reshaping?
      MatrixK fun_tmp;
      for (Size i = 0; i < S; ++i) {
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
    //! \frac{\partial \mathbf{F}_i}{\partial \mathbf{K}_i} \left(
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
    void implicit_jacobian(VectorN const &x_k, VectorK const &K_k, Real t_k, Real d_t, MatrixJ &jac) const
    {
      using Eigen::all;
      using Eigen::seqN;

      jac.setZero();

      // Loop through each equation of the system
      MatrixK K{K_k.reshaped(N, S)};
      Real t_i;
      VectorN x_i, x_dot_i;
      MatrixN JF_x, JF_x_dot;
      auto idx = seqN(0, N), jdx = seqN(0, N);
      for (Size i = 0; i < S; ++i) {
        t_i = t_k + this->m_tableau.c(i) * d_t;
        x_i = x_k + K * this->m_tableau.A.row(i).transpose();

        // Compute the Jacobians with respect to x and x_dot
        x_dot_i  = K.col(i) / d_t;
        JF_x     = this->m_system->JF_x(x_i, x_dot_i, t_i);
        JF_x_dot = this->m_system->JF_x_dot(x_i, x_dot_i, t_i);

        // Derivative of F(x_i, K(:,i)/d_t, t_i)
        idx = seqN(i*N, N);
        for (Size j = 0; j < S; ++j) {
          // Combine the Jacobians with respect to x and x_dot to obtain the
          // Jacobian with respect to K
          jdx = seqN(j*N, N);
          if (i == j) {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x + JF_x_dot / d_t;
          } else {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x;
          }
        }
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
    bool implicit_step(VectorN const &x_k, Real t_k, Real d_t, VectorN &x_o, Real &d_t_star)
    {
      // Check if the solver converged
      VectorK K;

      if (!this->m_newton.solve(
        [this, &x_k, t_k, d_t](VectorK const &K, VectorK &fun) {this->implicit_function(x_k, K, t_k, d_t, fun);},
        [this, &x_k, t_k, d_t](VectorK const &K, MatrixJ &jac) {this->implicit_jacobian(x_k, K, t_k, d_t, jac);},
        VectorK::Zero(), K)) {return false;}

      // Perform the step and obtain x_k+1
      x_o = x_k + K.reshaped(N, S) * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_tableau.is_embedded) {
        VectorN x_e{x_k + K.reshaped(N, S) * this->m_tableau.b_e};
        d_t_star = this->estimate_step(x_o, x_e, d_t);
      }
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
    bool step(VectorN const &x_k, Real t_k, Real d_t, VectorN &x_o, Real &d_t_star)
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
    bool advance(VectorN const &x_k, Real t_k, Real d_t, VectorN &x_new, Real &d_t_star)
    {
      #define CMD "Indigo.RungeKutta.advance(...): "

      // Check step size
      SANDALS_ASSERT(d_t > Real(0.0), CMD "in " << this->m_tableau.name <<
        "solver, d_t = "<< d_t << ", expected > 0.");

      // If the integration step failed, try again with substepping
      if (!this->step(x_k, t_k, d_t, x_new, d_t_star))
      {
        VectorN x_tmp   = x_k;
        Real t_tmp   = t_k;
        Real d_t_tmp = d_t / Real(2.0);

        // Substepping logic
        Size max_substeps = 10;
        Size max_k = max_substeps * max_substeps;
        Size k = 2;
        Real d_t_star_tmp;
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
                d_t_tmp = Real(2.0) * d_t_tmp;
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
            d_t_tmp = d_t_tmp / Real(2.0);
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
    bool solve(VectorD const &t, VectorN const &ics, Solution &sol)
    {
      using Eigen::last;

      // Instantiate output
      sol.resize(t.size());

      // Store first step
      sol.time(0)       = t(0);
      sol.states.col(0) = ics;

      // Update the current step
      Size s{0};
      VectorN x_s{ics};
      Real t_s{t(0)};
      Real d_t_s{t(1) - t(0)};
      Real d_t_tmp{d_t_s}, d_t_star;
      bool mesh_point_bool, saturation_bool;

      while (true) {
        // Integrate system
        if (!this->advance(sol.states.col(s), t_s, d_t_s, x_s, d_t_star)) {return false;}

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
          sol.time(s)       = t_s;
          sol.states.col(s) = x_s;

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
    bool adaptive_solve(VectorD const &t, VectorN const &ics, Solution &sol)
    {
      using Eigen::last;

      #define CMD "Sandals::RungeKutta::adaptive_solve(...): "

      // Instantiate output
      this->enable_adaptive_step();
      Real d_t{t(1) - t(0)}, d_t_star;
        Real scale{100.0};
        std::cout << "d_t = " << d_t << std::endl;

        Real t_min{std::max(this->m_d_t_min, d_t/scale)};
        Real t_max{scale*d_t};
        Size safety_length;
      if (this->m_tableau.is_embedded) {
        safety_length = std::ceil(Real(1.5)/this->m_min_sft_fac) * t.size();
        sol.resize(safety_length);
      } else {
        sol.resize(t.size());
      }

      // Store first step
      sol.time(0)      = t(0);
      sol.states.col(0) = ics;

      // Instantiate temporary variables
      Size s{0};
      VectorN x_s{ics};

      while (true) {
        // Integrate system
        this->advance(sol.states.col(s), sol.time(s), d_t, x_s, d_t_star);

        // Saturate the suggested timestep
        if (this->m_adaptive_step && this->m_tableau.is_embedded) {
          d_t = std::max(std::min(d_t_star, t_max), t_min);
        }
        std::cout << "d_t = " << d_t_star << std::endl;

        SANDALS_ASSERT(s < t.size(), CMD "safety length exceeded.");

        // Store solution
        sol.time(s+1)       = sol.time(s) + d_t;
        sol.states.col(s+1) = x_s;

        // Check if the current step is the last one
        if (sol.time(s+1) + d_t > t(last)) {break;}

        // Update steps counter
        s += 1;
      }

      // Resize the output
      sol.conservative_resize(s-1);

      return true;

      #undef CMD
    }

  }; // class RungeKutta

} // namespace Sandals

#endif // SANDALS_RUNGEKUTTA_HXX
