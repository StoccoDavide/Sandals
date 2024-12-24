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

  /*\
   |   _____     _     _
   |  |_   _|_ _| |__ | | ___  __ _ _   _
   |    | |/ _` | '_ \| |/ _ \/ _` | | | |
   |    | | (_| | |_) | |  __/ (_| | |_| |
   |    |_|\__,_|_.__/|_|\___|\__,_|\__,_|
   |
  \*/

  //! \brief Class container for the Butcher tableau of a Runge-Kutta methos.
  //!
  //! \includedoc docs/markdown/Tableau.md
  template <Size S>
  struct Tableau
  {
    using Type = enum class type : Size {ERK=0, IRK=1, DIRK=2}; //!< Runge-Kutta type enumeration.
    using Vector = Eigen::Vector<Real, S>;    //!< Templetized vector type.
    using Matrix = Eigen::Matrix<Real, S, S>; //!< Templetized matrix type.

    std::string name;        //!< Name of the Runge-Kutta method.
    Type        type;        //!< Runge-Kutta type.
    Size        order;       //!< Order of the Runge-Kutta method.
    bool        is_embedded; //!< Embedded method boolean.
    Matrix      A;           //!< Matrix \f$ \mathbf{A} \f$.
    Vector      b;           //!< Weights vector \f$ \mathbf{b} \f$.
    Vector      b_e;         //!< Embedded weights vector \f$ \hat{\mathbf{b}} \f$.
    Vector      c;           //!< Nodes vector \f$ \mathbf{c} \f$.
  };

  /*\
   |   ____        _       _   _
   |  / ___|  ___ | |_   _| |_(_) ___  _ __
   |  \___ \ / _ \| | | | | __| |/ _ \| '_ \
   |   ___) | (_) | | |_| | |_| | (_) | | | |
   |  |____/ \___/|_|\__,_|\__|_|\___/|_| |_|
   |
  \*/

  //! \brief Class container for the numerical solution of a system of ODEs/DAEs.
  //!
  //! Class container for the numerical solution of a system of Ordinary Differential Equations (ODEs)
  //! or Differential Algebraic Equations (DAEs). The solution is stored in a time vector and states
  //! matrix, where each row of the matrix represents the solution of a state at a given time. The
  //! solution can be accessed as standard library vectors or Eigen library vectors.
  template <Size N>
  struct Solution
  {
    using Vector = Eigen::Vector<Real, Eigen::Dynamic>;    //!< Templetized vector type.
    using Matrix = Eigen::Matrix<Real, N, Eigen::Dynamic>; //!< Templetized matrix type.

    Vector time;   //!< Time solution vector.
    Matrix states; //!< States solution matrix.

    //! Default class constructor for the Solution class.
    //!
    Solution() : time(0), states(Vector::Zero(N, 0)) {}

    //! Class constructor for the Solution class given the time and states sizes.
    //! \param[in] size The size of the time vector and states matrix.
    Solution(Size size)
      : time(size), states(Vector::Zero(N, size)) {}

    //! Resize the time vector and states matrix.
    //! \param[in] size The size of the time vector and states matrix.
    void resize(Size size) {
      this->time.resize(size);
      this->states.resize(N, size);
    }

    //! Resize conserving the data of the time vector and states matrix.
    //! \param[in] size The size of the time vector and states matrix.
    void conservative_resize(Size size) {
      this->time.conservativeResize(size);
      this->states.conservativeResize(Eigen::NoChange, size);
    }

    //! Clear the time vector and states matrix.
    //!
    void clear()
    {
      this->time.resize(0);
      this->states.resize(N, 0);
    }

    //! Check if the time vector and states matrix are empty.
    //! \return True if the time vector and states matrix are empty, false otherwise.
    bool is_empty() const {return this->time.size() == 0 && this->states.cols()== 0 ;}

    //! Get the size of the time vector and states matrix.
    //! \return The size of the time vector and states matrix.
    Size size() const {return this->time.size();}

    //! Get the time vector as a standard library vector.
    //! \return The time vector as a standard library vector.
    std::vector<Real> std_time() const
    {
      return std::vector<Real>(this->time.data(), this->time.data() + this->time.size());
    }

    //! Get the states i-th state vector as a standard library vector.
    //! \param[in] i The index of the state vector.
    //! \return The states i-th state vector as a standard library vector.
    std::vector<Real> std_state(Size i) const
    {
      Vector tmp(this->states.row(i));
      return std::vector<Real>(tmp.data(), tmp.data() + tmp.size());
    }

    //! Get a map of the states vectors as standard library vectors.
    //! \return A map of the states vectors as standard library vectors.
    std::map<Size, std::vector<Real>> std_states() const
    {
      std::map<Size, std::vector<Real>> states_map;
      for (Size i = 0; i < N; ++i) {states_map[i] = this->std_states(i);}
      return states_map;
    }
    //! Get a map of the states vectors as Eigen library vectors.
    //! \return A map of the states vectors as Eigen library vectors.
    std::map<Size, Vector> eig_states() const
    {
      std::map<Size, Vector> states_map;
      for (Size i = 0; i < N; ++i) {states_map[i] = this->states.row(i);}
      return states_map;
    }

    //! Get a map of the states vectors as standard vectors.
    //! \param[in] names The names of the states vectors.
    //! \return A map of the states vectors as standard vectors.
    std::map<std::string, std::vector<Real>> std_states(std::vector<std::string> names) const
    {
      std::map<std::string, std::vector<Real>> states_map;
      for (Size i = 0; i < N; ++i) {states_map[names[i]] = this->std_state(i);}
      return states_map;
    }

    //! Get a map of the states vectors as Eigen library vectors.
    //! \param[in] names The names of the states vectors.
    //! \return A map of the states vectors as Eigen library vectors.
    std::map<std::string, Vector> eig_states(std::vector<std::string> names) const
    {
      std::map<std::string, Vector> states_map;
      for (Size i = 0; i < N; ++i) {states_map[names[i]] = this->states.row(i);}
      return states_map;
    }
  };

  /*\
   |   ____                         _  __     _   _
   |  |  _ \ _   _ _ __   __ _  ___| |/ /   _| |_| |_ __ _
   |  | |_) | | | | '_ \ / _` |/ _ \ ' / | | | __| __/ _` |
   |  |  _ <| |_| | | | | (_| |  __/ . \ |_| | |_| || (_| |
   |  |_| \_\\__,_|_| |_|\__, |\___|_|\_\__,_|\__|\__\__,_|
   |                     |___/
  \*/

  //! \brief Class container for the generic \em implicit, \em explicit, and <em> diagonally implicit +
  //! </em> Runge-Kutta methods.
  //!
  //! \includedoc docs/markdown/RungeKutta.md
  template <Size S, Size N>
  class RungeKutta
  {
    using VectorK  = Eigen::Vector<Real, N*S>;      //!< Templetized vector type.
    using MatrixK  = Eigen::Matrix<Real, N, S>;     //!< Templetized matrix type.
    using MatrixJ  = Eigen::Matrix<Real, N*S, N*S>; //!< Templetized matrix type.
    using VectorS  = typename Tableau<S>::Vector;   //!< Templetized vector type.
    using MatrixS  = typename Tableau<S>::Matrix;   //!< Templetized matrix type.
    using VectorN  = typename Implicit<N>::Vector;  //!< Templetized vector type.
    using MatrixN  = typename Implicit<N>::Matrix;  //!< Templetized matrix type.
    using NLSolver = typename NonlinearSolver<N*S>; //!< Templetized nonlinear solver type.

  public:
    using System   = typename Implicit<N>::Pointer;       //!< Shared pointer to an implicit ODE/DAE system.
    using Type     = typename Tableau<S>::Type;           //!< Runge-Kutta type enumeration.
    using Time     = Eigen::Vector<Real, Eigen::Dynamic>; //!< Templetized vector type for the time.
    using Solution = Solution<N>;                         //!< Templetized structure type for solution.

  private:
    NLSolver    *m_nlsolver;                         //!< Nonlinear solver for implicit methods.
    Newton<N*S>  m_newton;                           //!< Newton solver for implicit methods.
    Broyden<N*S> m_broyden;                          //!< Broyden solver for implicit methods.
    Tableau<S>   m_tableau;                          //!< Butcher tableau of the Runge-Kutta method.
    System       m_system;                           //!< ODE/DAE system object pointer.
    Real         m_absolute_tolerance{EPSILON_HIGH}; //!< Absolute tolerance for adaptive step \f$ \epsilon_{\text{abs}} \f$.
    Real         m_relative_tolerance{EPSILON_HIGH}; //!< Relative tolerance for adaptive step \f$ \epsilon_{\text{rel}} \f$.
    Real         m_safety_factor{0.9};               //!< Safety factor for adaptive step \f$ f \f$.
    Real         m_min_safety_factor{0.2};           //!< Minimum safety factor for adaptive step \f$ f_{\max} \f$.
    Real         m_max_safety_factor{1.5};           //!< Maximum safety factor for adaptive step \f$ f_{\min} \f$.
    Real         m_min_step{EPSILON_HIGH};           //!< Minimum step for advancing \f$ h_{\min} \f$.
    bool         m_adaptive_step{false};             //!< Aadaptive step mode boolean.
    bool         m_verbose{false};                   //!< Verbose mode boolean.

  public:
    //! Class constructor for the Runge-Kutta method.
    //! \param[in] t_tableau The Tableau reference.
    RungeKutta(Tableau<S> const &t_tableau)
      : m_tableau(t_tableau) {
      this->m_nlsolver = &this->m_newton;
    }

    //! Class constructor for the Runge-Kutta method.
    //! \param[in] t_tableau The Tableau reference.
    //! \param[in] t_system The ODE/DAE system shared pointer.
    RungeKutta(Tableau<S> const &t_tableau, System t_system)
      : m_tableau(t_tableau), m_system(t_system) {
      this->m_nlsolver = &this->m_newton;
    }

    //! Get the enumeration type of the Runge-Kutta method.
    //! \return The enumeration type of the Runge-Kutta method.
    Type type() const {return this->m_tableau.type;}

    //! Check if the method is an explicit Runge-Kutta (ERK) method.
    //! \return True if the method is an explicit Runge-Kutta method, false otherwise.
    bool is_erk() const {return this->m_tableau.type == Type::ERK;}

    //! Check if the method is an implicit Runge-Kutta (IRK) method.
    //! \return True if the method is an implicit Runge-Kutta method, false otherwise.
    bool is_irk() const {return this->m_tableau.type == Type::IRK;}

    //! Check if the method is a diagonally implicit Runge-Kutta (DIRK) method.
    //! \return True if the method is a diagonally implicit Runge-Kutta method, false otherwise.
    bool is_dirk() const {return this->m_tableau.type == Type::DIRK;}

    //! Get the Butcher Tableau reference.
    //! \return The Tableau reference.
    Tableau<S> & tableau() {return this->m_tableau;}

    //! Get the Butcher Tableau const reference.
    //! \return The Tableau const reference.
    Tableau<S> const & tableau() const {return this->m_tableau;}

    //! Get the stages \f$ s \f$ number of the Runge-Kutta method.
    //! \return The stages \f$ s \f$ number of the Runge-Kutta method.
    Size stages() const {return S;}

    //! Get the name of the Runge-Kutta method.
    //! \return The name of the Runge-Kutta method.
    std::string name() const {return this->m_tableau.name;}

    //! Get the order \f$ p \f$ of the Runge-Kutta method.
    //! \return The order \f$ p \f$ of the Runge-Kutta method.
    Size order() const {return this->m_tableau.order;}

    //! Check if the Runge-Kutta method is embedded.
    //! \return True if the Runge-Kutta method is embedded, false otherwise.
    bool is_embedded() const {return this->m_tableau.is_embedded;}

    //! Get the Butcher tableau matrix \f$ \mathbf{A} \f$.
    //! \return The Butcher tableau matrix \f$ \mathbf{A} \f$.
    MatrixS A() const {return this->m_tableau.A;}

    //! Get the Butcher tableau weights vector \f$ \mathbf{b} \f$.
    //! \return The Butcher tableau weights vector \f$ \mathbf{b} \f$.
    VectorS b() const {return this->m_tableau.b;}

    //! Get the Butcher tableau embedded weights vector \f$ \hat{\mathbf{b}} \f$.
    //! \return The Butcher tableau embedded weights vector \f$ \hat{\mathbf{b}} \f$.
    VectorS b_embedded() const {return this->m_tableau.b_e;}

    //! Get the Butcher tableau nodes vector \f$ \mathbf{c} \f$.
    //! \return The Butcher tableau nodes vector \f$ \mathbf{c} \f$.
    VectorS c() const {return c;}

    //! Get the ODE/DAE system shared pointer.
    //! \return The ODE/DAE system shared pointer.
    System ode() {return this->m_system;}

    //! Set the ODE/DAE system shared pointer.
    //! \param[in] t_system The ODE/DAE system shared pointer.
    void ode(System t_system) {this->m_system = t_system;}

    //! Get the adaptive step absolute tolerance \f$ \epsilon_{\text{abs}} \f$.
    //! \return The adaptive step absolute tolerance \f$ \epsilon_{\text{abs}} \f$.
    Real absolute_tolerance() {return this->m_absolute_tolerance;}

    //! Get the adaptive step absolute tolerance reference \f$ \epsilon_{\text{abs}} \f$.
    //! \param[in] t_absolute_tolerance The adaptive step absolute tolerance \f$ \epsilon_{\text{abs}} \f$.
    void absolute_tolerance(Real t_absolute_tolerance) {this->m_absolute_tolerance = t_absolute_tolerance;}

    //! Get the adaptive step relative tolerance \f$ \epsilon_{\text{rel}} \f$.
    //! \return The adaptive step relative tolerance \f$ \epsilon_{\text{rel}} \f$.
    Real relative_tolerance() {return this->m_relative_tolerance;}

    //! Get the adaptive step relative tolerance reference \f$ \epsilon_{\text{rel}} \f$.
    //! \param[in] t_relative_tolerance The adaptive step relative tolerance \f$ \epsilon_{\text{rel}} \f$.
    void relative_tolerance(Real t_relative_tolerance) {this->m_relative_tolerance = t_relative_tolerance;}

    //! Get the safety factor for adaptive step \f$ f \f$.
    //! \return The safety factor for adaptive step \f$ f \f$.
    Real & safety_factor() {return this->m_safety_factor;}

    //! Set safety factor for adaptive step \f$ f \f$.
    //! \param[in] t_safety_factor The safety factor for adaptive step \f$ f \f$.
    void safety_factor(Real t_safety_factor) {this->m_safety_factor = t_safety_factor;}

    //! Get the minimum safety factor for adaptive step \f$ f_{\min} \f$.
    //! \return The minimum safety factor for adaptive step \f$ f_{\min} \f$.
    Real & min_safety_factor() {return this->m_min_safety_factor;}

    //! Set the minimum safety factor for adaptive step \f$ f_{\min} \f$.
    //! \param[in] t_min_safety_factor The minimum safety factor for adaptive step \f$ f_{\min} \f$.
    void min_safety_factor(Real t_min_safety_factor) {this->m_min_safety_factor = t_min_safety_factor;}

    //! Get the maximum safety factor for adaptive step \f$ f_{\max} \f$.
    //! \return The maximum safety factor for adaptive step \f$ f_{\max} \f$.
    Real & max_safety_factor() {return this->m_max_safety_factor;}

    //! Set the maximum safety factor for adaptive step \f$ f_{\max} \f$.
    //! \param[in] t_max_safety_factor The maximum safety factor for adaptive step \f$ f_{\max} \f$.
    void max_safety_factor(Real t_max_safety_factor) {this->m_max_safety_factor = t_max_safety_factor;}

    //! Get the minimum step for advancing \f$ h_{\min} \f$.
    //! \return The minimum step for advancing \f$ h_{\min} \f$.
    Real & min_step() {return this->m_min_step;}

    //! Set the minimum step for advancing \f$ h_{\min} \f$.
    //! \param[in] t_min_step The minimum step for advancing \f$ h_{\min} \f$.
    void min_step(Real t_min_step) {this->m_min_step = t_min_step;}

    //! Enable verbose mode.
    //!
    void enable_verbose_mode() {this->m_verbose = true;}

    //! Disable verbose mode.
    //!
    void disable_verbose_mode() {this->m_verbose = false;}

    //! Enable adaptive step mode.
    //!
    void enable_adaptive_step() {this->m_adaptive_step = true;}

    //! Disable adaptive step mode.
    //!
    void disable_adaptive_step() {this->m_adaptive_step = false;}

    //! Compute time step for the next advancing step according to the error control method. The
    //! error control method used is the local truncation error method, which is based on the
    //! following formula
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
    //! \f$ h_{k+1} \f$, The error is compared to \f$ 1 \f$ in order to find
    //! an optimal step size. From the error behaviour \f$ e \approx Ch^{q+1} \f$
    //! and from \f$ 1 \approx Ch_{opt}^{q+1} \f$ (where \f$ q = \min(p,\hat{p}) \f$)
    //! the optimal step size is obtained as
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
    //! h_{new} = h \min\left(f_{\max}, \max\left(f_{\min}, f \left(\frac{1}{e}\right)^{\frac{1}{q+1}}
    //! \right) \right)
    //! \f]
    //!
    //! for the new step size. Then, if \f$ e \leq 1 \f$, the computed step is
    //! accepted and the solution is advanced to \f$ \mathbf{x} \f$ and a new step
    //! is tried with \f$ h_{new} \f$ as step size. Else, the step is rejected
    //! and the computations are repeated with the new step size \f$ h_{new} \f$.
    //! Typially, \f$ f \f$ is set in the interval \f$ [0.8, 0.9] \f$,
    //! \f$ f_{\max} \f$ is set in the interval \f$ [1.5, 5] \f$, and \f$ f_{min} \f$
    //! is set in the interval \f$ [0.1, 0.2] \f$.
    //!
    //! \note The error control method is only available for adaptive step mode. Furthermore, the
    //! implementation of the error control method \em should be overridden in the derived class to
    //! provide a more accurate and efficient error control based on the specific Runge-Kutta method.
    //! \param[in] x States approximation \f$ \mathbf{x}_{k+1} \f$.
    //! \param[in] x_e States approximation \f$ \hat{\mathbf{x}}_{k+1} \f$ computed with the embedded
    //! weights vector \f$ \hat{\mathbf{b}} \f$.
    //! \param[in] h_k Actual advancing time step \f$ h_k \f$.
    //! \return The suggested time step for the next integration step \f$ h_{k+1}^\star \f$.
    Real estimate_step(VectorN const &x, VectorN const &x_e, Real h_k) const
    {
      return h_k * std::min(this->m_max_safety_factor, std::max(this->m_min_safety_factor,
        this->m_safety_factor * ((x - x_e) / (this->m_absolute_tolerance + this->m_relative_tolerance*std::max(x.cwiseAbs().maxCoeff(), x_e.cwiseAbs().maxCoeff())
        )).cwiseAbs().maxCoeff()));
    }

    //! Print the Runge-Kutta method information to the output stream.
    //! \param os The output stream.
    void info(std::ostream &os) const {
      os
        << "Runge-Kutta method:\t" << this->name() << std::endl
        << "\t- order:\t" << this->order() << std::endl
        << "\t- stages:\t" << this->stages() << std::endl
        << "\t- type:\t";
      switch (this->type()) {
        case Type::ERK: os << "explicit"; break;
        case Type::IRK: os << "implicit"; break;
        case Type::DIRK: os << "diagonally implicit"; break;
      }
      os
        << std::endl
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
    //!   \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, t_k + h_k c_i
    //!   \right), \qquad i = 1, 2, \ldots, s \text{,} \\
    //! \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i \text{,}
    //! \end{array}
    //! \f]
    //!
    //! Beacuse of the nature of the matrix \f$ \mathbf{A} \f$ (lower triangular)
    //! the \f$ s\f$ stages for a generic explicit Runge-Kutta method take the
    //! form:
    //!
    //! \f[
    //! \mathbf{K}_i = \mathbf{f} \left(
    //!   \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{i-1} a_{ij}
    //!   \mathbf{K}_j, t_k + h_k c_i
    //!   \right), \qquad i = 1, 2, \ldots, s \text{.}
    //! \f]
    //!
    //! Then the explicit Runge-Kutta method for an implicit system of the form
    //! \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{0} \f$ can be
    //! written as
    //!
    //! \f[
    //! \begin{array}{l}
    //! \mathbf{F}_i \left(
    //!   \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{i-1} a_{ij}
    //!     \mathbf{K}_j, \mathbf{K}_i, t_k + h_k c_i
    //! \right) = \mathbf{0}, \qquad i = 1, 2, \ldots, s \text{,} \\
    //! \mathbf{x}_{k+1} = \mathbf{x}_k + \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i.
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
    //! \mathbf{F}_1 \left(\mathbf{x}_k, \mathbf{K}_1, t_k + c_1 h_k\right) = \mathbf{0} \\
    //! \mathbf{F}_2 \left(\mathbf{x}_k + h_k a_{21} \mathbf{K}_1, \mathbf{K}_2, t_k + c_2 h_k
    //! \right) = \mathbf{0} \\
    //! ~~ \vdots \\
    //! \mathbf{F}_s \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, \mathbf{K}_s, t_k + c_s h_k \right) = \mathbf{0}
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
    //! The suggested time step for the next advancing step \f$ h_{k+1} \f$,
    //! is the same as the input time step \f$ h_k \f$ since in the explicit
    //! Runge-Kutta method the time step is not modified through any error control
    //! method.
    //!
    //! \param[in] x_k States at \f$ k \f$-th step \f$ \mathbf{x}_k \f$.
    //! \param[in] t_k Time \f$ t_k \f$.
    //! \param[in] t_step Advancing time step \f$ h_k\f$.
    //! \param[out] x_out States at \f$ k+1 \f$-th step \f$ \mathbf{x}_{k+1} \f$.
    //! \param[out] t_adap The suggested adaptive time step for the next advancing step
    //!               \f$ h_{k+1} \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool erk_step(VectorN const &x_k, Real t_k, Real t_step, VectorN &x_out, Real &t_adap) const
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
    bool dirk_step(VectorN const &x_k, Real t_k, Real t_step, VectorN &x_out, Real &t_adap) const
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
    //! \param[in] x_k States at \f$ k \f$-th step \f$ \mathbf{x}_k \f$.
    //! \param[in] K_k Variables \f$ \mathbf{K} \f$ of the system to be solved.
    //! \param[in] t_k Time \f$ t_k \f$.
    //! \param[in] d_t Advancing time step \f$ \Delta t \f$.
    //! \param[out] fun The residual of system to be solved.
    void irk_function(VectorN const &x_k, VectorK const &K_k, Real t_k, Real d_t, VectorK &fun) const
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
    //! \param[in] x_k States at \f$ k \f$-th step \f$ \mathbf{x}_k \f$.
    //! \param[in] K_k Variables \f$ \mathbf{K} \f$ of the system to be solved.
    //! \param[in] t_k Time \f$ t_k \f$.
    //! \param[in] d_t Advancing time step \f$ \Delta t \f$.
    //! \param[out] jac The Jacobian of the system of equations to be solved.
    void irk_jacobian(VectorN const &x_k, VectorK const &K_k, Real t_k, Real d_t, MatrixJ &jac) const
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
    //! \f$ \mathbf{x}^{\prime} = \mathbf{f}(\mathbf{x}, t) \f$:
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
    //! as
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
    //!   \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^s a_{1j}
    //!   \mathbf{K}_j, \, \mathbf{K}_1, \, t_k + h_k c_1
    //! \right) = \mathbf{0} \\
    //! \mathbf{F}_2 \left(
    //!   \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^s a_{2j}
    //!   \mathbf{K}_j, \, \mathbf{K}_2, \, t_k + h_k c_2
    //! \right) = \mathbf{0} \\
    //! ~~ \vdots \\
    //! \mathbf{F}_s \left(
    //!   \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^s a_{sj}
    //!   \mathbf{K}_j, \, \mathbf{K}_s, \, t_k + h_k c_s
    //! \right) = \mathbf{0}
    //! \end{array}\right.
    //! \f]
    //!
    //! The \f$ \mathbf{K} \f$ variables are computed using the Newton's method.
    //!
    //! The suggested time step for the next advancing step
    //! \f$ h_{k+1} \f$, is the same as the input time step
    //! \f$ \Delta t \f$ since in the implicit Runge-Kutta method the time step
    //! is not modified through any error control method.
    //!
    //! \param[in] x_k States at \f$ k \f$-th step \f$ \mathbf{x}_k \f$.
    //! \param[in] t_k Time \f$ t_k \f$.
    //! \param[in] d_t Advancing time step \f$ \Delta t\f$.
    //! \param[out] x_o States at \f$ k+1 \f$-th step \f$ \mathbf{x}_{k+1} \f$.
    //! \param[out] d_t_star Suggested time step for the next integration step \f$ h_{k+1}^\star \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool irk_step(VectorN const &x_k, Real t_k, Real d_t, VectorN &x_o, Real &d_t_star)
    {
      // Check if the solver converged
      VectorK K;

      if (!this->m_newton.solve(
        [this, &x_k, t_k, d_t](VectorK const &K, VectorK &fun) {this->irk_function(x_k, K, t_k, d_t, fun);},
        [this, &x_k, t_k, d_t](VectorK const &K, MatrixJ &jac) {this->irk_jacobian(x_k, K, t_k, d_t, jac);},
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
    //! form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{0}
    //! \f$. The step is based on the following formula:
    //!
    //! \f[ \mathbf{x}_{k+1} = \mathbf{x}_k + \mathcal{S}(\mathbf{x}_k, \mathbf{x}^{\prime}_k, t_k, h_k) \f]
    //!
    //! where \f$ \mathcal{S} \f$ is the generic advancing step of the solver.
    //!
    //! \param x_k States at \f$ k \f$-th step \f$ \mathbf{x}_k \f$.
    //! \param t_k Time \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ h_k\f$.
    //! \param x_o States at \f$ k+1 \f$-th step \f$ \mathbf{x}_{k+1} \f$.
    //! \param d_t_star Suggested time step for the next integration step \f$ h_{k+1}^\star \f$.
    //! \return True if the step is successfully computed, false otherwise.
    bool step(VectorN const &x_k, Real t_k, Real d_t, VectorN &x_o, Real &d_t_star)
    {
      if (this->is_erk() && this->m_system->is_explicit()) {
        return this->erk_step(x_k, t_k, d_t, x_o, d_t_star);
      //} else if (this->is_dirk()) {
      //  return this->dirk_step(x_k, t_k, d_t, x_o, d_t_star);
      } else {
        return this->irk_step(x_k, t_k, d_t, x_o, d_t_star);
      }
    }

    //! Advance using a generic integration method for a system of the form
    //! \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^{\prime}, t) = \mathbf{0} \f$.
    //! The step is based on the following formula:
    //!
    //! \f[
    //! \mathbf{x}_{k+1}(t_k + h_k) = \mathbf{x}_k(t_k) +
    //! \mathbf{S}(\mathbf{x}_k, \mathbf{x}^{\prime}_k, t_k, \Delta t)
    //! \f]
    //!
    //! where \f$ \mathbf{S} \f$ is the generic advancing step of the solver.
    //! In the advvancing step, the system solution is also projected on the
    //! manifold \f$ \mathbf{h}(\mathbf{x}, t) \f$. Substepping is
    //! also used to ensure that the solution is accurate.
    //!
    //! \param x_k States at \f$ k \f$-th step \f$ \mathbf{x}_k \f$.
    //! \param t_k Time \f$ t_k \f$.
    //! \param d_t Advancing time step \f$ \Delta t\f$.
    //! \param x_new States at \f$ k+1 \f$-th step \f$ \mathbf{x}_{k+1} \f$.
    //! \param d_t_star Suggested time step for the next integration step \f$ h_{k+1}^\star \f$.
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

    //! Solve the system and calculate the approximate solution over the mesh of time points.
    //! \param[in] t Time mesh points \f$ \mathbf{t} = \left[ t_1, t_2, \ldots, t_n \right]^T \f$.
    //! \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    //! \param[out] sol A matrix \f$ \left[(\mathrm{size}(\mathbf{x}) \times \mathrm{size}(\mathbf{t})
    //! \right] \f$ containing the approximated solution over the mesh of time points.
    //! \return True if the system is successfully solved, false otherwise.
    bool solve(VectorX const &t, VectorN const &ics, Solution &sol)
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

    //! Solve the system and calculate the approximate solution over the suggested mesh of time
    //! points with adaptive step size control.
    //! \param[in] t Time mesh points \f$ \mathbf{t} = \left[ t_1, t_2, \ldots, t_n \right]^T \f$.
    //! \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    //! \param[out] sol A matrix \f$ \left[(\mathrm{size}(\mathbf{x}) \times \mathrm{size}(\mathbf{t})
    //! \right] \f$ containing the approximated solution over the mesh of time points.
    //! \return True if the system is successfully solved, false otherwise.
    bool adaptive_solve(VectorX const &t, VectorN const &ics, Solution &sol)
    {
      using Eigen::last;

      #define CMD "Sandals::RungeKutta::adaptive_solve(...): "

      // Instantiate output
      this->enable_adaptive_step();
      Real d_t{t(1) - t(0)}, d_t_star;
        Real scale{100.0};
        std::cout << "d_t = " << d_t << std::endl;

        Real t_min{std::max(this->m_min_step, d_t/scale)};
        Real t_max{scale*d_t};
        Size safety_length;
      if (this->m_tableau.is_embedded) {
        safety_length = std::ceil(Real(1.5)/this->m_min_safety_factor) * t.size();
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
