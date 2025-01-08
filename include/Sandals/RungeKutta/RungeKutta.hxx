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


    //! Check Butcher tableau consistency for an explicit Runge-Kutta method.
    //! \param[in] tbl.A   Matrix \f$ \mathbf{A} \f$.
    //! \param[in] tbl.b   Weights vector \f$ \mathbf{b} \f$.
    //! \param[in] tbl.b_e [optional] Embedded weights vector \f$ \mathbf{b}_{e} \f$.
    //! \param[in] tbl.c   Nodes vector \f$ \mathbf{c} \f$.
    //! \return True if the Butcher tableau is consistent, false otherwise.
    //[out,order,e_order] = check_tableau( this, tbl )

    //! Check the order of a Runge-Kutta tableau according to the conditions taken from:
    //! *A family of embedded Runge-Kutta formulae*, J. R. Dormand and P. J. Prince,
    //! Journal of Computational and Applied Mathematics, volume 6(1), 1980.
    //! Doi: [10.1016/0771-0509(80)90013-3](https://doi.org/10.1016/0771-0509(80)90013-3)
    //! \param[in] tableau The Runge-Kutta tableau to be checked.
    //! \return The order of the Runge-Kutta tableau.
    /*Size tableau_order(matS const &A, vecS const &b, vecS const &c) const
    {
      #define CMD "Sandals::RungeKutta::tableau_order(...): "

      // Temporary variables initialization
      Real tol{std::pow(EPSILON, Real(2.0/3.0))};
      vecN one{vecN::Ones()};
      vecS Ac{A*c};
      matS bA{(b*A).transpose()};
      vecS err{A*one - c};
      Size order{0};
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
  template <Size N, Size M = 0>
  struct Solution
  {
    using Vector  = Eigen::Vector<Real, Eigen::Dynamic>;    //!< Templetized vector type.
    using MatrixN = Eigen::Matrix<Real, N, Eigen::Dynamic>; //!< Templetized matrix type.
    using MatrixM = Eigen::Matrix<Real, M, Eigen::Dynamic>; //!< Templetized matrix type.

    Vector time;       //!< Time solution vector.
    MatrixN states;     //!< States solution matrix.
    MatrixM invariants; //!< Invariants solution matrix.

    //! Default class constructor for the Solution class.
    //!
    Solution() : time(0), states(Vector::Zero(N, 0)), invariants(Vector::Zero(M, 0)) {}

    //! Class constructor for the Solution class given the time and states sizes.
    //! \param[in] size The size of the time vector and states matrix.
    Solution(Size size)
      : time(size), states(Vector::Zero(N, size)), invariants(Vector::Zero(M, size)) {}

    //! Resize the time vector and states matrix.
    //! \param[in] size The size of the time vector and states matrix.
    void resize(Size size) {
      this->time.resize(size);
      this->states.resize(Eigen::NoChange, size);
      this->invariants.resize(Eigen::NoChange, size);
    }

    //! Resize conserving the data of the time vector and states matrix.
    //! \param[in] size The size of the time vector and states matrix.
    void conservative_resize(Size size) {
      this->time.conservativeResize(size);
      this->states.conservativeResize(Eigen::NoChange, size);
      this->invariants.conservativeResize(Eigen::NoChange, size);
    }

    //! Clear the time vector and states matrix.
    //!
    void clear()
    {
      this->time.resize(0);
      this->states.resize(Eigen::NoChange, 0);
      this->invariants.resize(Eigen::NoChange, 0);
    }

    //! Check if the time vector and states matrix are empty.
    //! \return True if the time vector and states matrix are empty, false otherwise.
    bool is_empty() const
    {
      return this->time.size() == 0 && this->states.cols()== 0 && this->invariants.cols() == 0;
    }

    //! Get the size of the time vector and states matrix.
    //! \return The size of the time vector and states matrix.
    Size size() const {return this->time.size();}

    //! Get the time vector as a standard library vector.
    //! \return The time vector as a standard library vector.
    std::vector<Real> std_time() const
    {
      return std::vector<Real>(this->time.data(), this->time.data() + this->time.size());
    }

    //! Get a the time vector as Eigen library vector.
    //! \return The time vector as Eigen library vector.
    Vector eig_time() const {return this->time;}

    //! Get the i-th state vector as a standard library vector.
    //! \param[in] i The index of the state vector.
    //! \return The i-th state vector as a standard library vector.
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
      for (Size i{0}; i < N; ++i) {states_map[i] = this->std_states(i);}
      return states_map;
    }

    //! Get a map of the states vectors as Eigen library vectors.
    //! \return A map of the states vectors as Eigen library vectors.
    std::map<Size, Vector> eig_states() const
    {
      std::map<Size, Vector> states_map;
      for (Size i{0}; i < N; ++i) {states_map[i] = this->states.row(i);}
      return states_map;
    }

    //! Get a map of the states vectors as standard vectors.
    //! \param[in] names The names of the states vectors.
    //! \return A map of the states vectors as standard vectors.
    std::map<std::string, std::vector<Real>> std_states(std::vector<std::string> names) const
    {
      std::map<std::string, std::vector<Real>> states_map;
      for (Size i{0}; i < N; ++i) {states_map[names[i]] = this->std_state(i);}
      return states_map;
    }

    //! Get a map of the states vectors as Eigen library vectors.
    //! \param[in] names The names of the states vectors.
    //! \return A map of the states vectors as Eigen library vectors.
    std::map<std::string, Vector> eig_states(std::vector<std::string> names) const
    {
      std::map<std::string, Vector> states_map;
      for (Size i{0}; i < N; ++i) {states_map[names[i]] = this->states.row(i);}
      return states_map;
    }

    //! Get the i-th invariants vector as a standard library vector.
    //! \param[in] i The index of the invariant vector.
    //! \return The i-th invariants vector as a standard library vector.
    std::vector<Real> std_invariant(Size i) const
    {
      Vector tmp(this->invariants.row(i));
      return std::vector<Real>(tmp.data(), tmp.data() + tmp.size());
    }

    //! Get a map of the invariants vectors as standard library vectors.
    //! \return A map of the invariants vectors as standard library vectors.
    std::map<Size, std::vector<Real>> std_invariants() const
    {
      std::map<Size, std::vector<Real>> invariants_map;
      for (Size i{0}; i < M; ++i) {invariants_map[i] = this->std_invariants(i);}
      return invariants_map;
    }

    //! Get a map of the invariants vectors as Eigen library vectors.
    //! \return A map of the invariants vectors as Eigen library vectors.
    std::map<Size, Vector> eig_invariants() const
    {
      std::map<Size, Vector> invariants_map;
      for (Size i{0}; i < M; ++i) {invariants_map[i] = this->invariants.row(i);}
      return invariants_map;
    }

    //! Get a map of the invariants vectors as standard vectors.
    //! \param[in] names The names of the invariants vectors.
    //! \return A map of the invariants vectors as standard vectors.
    std::map<std::string, std::vector<Real>> std_invariants(std::vector<std::string> names) const
    {
      std::map<std::string, std::vector<Real>> invariants_map;
      for (Size i{0}; i < M; ++i) {invariants_map[names[i]] = this->std_invariant(i);}
      return invariants_map;
    }

    //! Get a map of the invariants vectors as Eigen library vectors.
    //! \param[in] names The names of the invariants vectors.
    //! \return A map of the invariants vectors as Eigen library vectors.
    std::map<std::string, Vector> eig_invariants(std::vector<std::string> names) const
    {
      std::map<std::string, Vector> invariants_map;
      for (Size i{0}; i < M; ++i) {invariants_map[names[i]] = this->invariants.row(i);}
      return invariants_map;
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
  template <Size S, Size N, Size M = 0>
  class RungeKutta
  {
    using VectorK = Eigen::Vector<Real, N*S>;         //!< Templetized vector type.
    using MatrixK = Eigen::Matrix<Real, N, S>;        //!< Templetized matrix type.
    using MatrixJ = Eigen::Matrix<Real, N*S, N*S>;    //!< Templetized matrix type.
    using SolverK = NonlinearSolver<N*S>;             //!< Templetized nonlinear solver type for IRK methods.
    using SolverN = NonlinearSolver<N>;               //!< Templetized nonlinear solver type for ERK and DIRK methods.
    using VectorS = typename Tableau<S>::Vector;      //!< Templetized vector type.
    using MatrixS = typename Tableau<S>::Matrix;      //!< Templetized matrix type.
    using VectorN = typename Implicit<N, M>::VectorN; //!< Templetized vector type.
    using MatrixN = typename Implicit<N, M>::MatrixN; //!< Templetized matrix type.
    using VectorM = typename Implicit<N, M>::VectorM; //!< Templetized vector type.
    using MatrixM = typename Implicit<N, M>::MatrixM; //!< Templetized matrix type.
    using VectorP = Eigen::Matrix<Real, N+M, 1>;      //!< Templetized vector type.
    using MatrixP = Eigen::Matrix<Real, N+M, N+M>;    //!< Templetized matrix type.

  public:
    using System   = typename Implicit<N>::Pointer;       //!< Shared pointer to an implicit ODE/DAE system.
    using Type     = typename Tableau<S>::Type;           //!< Runge-Kutta type enumeration.
    using Time     = Eigen::Vector<Real, Eigen::Dynamic>; //!< Templetized vector type for the time.
    using Solution = Solution<N>;                         //!< Templetized structure type for solution.

  private:
    SolverN     *m_solverN;                          //!< Nonlinear solver for ERK and DIRK methods.
    Newton<N>    m_newtonN;                          //!< Newton solver for ERK and DIRK methods.
    Broyden<N>   m_broydenN;                         //!< Broyden solver for ERK and DIRK methods.
    SolverK     *m_solverK;                          //!< Nonlinear solver for IRK methods.
    Newton<N*S>  m_newtonK;                          //!< Newton solver for IRK methods.
    Broyden<N*S> m_broydenK;                         //!< Broyden solver for IRK methods.
    Tableau<S>   m_tableau;                          //!< Butcher tableau of the Runge-Kutta method.
    System       m_system;                           //!< ODE/DAE system object pointer.
    Real         m_absolute_tolerance{EPSILON_HIGH}; //!< Absolute tolerance for adaptive step \f$ \epsilon_{\text{abs}} \f$.
    Real         m_relative_tolerance{EPSILON_HIGH}; //!< Relative tolerance for adaptive step \f$ \epsilon_{\text{rel}} \f$.
    Real         m_safety_factor{0.9};               //!< Safety factor for adaptive step \f$ f \f$.
    Real         m_min_safety_factor{0.2};           //!< Minimum safety factor for adaptive step \f$ f_{\max} \f$.
    Real         m_max_safety_factor{1.5};           //!< Maximum safety factor for adaptive step \f$ f_{\min} \f$.
    Real         m_min_step{EPSILON_HIGH};           //!< Minimum step for advancing \f$ h_{\min} \f$.
    Size         m_max_substeps{5};                  //!< Maximum number of substeps.
    bool         m_adaptive_step{false};             //!< Aadaptive step mode boolean.
    bool         m_verbose{false};                   //!< Verbose mode boolean.

    Eigen::FullPivLU<MatrixP> m_lu;                         //!< LU decomposition for the projection matrix.
    Real              m_projection_tolerance{EPSILON_HIGH}; //!< Projection tolerance \f$ \epsilon_{\text{proj}} \f$.
    Size              m_max_projection_iterations{5};       //!< Maximum number of projection steps.
    bool              m_projection{true};                   //!< Aadaptive step mode boolean.


  public:
    //! Class constructor for the Runge-Kutta method.
    //! \param[in] t_tableau The Tableau reference.
    RungeKutta(Tableau<S> const &t_tableau)
      : m_tableau(t_tableau) {
      this->m_solverN = &this->m_newtonN;
      this->m_solverK = &this->m_newtonK;
    }

    //! Class constructor for the Runge-Kutta method.
    //! \param[in] t_tableau The Tableau reference.
    //! \param[in] t_system The ODE/DAE system shared pointer.
    RungeKutta(Tableau<S> const &t_tableau, System t_system)
      : m_tableau(t_tableau), m_system(t_system) {
      this->m_solverN = &this->m_newtonN;
      this->m_solverK = &this->m_newtonK;
    }

    //! Enable the Newton nonlinear solver for implicit methods.
    void enable_newton() {
      this->m_solverN = &this->m_newtonN;
      this->m_solverK = &this->m_newtonK;
    }

    //! Enable the Broyden nonlinear solver for implicit methods.
    void enable_broyden() {
      this->m_solverN = &this->m_broydenN;
      this->m_solverK = &this->m_broydenK;
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

    //! Get the maximum number of substeps.
    //! \return The maximum number of substeps.
    Size & max_substeps() {return this->m_max_substeps;}

    //! Set the maximum number of substeps.
    //! \param[in] t_max_substeps The maximum number of substeps.
    void max_substeps(Size t_max_substeps) {this->m_max_substeps = t_max_substeps;}

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

    //! Get the projection tolerance.
    //! \return The projection tolerance.
    Real & projection_tolerance() {return this->m_projection_tolerance;}

    //! Set the projection tolerance.
    //! \param[in] t_projection_tolerance The projection tolerance.
    void projection_tolerance(Real t_projection_tolerance) {this->m_projection_tolerance = t_projection_tolerance;}

    //! Get the maximum number of projection iterations.
    //! \return The maximum number of projection iterations.
    Size & max_projection_iterations() {return this->m_max_projection_iterationsations;}

    //! Set the maximum number of projection iterations.
    //! \param[in] t_max_projection_iterations The maximum number of projection iterations.
    void max_projection_iterations(Size t_max_projection_iterations) {this->m_max_projection_iterationsations = t_max_projection_iterations;}

    //! Enable projection mode.
    //!
    void enable_projection() {this->m_projection = true;}

    //! Disable projection mode.
    //!
    void disable_projection() {this->m_projection = false;}

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
        this->m_safety_factor * ((x - x_e) / (this->m_absolute_tolerance + this->m_relative_tolerance
        * std::max(x.cwiseAbs().maxCoeff(), x_e.cwiseAbs().maxCoeff())
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

    /*\
     |   _____ ____  _  __
     |  | ____|  _ \| |/ /
     |  |  _| | |_) | ' /
     |  | |___|  _ <| . \
     |  |_____|_| \_\_|\_\
     |
    \*/

    //! Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    //! + h_k \f$ using an explicit Runge-Kutta method, given an explicit system of the form \f$
    //! \mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t) \f$. If the method is embedded, the time
    //! step for the next advancing step \f$ h_{k+1}^\star \f$ is also computed according to the
    //! error control method.
    //! \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    //! \param[in] t_old Time \f$ t_k \f$ at the \f$ k \f$-th step.
    //! \param[in] h_old Advancing time step \f$ h_k \f$ at the \f$ k \f$-th step.
    //! \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    //! \param[out] h_new The suggested time step \f$ h_{k+1}^\star \f$ for the next advancing step.
    //! \return True if the step is successfully computed, false otherwise.
    bool erk_explicit_step(VectorN const &x_old, Real t_old, Real h_old, VectorN &x_new, Real &h_new) const
    {
      using Eigen::all;
      using Eigen::seqN;

      // Compute the K variables in the case of an explicit method and explicit system
      VectorN x_node;
      MatrixK K;
      for (Size i{0}; i < S; ++i) {
        x_node = x_old + K(all, seqN(0, i)) * this->m_tableau.A(i, seqN(0, i)).transpose();
        K.col(i) = h_old * static_cast<Explicit<N> const *>(this->m_system.get())->f(x_node, t_old + h_old*this->m_tableau.c(i));
      }
      if (!K.allFinite()) {return false;}

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_tableau.is_embedded) {
        VectorN x_emb = x_old + K * this->m_tableau.b_e;
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    //! Compute the residual of system to be solved, which is given by the values of the system
    //!
    //! \f[
    //! \begin{array}{l}
    //!   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^{n-1} a_{nj}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    //! \end{array} \text{.}
    //! \f]
    //!
    //! where \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$.
    //! \note If \f$ h \gets 0 \f$, then the first guess to solve the nonlinear system
    //! is given by \f$ \tilde{\mathbf{K}} = \mathbf{0} \f$.
    //! \param[in] s Stage index \f$ s \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Time \f$ t \f$.
    //! \param[in] h Advancing time step \f$ h \f$.
    //! \param[in] K Variables \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ of the system to be solved.
    //! \param[out] fun The residual of system to be solved.
    void erk_implicit_function(Size s, VectorN const &x, Real t, Real h, MatrixK const &K, VectorN &fun) const
    {
      using Eigen::seqN;
      using Eigen::all;
      VectorN x_node(x + K(all, seqN(0, s)) * this->m_tableau.A(s, seqN(0, s)).transpose());
      fun = this->m_system->F(x_node, K.col(s)/h, t + h * this->m_tableau.c(s));
    }

    //! Compute the Jacobian of the system of equations
    //!
    //! \f[
    //! \begin{array}{l}
    //!   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^{n-1} a_{nj}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    //! \end{array} \text{.}
    //! \f]
    //!
    //! with respect to the \f$ \tilde{\mathbf{K}}_n = h \mathbf{K}_n \f$ variables. The Jacobian
    //! matrix of the \f$ n \f$-th equation with respect to the \f$ \tilde{\mathbf{K}}_n \f$ variables
    //! is given by
    //!
    //! \f[
    //!   \frac{\partial\mathbf{F}_n}{\partial\tilde{\mathbf{K}}_n}\left(\mathbf{x} + \displaystyle
    //!   \sum_{j=1}^{n-1} a_{nj}\tilde{\mathbf{K}}_j, \displaystyle\frac{\tilde{\mathbf{K}}_n}{h},
    //!   t + h c_n \right) = a_{nj} \mathbf{JF}_x + \displaystyle\frac{1}{h}\begin{cases}
    //!   \mathbf{JF}_{x^{\prime}} & n = j \\
    //!   0 & n \neq j
    //! \end{cases} \text{,}
    //! \f]
    //!
    //! for \f$ j = 1, 2, \ldots, s \f$.
    //! \note If \f$ h \rightarrow 0 \f$, then the first guess to solve the nonlinear system
    //! is given by \f$ \tilde{\mathbf{K}} = \mathbf{0} \f$.
    //! \param[in] s Stage index \f$ s \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Time \f$ t \f$.
    //! \param[in] h Advancing time step \f$ h \f$.
    //! \param[in] K Variables \f$ h \mathbf{K} \f$ of the system to be solved.
    //! \param[out] fun The Jacobian of system to be solved.
    void erk_implicit_jacobian(Size s, VectorN const &x, Real t, Real h, MatrixK const &K, MatrixN &jac) const
    {
      using Eigen::all;
      using Eigen::seqN;
      VectorN x_node(x + K(all, seqN(0, s)) * this->m_tableau.A(s, seqN(0, s)).transpose());
      jac = this->m_system->JF_x_dot(x_node, K.col(s)/h, t + h * this->m_tableau.c(s)) / h;
    }

    //! Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    //! + h_k \f$ using an explicit Runge-Kutta method, given an explicit system of the form \f$
    //! \mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t) \f$. If the method is embedded, the time
    //! step for the next advancing step \f$ h_{k+1}^\star \f$ is also computed according to the
    //! error control method.
    //! \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    //! \param[in] t_old Time \f$ t_k \f$ at the \f$ k \f$-th step.
    //! \param[in] h_old Advancing time step \f$ h_k \f$ at the \f$ k \f$-th step.
    //! \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    //! \param[out] h_new The suggested time step \f$ h_{k+1}^\star \f$ for the next advancing step.
    //! \return True if the step is successfully computed, false otherwise.
    bool erk_implicit_step(VectorN const &x_old, Real t_old, Real h_old, VectorN &x_new, Real &h_new) const
    {
      MatrixK K;
      VectorN K_sol;
      VectorN K_ini(VectorN::Zero());

      // Check if the solver converged
      for (Size s{0}; s < S; ++s) {
        if (this->m_solverN->solve(
            [this, s, &K, &x_old, t_old, h_old](VectorN const &K_fun, VectorN &fun)
              {K.col(s) = K_fun; this->erk_implicit_function(s, x_old, t_old, h_old, K, fun);},
            [this, s, &K, &x_old, t_old, h_old](VectorN const &K_jac, MatrixN &jac)
              {K.col(s) = K_jac; this->erk_implicit_jacobian(s, x_old, t_old, h_old, K, jac);},
            K_ini, K_sol)) {
          K.col(s) = K_sol;
        } else {
          return false;
        }
      }

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_tableau.is_embedded) {
        VectorN x_emb(x_old + K * this->m_tableau.b_e);
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    /*\
     |   ___ ____  _  __
     |  |_ _|  _ \| |/ /
     |   | || |_) | ' /
     |   | ||  _ <| . \
     |  |___|_| \_\_|\_\
     |
    \*/

    //! Compute the residual of system to be solved, which is given by the values of the system
    //!
    //! \f[
    //! \begin{cases}
    //!   \mathbf{F}_1 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{1j}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_1}{h}, t + h c_1\right) \\
    //!   \mathbf{F}_2 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{2j}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_2}{h}, t + h c_2\right) \\[-0.5em]
    //!   \vdots \\[-0.5em]
    //!   \mathbf{F}_s \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{sj}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_s}{h}, t + h c_s\right)
    //! \end{cases}
    //! \f]
    //!
    //! where \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$.
    //! \note If \f$ h \gets 0 \f$, then the first guess to solve the nonlinear system
    //! is given by \f$ \tilde{\mathbf{K}} = \mathbf{0} \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Time \f$ t \f$.
    //! \param[in] h Advancing time step \f$ h \f$.
    //! \param[in] K Variables \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ of the system to be solved.
    //! \param[out] fun The residual of system to be solved.
    void irk_function(VectorN const &x, Real t, Real h, VectorK const &K, VectorK &fun) const
    {
      VectorN x_node; // FIXME: can i avoid this temporary variables?
      MatrixK K_mat{K.reshaped(N, S)}; // FIXME: can i avoid reshaping?
      MatrixK fun_mat;
      for (Size i{0}; i < S; ++i) {
        x_node = x + K_mat * this->m_tableau.A.row(i).transpose();
        fun_mat.col(i) = this->m_system->F(x_node, K_mat.col(i)/h, t + h * this->m_tableau.c(i));
      }
      fun = fun_mat.reshaped(N*S, 1);
    }

    //! Compute the Jacobian of the system of equations
    //!
    //! \f[
    //! \begin{cases}
    //!   \mathbf{F}_1 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{1j}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_1}{h}, t + h c_1\right) \\
    //!   \mathbf{F}_2 \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{2j}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_2}{h}, t + h c_2\right) \\[-0.5em]
    //!   \vdots \\[-0.5em]
    //!   \mathbf{F}_s \left(\mathbf{x} + \displaystyle\sum_{j=1}^{s} a_{sj}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_s}{h}, t + h c_s\right)
    //! \end{cases} \text{.}
    //! \f]
    //!
    //! with respect to the \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ variables. The \f$ ij \f$-th
    //! Jacobian matrix entry is given by
    //!
    //! \f[
    //!   \frac{\partial\mathbf{F}_i}{\partial\tilde{\mathbf{K}}_j}\left(\mathbf{x} + \displaystyle
    //!   \sum_{j=1}^s a_{ij}\tilde{\mathbf{K}}_j, \displaystyle\frac{\tilde{\mathbf{K}}_i}{h}, t +
    //!   h c_i \right) = a_{ij} \mathbf{JF}_x + \displaystyle\frac{1}{h}\begin{cases}
    //!   \mathbf{JF}_{x^{\prime}} & i = j \\
    //!   0 & i \neq j
    //! \end{cases} \text{,}
    //!
    //! \f]
    //! for \f$ i = 1, 2, \ldots, s \f$ and \f$ j = 1, 2, \ldots, s \f$.
    //! \note If \f$ h \rightarrow 0 \f$, then the first guess to solve the nonlinear system
    //! is given by \f$ \tilde{\mathbf{K}} = \mathbf{0} \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Time \f$ t \f$.
    //! \param[in] h Advancing time step \f$ h \f$.
    //! \param[in] K Variables \f$ h \mathbf{K} \f$ of the system to be solved.
    //! \param[out] fun The Jacobian of system to be solved.
    void irk_jacobian(VectorN const &x, Real t, Real h, VectorK const &K, MatrixJ &jac) const
    {
      using Eigen::all;
      using Eigen::seqN;

      jac.setZero();

      // Loop through each equation of the system
      MatrixK K_mat{K.reshaped(N, S)};
      Real t_node;
      VectorN x_node, x_dot_node;
      MatrixN JF_x, JF_x_dot;
      auto idx = seqN(0, N), jdx = seqN(0, N);
      for (Size i{0}; i < S; ++i) {
        t_node = t + h * this->m_tableau.c(i);
        x_node = x + K_mat * this->m_tableau.A.row(i).transpose();

        // Compute the Jacobians with respect to x and x_dot
        x_dot_node = K_mat.col(i) / h;
        JF_x       = this->m_system->JF_x(x_node, x_dot_node, t_node);
        JF_x_dot   = this->m_system->JF_x_dot(x_node, x_dot_node, t_node);

        // Combine the Jacobians with respect to x and x_dot to obtain the Jacobian with respect to K
        idx = seqN(i*N, N);
        for (Size j{0}; j < S; ++j) {
          jdx = seqN(j*N, N);
          if (i == j) {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x + JF_x_dot / h;
          } else {
            jac(idx, jdx) = this->m_tableau.A(i,j) * JF_x;
          }
        }
      }
    }

    //! Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    //! + h_k \f$ using an implicit Runge-Kutta method, given an implicit system of the form \f$
    //! \mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0} \f$. If the method is embedded,
    //! the time step for the next advancing step \f$ h_{k+1}^\star \f$ is also computed according
    //! to the error control method.
    //! \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    //! \param[in] t_old Time \f$ t_k \f$ at the \f$ k \f$-th step.
    //! \param[in] h_old Advancing time step \f$ h_k \f$ at the \f$ k \f$-th step.
    //! \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    //! \param[out] h_new The suggested time step \f$ h_{k+1}^\star \f$ for the next advancing step.
    //! \return True if the step is successfully computed, false otherwise.
    bool irk_step(VectorN const &x_old, Real t_old, Real h_old, VectorN &x_new, Real &h_new)
    {
      VectorK K;

      // Check if the solver converged
      if (!this->m_solverK->solve(
          [this, &x_old, t_old, h_old](VectorK const &K_fun, VectorK &fun)
            {this->irk_function(x_old, t_old, h_old, K_fun, fun);},
          [this, &x_old, t_old, h_old](VectorK const &K_jac, MatrixJ &jac)
            {this->irk_jacobian(x_old, t_old, h_old, K_jac, jac);},
          VectorK::Zero(), K))
        {return false;}

      // Perform the step and obtain the next state
      x_new = x_old + K.reshaped(N, S) * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_tableau.is_embedded) {
        VectorN x_emb(x_old + K.reshaped(N, S) * this->m_tableau.b_e);
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    /*\
     |   ____ ___ ____  _  __
     |  |  _ \_ _|  _ \| |/ /
     |  | | | | || |_) | ' /
     |  | |_| | ||  _ <| . \
     |  |____/___|_| \_\_|\_\
     |
    \*/

    //! Compute the residual of system to be solved, which is given by the values of the system
    //!
    //! \f[
    //! \begin{array}{l}
    //!   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^n a_{nj}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    //! \end{array} \text{.}
    //! \f]
    //!
    //! where \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$.
    //! \note If \f$ h \gets 0 \f$, then the first guess to solve the nonlinear system
    //! is given by \f$ \tilde{\mathbf{K}} = \mathbf{0} \f$.
    //! \param[in] n Equation index \f$ n \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Time \f$ t \f$.
    //! \param[in] h Advancing time step \f$ h \f$.
    //! \param[in] K Variables \f$ \tilde{\mathbf{K}} = h \mathbf{K} \f$ of the system to be solved.
    //! \param[out] fun The residual of system to be solved.
    void dirk_function(Size n, VectorN const &x, Real t, Real h, MatrixK const &K, VectorN &fun) const
    {
      using Eigen::seqN;
      using Eigen::all;
      VectorN x_node(x + K(all, seqN(0, n+1)) * this->m_tableau.A(n, seqN(0, n+1)).transpose());
      fun = this->m_system->F(x_node, K.col(n)/h, t + h * this->m_tableau.c(n));
    }

    //! Compute the Jacobian of the system of equations
    //!
    //! \f[
    //! \begin{array}{l}
    //!   \mathbf{F}_n \left(\mathbf{x} + \displaystyle\sum_{j=1}^n a_{nj}\tilde{\mathbf{K}}_j,
    //!   \displaystyle\frac{\tilde{\mathbf{K}}_n}{h}, t + h c_n\right)
    //! \end{array} \text{.}
    //! \f]
    //!
    //! with respect to the \f$ \tilde{\mathbf{K}}_n = h \mathbf{K}_n \f$ variables. The Jacobian
    //! matrix of the \f$ n \f$-th equation with respect to the \f$ \tilde{\mathbf{K}}_n \f$ variables
    //! is given by
    //!
    //! \f[
    //!   \frac{\partial\mathbf{F}_n}{\partial\tilde{\mathbf{K}}_n}\left(\mathbf{x} + \displaystyle
    //!   \sum_{j=1}^n a_{nj}\tilde{\mathbf{K}}_j, \displaystyle\frac{\tilde{\mathbf{K}}_n}{h},
    //!   t + h c_n \right) = a_{nj} \mathbf{JF}_x + \displaystyle\frac{1}{h}\begin{cases}
    //!   \mathbf{JF}_{x^{\prime}} & n = j \\
    //!   0 & n \neq j
    //! \end{cases} \text{,}
    //! \f]
    //!
    //! for \f$ j = 1, 2, \ldots, s \f$.
    //! \note If \f$ h \rightarrow 0 \f$, then the first guess to solve the nonlinear system
    //! is given by \f$ \tilde{\mathbf{K}} = \mathbf{0} \f$.
    //! \param[in] n Equation index \f$ n \f$.
    //! \param[in] x States \f$ \mathbf{x} \f$.
    //! \param[in] t Time \f$ t \f$.
    //! \param[in] h Advancing time step \f$ h \f$.
    //! \param[in] K Variables \f$ h \mathbf{K} \f$ of the system to be solved.
    //! \param[out] fun The Jacobian of system to be solved.
    void dirk_jacobian(Size n, VectorN const &x, Real t, Real h, MatrixK const &K, MatrixN &jac) const
    {
      using Eigen::all;
      using Eigen::seqN;
      Real t_node{t + h * this->m_tableau.c(n)};
      VectorN x_node(x + K(all, seqN(0, n+1)) * this->m_tableau.A(n, seqN(0, n+1)).transpose());
      VectorN x_dot_node(K.col(n)/h);
      jac = this->m_tableau.A(n,n) * this->m_system->JF_x(x_node, x_dot_node, t_node) +
        this->m_system->JF_x_dot(x_node, x_dot_node, t_node) / h;
    }

    //! Compute the new states \f$ \mathbf{x}_{k+1} \f$ at the next advancing step \f$ t_{k+1} = t_k
    //! + h_k \f$ using an diagonally implicit Runge-Kutta method, given an implicit system of the
    //! form \f$ \mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0} \f$. If the method is
    //! embedded, the time step for the next advancing step \f$ h_{k+1}^\star \f$ is also computed
    //! according to the error control method.
    //! \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    //! \param[in] t_old Time \f$ t_k \f$ at the \f$ k \f$-th step.
    //! \param[in] h_old Advancing time step \f$ h_k \f$ at the \f$ k \f$-th step.
    //! \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    //! \param[out] h_new The suggested time step \f$ h_{k+1}^\star \f$ for the next advancing step.
    //! \return True if the step is successfully computed, false otherwise.
    bool dirk_step(VectorN const &x_old, Real t_old, Real h_old, VectorN &x_new, Real &h_new) const
    {
      MatrixK K;
      VectorN K_sol;
      VectorN K_ini(VectorN::Zero());

      // Check if the solver converged at each step
      for (Size n{0}; n < S; ++n) {
        if (this->m_solverN->solve(
            [this, n, &K, &x_old, t_old, h_old](VectorN const &K_fun, VectorN &fun)
              {K.col(n) = K_fun; this->dirk_function(n, x_old, t_old, h_old, K, fun);},
            [this, n, &K, &x_old, t_old, h_old](VectorN const &K_jac, MatrixN &jac)
              {K.col(n) = K_jac; this->dirk_jacobian(n, x_old, t_old, h_old, K, jac);},
            K_ini, K_sol)) {
            K.col(n) = K_sol;
          } else {
            return false;
          }
      }

      // Perform the step and obtain the next state
      x_new = x_old + K * this->m_tableau.b;

      // Adapt next time step
      if (this->m_adaptive_step && this->m_tableau.is_embedded) {
        VectorN x_emb(x_old + K * this->m_tableau.b_e);
        h_new = this->estimate_step(x_new, x_emb, h_old);
      }
      return true;
    }

    //! Compute a step using a generic integration method for a system of the form \f$ \mathbf{F}(
    //! \mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0} \f$. The step is automatically selected
    //! based on the system properties and the integration method properties.
    //! \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    //! \param[in] t_old Time \f$ t_k \f$ at the \f$ k \f$-th step.
    //! \param[in] h_old Advancing time step \f$ h_k \f$ at the \f$ k \f$-th step.
    //! \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    //! \param[out] h_new The suggested time step \f$ h_{k+1}^\star \f$ for the next advancing step.
    //! \return True if the step is successfully computed, false otherwise.
    bool step(VectorN const &x_old, Real t_old, Real h_old, VectorN &x_new, Real &h_new)
    {
      if (this->is_erk() && this->m_system->is_explicit()) {
        return this->erk_explicit_step(x_old, t_old, h_old, x_new, h_new);
      } else if (this->is_erk() && this->m_system->is_implicit()) {
        return this->erk_implicit_step(x_old, t_old, h_old, x_new, h_new);
      } else if (this->is_dirk()) {
        return this->dirk_step(x_old, t_old, h_old, x_new, h_new);
      } else {
        return this->irk_step(x_old, t_old, h_old, x_new, h_new);
      }
    }

    //! Advance using a generic integration method for a system of the form \mathbf{F}(\mathbf{x},
    //! \mathbf{x}^\prime, t) = \mathbf{0} \f$. The step is automatically selected based on the
    //! system properties and the integration method properties. In the advvancing step, the system
    //! solution is also projected on the manifold \f$ \mathbf{h}(\mathbf{x}, t) \f$. Substepping
    //! may be also used to if the integration step fails with the current step size.
    //! \param[in] x_old States \f$ \mathbf{x}_k \f$ at the \f$ k \f$-th step.
    //! \param[in] t_old Time \f$ t_k \f$ at the \f$ k \f$-th step.
    //! \param[in] h_old Advancing time step \f$ h_k \f$ at the \f$ k \f$-th step.
    //! \param[out] x_new Computed states \f$ \mathbf{x}_{k+1} \f$ at the \f$ (k+1) \f$-th step.
    //! \param[out] h_new The suggested time step \f$ h_{k+1}^\star \f$ for the next advancing step.
    //! \return True if the step is successfully computed, false otherwise.
    bool advance(VectorN const &x_old, Real t_old, Real h_old, VectorN &x_new, Real &h_new)
    {
      #define CMD "Sandals::RungeKutta::advance(...): "

      // Check step size
      SANDALS_ASSERT(h_old > Real(0.0), CMD "in " << this->m_tableau.name << "solver, d_t = "<<
        h_old << ", expected > 0.");

      // If the integration step failed, try again with substepping
      if (!this->step(x_old, t_old, h_old, x_new, h_new))
      {
        VectorN x_tmp(x_old);
        Real t_tmp{t_old}, h_tmp{h_old / Real(2.0)};

        // Substepping logic
        Size max_k{this->m_max_substeps * this->m_max_substeps}, k{2};
        Real h_new_tmp;
        while (k > 0) {
          // Calculate the next time step with substepping logic
          if (this->step(x_tmp, t_tmp, h_tmp, x_new, h_new_tmp)) {

            // Accept the step
            h_tmp = h_new_tmp;

            // If substepping is enabled, double the step size
            if (k > 0 && k < max_k) {
              k = k - 1;
              // If the substepping index is even, double the step size
              if (k % 2 == 0) {
                h_tmp = Real(2.0) * h_tmp;
                if (this->m_verbose) {
                  SANDALS_WARNING(CMD "in " << this->m_tableau.name << " solver, at t = " << t_tmp <<
                    ", integration succedded disable one substepping layer.");
                }
              }
            }

            // Check the infinity norm of the solution
            SANDALS_ASSERT(std::isfinite(x_tmp.maxCoeff()), CMD "in " << this->m_tableau.name <<
              " solver, at t = " << t_tmp << ", ||x||_inf = inf, computation interrupted.");

          } else {

            // If the substepping index is too high, abort the integration
            k = k + 2;
            SANDALS_ASSERT(k < max_k, CMD "in " << this->m_tableau.name << " solver, at t = " <<
              t_tmp << ", integration failed with d_t = " << h_tmp << ", aborting.");
            return false;

            // Otherwise, try again with a smaller step
            if (this->m_verbose) {
              SANDALS_WARNING(CMD "in " << this->m_tableau.name << " solver, at t = " <<
                t_tmp << ", integration failed, adding substepping layer.");
            }
            h_tmp = h_tmp / Real(2.0);
            continue;
          }

          // Store time solution
          t_tmp = t_tmp + h_tmp;
        }

        // Store output states substepping solutions
        x_new = x_tmp;
        h_new = h_tmp;
      }

      // Project intermediate solution on the invariants
      if (this->m_projection) {
        VectorN x_projected;
        if (this->project(x_new, t_old + h_new, x_projected)) {
          x_new = x_projected;
          return true;
        } else {
          return false;
        }
      } else {
        return true;
      }

      #undef CMD
    }

    //! Solve the system and calculate the approximate solution over the time mesh \f$ \mathbf{t} =
    //! \left[ t_1, t_2, \ldots, t_n \right]^T \f$. The step size is fixed and given by \f$ h =
    //! t_{k+1} - t_k \f$.
    //! \param[in] t_mesh Time mesh points \f$ \mathbf{t} \f$.
    //! \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    //! \param[out] sol The solution of the system over the mesh of time points.
    //! \return True if the system is successfully solved, false otherwise.
    bool solve(VectorX const &t_mesh, VectorN const &ics, Solution &sol)
    {
      using Eigen::last;

      // Instantiate output
      sol.resize(t_mesh.size());

      // Store first step
      sol.time(0)       = t_mesh(0);
      sol.states.col(0) = ics;

      // Update the current step
      Size step{0};
      VectorN x_step{ics};
      Real t_step{t_mesh(0)}, h_step{t_mesh(1) - t_mesh(0)}, h_tmp_step{h_step}, h_new_step;
      bool mesh_point_bool, saturation_bool;

      while (true) {
        // Integrate system
        if (!this->advance(sol.states.col(step), t_step, h_step, x_step, h_new_step)) {return false;}

        // Update the current step
        t_step = t_step + h_step;

        // Saturate the suggested timestep
        mesh_point_bool = std::abs(t_step - t_mesh(step+1)) < SQRT_EPSILON;
        saturation_bool = t_step + h_new_step > t_mesh(step+1) + SQRT_EPSILON;
        if (this->m_adaptive_step && !mesh_point_bool && saturation_bool) {
          h_tmp_step = h_new_step;
          h_step     = t_mesh(step+1) - t_step;
        } else {
          h_step = h_new_step;
        }

        // Store solution if the step is a mesh point
        if (!this->m_adaptive_step || mesh_point_bool) {
          // Update temporaries
          step += 1;
          h_step = h_tmp_step;

          // Update outputs
          sol.time(step)       = t_step;
          sol.states.col(step) = x_step;

          // Check if the current step is the last one
          if (std::abs(t_step - t_mesh(last)) < SQRT_EPSILON) {break;}
        }
      }
      return true;
    }

    //! Solve the system and calculate the approximate solution over the suggested time mesh \f$
    //! \mathbf{t} = \left[ t_1, t_2, \ldots, t_n \right]^T \f$, the step size is automatically
    //! computed based on the error control method.
    //! \param[in] t_mesh Time mesh points \f$ \mathbf{t} \f$.
    //! \param[in] ics Initial conditions \f$ \mathbf{x}(t = 0) \f$.
    //! \param[out] sol The solution of the system over the mesh of time points.
    //! \return True if the system is successfully solved, false otherwise.
    bool adaptive_solve(VectorX const &t_mesh, VectorN const &ics, Solution &sol)
    {
      using Eigen::last;

      #define CMD "Sandals::RungeKutta::adaptive_solve(...): "

      // Instantiate output
      this->enable_adaptive_step();
      Real h_step{t_mesh(1) - t_mesh(0)}, h_new_step, scale{100.0};
      Real h_min{std::max(this->m_min_step, h_step/scale)}, h_max{scale*h_step};
      if (this->m_tableau.is_embedded) {
        Size safety_length{static_cast<Size>(std::ceil(std::abs(t_mesh(last) - t_mesh(0))/h_min))};
        sol.resize(safety_length);
      } else {
        sol.resize(t_mesh.size());
      }

      // Store first step
      sol.time(0)      = t_mesh(0);
      sol.states.col(0) = ics;

      // Instantiate temporary variables
      Size step{0};
      VectorN x_step{ics};

      while (true) {
        // Integrate system
        this->advance(sol.states.col(step), sol.time(step), h_step, x_step, h_new_step);

        // Saturate the suggested timestep
        if (this->m_adaptive_step && this->m_tableau.is_embedded) {
          h_step = std::max(std::min(h_new_step, h_max), h_min);
        }

        SANDALS_ASSERT(step < t_mesh.size(), CMD "safety length exceeded.");

        // Store solution
        sol.time(step+1)       = sol.time(step) + h_step;
        sol.states.col(step+1) = x_step;

        // Check if the current step is the last one
        if (sol.time(step+1) + h_step > t_mesh(last)) {break;}

        // Update steps counter
        step += 1;
      }

      // Resize the output
      sol.conservative_resize(step-1);

      return true;

      #undef CMD
    }

    //! Project the system solution \f$ \mathbf{x} \f$ on the invariants \f$ \mathbf{h} (\mathbf{x},
    //! t) = \mathbf{0} \f$.
    //! \param[in] x The initial guess for the states \f$\widetilde{\mathbf{x}} \f$.
    //! \param[in] t The time \f$ t \f$ at which the states are evaluated.
    //! \param[out] x_projected The projected states \f$ \mathbf{x} \f$ closest to the invariants
    //! manifold \f$ \mathbf{h} (\mathbf{x}, t) = \mathbf{0} \f$.
    //! \return True if the solution is successfully projected, false otherwise.
    bool project(VectorN const &x, Real t, VectorN &x_projected) {

      #define CMD "Sandals::RungeKutta::project(...): "

      // Check if there are any constraints
      x_projected = x;
      if (M > Size(0)) {

        VectorM h;
        MatrixM Jh_x;
        VectorP b, x_step;
        MatrixP A;
        for (Size k{0}; k < this->m_max_projection_iterations; ++k) {

          // Standard projection method
          //     [A]         {x}    =        {b}
          // / I  Jh_x^T \ /   dx   \   / x_t - x_k \
          // |           | |        | = |           |
          // \ Jh_x    0 / \ lambda /   \    -h     /

          // Evaluate the invariants vector and its Jacobian
          h    = this->m_system->h(x, t);
          Jh_x = this->m_system->Jh_x(x, t);

          // Check if the solution is found
          if (h.norm() < this->m_projection_tolerance) {return true;}

          // Build the linear system
          A << MatrixN::Identity(), Jh_x.transpose(), Jh_x, MatrixM::Zero();
          b << x - x_projected, -h;

          // Compute the solution of the linear system
          this->m_lu.compute(A);
          SANDALS_ASSERT(this->m_lu.rank() == N+M, CMD "singular Jacobian detected.");
          x_step = this->m_lu.solve(b);

          // Check if the step is too small
          if (x_step.norm() < this->m_projection_tolerance * this->m_projection_tolerance) {return false;} // return false or true?

          // Update the solution
          x_projected += x_step(Eigen::seqN(0, N));
        }
        SANDALS_WARNING(CMD "maximum number of iterations reached.");
        return false;
      } else {
        return true;
      }

      #undef CMD
    }

    //! Project the system solution \f$ \mathbf{x} \f$ on the invariants \f$ \mathbf{h} (\mathbf{x},
    //! t) = \mathbf{0} \f$.
    //! \param[in] x The initial guess for the states \f$\widetilde{\mathbf{x}} \f$.
    //! \param[in] t The time \f$ t \f$ at which the states are evaluated.
    //! \param[in] projected_equations The indices of the states to be projected.
    //! \param[in] projected_invariants The indices of the invariants to be projected.
    //! \param[out] x_projected The projected states \f$ \mathbf{x} \f$ closest to the invariants
    //! manifold \f$ \mathbf{h} (\mathbf{x}, t) = \mathbf{0} \f$.
    //! \return True if the solution is successfully projected, false otherwise.
    bool project_ics(VectorN const &x, Real t, std::vector<Size> const &
      projected_equations, std::vector<Size> const & projected_invariants, VectorN &x_projected) const {

      #define CMD "Sandals::RungeKutta::project_ics(...): "

      Size X{static_cast<Size>(projected_equations.size())};
      Size H{static_cast<Size>(projected_invariants.size())};

      // Check if there are any constraints
      x_projected = x;
      if (H > Size(0)) {

        VectorM h;
        MatrixM Jh_x;
        VectorX b(X+H), x_step(X+H);
        MatrixX A(X+H, X+H);
        Eigen::FullPivLU<MatrixX> lu;
        for (Size k{0}; k < this->m_max_projection_iterations; ++k) {

          // Standard projection method
          //     [A]         {x}    =        {b}
          // / I  Jh_x^T \ /   dx   \   / x_t - x_k \
          // |           | |        | = |           |
          // \ Jh_x    0 / \ lambda /   \    -h     /

          // Evaluate the invariants vector and its Jacobian
          h    = this->m_system->h(x, t);
          Jh_x = this->m_system->Jh_x(x, t);

          // Select only the projected invariants
          h    = h(projected_invariants);
          Jh_x = Jh_x(projected_invariants, projected_equations);

          // Check if the solution is found
          if (h.norm() < this->m_projection_tolerance) {return true;}

          // Build the linear system
          A << MatrixX::Identity(X+H, X+H), Jh_x.trasposed(), Jh_x, MatrixX::Zero(X+H, X+H);
          b << x(projected_equations) - x_projected(projected_equations), -h;

          // Compute the solution of the linear system
          lu.compute(A);
          SANDALS_ASSERT(lu.rank() == X+H, CMD "singular Jacobian detected.");
          x_step = this->m_lu.solve(b);

          // Check if the step is too small
          if (x_step.norm() < this->m_projection_tolerance * this->m_projection_tolerance) {return false;} // return false or true?

          // Update the solution
          x_projected(projected_equations) += x_step;
        }
        SANDALS_WARNING(CMD "maximum number of iterations reached.");
        return false;
      } else {
        return true;
      }

      #undef CMD
    }

  }; // class RungeKutta

} // namespace Sandals

#endif // SANDALS_RUNGEKUTTA_HXX
