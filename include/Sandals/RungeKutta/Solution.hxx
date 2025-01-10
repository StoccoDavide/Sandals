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

#ifndef SANDALS_SOLUTION_HXX
#define SANDALS_SOLUTION_HXX

namespace Sandals {

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
  //! or Differential Algebraic Equations (DAEs). The solution is stored in a independent variable
  //! (or time) vector and states matrix, where each row of the matrix represents the solution of a
  //! state at a given independent variable mesh. The solution can be accessed as standard library
  //! vectors or Eigen library vectors.
  template <Size N, Size M = 0>
  struct Solution
  {
    using Vector  = Eigen::Vector<Real, Eigen::Dynamic>;    //!< Templetized vector type.
    using MatrixN = Eigen::Matrix<Real, N, Eigen::Dynamic>; //!< Templetized matrix type.
    using MatrixM = Eigen::Matrix<Real, M, Eigen::Dynamic>; //!< Templetized matrix type.

    Vector  time;       //!< Independent variable (or time) solution vector.
    MatrixN states;     //!< States solution matrix.
    MatrixM invariants; //!< Invariants solution matrix.

    //! Default class constructor for the Solution class.
    //!
    Solution() : time(0), states(Vector::Zero(N, 0)), invariants(Vector::Zero(M, 0)) {}

    //! Class constructor for the Solution class given the independent variable (or time) and states sizes.
    //! \param[in] size The size of the independent variable (or time) vector and states matrix.
    Solution(Size size)
      : time(size), states(Vector::Zero(N, size)), invariants(Vector::Zero(M, size)) {}

    //! Resize the independent variable (or time) vector and states matrix.
    //! \param[in] size The size of the independent variable (or time) vector and states matrix.
    void resize(Size size) {
      this->time.resize(size);
      this->states.resize(Eigen::NoChange, size);
      this->invariants.resize(Eigen::NoChange, size);
    }

    //! Resize conserving the data of the independent variable (or time) vector and states matrix.
    //! \param[in] size The size of the independent variable (or time) vector and states matrix.
    void conservative_resize(Size size) {
      this->time.conservativeResize(size);
      this->states.conservativeResize(Eigen::NoChange, size);
      this->invariants.conservativeResize(Eigen::NoChange, size);
    }

    //! Clear the independent variable (or time) vector and states matrix.
    //!
    void clear()
    {
      this->time.resize(0);
      this->states.resize(Eigen::NoChange, 0);
      this->invariants.resize(Eigen::NoChange, 0);
    }

    //! Check if the independent variable (or time) vector and states matrix are empty.
    //! \return True if the independent variable (or time) vector and states matrix are empty, false otherwise.
    bool is_empty() const
    {
      return this->time.size() == 0 && this->states.cols()== 0 && this->invariants.cols() == 0;
    }

    //! Get the size of the independent variable (or time) vector and states matrix.
    //! \return The size of the independent variable (or time) vector and states matrix.
    Size size() const {return this->time.size();}

    //! Get the independent variable (or time) vector as a standard library vector.
    //! \return The independent variable (or time) vector as a standard library vector.
    std::vector<Real> std_time() const
    {
      return std::vector<Real>(this->time.data(), this->time.data() + this->time.size());
    }

    //! Get the independent variable (or time) vector as Eigen library vector.
    //! \return The independent variable (or time) vector as Eigen library vector.
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

} // namespace Sandals

#endif // SANDALS_SOLUTION_HXX
