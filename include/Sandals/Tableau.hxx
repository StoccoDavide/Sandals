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

#ifndef SANDALS_TABLEAU_HXX
#define SANDALS_TABLEAU_HXX

namespace Sandals {

  /*\
   |   _____     _     _
   |  |_   _|_ _| |__ | | ___  __ _ _   _
   |    | |/ _` | '_ \| |/ _ \/ _` | | | |
   |    | | (_| | |_) | |  __/ (_| | |_| |
   |    |_|\__,_|_.__/|_|\___|\__,_|\__,_|
   |
  \*/

  /**
  * \brief Struct container for the Butcher tableau of a Runge-Kutta method.
  *
  * \includedoc docs/markdown/Tableau.md
  *
  * \tparam S The size of the Butcher tableau.
  */
  template <Integer S>
  struct Tableau
  {
    using Type = enum class type : Integer {ERK=0, IRK=1, DIRK=2}; /**< Runge-Kutta type enumeration. */
    using Vector = Eigen::Vector<Real, S>;    /**< Templetized vector type. */
    using Matrix = Eigen::Matrix<Real, S, S>; /**< Templetized matrix type. */

    std::string name;                /**< Name of the Runge-Kutta method. */
    Type        type;                /**< Runge-Kutta type. */
    Integer     order;               /**< Order of the Runge-Kutta method. */
    Integer     order_e{-1};         /**< Order of the Runge-Kutta embedded method. */
    Matrix      A;                   /**< Matrix \f$ \mathbf{A} \f$. */
    Vector      b;                   /**< Weights vector \f$ \mathbf{b} \f$. */
    Vector      b_e{Vector::Zero()}; /**< Embedded weights vector \f$ \hat{\mathbf{b}} \f$. */
    Vector      c;                   /**< Nodes vector \f$ \mathbf{c} \f$. */
    bool        is_embedded{false};  /**< Embedded method boolean. */

    /**
    * Check the Butcher tableau consistency for a generic Runge-Kutta method.
    * \param[in] verbose Verbosity flag.
    * \return True if the Butcher tableau is consistent and its order is verified, false otherwise.
    */
    bool check(bool verbose = false) const {

      #define CMD "Sandals::" << this->name << "::check(...): "

      // Check the occupancy of the matrix A
      if (type == Type::ERK && (!this->A.isLowerTriangular() || !this->A.diagonal().isZero())) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "A matrix occupancy not consistent with an ERK method.");
        return false;
      } else if (type == Type::DIRK && (!this->A.isLowerTriangular() || this->A.diagonal().isZero())) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "A matrix occupancy not consistent with a DIRK method.");
        return false;
      } else if (type == Type::IRK && this->A.isZero()) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "A matrix occupancy not consistent with an IRK method.");
        return false;
      }

      // Check the order of the method
      Integer computed_order{compute_order(this->A, this->b, this->c, verbose)};
      if (this->order != computed_order) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "method order check failed, " << computed_order << " ≠ "
          << this->order << ".");
        return false;
      }

      // Check the embedded method consistency
      if (!this->is_embedded && (this->order_e != -1 || !this->b_e.isZero())) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "embedded method inconsistency.");
        return false;
      } else {
        return true;
      }

      // Check the embedded method order
      Integer computed_order_e{compute_order(this->A, this->b_e, this->c, verbose)};
      if (this->order_e != computed_order_e) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "embedded method order check failed, " << computed_order_e
          << " ≠ " << this->order_e << ".");
        return false;
      } else {
        return true;
      }

      #undef CMD
    }

  private:

    /**
    * Check the order of a Runge-Kutta method according to the conditions taken from: *A family of
    * embedded Runge-Kutta formulae*, J. R. Dormand and P. J. Prince, Journal of Computational and
    * Applied Mathematics, volume 6(1), 1980. DOI 10.1016/0771-0509(80)90013-3.
    * \param[in] A Matrix \f$ \mathbf{A} \f$.
    * \param[in] b Weights vector \f$ \mathbf{b} \f$ or embedded weights vector \f$ \hat{\mathbf{b}} \f$.
    * \param[in] c Nodes vector \f$ \mathbf{c} \f$.
    * \param[in] verbose Verbosity flag.
    * \return The calculated order of the Runge-Kutta method.
    */
    Integer compute_order(Matrix const &A, Vector const &b, Vector const &c, bool verbose = false) const {

      #define CMD "Sandals::" << this->name << "::tableau_order(...): "

      // Temporary variables initialization
      Integer order{0};
      Real tolerance{std::pow(EPSILON, 2.0/3.0)};

      // Precheck consistency
      Real tmp{(A*Vector::Ones() - c).norm()};
      if (tmp > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "precheck failed, ||A*I - c|| = " << tmp << " ≠ 0.");
        return order;
      }

      // Check order 1
      tmp = b.sum();
      if (std::abs(tmp - 1.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 1 check failed, a_1 = " << tmp << " ≠ 1.");
        return order;
      }

      order = 1; // Order 1 is the highest order that can be checked

      // Check order 2
      tmp = b.transpose() * c;
      if (std::abs(tmp - 1.0/2.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 2 check failed, a_1 = " << tmp << " ≠ 1/2.");
        return order;
      }

      order = 2; // Order 2 is the highest order that can be checked

      // Check order 3
      Vector c_2{c.array().pow(2)};
      tmp = b.transpose() * c_2;
      if (std::abs(tmp - 1.0/3.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 3 check failed, a_1 = " << tmp << " ≠ 1/3.");
        return order;
      }

      Vector Ac{A * c};
      tmp = b.transpose() * Ac;
      if (std::abs(tmp - 1.0/6.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 3 check failed, a_2 = " << tmp << " ≠ 1/6.");
        return order;
      }

      order = 3; // Order 3 is the highest order that can be checked

      // Check order 4
      Vector c_3{c.array().pow(3)};
      tmp = b.transpose() * c_3;
      if (std::abs(tmp - 1.0/4.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 4 check failed, a_1 = " << tmp << " ≠ 1/4.");
        return order;
      }

      Vector cAc{c.array() * Ac.array()};
      tmp = b.transpose() * cAc;
      if (std::abs(tmp - 1.0/8.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 4 check failed, a_2 = " << tmp << " ≠ 1/8.");
        return order;
      }

      Vector bA{(b.transpose() * A).transpose()};
      tmp = bA.transpose() * c_2;
      if (std::abs(tmp - 1.0/12.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 4 check failed, a_3 = " << tmp << " ≠ 1/12.");
        return order;
      }

      tmp = b.transpose() * A * A * c;
      if (std::abs(tmp - 1.0/24.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 4 check failed, a_4 = " << tmp << " ≠ 1/24.");
        return order;
      }

      order = 4; // Order 4 is the highest order that can be checked

      // Check order 5
      Vector c_4{c.array().pow(4)};
      tmp = b.transpose() * c_4;
      if (std::abs(tmp - 1.0/5.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_1 = " << tmp << " ≠ 1/5.");
        return order;
      }

      Vector b_c2{b.array() * c_2.array()};
      tmp = b_c2.transpose() * Ac;
      if (std::abs(tmp - 1.0/10.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_2 = " << tmp << " ≠ 1/10.");
        return order;
      }

      Vector b_Ac{b.array() * Ac.array()};
      tmp = b_Ac.transpose() * Ac;
      if (std::abs(tmp - 1.0/20.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_3 = " << tmp << " ≠ 1/20.");
        return order;
      }

      Vector b_c{b.array() * c.array()};
      Vector Ac2{A * c_2};
      tmp = b_c.transpose() * Ac2;
      if (std::abs(tmp - 1.0/15.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_4 = " << tmp << " ≠ 1/15.");
        return order;
      }

      Vector Ac3{A * c_3};
      tmp = b.transpose() * Ac3;
      if (std::abs(tmp - 1.0/20.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_5 = " << tmp << " ≠ 1/20.");
        return order;
      }

      Vector AAc{A * Ac};
      tmp = b_c.transpose() * AAc;
      if (std::abs(tmp - 1.0/30.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_6 = " << tmp << " ≠ 1/30.");
        return order;
      }

      tmp = bA.transpose() * (c.array() * Ac.array()).matrix();
      if (std::abs(tmp - 1.0/40.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_7 = " << tmp << " ≠ 1/40.");
        return order;
      }

      tmp = bA.transpose() * Ac2;
      if (std::abs(tmp - 1.0/60.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_8 = " << tmp << " ≠ 1/60.");
        return order;
      }

      tmp = bA.transpose() * AAc;
      if (std::abs(tmp - 1.0/120.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 5 check failed, a_9 = " << tmp << " ≠ 1/120.");
        return order;
      }

      order = 5; // Order 5 is the highest order that can be checked

      // Check order 6
      Vector c_5{c.array().pow(5)};
      tmp = b.transpose() * c_5;
      if (std::abs(tmp - 1.0/6.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_1 = " << tmp << " ≠ 1/6.");
        return order;
      }

      tmp = (b.array() * c_3.array()).matrix().transpose() * Ac;
      if (std::abs(tmp - 1.0/12.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_2 = " << tmp << " ≠ 1/12.");
        return order;
      }

      Vector Ac_2{Ac.array().pow(2)};
      tmp = b_c.transpose() * Ac_2;
      if (std::abs(tmp - 1.0/24.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_3 = " << tmp << " ≠ 1/24.");
        return order;
      }

      tmp = b_c2.transpose() * Ac2;
      if (std::abs(tmp - 1.0/18.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_4 = " << tmp << " ≠ 1/18.");
        return order;
      }

      tmp = b.transpose() * (Ac2.array() * Ac.array()).matrix();
      if (std::abs(tmp - 1.0/36.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_5 = " << tmp << " ≠ 1/36.");
        return order;
      }

      tmp = b_c.transpose() * Ac3;
      if (std::abs(tmp - 1.0/24.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_6 = " << tmp << " ≠ 1/24.");
        return order;
      }

      Vector Ac4{A * c_4};
      tmp = b.transpose() * Ac4;
      if (std::abs(tmp - 1.0/30.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_7 = " << tmp << " ≠ 1/30.");
        return order;
      }

      Vector bc2A{A.transpose() * b_c2};
      tmp = bc2A.transpose() * Ac;
      if (std::abs(tmp - 1.0/36.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_8 = " << tmp << " ≠ 1/36.");
        return order;
      }

      tmp = b_Ac.transpose() * AAc;
      if (std::abs(tmp - 1.0/72.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_9 = " << tmp << " ≠ 1/72.");
        return order;
      }

      Vector bcA{b_c.transpose() * A};
      tmp = bcA.transpose() * cAc;
      if (std::abs(tmp - 1.0/48.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_10 = " << tmp << " ≠ 1/48.");
        return order;
      }

      tmp = (bA.array() * c_2.array()).matrix().transpose() * Ac;
      if (std::abs(tmp - 1.0/60.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_11 = " << tmp << " ≠ 1/60.");
        return order;
      }

      tmp = bA.transpose() * Ac_2;
      if (std::abs(tmp - 1.0/120.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_12 = " << tmp << " ≠ 1/120.");
        return order;
      }

      tmp = bcA.transpose() * Ac2;
      if (std::abs(tmp - 1.0/72.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_13 = " << tmp << " ≠ 1/72.");
        return order;
      }

      Vector bA_c{bA.array() * c.array()};
      tmp = bA_c.transpose() * Ac2;
      if (std::abs(tmp - 1.0/90.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_14 = " << tmp << " ≠ 1/90.");
        return order;
      }

      tmp = bA.transpose() * Ac3;
      if (std::abs(tmp - 1.0/120.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_15 = " << tmp << " ≠ 1/120.");
        return order;
      }

      tmp = bcA.transpose() * AAc;
      if (std::abs(tmp - 1.0/144.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_16 = " << tmp << " ≠ 1/144.");
        return order;
      }

      tmp = bA_c.transpose() * AAc;
      if (std::abs(tmp - 1.0/180.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_17 = " << tmp << " ≠ 1/180.");
        return order;
      }

      Vector AcAc{A * cAc};
      tmp = bA.transpose() * AcAc;
      if (std::abs(tmp - 1.0/240.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_18 = " << tmp << " ≠ 1/240.");
        return order;
      }

      tmp = bA.transpose() * A * A * c_2;
      if (std::abs(tmp - 1.0/360.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_19 = " << tmp << " ≠ 1/360.");
        return order;
      }

      tmp = bA.transpose() * A * A * Ac;
      if (std::abs(tmp - 1.0/720.0) > tolerance) {
        SANDALS_ASSERT_WARNING(!verbose, CMD "order 6 check failed, a_20 = " << tmp << " ≠ 1/720.");
        return order;
      }

      order = 6; // Order 6 is the highest order that can be checked with this method

      return order;

      #undef CMD
    }
  }; // struct Tableau

} // namespace Sandals

#endif // SANDALS_TABLEAU_HXX
