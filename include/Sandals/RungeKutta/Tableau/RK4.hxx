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

#ifndef SANDALS_RK4_HXX
#define SANDALS_RK4_HXX

namespace Sandals
{

  //! Butcher tableau for the Runge-Kutta 4 method:
  //!
  //! \f[
  //! \begin{array}{c|cccc}
  //!     0 &   0 &   0 &   0 &   0 \\
  //!   1/2 & 1/2 &   0 &   0 &   0 \\
  //!   1/2 &   0 & 1/2 &   0 &   0 \\
  //!     1 &   0 &   0 &   1 &   0 \\
  //!   \hline
  //!       & 1/6 & 1/3 & 1/3 & 1/6
  //! \end{array} \text{.}
  //! \f]
  class RK4Tableau : public Tableau<4>
  {
    public:
    using Tableau<4>::Type;   //!< Type enumeration.
    using Tableau<4>::Vector; //!< Vector type.
    using Tableau<4>::Matrix; //!< Matrix type.

    //! Class constructor for the RK4 method.
    RK4Tableau() {
      this->name        = "RK4";
      this->rk_type     = Type::ERK;
      this->order       = 4;
      this->is_embedded = false;
      this->A << 0.0,     0.0,     0.0,     0.0,
                 1.0/2.0, 0.0,     0.0,     0.0,
                 0.0,     1.0/2.0, 0.0,     0.0,
                 0.0,     0.0,     1.0,     0.0;
      this->b << 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0;
      this->c << 0.0,     1.0/2.0, 1.0/2.0, 1.0;
      this->b_e << NAN_VEC4;
    }
  }; // class RK4Tableau

  //! Class container for the Runge-Kutta 4 method.
  template <Size N>
  class RK4 : public RungeKutta<4, N>
  {
  public:
    using System = typename Implicit<N>::Pointer; //!< System type.

    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    //! \param t_tableau The Tableau reference.
    RK4() : RungeKutta<4, N>(RK4Tableau()) {}

    //! Class constructor for a Runge-Kutta solver given a Tableau reference.
    //! \param t_tableau The Tableau reference.
    //! \param t_system The system reference.
    RK4(System t_system) : RungeKutta<4, N>(RK4Tableau(), t_system) {}

  }; // class RK4

} // namespace Sandals

#endif // SANDALS_RK4_HXX


    //! Check Butcher tableau consistency for an explicit Runge-Kutta method.
    //! \param tbl.A   Matrix \f$ \mathbf{A} \f$.
    //! \param tbl.b   Weights vector \f$ \mathbf{b} \f$.
    //! \param tbl.b_e [optional] Embedded weights vector \f$ \mathbf{b}_{e} \f$.
    //! \param tbl.c   Nodes vector \f$ \mathbf{c} \f$.
    //! \return True if the Butcher tableau is consistent, false otherwise.
    //[out,order,e_order] = check_tableau( this, tbl )

    //! Check the order of a Runge-Kutta tableau according to the conditions taken from:
    //! *A family of embedded Runge-Kutta formulae*, J. R. Dormand and P. J. Prince,
    //! Journal of Computational and Applied Mathematics, volume 6(1), 1980.
    //! Doi: [10.1016/0771-0509(80)90013-3](https://doi.org/10.1016/0771-0509(80)90013-3)
    //! \param tableau The Runge-Kutta tableau to be checked.
    //! \return The order of the Runge-Kutta tableau.
    /*Size tableau_order(matS const &A, vecS const &b, vecS const &c) const
    {
      #define CMD "Indigo.RungeKutta.tableau_order(...): "

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
