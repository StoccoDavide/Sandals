/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *
 *                                                                                               *
 * The Sandals project is distributed under the BSD 2-Clause License.                            *
 *                                                                                               *
 * Davide Stocco                                                               Enrico Bertolazzi *
 * University of Trento                                                     University of Trento *
 * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TESTS_PROBLEMS_SHAMPINE2_EXPLICIT_HH
#define TESTS_PROBLEMS_SHAMPINE2_EXPLICIT_HH

#include "Sandals.hh"
#include "Sandals/Problem.hh"
#include "Sandals/System/Explicit.hh"

using namespace Sandals;

template<typename Real = double>
class Shampine2Explicit : public Explicit<Real, 2, 0>
{
public:
  using typename Explicit<Real, 2, 0>::VectorF;
  using typename Explicit<Real, 2, 0>::MatrixJF;
  using typename Explicit<Real, 2, 0>::VectorH;
  using typename Explicit<Real, 2, 0>::MatrixJH;

private:
  Real m_p{1.0e-1}; // Parameter

public:
  Shampine2Explicit() : Explicit<Real, 2, 0>("Shampine2Explicit") {}

  ~Shampine2Explicit() {}

  void p(Real const p) {this->m_p = p;}

  Real p() const {return this->m_p;}

  VectorF f(VectorF const & x, Real const t) const override
  {
    VectorF f;
    f << x(1), -3.0*this->m_p*x(0) / std::pow(this->m_p + t*t, 2.0);
    return f;
  }

  MatrixJF Jf_x(VectorF const & /*x*/, Real const t) const override {
    MatrixJF Jf_x(MatrixJF::Zero());
    Jf_x(0, 1) = 1.0;
    Jf_x(1, 0) = -3.0*this->m_p / ((this->m_p + t*t)*(this->m_p + t*t));
    return Jf_x;
  }

  VectorH h(VectorF const & /*x*/, Real const /*t*/) const override {return VectorH::Zero();}

  MatrixJH Jh_x(VectorF const & /*x*/, Real const /*t*/) const override {return MatrixJH::Zero();}

  bool in_domain(VectorF const & /*x*/, Real const /*t*/) const override {return true;}

};

template<typename Real, typename Integrator>
class Shampine2ExplicitProblem : public Problem<Real, 2, 0, Integrator>
{
public:
  using typename Problem<Real, 2, 0, Integrator>::SystemPtr;
  using typename Problem<Real, 2, 0, Integrator>::IntegratorPtr;
  using typename Problem<Real, 2, 0, Integrator>::SolutionPtr;
  using typename Problem<Real, 2, 0, Integrator>::VectorF;
  using typename Problem<Real, 2, 0, Integrator>::MatrixJF;
  using VectorX = Eigen::Vector<Real, Eigen::Dynamic>;
  using MatrixX = Eigen::Matrix<Real, 2, Eigen::Dynamic>;

  Shampine2ExplicitProblem(IntegratorPtr rk)
    : Problem<Real, 2, 0, Integrator>("Shampine2ExplicitProblem", std::make_shared<Shampine2Explicit<Real>>(), rk)
  {
    rk->system(this->system());
  }

  ~Shampine2ExplicitProblem() {}

  void p(Real const p) {static_cast<Shampine2Explicit<Real>*>(this->integrator()->system().get())->p(p);}

  Real p() const {return static_cast<const Shampine2Explicit<Real>*>(this->integrator()->system().get())->p();}

  VectorF b(VectorF const & x_ini, VectorF const & x_end) const override
  {
    VectorF b;
    b(0) = x_ini(0) + 0.1 / std::sqrt(this->p()+0.01);
    b(1) = x_end(0) - 0.1 / std::sqrt(this->p()+0.01);
    return b;
  }

  MatrixJF Jb_x_ini(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
  {
    MatrixJF Jb_x_ini(MatrixJF::Zero());
    Jb_x_ini(0, 0) = 1.0;
    return Jb_x_ini;
  }

  MatrixJF Jb_x_end(VectorF const & /*x_ini*/, VectorF const & /*x_end*/) const override
  {
    MatrixJF Jb_x_end(MatrixJF::Zero());
    Jb_x_end(1, 0) = 1.0;
    return Jb_x_end;
  }

  VectorF analytical_solution(Real const t) const {
    VectorF x;
    x <<
      this->p() / std::pow(this->p() + t*t, 1.5),
      t / std::sqrt(this->p() + t*t);
    return x;
  }

  MatrixX analytical_solution(VectorX const & t) const {
    MatrixX x(2, t.size());
    for (Integer i{0}; i < t.size(); ++i) {x.col(i) = this->analytical_solution(t(i));}
    return x;
  }
};

template<typename Integrator>
Shampine2ExplicitProblem(std::shared_ptr<Integrator>)
    -> Shampine2ExplicitProblem<typename Integrator::real_type, Integrator>;

#endif // TESTS_PROBLEMS_SHAMPINE2_EXPLICIT_HH
