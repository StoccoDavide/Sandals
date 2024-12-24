# Nonlinear solvers

This section describes the nonlinear solvers implemented in `Sandals`. The available solvers are *Newton*'s method and *Broyden*'s method, both designed to solve nonlinear systems of equations. Newton's method is a well-known classical approach, while Broyden's method is a quasi-Newton technique that approximates the Jacobian matrix for improved efficiency in certain scenarios.

The nonlinear solvers are implemented for solving problems of the form
\f[
  \mathbf{F}(\mathbf{x}) = \mathbf{0} \text{,} \quad \mathbf{F}: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n} \text{.}
\f]
and finding the root of the function \f$\mathbf{F}\f$ by iteratively updating the current iterate \f$\mathbf{x}_k\f$ until convergence is achieved. The solvers require the function \f$\mathbf{F}\f$ and its Jacobian matrix \f$\mathbf{JF}_{\mathbf{x}}\f$ to be provided by the user. The Jacobian matrix can be computed analytically or numerically, depending on the problem's complexity and the user's preference. Nontheless, quasi-Newton methods, such as Broyden's method, are designed to approximate the Jacobian matrix iteratively, avoiding the need to compute it explicitly.
