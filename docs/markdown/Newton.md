# Newton's method

Newton's method, including its damped variant with an affine-invariant step, is implemented for solving problems of the form
\f[
  \mathbf{F}(\mathbf{x}) = \mathbf{0} \text{,} \quad \mathbf{F}: \mathbb{R}^{n} \rightarrow \mathbb{R}^{n} \text{.}
\f]
The method is based on the linearization of the function \f$\mathbf{F}\f$ around the current iterate \f$\mathbf{x}_k\f$, which leads to the linear system
\f[
  \mathbf{JF}_{\mathbf{x}}(\mathbf{x}_k) \mathbf{h}_k = -\mathbf{F}(\mathbf{x}_k) \text{.}
\f]
The advancing step \f$\mathbf{h}_k\f$ is then computed as
\f[
  \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{h}_k \text{,}
\f]
where \f$\alpha_k\f$ is a damping coefficient that ensures convergence by satisfying
\f[
  \left\|\mathbf{JF}_{\mathbf{x}}(\mathbf{x}_k)^{-1} \mathbf{F}(\mathbf{x}_{k+1})\right\|
  \leq \left(1 - \displaystyle\frac{\alpha_k}{2}\right) \left\|\mathbf{JF}_{\mathbf{x}}(\mathbf{x}_k)^{-1}
  \mathbf{F}(\mathbf{x}_k)\right\| = \left(1 - \displaystyle\frac{\alpha_k}{2} \right)
  \left\|\mathbf{h}_k\right\| \text{.}
\f]
For detailed information on the affine-invariant Newton's method, refer to [this link](https://www.zib.de/deuflhard/research/algorithm/ainewton.en.html).
