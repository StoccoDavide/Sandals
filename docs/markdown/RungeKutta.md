# Runge-Kutta methods and Butcher tableaus

In this section, we introduce Runge-Kutta methods for solving ordinary differential equations (ODEs) and differential-algebraic equations (DAEs). Such equation systems describe the dynamics of systems and are used to model a wide range of physical phenomena.

## Dynamic systems representation

ODEs and DAEs can be expressed in several forms. While these forms may be mathematically equivalent and convertible into one another, their numerical properties differ, making some more practical for specific computations. Below are the most common representations:

- *Explicit* form: \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$.
- *Implicit* form: \f$\mathbf{F}(t\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$.
- *Linear* form: \f$\mathbf{E}(t)\mathbf{x}^\prime = \mathbf{A}(t)\mathbf{x} + \mathbf{f}(t)\f$.
- *SemiExplicit* form: \f$\mathbf{M}(\mathbf{x}, t)\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$.

The implicit form is the most general and is frequently used to represent DAEs. However, the explicit and semi-explicit forms are more commonly employed for ODEs. Note that the semi-explicit and linear forms are subsets of the explicit form. Consequently, the Runge-Kutta methods we implement will be specialized for the explicit and implicit forms.

## Explicit Runge-Kutta methods

Explicit Runge-Kutta (ERK) methods are among the simplest numerical techniques for solving ordinary differential equations. These methods are termed "explicit" because they typically avoid solving nonlinear systems at each step. Instead, the solution is computed directly using function evaluations. However, this holds only if the derivatives of the dynamic system's states, \f$\mathbf{x}^\prime\f$, can be computed directly or via matrix inversion. This is possible when the dynamic system can be expressed in an explicit, semi-explicit, or linear form.

### ERK methods for explicit dynamic systems

For an explicit system of the form \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$, the ERK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{K}_i = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
where \f$s\f$ is the number of stages, \f$h_k\f$ is the step size, \f$\mathbf{x}_k\f$ is the state at time \f$t_k\f$, and \f$\mathbf{K}_i\f$ are the intermediate variables.

For explicit methods, the Runge-Kutta coefficient matrix \f$\mathbf{A}\f$ is strictly lower triangular. Thus, the stages can be computed sequentially as
\f[
  \mathbf{K}_i = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{i-1} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \text{,}
\f]
and by unrolling the stages, we obtain
\f[
  \begin{array}{l}
    \mathbf{K}_1 = \mathbf{f} \left(\mathbf{x}_k, t_k\right) \\
    \mathbf{K}_2 = \mathbf{f} \left(\mathbf{x}_k + h_k a_{21} \mathbf{K}_1, t_k + h_k c_2\right) \\
    \vdots \\[-0.5em]
    \mathbf{K}_s = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
  \end{array} \text{.}
\f]

### ERK methods for implicit dynamic systems

For an implicit dynamical system of the form  \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$, the ERK method becomes
\f[
  \begin{array}{l}
    \mathbf{F}_i \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{i-1} a_{ij} \mathbf{K}_j, \mathbf{K}_i, t_k + h_k c_i\right) = \mathbf{0} \\
      \mathbf{x}_{k+1} = \mathbf{x}_k + \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
Here, the \f$s\f$-stage system of equations \f$\mathbf{F}_i\f$ forms a lower triangular system, typically nonlinear in the intermediate variables \f$\mathbf{K}_i\f$. The solution is obtained using forward substitution, where each stage is solved sequentially, *i.e.,*
\f[
  \begin{array}{l}
    \mathbf{F}_1 \left(\mathbf{x}_k, \mathbf{K}_1, t_k + h_k c_1\right) = \mathbf{0} \\
    \mathbf{F}_2 \left(\mathbf{x}_k + h_k a_{21} \mathbf{K}_1, \mathbf{K}_2, t_k + h_k c_2\right) = \mathbf{0} \\
    \vdots \\[-0.5em]
    \mathbf{F}_s \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, \mathbf{K}_s, t_k + h_k c_s\right) = \mathbf{0}
  \end{array} \text{.}
\f]
Notice that the intermediate variables \f$\mathbf{K}_i\f$ are typically computed using Newton's method.

## Implicit Runge-Kutta methods

Implicit Runge-Kutta (IRK) methods are more general than ERK methods, as they can handle implicit systems of the form \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$. IRK methods are more stable but computationally expensive because they require solving a nonlinear system at each step. A key characteristic of IRK methods is that the Runge-Kutta coefficient matrix \f$\mathbf{A}\f$ is fully populated, meaning that the stages are coupled and must be solved simultaneously.

### IRK methods for explicit dynamic systems

For an explicit dynamical system of the form \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$, the IRK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{K}_i = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{,}
\f]
Unrolling the stages gives
\f[
  \begin{array}{l}
    \mathbf{K}_1 = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, t_k + h_k c_1\right) \\
    \mathbf{K}_2 = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, t_k + h_k c_2\right) \\[-0.5em]
    \vdots \\[-0.5em]
    \mathbf{K}_s = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
  \end{array} \text{.}
\f]

### IRK methods for implicit dynamic systems

For an implicit dynamic system of the form \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$, the IRK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{F}_i \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, \mathbf{K}_i, t_k + h_k c_i\right) = \mathbf{0} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
Unrolling the stages, we have
\f[
  \begin{array}{l}
    \mathbf{F}_1 \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, \mathbf{K}_1, t_k + h_k c_1\right) = \mathbf{0} \\
    \mathbf{F}_2 \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, \mathbf{K}_2, t_k + h_k c_2\right) = \mathbf{0} \\[-0.5em]
    \vdots \\[-0.5em]
    \mathbf{F}_s \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \mathbf{K}_s, t_k + h_k c_s\right) = \mathbf{0}
  \end{array} \text{.}
\f]

> **NOTE:** In `Sandals`, the implementation of IRK methods is unified for both explicit and implicit dynamical systems. IRK methods inherently require solving nonlinear systems, irrespective of the system's form. Creating a specialized version for explicit systems would not yield significant computational advantages, as the complexity of solving the nonlinear system remains. Consequently, only the IRK method for implicit dynamical systems is implemented.

## Diagonally implicit Runge-Kutta methods

Diagonally implicit Runge-Kutta (DIRK) methods are a specialized subset of IRK methods. They are characterized by a lower triangular \f$\mathbf{A}\f$-matrix with at least one nonzero diagonal entry. This structure allows each stage to be solved sequentially rather than simultaneously, making DIRK methods more computationally efficient compared to general IRK methods.

### DIRK methods for explicit dynamic systems

For an explicit dynamic system of the form \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$, the DIRK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{K}_i = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{i} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
When the stages are unrolled, the equations become
\f[
  \begin{array}{l}
    \mathbf{K}_1 = \mathbf{f} \left(\mathbf{x}_k + h_k a_{11} \mathbf{K}_1, t_k + h_k c_1\right) \\
    \mathbf{K}_2 = \mathbf{f} \left(\mathbf{x}_k + h_k a_{21} \mathbf{K}_1 + h_k a_{22} \mathbf{K}_2, t_k + h_k c_2\right) \\
    \vdots \\[-0.5em]
    \mathbf{K}_s = \mathbf{f} \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
  \end{array} \text{.}
\f]

### DIRK methods for implicit dynamic systems

For an implicit dynamic system of the form \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$, the DIRK method is written as
\f[
  \begin{array}{l}
    \mathbf{F}_i \left(
      \mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{i} a_{ij} \mathbf{K}_j, \mathbf{K}_i, t_k + h_k c_i
    \right) = \mathbf{0} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
If the stages are unrolled, we have
\f[
  \begin{array}{l}
    \mathbf{F}_1 \left(\mathbf{x}_k + h_k a_{11} \mathbf{K}_1, \mathbf{K}_1, t_k + h_k c_1\right) = \mathbf{0} \\
    \mathbf{F}_2 \left(\mathbf{x}_k + h_k a_{21} \mathbf{K}_1 + h_k a_{22} \mathbf{K}_2, \mathbf{K}_2, t_k + h_k c_2\right) = \mathbf{0} \\
    \vdots \\[-0.5em]
    \mathbf{F}_s \left(\mathbf{x}_k + h_k \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \mathbf{K}_s, t_k + h_k c_s\right) = \mathbf{0}
  \end{array} \text{.}
\f]

> **NOTE:** As with IRK methods, the implementation of DIRK methods for explicit and implicit dynamic systems in `Sandals` is unified. Since solving the nonlinear system is inherent to DIRK methods, a specialized version for explicit systems would not significantly reduce computational complexity. Consequently, only the DIRK method for implicit dynamic systems is implemented.
