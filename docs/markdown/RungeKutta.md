# Runge-Kutta methods and Butcher tableaus

In this section, we introduce Runge-Kutta methods for solving ordinary differential equations (ODEs) and differential-algebraic equations (DAEs). Such equation systems describe the dynamics of systems and are used to model a wide range of physical phenomena.

## Dynamic systems representation

ODEs and DAEs can be expressed in several forms. While these forms may be mathematically equivalent and convertible into one another, their numerical properties differ, making some more practical for specific computations. Below are the most common representations:

- *Explicit* form: \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$.
- *Implicit* form: \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$.
- *SemiExplicit* form: \f$\mathbf{A}(\mathbf{x}, t)\mathbf{x}^\prime = \mathbf{b}(\mathbf{x}, t)\f$.
- *Linear* form: \f$\mathbf{E}(t)\mathbf{x}^\prime = \mathbf{A}(t)\mathbf{x} + \mathbf{b}(t)\f$.

The implicit form is the most general and is frequently used to represent DAEs. However, the explicit forms is more commonly employed for ODEs. However, the explicit, semi-explicit, and linear forms are more commonly employed for ODEs. Note that the semi-explicit and the linear forms are a subset of the explicit form. Consequently, the Runge-Kutta methods we implement will be specialized only for the explicit and implicit forms.

## Explicit Runge-Kutta methods

Explicit Runge-Kutta (ERK) methods are among the simplest numerical techniques for solving ordinary differential equations. These methods are termed "explicit" because they typically avoid solving nonlinear systems at each step. Instead, the solution is computed directly using function evaluations. However, this holds only if the derivatives of the dynamic system's states, \f$\mathbf{x}^\prime\f$, can be computed directly or via matrix inversion. This is possible when the dynamic system can be expressed in an explicit, semi-explicit, or linear form.

### ERK methods for explicit dynamic systems

For an explicit system of the form \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$, the ERK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{K}_i = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
where \f$s\f$ is the number of stages, \f$h_k\f$ is the step size, \f$\mathbf{x}_k\f$ is the state at time \f$t_k\f$, and \f$\mathbf{K}_i\f$ are the intermediate variables.

For explicit methods, the Runge-Kutta coefficient matrix \f$\mathbf{A}\f$ is strictly lower triangular. Thus, the stages can be computed sequentially as
\f[
  \mathbf{K}_i = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{i-1} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \text{,}
\f]
and by unrolling the stages, we obtain
\f[
  \begin{array}{l}
    \begin{cases}
      \mathbf{K}_1 = h_k \mathbf{f} \left(\mathbf{x}_k, t_k\right) \\
      \mathbf{K}_2 = h_k \mathbf{f} \left(\mathbf{x}_k + a_{21} \mathbf{K}_1, t_k + h_k c_2\right) \\
      \mathbf{K}_3 = h_k \mathbf{f} \left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, t_k + h_k c_3\right) \\[-0.5em]
      \vdots \\[-0.5em]
      \mathbf{K}_s = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
    \end{cases} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    \mathbf{I} & \mathbf{0} & \cdots & \cdots & \mathbf{0} \\
    a_{21} & \mathbf{I} & \ddots & \mathbf{0} & \vdots \\
    a_{31} & a_{32} & \mathbf{I} & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & \mathbf{0} \\
    a_{s1} & a_{s2} & \cdots & a_{s,s-1} & \mathbf{I}
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{K}_1 \\
    \mathbf{K}_2 \\
    \mathbf{K}_3 \\
    \vdots \\
    \mathbf{K}_s
  \end{bmatrix} =
  h_k
  \begin{bmatrix}
    \mathbf{f} \left(\mathbf{x}_k, t_k\right) \\
    \mathbf{f} \left(\mathbf{x}_k + a_{21} \mathbf{K}_1, t_k + h_k c_2\right) \\
    \mathbf{f} \left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, t_k + h_k c_3\right) \\
    \vdots \\
    \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{K}\f$ with respect to \f$\mathbf{x}_k\f$

For explicit systems, the derivatives of the intermediate variables \f$\mathbf{K}\f$ with respect to the states \f$\mathbf{x}\f$ can be computed by the differentiation of the stage equations, which yields
\f[
  \begin{array}{l}
    \begin{cases}
      \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}(\mathbf{x}_k, t_k) \\
      \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1, t_k + h_k c_2\right) \left(\mathbf{I} + a_{21} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k}\right) \\
      \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, t_k + h_k c_3\right) \left(\mathbf{I} + a_{31} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} + a_{32} \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k}\right) \\[-0.5em]
      \vdots \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s-1} a_{sj} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right)
    \end{cases}
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    \mathbf{I} & \mathbf{0} & \cdots & \cdots & \mathbf{0} \\
    -h_k a_{21} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} & \ddots & \mathbf{0} & \vdots \\
    -h_k a_{31} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{32} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & \mathbf{0} \\
    -h_k a_{s1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{s2} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & -h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I}
  \end{bmatrix}
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k}
  \end{bmatrix} =
  h_k
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}(\mathbf{x}_k, t_k) \\
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1, t_k + h_k c_2\right) \\
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, t_k + h_k c_3\right) \\
    \vdots\\
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{x}_{k+1}\f$ with respect to \f$\mathbf{x}_k\f$

The derivatives of the next state \f$\mathbf{x}_{k+1}\f$ with respect to the states \f$\mathbf{x}\f$ can be computed as
\f[
  \displaystyle\frac{\partial\mathbf{x}_{k+1}}{\partial\mathbf{x}_k} =
  \mathbf{I} + h_k \displaystyle\sum_{i=1}^s b_i \frac{\partial\mathbf{K}_i}{\partial\mathbf{x}_k} \text{.}
\f]
where the derivatives of the intermediate variables \f$\mathbf{K}_i\f$ with respect to the states \f$\mathbf{x}_k\f$ are computed as described above.

### ERK methods for implicit dynamic systems

For an implicit dynamical system of the form  \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$, the ERK method becomes
\f[
  \begin{array}{l}
    \mathbf{F}_i \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{i-1} a_{ij} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_i}{h_k}, t_k + h_k c_i\right) = \mathbf{0} \\
      \mathbf{x}_{k+1} = \mathbf{x}_k + \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
Here, the \f$s\f$-stage system of equations \f$\mathbf{F}_i\f$ forms a lower triangular system, typically nonlinear in the intermediate variables \f$\mathbf{K}_i\f$. The solution is obtained using forward substitution, where each stage is solved sequentially, *i.e.,*
\f[
  \begin{array}{l}
    \begin{cases}
      \mathbf{F}_1 \left(\mathbf{x}_k, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) = \mathbf{0} \\
      \mathbf{F}_2 \left(\mathbf{x}_k + a_{21} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) = \mathbf{0} \\
      \mathbf{F}_3 \left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) = \mathbf{0} \\[-0.5em]
      \vdots \\[-0.5em]
      \mathbf{F}_s \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) = \mathbf{0}
    \end{cases} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
Notice that the intermediate variables \f$\mathbf{K}_i\f$ are typically computed using Newton's method.

#### Derivatives of \f$\mathbf{K}\f$ with respect to \f$\mathbf{x}_k\f$

For implicit systems, the derivatives of the intermediate variables \f$\mathbf{K}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed by differentiating the stage equations, which yields
\f[
  \begin{array}{l}
    \begin{cases}
      \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\mathbf{x}_k, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{K}_1}\left(\mathbf{x}_k, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} = \mathbf{0} \\
      \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \left(\mathbf{I} + a_{21} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{K}_2}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} = \mathbf{0} \\
      \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \left(\mathbf{I} + a_{31} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} + a_{32} \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{K}_3}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} = \mathbf{0} \\[-0.5em]
      \vdots \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s-1} a_{sj} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{K}_s}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k} = \mathbf{0}
    \end{cases}
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{K}_1}\left(\cdot\right) & \mathbf{0} & \cdots & \cdots & \mathbf{0} \\
    h_k a_{21} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{K}_2}\left(\cdot\right) & \ddots & \mathbf{0} & \vdots \\
    h_k a_{31} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{32} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{K}_3}\left(\cdot\right) & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & \mathbf{0} \\
    h_k a_{s1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{s2} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{K}_s}\left(\cdot\right)
  \end{bmatrix}
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k}
  \end{bmatrix} =
  -h_k
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\mathbf{x}_k, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \\
    \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \\
    \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \\
    \vdots\\
    \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s-1} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{x}_{k+1}\f$ with respect to \f$\mathbf{x}_k\f$

> **NOTE:** The derivatives of the next state \f$\mathbf{x}_{k+1}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed as previously described for explicit systems, as the structure of the ERK method remains the same.

## Implicit Runge-Kutta methods

Implicit Runge-Kutta (IRK) methods are more general than ERK methods, as they can handle implicit systems of the form \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$. IRK methods are more stable but computationally expensive because they require solving a nonlinear system at each step. A key characteristic of IRK methods is that the Runge-Kutta coefficient matrix \f$\mathbf{A}\f$ is fully populated, meaning that the stages are coupled and must be solved simultaneously.

### IRK methods for explicit dynamic systems

For an explicit dynamical system of the form \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$, the IRK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{K}_i = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{,}
\f]
Unrolling the stages gives
\f[
  \begin{array}{l}
    \begin{cases}
      \mathbf{K}_1 = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, t_k + h_k c_1\right) \\
      \mathbf{K}_2 = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, t_k + h_k c_2\right) \\
      \mathbf{K}_3 = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, t_k + h_k c_3\right) \\[-0.5em]
      \vdots \\[-0.5em]
      \mathbf{K}_s = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
    \end{cases} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]

#### Derivatives of \f$\mathbf{K}\f$ with respect to \f$\mathbf{x}_k\f$

For explicit systems, the derivatives of the intermediate variables \f$\mathbf{K}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed by differentiating the stage equations, which yields
\f[
  \begin{array}{l}
    \begin{cases}
      \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, t_k + h_k c_1\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{1j} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) \\
      \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, t_k + h_k c_2\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{2j} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) \\
      \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, t_k + h_k c_3\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{3j} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) \\[-0.5em]
      \vdots \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{sj} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right)
    \end{cases}
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    \mathbf{I} - h_k a_{11} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{12} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{13} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & -h_k a_{1s} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    -h_k a_{21} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} - h_k a_{22} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{23} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & \vdots \\
    -h_k a_{31} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{32} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} - h_k a_{33} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & -h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    -h_k a_{s1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{s2} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & -h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} - h_k a_{ss} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right)
  \end{bmatrix}
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k}
  \end{bmatrix} =
  h_k
  \begin{bmatrix}
  \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, t_k + h_k c_1\right) \\
  \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, t_k + h_k c_2\right) \\
  \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, t_k + h_k c_3\right) \\
  \vdots \\
  \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{x}_{k+1}\f$ with respect to \f$\mathbf{x}_k\f$

The derivatives of the next state \f$\mathbf{x}_{k+1}\f$ with respect to the states \f$\mathbf{x}\f$ can be computed as
\f[
  \displaystyle\frac{\partial\mathbf{x}_{k+1}}{\partial\mathbf{x}_k} =
  \mathbf{I} + h_k \displaystyle\sum_{i=1}^s b_i \frac{\partial\mathbf{K}_i}{\partial\mathbf{x}_k} \text{.}
\f]
where the derivatives of the intermediate variables \f$\mathbf{K}_i\f$ with respect to the states \f$\mathbf{x}_k\f$ are computed as described above.

### IRK methods for implicit dynamic systems

For an implicit dynamic system of the form \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$, the IRK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{F}_i \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{ij} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_i}{h_k}, t_k + h_k c_i\right) = \mathbf{0} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
Unrolling the stages, we have
\f[
  \begin{array}{l}
    \begin{cases}
      \mathbf{F}_1 \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) = \mathbf{0} \\
      \mathbf{F}_2 \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) = \mathbf{0} \\
      \mathbf{F}_3 \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) = \mathbf{0} \\[-0.5em]
      \vdots \\[-0.5em]
      \mathbf{F}_s \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) = \mathbf{0}
    \end{cases} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]

> **NOTE:** In `Sandals`, the implementation of IRK methods is unified for both explicit and implicit dynamical systems. IRK methods inherently require solving nonlinear systems, irrespective of the system's form. Creating a specialized version for explicit systems would not yield significant computational advantages, as the complexity of solving the nonlinear system remains. Consequently, only the IRK method for implicit dynamical systems is implemented.

#### Derivatives of \f$\mathbf{K}\f$ with respect to \f$\mathbf{x}_k\f$

For implicit systems, the derivatives of the intermediate variables \f$\mathbf{K}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed by differentiating the stage equations, which yields
\f[
  \begin{array}{l}
    \begin{cases}
      \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{1j} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{K}_1}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} = \mathbf{0} \\
      \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{2j} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{K}_2}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} = \mathbf{0} \\
      \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{3j} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{K}_3}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} = \mathbf{0} \\[-0.5em]
      \vdots \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{sj} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{K}_s}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k} = \mathbf{0}
    \end{cases}
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{K}_1}\left(\cdot\right) & h_k a_{12} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{13} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & h_k a_{1s} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    h_k a_{21} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{K}_2}\left(\cdot\right) & h_k a_{23} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & \vdots \\
    h_k a_{31} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{32} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{K}_3}\left(\cdot\right) & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    h_k a_{s1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{s2} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{K}_s}\left(\cdot\right)
  \end{bmatrix}
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k}
  \end{bmatrix} =
  -h_k
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{1j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \\
    \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{2j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \\
    \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{3j} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{x}_{k+1}\f$ with respect to \f$\mathbf{x}_k\f$

> **NOTE:** The derivatives of the next state \f$\mathbf{x}_{k+1}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed as previously described for explicit systems, as the structure of the IRK method remains the same.

## Diagonally implicit Runge-Kutta methods

Diagonally implicit Runge-Kutta (DIRK) methods are a specialized subset of IRK methods. They are characterized by a lower triangular \f$\mathbf{A}\f$-matrix with at least one nonzero diagonal entry. This structure allows each stage to be solved sequentially rather than simultaneously, making DIRK methods more computationally efficient compared to general IRK methods.

### DIRK methods for explicit dynamic systems

For an explicit dynamic system of the form \f$\mathbf{x}^\prime = \mathbf{f}(\mathbf{x}, t)\f$, the DIRK method is expressed as
\f[
  \begin{array}{l}
    \mathbf{K}_i = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{i} a_{ij} \mathbf{K}_j, t_k + h_k c_i\right) \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
When the stages are unrolled, the equations become
\f[
  \begin{array}{l}
    \begin{cases}
      \mathbf{K}_1 = h_k \mathbf{f} \left(\mathbf{x}_k + a_{11} \mathbf{K}_1, t_k + h_k c_1\right) \\
      \mathbf{K}_2 = h_k \mathbf{f} \left(\mathbf{x}_k + a_{21} \mathbf{K}_1 + h_k a_{22} \mathbf{K}_2, t_k + h_k c_2\right) \\
      \vdots \\[-0.5em]
      \mathbf{K}_s = h_k \mathbf{f} \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right)
    \end{cases} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]

#### Derivatives of \f$\mathbf{K}\f$ with respect to \f$\mathbf{x}_k\f$

For explicit systems, the derivatives of the intermediate variables \f$\mathbf{K}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed by differentiating the stage equations, which yields
\f[
  \begin{array}{l}
    \begin{cases}
      \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{11} \mathbf{K}_1, t_k + h_k c_1\right) \left(\mathbf{I} + a_{11} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k}\right) \\
      \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1 + a_{22} \mathbf{K}_2, t_k + h_k c_2\right) \left(\mathbf{I} + a_{21} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} + a_{22} \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k}\right) \\
      \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2 + a_{33} \mathbf{K}_3, t_k + h_k c_3\right) \left(\mathbf{I} + a_{31} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} + a_{32} \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} + a_{33} \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k}\right) \\[-0.5em]
      \vdots \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k} = h_k \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, t_k + h_k c_s\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{sj} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right)
    \end{cases}
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    \mathbf{I} - h_k a_{11} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \\
    -h_k a_{21} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} - h_k a_{22} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \ddots & \mathbf{0} & \vdots \\
    -h_k a_{31} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{32} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} - h_k a_{33} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & \mathbf{0} \\
    -h_k a_{s1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{s2} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & -h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{I} - h_k a_{ss} \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right)
  \end{bmatrix}
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k}
  \end{bmatrix} =
  h_k
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right) \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{f}}{\partial\mathbf{x}_k}\left(\cdot\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{x}_{k+1}\f$ with respect to \f$\mathbf{x}_k\f$

The derivatives of the next state \f$\mathbf{x}_{k+1}\f$ with respect to the states \f$\mathbf{x}\f$ can be computed as
\f[
  \displaystyle\frac{\partial\mathbf{x}_{k+1}}{\partial\mathbf{x}_k} =
  \mathbf{I} + h_k \displaystyle\sum_{i=1}^s b_i \frac{\partial\mathbf{K}_i}{\partial\mathbf{x}_k} \text{.}
\f]
where the derivatives of the intermediate variables \f$\mathbf{K}_i\f$ with respect to the states \f$\mathbf{x}_k\f$ are computed as described above.

### DIRK methods for implicit dynamic systems

For an implicit dynamic system of the form \f$\mathbf{F}(\mathbf{x}, \mathbf{x}^\prime, t) = \mathbf{0}\f$, the DIRK method is written as
\f[
  \begin{array}{l}
    \mathbf{F}_i \left(
      \mathbf{x}_k + \displaystyle\sum_{j=1}^{i} a_{ij} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_i}{h_k}, t_k + h_k c_i
    \right) = \mathbf{0} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]
If the stages are unrolled, we have
\f[
  \begin{array}{l}
    \begin{cases}
      \mathbf{F}_1 \left(\mathbf{x}_k + a_{11} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) = \mathbf{0} \\
      \mathbf{F}_2 \left(\mathbf{x}_k + a_{21} \mathbf{K}_1 + a_{22} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_2\right) = \mathbf{0} \\
      \mathbf{F}_3 \left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2 + a_{33} \mathbf{K}_3, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) = \mathbf{0} \\[-0.5em]
      \vdots \\[-0.5em]
      \mathbf{F}_s \left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) = \mathbf{0}
    \end{cases} \\
    \mathbf{x}_{k+1} = \mathbf{x}_k + h_k \displaystyle\sum_{i=1}^s b_i \mathbf{K}_i
  \end{array} \text{.}
\f]

> **NOTE:** As with IRK methods, the implementation of DIRK methods for explicit and implicit dynamic systems in `Sandals` is unified. Since solving the nonlinear system is inherent to DIRK methods, a specialized version for explicit systems would not significantly reduce computational complexity. Consequently, only the DIRK method for implicit dynamic systems is implemented.

#### Derivatives of \f$\mathbf{K}\f$ with respect to \f$\mathbf{x}_k\f$

For implicit systems, the derivatives of the intermediate variables \f$\mathbf{K}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed by differentiating the stage equations, which yields
\f[
  \begin{array}{l}
    \begin{cases}
      \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{11} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \left(\mathbf{I} + a_{11} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{K}_1}\left(\mathbf{x}_k + a_{11} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} = \mathbf{0} \\
      \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1 + a_{22} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \left(\mathbf{I} + a_{21} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} + a_{22} \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{K}_2}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1 + a_{22} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} = \mathbf{0} \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2 + a_{33} \mathbf{K}_3, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \left(\mathbf{I} + a_{31} \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} + a_{32} \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} + a_{33} \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{K}_3}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2 + a_{33} \mathbf{K}_3, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} = \mathbf{0} \\[-0.5em]
      \vdots \\[-0.5em]
      \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) \left(\mathbf{I} + \displaystyle\sum_{j=1}^{s} a_{sj} \displaystyle\frac{\partial\mathbf{K}_j}{\partial\mathbf{x}_k}\right) + \frac{1}{h_k} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{K}_s}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right) \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k} = \mathbf{0}
    \end{cases}
  \end{array} \text{.}
\f]
In matrix form, this can be expressed as
\f[
  \begin{bmatrix}
    h_k a_{11} \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\cdot\right) & \mathbf{0} & \mathbf{0} & \cdots & \mathbf{0} \\
    -h_k a_{21} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{22} \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\cdot\right) & \ddots & \mathbf{0} & \vdots \\
    -h_k a_{31} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{32} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & h_k a_{33} \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\cdot\right) & \ddots & \vdots \\
    \vdots & \vdots & \ddots & \ddots & \mathbf{0} \\
    -h_k a_{s1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & -h_k a_{s2} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & \cdots & -h_k a_{s,s-1} \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\cdot\right) & \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{K}_s}\left(\cdot\right)
  \end{bmatrix}
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{K}_1}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_2}{\partial\mathbf{x}_k} \\
    \displaystyle\frac{\partial\mathbf{K}_3}{\partial\mathbf{x}_k} \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{K}_s}{\partial\mathbf{x}_k}
  \end{bmatrix} =
  -h_k
  \begin{bmatrix}
    \displaystyle\frac{\partial\mathbf{F}_1}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{11} \mathbf{K}_1, \displaystyle\frac{\mathbf{K}_1}{h_k}, t_k + h_k c_1\right) \\
    \displaystyle\frac{\partial\mathbf{F}_2}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{21} \mathbf{K}_1 + a_{22} \mathbf{K}_2, \displaystyle\frac{\mathbf{K}_2}{h_k}, t_k + h_k c_2\right) \\
    \displaystyle\frac{\partial\mathbf{F}_3}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + a_{31} \mathbf{K}_1 + a_{32} \mathbf{K}_2 + a_{33} \mathbf{K}_3, \displaystyle\frac{\mathbf{K}_3}{h_k}, t_k + h_k c_3\right) \\
    \vdots \\
    \displaystyle\frac{\partial\mathbf{F}_s}{\partial\mathbf{x}_k}\left(\mathbf{x}_k + \displaystyle\sum_{j=1}^{s} a_{sj} \mathbf{K}_j, \displaystyle\frac{\mathbf{K}_s}{h_k}, t_k + h_k c_s\right)
  \end{bmatrix} \text{.}
\f]

#### Derivatives of \f$\mathbf{x}_{k+1}\f$ with respect to \f$\mathbf{x}_k\f$

> **NOTE:** The derivatives of the next state \f$\mathbf{x}_{k+1}\f$ with respect to the states \f$\mathbf{x}_k\f$ can be computed as previously described for explicit systems, as the structure of the IRK method remains the same.

# Projection on the invariants manifold

In many applications, the dynamic system's states must satisfy certain constraints, which are typically expressed as invariants. For example, in mechanical systems, the constraints may represent the conservation of energy or momentum. In such cases, the states must be projected onto the manifold defined by the invariants to ensure that the system's dynamics are physically consistent. In `Sandals`, the projection is performed using a simple iterative method, where the states are updated until the invariants are satisfied within a specified tolerance. In other words, at each step, the states are adjusted to ensure that the invariants are preserved.

## Constrained minimization problem

Consider the constrained minimization problem
\f[
  \underset{\mathbf{x}}{\text{minimize}} \quad \frac{1}{2}\left(\mathbf{x} - \tilde{\mathbf{x}}\right)^2 \quad \text{subject to} \quad \mathbf{h}(\mathbf{x}) = \mathbf{0} \text{,}
\f]
whose Lagrangian is given by
\f[
  \mathcal{L}(\mathbf{x}, \boldsymbol{\lambda}) = \frac{1}{2}\left(\mathbf{x} - \tilde{\mathbf{x}}\right)^2 + \boldsymbol{\lambda} \cdot \mathbf{h}(\mathbf{x})  \text{.}
\f]
The first-order Karush-Kuhn-Tucker conditions for this problem are then
\f[
  \begin{cases}
    \mathbf{x} + \mathbf{Jh}_{\mathbf{x}}^\top \boldsymbol{\lambda} = \tilde{\mathbf{x}} \\
    \mathbf{h}(\mathbf{x}) = \mathbf{0}
  \end{cases} \text{.}
\f]
This system of equations can be solved though the iterative solution of a linear system derived from the Taylor expansion
\f[
  \begin{cases}
    \mathbf{x} + \delta\mathbf{x} + \mathbf{Jh}_{\mathbf{x}}^\top(\mathbf{x} + \delta\mathbf{x}, t) \boldsymbol{\lambda} = \tilde{\mathbf{x}} \\
    \mathbf{h}(\mathbf{x}) + \mathbf{Jh}_{\mathbf{x}}(\mathbf{x}, t) \delta\mathbf{x} + \mathcal{O}\left(\| \delta\mathbf{x} \|^2\right) = \mathbf{0}
  \end{cases} \text{,}
\f]
where \f$\mathbf{Jh}_{\mathbf{x}}\f$ is the Jacobian of the constraint function \f$\mathbf{h}(\mathbf{x})\f$ with respect to the states \f$\mathbf{x}\f$. The linear system can be written in matrix form as
\f[
  \begin{bmatrix}
    \mathbf{I}              & \mathbf{Jh}_{\mathbf{x}}^\top \\
    \mathbf{Jh}_{\mathbf{x}} & \mathbf{0}
  \end{bmatrix}
  \begin{bmatrix}
    \delta\mathbf{x} \\
    \boldsymbol{\lambda}
  \end{bmatrix} = \begin{bmatrix}
    \tilde{\mathbf{x}} - \mathbf{x} \\
    -\mathbf{h}(\mathbf{x})
  \end{bmatrix} \text{,}
\f]
and the update step for the states is then \f$\mathbf{x} = \tilde{\mathbf{x}} + \delta\mathbf{x}\f$.
