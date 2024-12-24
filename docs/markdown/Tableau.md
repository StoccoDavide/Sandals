# Butcher tableau

Runge-Kutta methods are a family of numerical techniques for solving ODEs and DAEs. These methods approximate the solution by combining weighted evaluations of the system's functions at specific points called *nodes* within the integration domain.

A Runge-Kutta method is characterized by its Butcher tableau, a tabular representation of the method's parameters. For a method of order \f$p\f$ and \f$s\f$ stages, the Butcher tableau is defined as
\f[
  \begin{array}{c|ccc}
    c_1 & a_{11} & a_{12} & \dots & a_{1s} \\
    c_2 & a_{21} & a_{22} & \dots & a_{2s} \\
        \vdots & \vdots & \vdots & \ddots & \vdots \\
    c_s    & a_{s1} & a_{s2} & \dots & a_{ss} \\ \hline
    & b_1 & b_2 & \dots & b_s
  \end{array} \text{,}
\f]
where \f$\mathbf{A}\f$ is the Runge-Kutta coefficient matrix, \f$\mathbf{b}\f$ is the Runge-Kutta weights vector, \f$\hat{\mathbf{b}}\f$ is the (optional) embedded Runge-Kutta weights vector, and \f$\mathbf{c}\f$ is the Runge-Kutta nodes vector. These parameters define the method's structure and determine its accuracy and stability.
