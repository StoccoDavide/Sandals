
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Chemical Akzo Nobel problem
# https://archimede.uniba.it/~testset/problems/chemakzo.php

# Model
NAME  := "ChemAkzo":
TYPE  := "Implicit":
INDEX := 1:
TIME  := [0.0, 180.0]:
DIMS  := [5, 1]:

# Variables
VARS := [seq(y_||i(t), i = 1..add(DIMS))];

# Mass matrix
MASS_MATRIX := Matrix(add(DIMS)):
for i to DIMS[1] do
  MASS_MATRIX[i,i] := 1
end do:

# Auxiliar variables
r_1  := k_1 * y_1(t)^4 * y_2(t)^(1/2):
r_2  := k_2 * y_3(t) * y_4(t):
r_3  := k_2/K * y_1(t) * y_5(t):
r_4  := k_3 * y_1(t) * y_4(t)^2:
r_5  := k_4 * y_6(t)^2 * y_2(t)^(1/2):
F_in := klA * (p_CO2/H - y_2(t)):

# Right-hand side
RHS_VECTOR := <
  -2 * r_1 + r_2 - r_3 - r_4,
  -r_1 / 2 - r_4 - r_5 / 2 + F_in,
  r_1 - r_2 + r_3,
  -r_2 + r_3 - 2 * r_4,
  r_2 - r_3 + r_5,
  K_s * y_1(t) * y_4(t) - y_6(t)
>;

# Create system of DAEs (Implicit)
EQNS := convert(MASS_MATRIX.diff(convert(VARS, Vector), t) =~ RHS_VECTOR, list):
EQNS := map(x -> `if`(type(denom(x), numeric), numer(x), x), lhs~(EQNS) - rhs~(EQNS)) =~0:

# Parameters
DATA := [
  k_1   = 18.7,
  k_2   = 0.58,
  k_3   = 0.09,
  k_4   = 0.42,
  K     = 34.4,
  klA   = 3.3,
  K_s   = 115.83,
  p_CO2 = 0.9,
  H     = 737.0
];

# Initial conditions
ICS := [
  0.444,
  0.00123,
  0.0,
  0.007,
  0.0,
  K_s * 0.444 * 0.007 # K_s * ICS[1] * ICS[4]
];
