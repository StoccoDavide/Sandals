# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

NUM_TESTS := 33:
NAME  := "BVP":
TIME  := [seq([], i=1..NUM_TESTS)]:
DIMS  := 2:

# Variables
VARS := [seq(y_||i(t), i = 1..add(DIMS))];

# Mass matrix
MASS_MATRIX := Matrix(add(DIMS)):
for i to DIMS[1] do
  MASS_MATRIX[i, i] := 1
end do:

# Right-hand side of the dynamic system
RHS_VECTOR := (func) -> <y_2(t), func(y_1(t), y_1(t))/lambda>;

# Create system of dynamic system
EQNS := convert(MASS_MATRIX.diff(convert(VARS, Vector), t) =~ RHS_VECTOR, list):
EQNS := map(x -> `if`(type(denom(x), numeric), numer(x), x), lhs~(EQNS) - rhs~(EQNS)) =~0:

# Initial conditions
ICS := [0.0, 0.0];

# Parameters
DATA := [lambda = 1.0e-3];

# Create an Sandals object
Sandals_obj := Object(Sandals);

# Settings
Sandals_obj:-SetTimeLimit(Sandals_obj, 2.0);
Sandals_obj:-EnableVerboseMode(Sandals_obj);
Sandals_obj:-EnableWarningMode(Sandals_obj);

# Retrieve directories contents
base_dir := "./GeneratedTestSet/";

folders := FileTools:-ListDirectory(base_dir);
for folder in folders do
  files := FileTools:-ListDirectory(cat(base_dir, folder));
  for file in files do
    if StringTools:-IsSuffix(".mpl", file) then

      # Print the file name
      printf("PROCESSING FILE: ""%s""\n", file);

      # Load the file
      try
        read(cat(base_dir, folder, "/", file));
      catch:
        printf("\tERROR reading file: %s\n", file);
      end try;

      # Reset the object
      Sandals_obj:-Reset(Sandals_obj);

      # Load the equations
      Sandals_obj:-LoadSystem(Sandals_obj, VARS, EQNS, parse("type") = TYPE);

      # Code generation settings
      Sandals_obj:-SetName(Sandals_obj, cat(NAME, i));
      #Sandals_obj:-SetInfo(Sandals_obj, ...);
      Sandals_obj:-SetPath(Sandals_obj, cat(base_dir, folder));
      #Sandals_obj:-SetVarsInfo(Sandals_obj, ...);
      Sandals_obj:-SetData(Sandals_obj, DATA);
      #Sandals_obj:-SetDataInfo(Sandals_obj, ...);
      #Sandals_obj:-SetDomain(Sandals_obj, ...);
      #Sandals_obj:-SetTime(Sandals_obj, TIME);
      #Sandals_obj:-SetCompSequence(Sandals_obj, ...);
      #Sandals_obj:-SetUserFunction(Sandals_obj, ...);
      #Sandals_obj:-SetIcs(Sandals_obj, ICS);
      #Sandals_obj:-SetGuess(Sandals_obj, ICS);

      # Generate the code
      Sandals_obj:-GenerateCode(Sandals_obj);

      # Print the content of the file
      printf("DONE: %s\n", file);
    end if;
  end do;
end do;
