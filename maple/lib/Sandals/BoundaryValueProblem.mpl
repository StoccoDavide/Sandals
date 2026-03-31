# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2026, Davide Stocco and Enrico Bertolazzi.                    #
#                                                                             #
# The Sandals project is distributed under the BSD 2-Clause License.          #
#                                                                             #
# Davide Stocco                                             Enrico Bertolazzi #
# University of Trento                                   University of Trento #
# davide.stocco@unitn.it                           enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Boundary Value Problem handling
local m_b_ini := [];
local m_b_end := [];

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ResetBoundaryValueProblem::static := proc(
  _self::Sandals,
  $)

  description "Reset the BVP to an empty state.";

  _self:-m_b_ini := [];
  _self:-m_b_end := [];
  return NULL;
end proc: # ResetBoundaryValueProblem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CopyBoundaryValueProblem::static := proc(
  _self::Sandals,
  proto::Sandals,
  $)

  description "Copy the BVP from <proto> into <self>.";

  _self:-m_b           := proto:-m_b;
  _self:-m_h           := proto:-m_h;
  return NULL;
end proc: # CopyBoundaryValueProblem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetBini::static := proc(
  _self::Sandals,
  $)::Vector;

  description "Get the initial boundary conditions vector.";

  return _self:-m_b_ini;
end proc: # GetBini

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetBend::static := proc(
  _self::Sandals,
  $)::Vector;

  description "Get the final boundary conditions vector.";

  return _self:-m_b_end;
end proc: # GetBend

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetB::static := proc(
  _self::Sandals,
  $)::Vector;

  description "Get the boundary conditions vector.";

  return <_self:-m_b_ini; _self:-m_b_end>;
end proc: # GetB

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
