# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# System handling
local m_system_type   := "Empty";
local m_vars          := [];
local m_F             := [];
local m_f             := [];
local m_A             := [];
local m_b             := [];
local m_h             := [];

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ResetSystem::static := proc(
  _self::Sandals,
  $)

  description "Reset the system to an empty state.";

  _self:-m_system_type := "Empty";
  _self:-m_vars        := [];
  _self:-m_F           := [];
  _self:-m_f           := [];
  _self:-m_A           := [];
  _self:-m_b           := [];
  _self:-m_h           := [];
  return NULL;
end proc: # ResetSystem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CopySystem::static := proc(
  _self::Sandals,
  proto::Sandals,
  $)

  description "Copy the system from <proto> into <self>.";

  _self:-m_system_type := proto:-m_system_type;
  _self:-m_vars        := proto:-m_vars;
  _self:-m_F           := proto:-m_F;
  _self:-m_f           := proto:-m_f;
  _self:-m_A           := proto:-m_A;
  _self:-m_b           := proto:-m_b;
  _self:-m_h           := proto:-m_h;
  return NULL;
end proc: # CopySystem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetSystemType::static := proc(
  _self::Sandals,
  $)::symbol;

  description "Return the type of the system.";

  return _self:-m_system_type;
end proc: # GetSystemType

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetVars::static := proc(
  _self::Sandals,
  $)::list(function);

  description "Get the vector of state variables.";

  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") or (_self:-m_system_type = "Explicit") or
    (_self:-m_system_type = "SemiExplicit") then
    return _self:-m_vars;
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
end proc: # GetVars

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SetVars::static := proc(
  _self::Sandals,
  vars::list(function),
  {
  desc::list(string) := [seq(cat("State variable ", i), i = 1..nops(vars))]
  }, $)

  description "Set the state variables list to <vars>. Variables' description will be updated to "
    "<desc>.";

  _self:-m_vars := vars;
  _self:-SetVarsInfo(_self, desc);
  return NULL;
end proc: # SetVars

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetF::static := proc(
  _self::Sandals,
  $)::Vector;

  description "Get the vector F(x,x',t) of the implicit system F(x,x',t) = 0.";

  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Explicit") then
    return _self:-m_f;
  elif (_self:-m_system_type = "Implicit") or (_self:-m_system_type = "SemiExplicit") then
    return _self:-m_F;
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
end proc: # GetF

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetA::static := proc(
  _self::Sandals,
  $)::Matrix;

  description "Get the mass matrix of the semi-explicit system.";

  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") or (_self:-m_system_type = "Explicit") or
    (_self:-m_system_type = "SemiExplicit") then
    return _self:-m_A;
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
end proc: # GetA

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetB::static := proc(
  _self::Sandals,
  $)::Vector;

  description "Get the right-hand side vector of the semi-explicit system.";

  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") or (_self:-m_system_type = "Explicit") or
    (_self:-m_system_type = "SemiExplicit") then
    return _self:-m_b;
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
end proc: # GetB

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetH::static := proc(
  _self::Sandals,
  $)::Vector;

  description "Get the invariants of the system.";

  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") or (_self:-m_system_type = "Explicit") or
    (_self:-m_system_type = "SemiExplicit") then
    return _self:-m_h;
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
end proc: # GetH

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetDifferentalEquations::static := proc(
  _self::Sandals,
  {
  invert::boolean := false
  }, $)::list;

  description "Get the differential equations of the system as an implicit system F(x,x',t) = 0, "
    "an explicit system x' = f(x,t), or a semi-explicit system A(x,t)x' = b(x,t). If <invert> is "
    "activated, than the system x' = A(x,t)^{-1}b(x,t) is returned.";

  local vars, out;

  # Store the differential equations
  vars := convert(_self:-m_vars, Vector);
  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") then
    out := _self:-m_F;
  elif (_self:-m_system_type = "Explicit") then
    out := diff(vars, t) - _self:-m_f;
  elif (_self:-m_system_type = "SemiExplicit") and not invert then
    out := _self:-m_A.diff(vars, t) - _self:-m_b;
  elif (_self:-m_system_type = "SemiExplicit") and invert then
    out := diff(vars, t) - LinearAlgebra:-LinearSolve(_self:-m_A, _self:-m_b);
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;

  # Try to simplify
  try
    out := timelimit(_self:-m_time_limit, simplify(out));
  catch "time expired":
    WARNING("time expired, simplify(ODEs) interrupted.");
  end try;

  # Return the results
  return convert(out, list);
end proc: # GetDifferentialEquations

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetInvariants::static := proc(
  _self::Sandals,
  $)::list;

  description "Get the invariants of the system.";

  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") or (_self:-m_system_type = "Explicit") or
    (_self:-m_system_type = "SemiExplicit") then
    return convert(_self:-m_h, list);
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
end proc: # GetInvariants

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetEquations::static := proc(
  _self::Sandals,
  {
  invert::boolean := false
  }, $)::list;

  description "Get the differential equations of the system as an implicit system F(x,x',t) = 0, "
    "an explicit system x' = f(x,t), or a semi-explicit system A(x,t)x' = b(x,t), plus the invariants "
    "h(x,t) = 0. If <invert> is activated, than the system x' = A(x,t)^{-1}b(x,t) is returned.";

  return [
    op(_self:-GetDifferentalEquations(_self, parse("invert") = invert)),
    op(_self:-GetInvariants(_self))
  ];
end proc: # GetEquations

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DiffOrder::static := proc(
  eqns::list,
  vars::list,
  $)::nonnegint;

  description "Return the differential order of the system.";

  if (nops(eqns) = 0) then
    return 0;
  else
    return max(map(y -> op(map(x -> PDEtools:-difforder(x), eqns)), vars));
  end if;
end proc: # SystemDegree

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LoadSystem::static := proc(
  _self::Sandals,
  vars::{list(function), Vector(function)},
  eqns::{list({`=`, algebraic}), Vector({`=`, algebraic})},
  invs::{list({`=`, algebraic}), Vector({`=`, algebraic})} := [],
  {
  type::string := "Implicit"
  }, $)

  description "Load the equations from a system of type <type> with equations <eqns>, states "
    "<vars>, equations <eqns>, and optional invariants <invs>.";

  if (type = "Implicit") then
    return _self:-LoadImplicitSystem(_self, vars, eqns, invs);
  elif (type = "Explicit") then
    return _self:-LoadExplicitSystem(_self, vars, eqns, invs, parse("invert") = true);
  elif (type = "SemiExplicit") then
    return _self:-LoadExplicitSystem(_self, vars, eqns, invs, parse("invert") = false);
  else
    error("unknown system type ""%1"".", type);
  end if;

  return NULL;
end proc: # LoadSystem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LoadImplicitSystem::static := proc(
  _self::Sandals,
  vars::{list(function), Vector(function)},
  eqns::{list({`=`, algebraic}), Vector({`=`, algebraic})},
  invs::{list({`=`, algebraic}), Vector({`=`, algebraic})} := [],
  $)

  description "Load the equations from an implicit system F(x,x',t) = 0 with equations <eqns>, "
    "states <vars>, equations <eqns>, and optional invariants <invs>.";

  local vars_tmp, eqns_tmp, invs_tmp;

  # Check if the system differential order is 1 on the equations
  vars_tmp := `if`(not type(vars, list), convert(vars, list), vars);
  eqns_tmp := `if`(not type(eqns, list), convert(eqns, list), eqns);
  eqns_tmp := map(x -> `if`(type(x, `=`), lhs(x) - rhs(x), x), eqns_tmp);
  if (_self:-DiffOrder(eqns_tmp, vars_tmp) <> 1) then
    error("equations' differential order must be 1 but got %1.", _self:-DiffOrder(eqns_tmp, vars_tmp));
  end if;

  # Check if the system differential order is 0 on the invariants
  invs_tmp := `if`(not type(invs, list), convert(invs, list), invs);
  invs_tmp :=  map(x -> `if`(type(x, `=`), lhs(x) - rhs(x), x), invs_tmp);
  if (_self:-DiffOrder(invs_tmp, vars_tmp) <> 0) then
    error("invariants' differential order must be 0 but got %1.", _self:-DiffOrder(invs_tmp, vars_tmp));
  end if;

  # Reset the system
  _self:-ResetSystem(_self);

  # Set the state variables
  _self:-SetVars(_self, vars_tmp);

  # Store the system data
  _self:-m_F := `if`(not type(eqns_tmp, Vector), convert(eqns_tmp, Vector), eqns_tmp);
  _self:-m_h := convert(invs_tmp, Vector);;
  _self:-m_system_type := "Implicit";

  # Check if the differential order of h is null
  if has(_self:-m_h, diff) then
    error("invalid invariants detected.");
  end if;

  # Simplify the vector F
  try
    _self:-m_F := timelimit(_self:-m_time_limit, simplify(_self:-m_F));
  catch "time expired":
    WARNING("time expired, simplify(F) interrupted.");
  end try;

  # Simplify the vector h
  try
    _self:-m_h := timelimit(_self:-m_time_limit, simplify(_self:-m_h));
  catch "time expired":
    WARNING("time expired, simplify(h) interrupted.");
  end try;

  return NULL;
end proc: # LoadImplicitSystem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LoadExplicitSystem::static := proc(
  _self::Sandals,
  vars::{list(function), Vector(function)},
  eqns::{list({`=`, algebraic}), Vector({`=`, algebraic})},
  invs::{list({`=`, algebraic}), Vector({`=`, algebraic})} := [],
  {
  invert::boolean := true
  }, $)

  description "Load the equations from an explicit system x' = f(x,t) or semi-explicit system "
    "M(x,t)x' = b(x,t) with equations <eqns>, states <vars>, equations <eqns>, and optional "
    "invariants <invs>. If <invert> is activated, than the system is transformed into an explicit "
    "form x' = A(x,t)^{-1}b(x,t).";

  local vars_tmp, eqns_tmp, invs_tmp;

  # Check if the system differential order is 1 on the equations
  vars_tmp := `if`(not type(vars, list), convert(vars, list), vars);
  eqns_tmp := `if`(not type(eqns, list), convert(eqns, list), eqns);
  eqns_tmp := map(x -> `if`(type(x, `=`), lhs(x) - rhs(x), x), eqns_tmp);
  if (_self:-DiffOrder(eqns_tmp, vars_tmp) <> 1) then
    error("equations' differential order must be 1 but got %1.", _self:-DiffOrder(eqns_tmp, vars_tmp));
  end if;

  # Check if the system differential order is 0 on the invariants
  invs_tmp := `if`(not type(invs, list), convert(invs, list), invs);
  invs_tmp :=  map(x -> `if`(type(x, `=`), lhs(x) - rhs(x), x), invs_tmp);
  if (_self:-DiffOrder(invs_tmp, vars_tmp) <> 0) then
    error("invariants' differential order must be 0 but got %1.", _self:-DiffOrder(invs_tmp, vars_tmp));
  end if;

  # Reset the system
  _self:-ResetSystem(_self);

  # Set the state variables
  _self:-SetVars(_self, vars_tmp);

  # Store the system data
  _self:-m_h := convert(invs_tmp, Vector);
  _self:-m_A, _self:-m_b := LinearAlgebra:-GenerateMatrix(eqns_tmp, diff(vars_tmp, t));

  # Check if the system can be transformed into an explicit form
  if (LinearAlgebra:-RowDimension(_self:-m_A) <> LinearAlgebra:-ColumnDimension(_self:-m_A)) then
    _self:-ResetSystem(_self);
    error("the system cannot be transformed into an explicit form, the "
      "matrix M is not square.");
    return NULL;
  end if;

  if (LinearAlgebra:-Rank(_self:-m_A) <> LinearAlgebra:-ColumnDimension(_self:-m_A)) then
    _self:-ResetSystem(_self);
    error("the system cannot be transformed into an explicit form, the "
      "matrix M is not full rank.");
    return NULL;
  end if;

  # Store the system data
  if invert then
    _self:-m_f := LinearAlgebra:-LinearSolve(_self:-m_A, _self:-m_b);
    _self:-m_system_type := "Explicit";
    _self:-m_A := [];
    _self:-m_b := [];
  else
    _self:-m_system_type := "SemiExplicit";
  end if;

  # Check if the differential order of h is null
  if has(_self:-m_h, diff) then
    error("invalid invariants detected.");
  end if;

  # Simplify the vector F
  try
    _self:-m_F := timelimit(_self:-m_time_limit, simplify(_self:-m_F));
  catch "time expired":
    WARNING("time expired, simplify(F) interrupted.");
  end try;

  # Simplify the vector f
  try
    _self:-m_f := timelimit(_self:-m_time_limit, simplify(_self:-m_f));
  catch "time expired":
    WARNING("time expired, simplify(f) interrupted.");
  end try;

  # Simplify the matrix A
  try
    _self:-m_A := timelimit(_self:-m_time_limit, simplify(_self:-m_A));
  catch "time expired":
    WARNING("time expired, simplify(A) interrupted.");
  end try;

  # Simplify the vector b
  try
    _self:-m_b := timelimit(_self:-m_time_limit, simplify(_self:-m_b));
  catch "time expired":
    WARNING("time expired, simplify(b) interrupted.");
  end try;

  # Simplify the vector h
  try
    _self:-m_h := timelimit(_self:-m_time_limit, simplify(_self:-m_h));
  catch "time expired":
    WARNING("time expired, simplify(h) interrupted.");
  end try;

  return NULL;
end proc: # LoadExplicitSystem

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ToImplicit::static := proc(
  _self::Sandals,
  $)

  description "Convert the system to an implicit form F(x,x',t) = 0.";

  local vars, eqns;

  vars := convert(_self:-m_vars, Vector);
  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") then
    # Nothing to do
  elif (_self:-m_system_type = "Explicit") then
    eqns := diff(vars, t) - _self:-m_f;
  elif (_self:-m_system_type = "SemiExplicit") then
    eqns := _self:-m_A.diff(vars, t) - _self:-m_b;
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
  _self:-LoadImplicitSystem(_self, vars, eqns, _self:-m_h);

  return NULL;
end proc: # ToImplicit

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ToExplicit::static := proc(
  _self::Sandals,
  $)

  description "Convert the system to an explicit form x' = f(x,t).";

  local vars, eqns;

  vars := convert(_self:-m_vars, Vector);
  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") then
    eqns := _self:-m_F;
  elif (_self:-m_system_type = "Explicit") then
    # Nothing to do
  elif (_self:-m_system_type = "SemiExplicit") then
    eqns := _self:-m_A.diff(vars, t) - _self:-m_b
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
  _self:-LoadExplicitSystem(_self, vars, eqns, _self:-m_h, parse("invert") = true);

  return NULL;
end proc: # ToExplicit

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ToSemiExplicit::static := proc(
  _self::Sandals,
  $)

  description "Convert the system to a semi-explicit form A(x,t)x' = b(x,t).";

  local vars, eqns;

  vars := convert(_self:-m_vars, Vector);
  if (_self:-m_system_type = "Empty") then
    error("system not yet available, please load the system first.");
  elif (_self:-m_system_type = "Implicit") then
    eqns := _self:-m_F;
  elif (_self:-m_system_type = "Explicit") then
    eqns := diff(vars, t) - _self:-m_f;
  elif (_self:-m_system_type = "SemiExplicit") then
    # Nothing to do
  else
    error("unknown system type ""%1"".", _self:-m_system_type);
  end if;
  _self:-LoadExplicitSystem(_self, vars, eqns, _self:-m_h, parse("invert") = false);

  return NULL;
end proc: # ToSemiExplicit

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
