# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

unprotect('Sandals');
module Sandals()

  description "'Sandals' module.";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local m_VerboseMode := false;
  local m_WarningMode := true;
  local m_TimeLimit   := 0.1;
  local m_SystemType  := "Empty";
  local m_x           := [];
  local m_F           := [];
  local m_f           := [];
  local m_A           := [];
  local m_b           := [];
  local m_h           := [];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print 'Sandals' module information.";

    printf(
      "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n"
      "* Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *\n"
      "*                                                                                               *\n"
      "* The Sandals project is distributed under the BSD 2-Clause License.                            *\n"
      "*                                                                                               *\n"
      "* Davide Stocco                                                               Enrico Bertolazzi *\n"
      "* University of Trento                                                     University of Trento *\n"
      "* e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *\n"
      "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'Sandals' module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("Sandals", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "cannot find 'Sandals' module.";
    end if;
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()

    description "'Sandals' module unload procedure.";

    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleCopy::static := proc(
    _self::Sandals,
    proto::Sandals,
    $)

    description "Copy the objects <proto> into <self>.";

    _self:-m_VerboseMode := proto:-m_VerboseMode;
    _self:-m_WarningMode := proto:-m_WarningMode;
    _self:-m_TimeLimit   := proto:-m_TimeLimit;
    _self:-m_SystemType  := proto:-m_SystemType;
    _self:-m_x           := proto:-m_x;
    _self:-m_F           := proto:-m_F;
    _self:-m_f           := proto:-m_f;
    _self:-m_A           := proto:-m_A;
    _self:-m_b           := proto:-m_b;
    _self:-m_h           := proto:-m_h;
  end proc: # ModuleCopy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode::static := proc(
    _self::Sandals,
    $)

    description "Enable the verbose mode.";

    _self:-m_VerboseMode := true;
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode::static := proc(
    _self::Sandals,
    $)

    description "Disable the verbose mode.";

    _self:-m_VerboseMode := false;
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerboseMode::static := proc(
    _self::Sandals,
    mode::boolean,
    $)

    description "Set the verbosity of the module to <mode>.";

    _self:-m_VerboseMode := mode;
    return NULL;
  end proc: # SetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode::static := proc(
    _self::Sandals,
    $)

    description "Enable the warning mode of the module.";

    _self:-m_WarningMode := true;
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode::static := proc(
    _self::Sandals,
    $)

    description "Disable the warning mode of the module.";

    _self:-m_WarningMode := false;
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWarningMode::static := proc(
    _self::Sandals,
    mode::boolean,
    $)

    description "Set the warning mode of the module to <mode>.";

    _self:-m_WarningMode := mode;
    return NULL;
  end proc: # SetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTimeLimit::static := proc(
    _self::Sandals,
    x::numeric,
    $)

    description "Set the time limit of the module to <x>.";

    if (x < 0) then
      error("time limit must be a non-negative number.");
    end if;

    _self:-m_TimeLimit := x;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetTimeLimit::static := proc(
    _self::Sandals,
    $)::numeric;

    description "Get the time limit of the module.";

    return _self:-m_TimeLimit;
  end proc: # GetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Reset::static := proc(
    _self::Sandals,
    label::{symbol, string} := NULL,
    $)

    description "Reset the library.";

    _self:-m_SystemType := "Empty";
    _self:-m_x := [];
    _self:-m_F := [];
    _self:-m_f := [];
    _self:-m_A := [];
    _self:-m_b := [];
    _self:-m_h := [];
    return NULL;
  end proc: # Reset

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetSystemType::static := proc(
    _self::Sandals,
    $)::symbol;

    description "Return the type of the system.";

    return _self:-m_SystemType;
  end proc: # GetSystemType

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetX::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the vector of state variables.";

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") or (_self:-m_SystemType = "Explicit") or
      (_self:-m_SystemType = "SemiExplicit") then
      return _self:-m_x;
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
  end proc: # GetX

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetF::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the vector F(x,x',t) of the implicit system F(x,x',t) = 0.";

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Explicit") then
      return _self:-m_f;
    elif (_self:-m_SystemType = "Implicit") or (_self:-m_SystemType = "SemiExplicit") then
      return _self:-m_F;
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
  end proc: # GetF

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetA::static := proc(
    _self::Sandals,
    $)::Matrix;

    description "Get the mass matrix of the semi-explicit system.";

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") or (_self:-m_SystemType = "Explicit") or
      (_self:-m_SystemType = "SemiExplicit") then
      return _self:-m_A;
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
  end proc: # GetA

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetB::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the right-hand side vector of the semi-explicit system.";

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") or (_self:-m_SystemType = "Explicit") or
      (_self:-m_SystemType = "SemiExplicit") then
      return _self:-m_b;
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
  end proc: # GetB

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetH::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the invariants of the system.";

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") or (_self:-m_SystemType = "Explicit") or
      (_self:-m_SystemType = "SemiExplicit") then
      return _self:-m_h;
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
  end proc: # GetH

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetDifferentalEquations::static := proc(
    _self::Sandals,
    {
    invert::boolean := false
    }, $)::list;

    description "Get the differential equations of the system as an implicit system F(x,x',t) = 0, "
      "an explicit system x' = f(x,t), or a semi-explicit system A(x,t)x' = b(x,t). If <invert> is "
      "activated, than the system x' = A(x,t)^{-1}b(x,t) is returned.";

    local out;

    # Store the differential equations
    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") then
      out := _self:-m_F;
    elif (_self:-m_SystemType = "Explicit") then
      out := diff(_self:-m_x, t) - _self:-m_f;
    elif (_self:-m_SystemType = "SemiExplicit") and not invert then
      out := _self:-m_A.diff(_self:-m_x, t) - _self:-m_b;
    elif (_self:-m_SystemType = "SemiExplicit") and invert then
      out := diff(_self:-m_x, t) - LinearAlgebra:-LinearSolve(_self:-m_A, _self:-m_b);
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;

    # Try to simplify
    try
      out := timelimit(_self:-m_TimeLimit, simplify(out));
    catch "time expired":
      WARNING("time expired, simplify(ODEs) interrupted.");
    end try;

    # Return the results
    return convert(out, list);
  end proc: # GetDifferentialEquations

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetInvariants::static := proc(
    _self::Sandals,
    $)::list;

    description "Get the invariants of the system.";

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") or (_self:-m_SystemType = "Explicit") or
      (_self:-m_SystemType = "SemiExplicit") then
      return convert(_self:-m_h, list);
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
  end proc: # GetInvariants

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DiffOrder := proc(
    _self::Sandals,
    eqns::list,
    vars::list,
    $)::nonnegint;

    description "Return the differential order of the system.";

    return max(map(y -> op(map(x -> PDEtools:-difforder(x), eqns)), vars))
  end proc: # SystemDegree

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LoadMatrices::static := proc(
    _self::Sandals,
    vars::{list, Vector},
    eqns::{list, Vector},
    invs::{list, Vector} := [],
    {
    type::{string} := "Implicit"
    }, $)

    description "Load the equations from a system of type <type> with equations <eqns>, states "
      "<vars>, equations <eqns>, and optional invariants <invs>.";

    # Check if the system differential order is 1
    if (_self:-DiffOrder(_self, eqns, vars) <> 1) then
      error("differential order must be 1 but got %1.", _self:-DiffOrder(_self, eqns, vars));
    end if;

    if (type = "Implicit") then
      return _self:-LoadImplicitSystem(_self, vars, eqns, invs);
    elif (type = "Explicit") then
      return _self:-LoadExplicitSystem(_self, vars, eqns, invs, parse("invert") = true);
    elif (type = "SemiExplicit") then
      return _self:-LoadExplicitSystem(_self, vars, eqns, invs, parse("invert") = false);
    else
      error("unknown system type '%1'.", type);
    end if;

    return NULL;
  end proc: # LoadSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LoadSystem::static := proc(
    _self::Sandals,
    vars::{list, Vector},
    eqns::{list, Vector},
    invs::{list, Vector} := [],
    {
    type::{string} := "Implicit"
    }, $)

    description "Load the equations from a system of type <type> with equations <eqns>, states "
      "<vars>, equations <eqns>, and optional invariants <invs>.";

    # Check if the system differential order is 1
    if (_self:-DiffOrder(_self, eqns, vars) <> 1) then
      error("differential order must be 1 but got %1.", _self:-DiffOrder(_self, eqns, vars));
    end if;

    if (type = "Implicit") then
      return _self:-LoadImplicitSystem(_self, vars, eqns, invs);
    elif (type = "Explicit") then
      return _self:-LoadExplicitSystem(_self, vars, eqns, invs, parse("invert") = true);
    elif (type = "SemiExplicit") then
      return _self:-LoadExplicitSystem(_self, vars, eqns, invs, parse("invert") = false);
    else
      error("unknown system type '%1'.", type);
    end if;

    return NULL;
  end proc: # LoadSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LoadImplicitSystem::static := proc(
    _self::Sandals,
    vars::{list, Vector},
    eqns::{list, Vector},
    invs::{list, Vector} := [],
    $)

    description "Load the equations from an implicit system F(x,x',t) = 0 with equations <eqns>, "
      "states <vars>, equations <eqns>, and optional invariants <invs>.";

      # Reset the system
      _self:-Reset(_self);

    # Store the system data
    _self:-m_x := `if`(not type(vars, Vector), convert(vars, Vector), vars);
    _self:-m_F := `if`(not type(eqns, Vector), convert(eqns, Vector), eqns);
    _self:-m_h := `if`(not type(invs, list), convert(invs, list), invs);
    _self:-m_SystemType := "Implicit";

    # Check if the differential order of h is null
    if has(_self:-m_h, diff) then
      error("invalid invariants detected.");
    end if;

    # Simplify the vector F
    try
      _self:-m_F := timelimit(_self:-m_TimeLimit, simplify(_self:-m_F));
    catch "time expired":
      WARNING("time expired, simplify(F) interrupted.");
    end try;

    # Simplify the vector h
    try
      _self:-m_h := timelimit(_self:-m_TimeLimit, simplify(_self:-m_h));
    catch "time expired":
      WARNING("time expired, simplify(h) interrupted.");
    end try;

    return NULL;
  end proc: # LoadImplicitSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LoadExplicitSystem::static := proc(
    _self::Sandals,
    vars::{list, Vector},
    eqns::{list, Vector},
    invs::{list, Vector} := [],
    {
    invert::boolean := true
    }, $)

    description "Load the equations from an explicit system x' = f(x,t) or semi-explicit system "
      "M(x,t)x' = b(x,t) with equations <eqns>, states <vars>, equations <eqns>, and optional "
      "invariants <invs>. If <invert> is activated, than the system is transformed into an explicit "
      "form x' = A(x,t)^{-1}b(x,t).";

    # Reset the system
    _self:-Reset(_self);

    # Store the system data
    _self:-m_x := `if`(not type(vars, Vector), convert(vars, Vector), vars);
    _self:-m_h := `if`(not type(invs, list), convert(invs, list), invs);
    _self:-m_A, _self:-m_b := LinearAlgebra:-GenerateMatrix(
      `if`(not type(eqns, list), convert(eqns, list), eqns),
      diff(`if`(not type(vars, list), convert(vars, list), vars), t)
    );

    # Check if the system can be transformed into an explicit form
    if (LinearAlgebra:-RowDimension(_self:-m_A) <> LinearAlgebra:-ColumnDimension(_self:-m_A)) then
      _self:-Reset(_self);
      error("the system cannot be transformed into an explicit form, the "
        "matrix M is not square.");
      return NULL;
    end if;

    if (LinearAlgebra:-Rank(_self:-m_A) <> LinearAlgebra:-ColumnDimension(_self:-m_A)) then
      _self:-Reset(_self);
      error("the system cannot be transformed into an explicit form, the "
        "matrix M is not full rank.");
      return NULL;
    end if;

    # Store the system data
    if invert then
      _self:-m_f := LinearAlgebra:-LinearSolve(_self:-m_A, _self:-m_b);
      _self:-m_SystemType := "Explicit";
      _self:-m_A := [];
      _self:-m_b := [];
    else
      _self:-m_SystemType := "SemiExplicit";
    end if;

    # Check if the differential order of h is null
    if has(_self:-m_h, diff) then
      error("invalid invariants detected.");
    end if;

    # Simplify the vector F
    try
      _self:-m_F := timelimit(_self:-m_TimeLimit, simplify(_self:-m_F));
    catch "time expired":
      WARNING("time expired, simplify(F) interrupted.");
    end try;

    # Simplify the vector f
    try
      _self:-m_f := timelimit(_self:-m_TimeLimit, simplify(_self:-m_f));
    catch "time expired":
      WARNING("time expired, simplify(f) interrupted.");
    end try;

    # Simplify the matrix A
    try
      _self:-m_A := timelimit(_self:-m_TimeLimit, simplify(_self:-m_A));
    catch "time expired":
      WARNING("time expired, simplify(A) interrupted.");
    end try;

    # Simplify the vector b
    try
      _self:-m_b := timelimit(_self:-m_TimeLimit, simplify(_self:-m_b));
    catch "time expired":
      WARNING("time expired, simplify(b) interrupted.");
    end try;

    # Simplify the vector h
    try
      _self:-m_h := timelimit(_self:-m_TimeLimit, simplify(_self:-m_h));
    catch "time expired":
      WARNING("time expired, simplify(h) interrupted.");
    end try;

    return NULL;
  end proc: # LoadExplicitSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ToImplicit::static := proc(
    _self::Sandals,
    $)

    description "Convert the system to an implicit form F(x,x',t) = 0.";

    local eqns;

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") then
      # Nothing to do
    elif (_self:-m_SystemType = "Explicit") then
      eqns := diff(_self:-m_x, t) - _self:-m_f;
    elif (_self:-m_SystemType = "SemiExplicit") then
      eqns := _self:-m_A.diff(_self:-m_x, t) - _self:-m_b;
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
    _self:-LoadImplicitSystem(_self, _self:-m_x, eqns, _self:-m_h);

    return NULL;
  end proc: # ToImplicit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ToExplicit::static := proc(
    _self::Sandals,
    $)

    description "Convert the system to an explicit form x' = f(x,t).";

    local eqns;

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") then
      eqns := _self:-m_F;
    elif (_self:-m_SystemType = "Explicit") then
      # Nothing to do
    elif (_self:-m_SystemType = "SemiExplicit") then
      eqns := _self:-m_A.diff(_self:-m_x, t) - _self:-m_b
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
    _self:-LoadExplicitSystem(_self, _self:-m_x, eqns, _self:-m_h, parse("invert") = true);

    return NULL;
  end proc: # ToExplicit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ToSemiExplicit::static := proc(
    _self::Sandals,
    $)

    description "Convert the system to a semi-explicit form A(x,t)x' = b(x,t).";

    local eqns;

    if (_self:-m_SystemType = "Empty") then
      error("system not yet available, please load the system first.");
    elif (_self:-m_SystemType = "Implicit") then
      eqns := _self:-m_F;
    elif (_self:-m_SystemType = "Explicit") then
      eqns := diff(_self:-m_x, t) - _self:-m_f;
    elif (_self:-m_SystemType = "SemiExplicit") then
      # Nothing to do
    else
      error("unknown system type '%1'.", _self:-m_SystemType);
    end if;
    _self:-LoadExplicitSystem(_self, _self:-m_x, eqns, _self:-m_h, parse("invert") = false);

    return NULL;
  end proc: # ToSemiExplicit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TranslateToCpp::static := proc(
    _self::Sandals,
    name::string,
    {
    data::list(symbol = algebraic)               := [],
    time::{list(numeric), range(numeric)}        := [],
    ics::{Vector(algebraic), list(algebraic)}    := [],
    domain::{Vector(algebraic), list(algebraic)} := [],
    info::string                                 := "No class description provided."
    }, $)::string;

    description "Generate C++ code for the loaded system with name <name>, Sandals class <type>, "
      "output file './<fname>.m', optional internal class data <data>, integration time range "
      "<time>, initial conditions <ics>, domain <domain> d(x,t) > 0, and class information string "
      "<info>.";

    local time_lst, ics_vec;

    # Check if the system is loaded
    if (_self:-m_SystemType = "Empty") then
      error("no system loaded yet.");
    end if;

    # Check if time range is valid
    if not type(time, list) then
      time_lst := convert(time, list)
    else
      time_lst := time;
    end if;

    if (nops(time_lst) > 0) then
      if (time_lst[1] = time_lst[2]) then
        error("degenerate time range detected.");
      elif (time_lst[2] < time_lst[1]) then
        WARNING("decreasing time range detected, consider using the reverse time mode.");
      end if;
    end if;

    # Check if the initial conditions are valid
    if not type(ics, Vector) then
      ics_vec := convert(ics, Vector)
    else
      ics_vec := ics;
    end if;

    if (LinearAlgebra:-Dimension(ics_vec) > 0) then
      if (LinearAlgebra:-Dimension(ics_vec) <> LinearAlgebra:-Dimension(_self:-m_x)) then
        error("invalid initial conditions detected.");
      elif (nops(indets(ics_vec) minus convert(lhs~(data), set)) > 0) then
        WARNING("initial conditions may contain unknown variables.");
      end if;
    end if;

    # Check if the domain is valid
    # FIXME:if nops(domain) > 0 then
    #   if (nops(indets(domain) minus lhs~(data) minus {t, op(_self:-m_x)}) > 0) then
    #     WARNING("domain may contain unknown variables.");
    #   end if;
    # end if;

    # Generate class body string
    if (_self:-m_SystemType = "Implicit") then
      return SandalsCodegen:-ImplicitToCpp(name, _self:-m_x, _self:-m_F, _self:-m_h,
        parse("data") = data, parse("time") = time, parse("ics") = ics, parse("domain") = domain,
        parse("info") = info);
    elif (_self:-m_SystemType = "Explicit") then
      return SandalsCodegen:-ExplicitToCpp(name, _self:-m_x, _self:-m_f, _self:-m_h,
        parse("data") = data, parse("time") = time, parse("ics") = ics, parse("domain") = domain,
        parse("info") = info);
    elif (_self:-m_SystemType = "SemiExplicit") then
      return SandalsCodegen:-SemiExplicitToCpp(name, _self:-m_x, _self:-m_A, _self:-m_b, _self:-m_h,
        parse("data") = data, parse("time") = time, parse("ics") = ics, parse("domain") = domain,
        parse("info") = info);
    else
      error("unknown Sandals class type '%1'.", _self:-m_SystemType);
    end if;
    return "";
  end proc: # TranslateToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateCppCode::static := proc(
    _self::Sandals,
    name::string,
    {
    type::string                                 := _self:-m_SystemType,
    path::string                                 := "./",
    data::list(symbol = algebraic)               := [],
    time::{list(numeric), range(numeric)}        := [],
    ics::{Vector(algebraic), list(algebraic)}    := [],
    domain::{Vector(algebraic), list(algebraic)} := [],
    info::string                                 := "No class description provided."
    }, $)

    description "Generate C++ code for the loaded system with name <name>, Sandals class <type> "
      "(""Implicit"", ""Explicit"" ,""SemiExplicit""), output file '<path>/<name>.hh', optional "
      "internal class data <data>, integration time range <time>, initial conditions <ics>, domain "
      "<domain> d(x,t) > 0, and class information string <info>. Notice that if the type is not "
      "specified, the system type is used. Instead, if the type is specified, and it is different "
      "from the system type, the system is converted to the specified type.";

    if (type <> _self:-m_SystemType) then
      if (type = "Implicit") then
        _self:-ToImplicit(_self);
      elif (type = "Explicit") then
        _self:-ToExplicit(_self);
      elif (type = "SemiExplicit") then
        _self:-ToSemiExplicit(_self);
      else
        error("unknown Sandals class type '%1'.", type);
      end if;
        WARNING("system converted to '%1' class.", type);
    end if;

    SandalsCodegen:-GenerateFile(cat(path, "/", name, ".hh"), _self:-TranslateToCpp(_self, name,
      parse("data") = data, parse("time") = time, parse("ics") = ics, parse("info") = info));
    return NULL;
  end proc: # GenerateCppCode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # Sandals
