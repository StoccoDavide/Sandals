# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                    ____                  _       _                          #
#                   / ___|  __ _ _ __   __| | __ _| |___                      #
#                   \___ \ / _` | '_ \ / _` |/ _` | / __|                     #
#                    ___) | (_| | | | | (_| | (_| | \__ \                     #
#                   |____/ \__,_|_| |_|\__,_|\__,_|_|___/                     #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Davide Stocco     (University of Trento)
#   Enrico Bertolazzi (University of Trento)
#
# License: BSD 2-Clause License
#
# This is a module for the 'Sandals' module.

unprotect('Sandals');
module Sandals()

  description "'Sandals' module.";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local m_VerboseMode := false;
  local m_WarningMode := true;
  local m_TimeLimit   := 0.1;
  local m_SystemType  := 'Empty';
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
      "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n"
      "* Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *\n"
      "*                                                                           *\n"
      "* The Sandals project is distributed under the BSD 2-Clause License.        *\n"
      "*                                                                           *\n"
      "* Davide Stocco                                           Enrico Bertolazzi *\n"
      "* University of Trento                                 University of Trento *\n"
      "* e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *\n"
      "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad::static := proc()

    description "'Sandals' module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if evalb(StringTools:-Search("Sandals", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if evalb(lib_base_path = NULL) then
      error "cannot find 'Sandals' module.";
    end if;
    protect('Empty', 'Implicit', 'Explicit', 'SemiExplicit');
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()

    description "'Sandals' module unload procedure.";

    unprotect('Empty', 'Implicit', 'Explicit', 'SemiExplicit');
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

    if evalb(x < 0) then
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

    _self:-m_SystemType := 'Empty';
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

    if evalb(_self:-m_x = []) then
      error("system not yet available, please load the system first.");
    end if;
    return _self:-m_x;
  end proc: # GetX

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetImplicitF::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the vector F(x,x',t) of the implicit system F(x,x',t) = 0.";

    if evalb(_self:-m_F = []) then
      error("system not yet available, please load the system first.");
    end if;
    return _self:-m_F;
  end proc: # GetImplicitF

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetExplicitF::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the right-hand side vector f(x,t) of the explicit system x' = f(x,t).";

    if evalb(_self:-m_F = []) then
      error("system not yet available, please load the system first.");
    end if;
    return _self:-m_F;
  end proc: # GetExplicitF

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetA::static := proc(
    _self::Sandals,
    $)::Matrix;

    description "Get the mass matrix of the semi-explicit system.";

    if evalb(_self:-m_A = []) then
      error("system not yet available, please load the system first.");
    end if;
    return _self:-m_A;
  end proc: # GetA

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetB::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the force vector of the semi-explicit system.";

    if evalb(_self:-m_b = []) then
      error("system not yet available, please load the system first.");
    end if;
    return _self:-m_b;
  end proc: # GetB

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetSemiExplicitF::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the right-hand side of the semi-explicit system x' = A(x,t)^{-1} b(x,t) = f(x,t) "
      "(no simplification is performed).";

    if evalb(_self:-m_A = []) or evalb(_self:-m_b = []) then
      error("system not yet available, please load the system first.");
    end if;
    return LinearAlgebra:-LinearSolve(_self:-m_A, _self:-m_b);
  end proc: # GetSemiExplicitF

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetH::static := proc(
    _self::Sandals,
    $)::Vector;

    description "Get the invariants of the system.";

    if evalb(_self:-m_h = []) then
      error("system not yet available, please load the system first.");
    end if;
    return _self:-m_h;
  end proc: # GetH

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetDifferentialEquations::static := proc(
    _self::Sandals,
    {
      explicit::boolean := false
    }, $)::list;

    description "Get the differential equations of the system as an implicit system F(x,x',t) = 0, "
    "an explicit system x' = f(x,t), or a semi-explicit system A(x,t)x' = b(x,t) or x' = A(x,t)^{-1} "
    "b(x,t) = f(x,t).";

    local out;

    # Store the differential equations
    if evalb(explicit = true) then
      out := _self:-m_A.diff(_self:-m_x, t) - _self:-m_b;
    else
      out := diff(_self:-m_x, t) - LinearAlgebra:-LinearSolve(_self:-m_A, _self:-m_b);
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

    if evalb(_self:-m_h = []) then
      error("system not yet available, please load the system first.");
    end if;
    return convert(_self:-m_h, list);
  end proc: # GetInvariants

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetSystemEquations::static := proc(
    _self::Sandals,
    {
      explicit::boolean := false
    }, $)::list;

    description "Get the differential equations of the system as an implicit system F(x,x',t) = 0, "
    "an explicit system x' = f(x,t), or a semi-explicit system A(x,t)x' = b(x,t) or x' = A(x,t)^{-1} "
    "b(x,t) = f(x,t), plus the invariants h(x,t) = 0.";

    return [
      op(_self:-GetDifferentialEquations(_self, parse("explicit") = explicit)),
      op(_self:-GetInvariants(_self))
    ];
  end proc: # GetSystemEquations

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export RemoveTimeStates::static := proc(
    _self::Sandals,
    expr::{algebraic, list(algebraic),
           algebraic = algebraic, list(algebraic = algebraic)},
    $)::{algebraic, list(algebraic),
         algebraic = algebraic, list(algebraic = algebraic)};

    description "Remove the time states from the system.";

    return subs(
      op(map(x -> diff(x, t) = cat(op(0, x), __dot), _self:-m_x)),
      op(map(x -> x = op(0, x), _self:-m_x)),
      expr
    );
  end proc: # RemoveTimeStates

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DiffOrder := proc(
    _self::Sandals,
    vars::{list, Vector},
    eqns::{list, Vector},
    $)::nonnegint;

    description "Return the differential order of the system.";

    local vars_tmp, eqns_tmp;

    vars_tmp := `if`(not type(vars, list), convert(vars, list), vars);
    eqns_tmp := `if`(not type(eqns, list), convert(eqns, list), eqns);

    return max(map(y -> op(map(x -> PDEtools:-difforder(x), eqns_tmp)), vars_tmp));
  end proc: # DiffOrder

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LoadMatrices::static := proc(
    _self::Sandals
    # _passed::{list, Vector},
    )

    description "Load the matrices from an explicit system x' = f(x,t) with "
      "and state vector <x> and right-hand side <f>, or a semi-explicit system "
      "A(x,t)x' = b(x,t) with mass matrix <A>, forcing term <b> and state "
      "vector <x>. Lastly, system invariants <h> must be specified.";

    # Check if the system is already loaded
    if evalb(_self:-m_SystemType <> 'Empty') then
      if _self:-m_VerboseMode then
        WARNING("Sandals::LoadMatrices(...): a system of ODEs is already loaded, internal data will "
        "be overwritten.");
      end if;
      _self:-Reset(_self);
    end if;

    # Set system type
    _self:-m_SystemType := 'Loaded';

    # Collect input arguments
    if evalb(_npassed < 4) then
      error("no system to load.");
    elif evalb(_npassed = 4) then
      _self:-m_x := `if`(not type(_passed[2], Vector), convert(_passed[2], Vector), _passed[2]);
      _self:-m_A := Matrix(LinearAlgebra:-Dimension(_self:-m_x), shape = identity);
      _self:-m_b := `if`(not type(_passed[3], Vector), convert(_passed[3], Vector), _passed[3]);
      _self:-m_h := `if`(not type(_passed[4], Vector), convert(_passed[4], Vector), _passed[4]);
    elif evalb(_npassed = 5) then
      _self:-m_x := `if`(not type(_passed[2], Vector), convert(_passed[2], Vector), _passed[2]);
      _self:-m_A := `if`(not type(_passed[3], Matrix), convert(_passed[3], Matrix), _passed[3]);
      _self:-m_b := `if`(not type(_passed[4], Vector), convert(_passed[4], Vector), _passed[4]);
      _self:-m_h := `if`(not type(_passed[5], Vector), convert(_passed[5], Vector), _passed[5]);
    else
      error("invalid input arguments.");
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
  end proc: # LoadMatrices

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export LoadEquations::static := proc(
    _self::Sandals,
    vars::{list, Vector},
    eqns::{list, Vector},
    invs::{list, Vector} := [],
    $)

    description "Load the equations from an explicit system x' = f(x,t) with "
      "and states <vars>, equations <eqns>, and optional invariants <invs>.";

    local x, A, b, h;

    x    := `if`(not type(vars, Vector), convert(vars, Vector), vars);
    A, b := LinearAlgebra:-GenerateMatrix(
      `if`(not type(eqns, list), convert(eqns, list), eqns),
      diff(`if`(not type(vars, list), convert(vars, list), vars), t)
    );
    h := `if`(not type(invs, list), convert(invs, list), invs);
    _self:-LoadMatrices(_self, x, A, b, h);
    return NULL;
  end proc: # LoadEquations

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TranslateToCpp::static := proc(
    _self::Sandals,
    name::string,
    {
    type::string                   := "Explicit",
    data::list(symbol = algebraic) := [],
    info::string                   := "No class description provided."
    }, $)::string;

    description "Generate C++ code for the loaded system with name <name>, Sandals class <type>, "
      "output file './<fname>.m', optional internal class data <data>, and class information string "
      "<info>.";

    local x, A, b, f, h;

    # Get system data
    if evalb(_self:-m_SystemType = 'Empty') then
      error("no system loaded yet.");
    end if;

    # Get system data
    x := _self:-m_x;
    A := _self:-m_A;
    b := _self:-m_b;
    h := _self:-m_h;

    # Generate class body string
    if evalb(type = "Explicit") then

      # Check if the system can be transformed into an explicit form
      if evalb(LinearAlgebra:-RowDimension(A) <> LinearAlgebra:-ColumnDimension(A) ) then
        error("the system cannot be transformed into an explicit form, the "
          "matrix A is not square.");
        return NULL;
      end if;

      if evalb(LinearAlgebra:-Rank(A) <> LinearAlgebra:-ColumnDimension(A)) then
        error("the system cannot be transformed into an explicit form, the "
          "matrix A is not full rank.");
        return NULL;
      end if;

      # Try to solve the system
      f := LinearAlgebra:-LinearSolve(A, b);
      try
        f := timelimit(_self:-m_TimeLimit, simplify(f));
      catch "time expired":
        WARNING("time expired, simplify(f) interrupted.");
      end try;

      # Generate class body string
      return SandalsCodegen:-ExplicitSystemToCpp(
        name, x, f, h, parse("data") = data, parse("info") = info
      );

    elif evalb(type = "Implicit") then

      # Generate class body string
      return SandalsCodegen:-ImplicitSystemToCpp(
        name, x, A, b, h, parse("data") = data, parse("info") = info
      );

    else
      error("unknown Sandals class type '%1'.", type);
    end if;
    return "";
  end proc: # TranslateToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateCppCode::static := proc(
    _self::Sandals,
    name::string,
    {
    type::string                   := "Explicit",
    path::string                   := "./",
    data::list(symbol = algebraic) := [],
    info::string                   := "No class description provided."
    }, $)

    description "Generate C++ code for the loaded system with name <name>, Sandals class <type>, "
    "output file '<path>/<name>.hh', optional internal class data <data>, and class information "
    "string <info>.";

    SandalsCodegen:-GenerateFile(
      cat(path, "/", name, ".hh"),
      _self:-TranslateToCpp(
        _self, name,
        parse("type") = type,
        parse("data") = data,
        parse("info") = info
      )
    );
    return NULL;
  end proc: # GenerateCppCode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # Sandals

# That's all folks!
