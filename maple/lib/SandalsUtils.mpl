# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#              ____                  _       _ _   _ _   _ _                  #
#             / ___|  __ _ _ __   __| | __ _| | | | | |_(_) |___              #
#             \___ \ / _` | '_ \ / _` |/ _` | | | | | __| | / __|             #
#              ___) | (_| | | | | (_| | (_| | | |_| | |_| | \__ \             #
#             |____/ \__,_|_| |_|\__,_|\__,_|_|\___/ \__|_|_|___/             #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors: Davide Stocco (University of Trento)
#          Enrico Bertolazzi (University of Trento)
#
# License: BSD 2-Clause License
#
# This is a module for the 'SandalsUtils' module. It contains the auxiliary functions to be used in
# the 'Sandals' module.

SandalsUtils := module()

  description "Utilities for the 'Sandals' module.";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ModuleLoad := proc()
    description "'SandalsUtils' module load procedure.";
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ModuleUnload := proc()
    description "'SandalsUtils' module unload procedure.";
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   ____  _  __  __                     _   _       _   _
  #  |  _ \(_)/ _|/ _| ___ _ __ ___ _ __ | |_(_) __ _| |_(_) ___  _ __
  #  | | | | | |_| |_ / _ \ '__/ _ \ '_ \| __| |/ _` | __| |/ _ \| '_ \
  #  | |_| | |  _|  _|  __/ | |  __/ | | | |_| | (_| | |_| | (_) | | | |
  #  |____/|_|_| |_|  \___|_|  \___|_| |_|\__|_|\__,_|\__|_|\___/|_| |_|
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DoDiff := proc(
    # _passed
    )::anything;

    description "Differentiate an expression with respect to a function.";

    local i, tmp1, tmp2, out;

    tmp2 := diff(subs(f = tmp1, out), tmp1);
    subs(tmp1 = f, convert(tmp2, D));
    out := _passed[1];
    for i from 2 to _npassed do
      tmp2 := diff(subs(_passed[i] = tmp1, out), tmp1);
      out  := subs(tmp1 = _passed[i], convert(tmp2, D));
    end do;
    return out;
  end proc: # DoDiff

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DoGradient := proc(
    fnc::algebraic,
    lst::{Vector, list},
    $)::Vector;

    description "Differentiate a scalar expression <fnc> with respect to a list of functions <lst>.";

    local i, n, out;

    # Extract dimensions
    if type(lst, Vector) then
      n := LinearAlgebra:-Dimension(lst);
    else
      n := nops(lst);
    end if;

    # Differentiate
    out := Vector(n);
    for i from 1 to n do
      out[i] := SandalsUtils:-DoDiff(fnc, lst[i]):
    end do;
    return out;
  end proc: # DoGradient

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DoHessian := proc(
    fnc::algebraic,
    lst::{Vector, list},
    $)::Matrix;

    description "Differentiate a vector of expressions (gradient) <fnc> with respect to a list <lst> "
      "of functions.";

    local i, j, n, out;

    # Extract dimensions
    if type(lst, Vector) then
      n := LinearAlgebra:-Dimension(lst);
    else
      n := nops(lst);
    end if;

    # Differentiate
    out := Matrix(n, n);
    for i from 1 to n do
      for j from i to n do
        out[i, j] := SandalsUtils:-DoDiff(fnc, lst[i], lst[j]);
        out[j, i] := out[i, j];
      end do;
    end do;
    return out;
  end proc: # DoHessian

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DoJacobian := proc(
    fnc::Vector,
    lst::{Vector, list},
    $)::Matrix;

    description "Differentiate a vector of expressions <fnc> with respect to a list <lst> of "
      "functions.";

    local i, j, m, n, out;

    # Extract dimensions
    m := LinearAlgebra:-Dimension(fnc);
    if type(lst, Vector) then
      n := LinearAlgebra:-Dimension(lst);
    else
      n := nops(lst);
    end if;

    # Differentiate
    out := Matrix(m, n);
    for i from 1 to m do
      for j from 1 to n do
        out[i, j] := SandalsUtils:-DoDiff(fnc[i], lst[j]);
      end do;
    end do;
    return out;
  end proc: # DoJacobian

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # SandalsUtils

# That's all folks!
