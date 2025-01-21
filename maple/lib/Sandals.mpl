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

  description """Sandals"" module.";

  option object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Miscellaneous options
  local m_verbose_mode  := false;
  local m_warning_mode  := true;
  local m_time_limit    := 0.1;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info::static := proc()

    description "Print ""Sandals"" module information.";

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

    description """Sandals"" module load procedure.";

    local i, lib_base_path;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("Sandals", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "cannot find ""Sandals"" module.";
    end if;
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload::static := proc()
    description """Sandals"" module unload procedure.";
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleCopy := proc(
    _self::Sandals,
    proto::Sandals,
    $)

    description "Copy the objects <proto> into <self>.";

    # Miscellaneous options
    _self:-m_verbose_mode := proto:-m_verbose_mode;
    _self:-m_warning_mode := proto:-m_warning_mode;
    _self:-m_time_limit   := proto:-m_time_limit;

    # Working directory
    _self:-CopyWorkingDirectory(_self, proto);

    # System
    _self:-CopySystem(_self, proto);

    # Code generation
    _self:-CopyCodeGeneration(_self, proto);
  end proc: # ModuleCopy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Reset := proc(
    _self::Sandals,
    label::{symbol, string} := NULL,
    $)

    description "Reset the library to the default.";

    # Miscellaneous options
    _self:-m_verbose_mode  := false;
    _self:-m_warning_mode  := true;
    _self:-m_time_limit    := 0.1;

    # Working directory
    _self:-ResetWorkingDirectory(_self);

    # System
    _self:-ResetSystem(_self);

    # Code generation
    _self:-ResetCodeGeneration(_self);
    return NULL;
  end proc: # Reset

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode := proc(
    _self::Sandals,
    $)

    description "Enable the verbose mode.";

    _self:-m_verbose_mode := true;
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode := proc(
    _self::Sandals,
    $)

    description "Disable the verbose mode.";

    _self:-m_verbose_mode := false;
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerboseMode := proc(
    _self::Sandals,
    mode::boolean,
    $)

    description "Set the verbosity of the module to <mode>.";

    _self:-m_verbose_mode := mode;
    return NULL;
  end proc: # SetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode := proc(
    _self::Sandals,
    $)

    description "Enable the warning mode of the module.";

    _self:-m_warning_mode := true;
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode := proc(
    _self::Sandals,
    $)

    description "Disable the warning mode of the module.";

    _self:-m_warning_mode := false;
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWarningMode := proc(
    _self::Sandals,
    mode::boolean,
    $)

    description "Set the warning mode of the module to <mode>.";

    _self:-m_warning_mode := mode;
    return NULL;
  end proc: # SetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTimeLimit := proc(
    _self::Sandals,
    x::numeric,
    $)

    description "Set the time limit of the module to <x>.";

    if (x < 0) then
      error("time limit must be a non-negative number.");
    end if;

    _self:-m_time_limit := x;
    return NULL;
  end proc: # SetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetTimeLimit := proc(
    _self::Sandals,
    $)::numeric;

    description "Get the time limit of the module.";

    return _self:-m_time_limit;
  end proc: # GetTimeLimit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/Sandals/WorkingDirectory.mpl"
$include "./lib/Sandals/Differentiation.mpl"
$include "./lib/Sandals/System.mpl"
$include "./lib/Sandals/CodeGeneration.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end module: # Sandals
