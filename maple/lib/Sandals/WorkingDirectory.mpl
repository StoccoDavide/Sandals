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

# Working directory
local m_working_directory := "";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ResetWorkingDirectory := proc(
  _self::Sandals,
  $)

  description "Reset the working directory.";

  _self:-m_working_directory := "";
  return NULL;
end proc: # ResetWorkingDirectory

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CopyWorkingDirectory := proc(
  _self::Sandals,
  proto::Sandals,
  $)

  description "Copy the working directory from <proto> into <self>.";

  _self:-m_working_directory := proto:-m_working_directory;
  return NULL;
end proc: # CopyWorkingDirectory

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetWorkingDirectory := proc(
  _self::Sandals,
  $)

  description "Return the working directory.";

  return _self:-m_working_directory;
end proc: # GetWorkingDirectory

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SetWorkingDirectory := proc(
  _self::Sandals,
  dname::string,
  $)

  description "Set the working directory to <dname> for the code generation.";

  _self:-m_working_directory := cat(currentdir(), "/", dname, "/");
  try
    mkdir(_self:-m_working_directory);
  catch:
    if (_self:-m_verbose_mode) then
      printf("working directory already exists.\n");
    end if;
  end try;
  if (_self:-m_verbose_mode) then
    printf("working directory: %s\n", _self:-m_working_directory);
  end if;
  return NULL;
end proc: # SetWorkingDirectory

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
