# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

unprotect('SandalsCodegen');
SandalsCodegen := module()

  description "Code Generation for the ""SandalsCodegen"" module.";

  option object;

  local m_working_directory := "";
  local m_verbose_mode      := true;
  local m_warning_mode      := true;
  local m_indent            := "  ";
  local m_name              := [];
  local m_info              := "";
  local m_vars              := [];
  local m_vars_info         := [];
  local m_data              := [];
  local m_data_info         := [];
  local m_time              := [];
  local m_ics               := [];
  local m_domain            := [];
  local m_user_function     := [];
  local m_comp_sequence     := [];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ModuleLoad::static := proc()

    description """SandalsCodegen"" module load procedure.";

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

  local ModuleUnload::static := proc()
    description """SandalsCodegen"" module unload procedure.";
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleCopy::static := proc(
    _self::SandalsCodegen,
    proto::SandalsCodegen,
    $)

    description "Copy the objects <proto> into <self>.";

    _self:-m_working_directory := proto:-m_working_directory;
    _self:-m_warning_mode      := proto:-m_warning_mode;
  end proc: # ModuleCopy

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateFile::static := proc(
    fname::string,
    str::string,
    $)

    description "Generate a file named <fname> with the content <str>.";

    local file;

    file := FileTools:-Text:-Open(fname, create = true, overwrite = true);
    FileTools:-Text:-WriteString(file, str);
    FileTools:-Text:-Close(file);
    return NULL;
  end proc: # GenerateFile

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearFile::static := proc(
    fname::string,
    $)

    description "Procedure that creates an empty file <fname>.";

    local file;

    file := FileTools:-Text:-Open(fname, create = true, overwrite = true);
    FileTools:-Text:-WriteString(file, "\n");
    FileTools:-Text:-Close(file);
    return NULL;
  end proc: # ClearFile

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ApplyPrefix::static := proc(
    pre::string,
    str::string,
    $)::string;

    description "Apply prefix <pre> to each line of the string <str>.";

    local out;

    # Apply initial prefix
    out := cat(pre, str);

    # Apply prefix to each line
    out := StringTools:-SubstituteAll(out, "\n", cat("\n", pre));

    # Remove prefix from last line
    out := StringTools:-Take(out, length(out) - length(pre));

    # Remove prefix empty lines or lines with only spaces
    out := StringTools:-SubstituteAll(out, "\n" || pre || "\n", "\n\n");

    return out
  end proc: # ApplyPrefix

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetWorkingDirectory := proc(
    _self::SandalsCodegen,
    $)

    description "Return the working directory.";

    return _self:-m_working_directory;
  end proc: # GetWorkingDirectory

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWorkingDirectory := proc(
    _self::SandalsCodegen,
    dname::string,
    $)

    description "Set the working directory to <dname> for the code generation.";

    _self:-m_working_directory := cat(currentdir(), "/", dname, "/");
    try
      mkdir(_self:-m_working_directory);
    catch:
      printf("working directory already exists.\n");
    end try;
    printf("working directory: %s\n", _self:-m_working_directory);
    return NULL;
  end proc: # SetWorkingDirectory

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode::static := proc(
    _self::SandalsCodegen,
    $)

    description "Enable the verbose mode.";

    _self:-m_verbose_mode := true;
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode::static := proc(
    _self::SandalsCodegen,
    $)

    description "Disable the verbose mode.";

    _self:-m_verbose_mode := false;
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVerboseMode::static := proc(
    _self::SandalsCodegen,
    mode::boolean,
    $)

    description "Set the verbosity of the module to <mode>.";

    _self:-m_verbose_mode := mode;
    return NULL;
  end proc: # SetVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode::static := proc(
    _self::SandalsCodegen,
    $)

    description "Enable the warning mode of the module.";

    _self:-m_warning_mode := true;
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode::static := proc(
    _self::SandalsCodegen,
    $)

    description "Disable the warning mode of the module.";

    _self:-m_warning_mode := false;
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWarningMode::static := proc(
    _self::SandalsCodegen,
    mode::boolean,
    $)

    description "Set the warning mode of the module to <mode>.";

    _self:-m_warning_mode := mode;
    return NULL;
  end proc: # SetWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ApplyIndent := proc(
    _self::SandalsCodegen,
    str::string,
    $)::string;

    description "Apply indentation at the beginning of each line of the string <str>.";

    local out, ind;

    # Retrieve indentation
    ind := _self:-m_indent;

    # Apply initial indentation
    out := cat(ind, str);

    # Apply indentation to each line
    out := StringTools:-SubstituteAll(out, "\n", cat("\n", ind));

    # Remove indentation from last line
    out := StringTools:-Take(out, length(out) - length(ind));

    # Remove indentation empty lines or lines with only spaces
    out := StringTools:-SubstituteAll(out, "\n" || ind || "\n", "\n\n");

    return out
  end proc: # ApplyIndent

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetIndent := proc(
    _self::SandalsCodegen,
    $)

    description "Return the indentation string.";

    return _self:-m_indent;
  end proc: # GetIndent

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetIndent := proc(
    _self::SandalsCodegen,
    indent::string,
    $)

    description "Set the indentation string to <indent>.";

    _self:-m_indent := indent;
    return NULL;
  end proc: # SetIndent

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetInfo := proc(
    _self::SandalsCodegen,
    $)

    description "Return the information string.";

    return _self:-m_indent;
  end proc: # GetInfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetInfo := proc(
    _self::SandalsCodegen,
    info::string,
    $)

    description "Set the information string to <info>.";

    _self:-m_info := info;
    return NULL;
  end proc: # SetInfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetVars := proc(
    _self::SandalsCodegen,
    $)::list(function);

    description "Return the variables list.";

    return _self:-m_vars;
  end proc: # GetVars

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVars := proc(
    _self::SandalsCodegen,
    vars::list(function),
    $)

    description "Set the variables list to <vars>.";

    _self:-m_vars      := vars;
    _self:-m_vars_info := _self:-SetVarsInfo(_self,
      [seq(cat("State variable ", i), i = 1..nops(_self:-m_vars))]);
    return NULL;
  end proc: # SetVars

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetVarsInfo := proc(
    _self::SandalsCodegen,
    $)::list(string);

    description "Return the variables information list.";

    return _self:-m_vars_info;
  end proc: # GetVarsInfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetVarsInfo := proc(
    _self::SandalsCodegen,
    vars_info::list(string),
    $)

    description "Set the variables information list to <vars_info>.";

    if not (nops(_self:-m_vars) > 0) and not (nops(vars_info) <> nops(_self:-m_vars)) then
      error("the number of variables' information must be equal to the number of variables.");
    end if;

    _self:-m_vars_info := vars_info;
    return NULL;
  end proc: # SetVarsInfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetData := proc(
    _self::SandalsCodegen,
    $)::list(symbol = algebraic);

    description "Return the data variables list.";

    return _self:-m_data;
  end proc: # GetData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetData := proc(
    _self::SandalsCodegen,
    data::list(symbol = algebraic),
    $)

    description "Set the data variables list to <data>.";

    _self:-m_data      := data;
    _self:-m_data_info := _self:-SetDataInfo(_self,
      [seq(cat("Data variable ", i), i = 1..nops(_self:-m_data))]);
    return NULL;
  end proc: # SetData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetDataInfo := proc(
    _self::SandalsCodegen,
    $)::list(string);

    description "Return the data variables information list.";

    return _self:-m_data_info;
  end proc: # GetDataInfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetDataInfo := proc(
    _self::SandalsCodegen,
    data_info::list(string),
    $)

    description "Set the data variables information list to <data_info>.";

    if not (nops(_self:-m_data) > 0) and not (nops(data_info) <> nops(_self:-m_data)) then
      error("the number of data' information must be equal to the number of data variables.");
    end if;

    _self:-m_data_info := data_info;
    return NULL;
  end proc: # SetDataInfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetTime := proc(
    _self::SandalsCodegen,
    $)::list(numeric);

    description "Return the integration time.";

    return _self:-m_time;
  end proc: # GetTime

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetTime := proc(
    _self::SandalsCodegen,
    time::{list(numeric), range(numeric)},
    $)

    description "Set the integration time to <time>.";

    local time_lst;

    if not type(time, list) then
      time_lst := convert(time, list)
    else
      time_lst := time;
    end if;

    if (nops(time_lst) > 0) then
      if (time_lst[1] = time_lst[2]) then
        error("degenerate time range detected.");
      elif (time_lst[2] < time_lst[1]) then
        if _self:-m_warning_mode then
          WARNING("decreasing time range detected, consider using the reverse time mode.");
        end if;
      end if;
    end if;

    _self:-m_time := time;
    return NULL;
  end proc: # SetTime

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetIcs := proc(
    _self::SandalsCodegen,
    $)::list(algebraic);

    description "Return the integration initial conditions.";

    return _self:-m_ics;
  end proc: # GetIcs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetIcs := proc(
    _self::SandalsCodegen,
    ics::list(algebraic),
    $)

    description "Set the integration initial conditions to <ics>.";

    local tmp;

    tmp := indets(ics) minus convert(lhs~(_self:-m_data), set);

    if (nops(ics) > 0) then
      if (nops(ics) <> nops(_self:-m_vars)) then
        error("invalid initial conditions detected.");
      elif (nops(tmp) > 0) then
        if _self:-m_warning_mode then
          WARNING("initial conditions contain the unknown variables {%1}.", tmp);
        end if;
      end if;
    end if;

    _self:-m_ics := ics;
    return NULL;
  end proc: # SetIcs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetUserFunction := proc(
    _self::SandalsCodegen,
    user_function::list(string),
    $)

    description "Set the user function to <user_function>.";

    _self:-m_user_function := user_function;
    return NULL;
  end proc: # SetUserFunction

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetUserFunction := proc(
    _self::SandalsCodegen,
    $)::list(string);

    description "Return the user function.";

    return _self:-m_user_function;
  end proc: # GetUserFunction

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetCompSequence := proc(
    _self::SandalsCodegen,
    comp_sequence::list(symbol = algebraic),
    $)

    description "Set the computation sequence to <comp_sequence>.";

    _self:-m_comp_sequence := comp_sequence;
    return NULL;
  end proc: # SetCompSequence

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetCompSequence := proc(
    _self::SandalsCodegen,
    $)::list(symbol = algebraic);

    description "Return the computation sequence.";

    return _self:-m_comp_sequence;
  end proc: # GetCompSequence

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetDomain := proc(
    _self::SandalsCodegen,
    domain::list({`<`, `<=`}),
    $)

    description "Set the domain to <domain>.";

    _self:-m_domain := domain;
    return NULL;
  end proc: # SetDomain

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetDomain := proc(
    _self::SandalsCodegen,
    $)::list({`<`, `<=`});

    description "Return the domain.";

    return _self:-m_domain;
  end proc: # GetDomain

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TranslateToCpp := proc(
    _self::SandalsCodegen,
    expr_list::list,
    $)

    description "Convert a list of expressions <expr_list> into C++ code.";

    local i, fun;

    # Define language functions
    fun := ["floor", "ceil", "round", "trunc", "erf", op(_self:-m_user_function)];

    # Define new language
    CodeGeneration:-LanguageDefinition:-Define("NewC", extend = "C",
      seq(AddFunction(fun[i], anything::numeric, fun[i]), i = 1..nops(fun)));

    # Generate code
    return _self:-ApplyPrefix("Real ", CodeGeneration:-Translate(
      expr_list,
      language          = "NewC",
      optimize          = true,
      digits            = 30,
      deducereturn      = false,
      coercetypes       = false,
      deducetypes       = false,
      reduceanalysis    = false,
      defaulttype       = numeric,
      functionprecision = double,
      output            = string
    ));
  end proc:

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateData := proc(
    _self::SandalsCodegen,
    data_bool::list(boolean),
    $)::string;

    description "Generate data extraction code, boolean flags <data_bool> are used to select the "
      "data variables to extract.";

    local data, i, out;

    data := lhs~(_self:-m_data);

    if has(data_bool, true) then
      out := "";
      for i from 1 to nops(data) do
        if data_bool[i] then
          convert(data[i], string);
          out := cat(out, _self:-m_indent, "Real const & ", %, " = this->m_", %, ";\n");
        end if;
      end do;
    else
      out := cat(_self:-m_indent, "// No data\n");
    end if;
    return out;
  end proc: # GenerateData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateInputs := proc(
    _self::SandalsCodegen,
    vars::list(list(symbol)),
    vars_bool::list(list(boolean)),
    vars_info::list(list(string)),
    {
    label::string := "in"
    }, $)::string;

    description "Generate inputs extraction code, boolean flags <vars_bool> are used to select the "
      "variables to extract. Variables information <vars_info> are used to provide a description of "
      "them. A label <label> is used to prefix the variables.";

    local i, j, out;

    if (nops(vars_bool) <> nops(vars)) then
      error("the number of variables' boolean must be equal to the number of variables.");
    end if;

    if has(vars_bool, true) then
      out := "";
      for j from 1 to nops(vars) do
        if (nops(vars_bool[j]) <> nops(vars[j])) then
          error("the number of variables' boolean must be equal to the number of variables.");
        end if;
        if (nops(vars[j]) > 1) then
          for i from 1 to nops(vars[j]) do
            if vars_bool[j][i] then
              convert(vars[j][i], string);
              out := cat(out, _self:-m_indent, "Real const & ", %, " = ", label, "_", j-1, "(", i-1,
                "); // ", vars_info[j][i], "\n");
            end if;
          end do;
        elif (nops(vars[j]) = 1) and vars_bool[j][1] then
          convert(vars[j][1], string);
          out := cat(out, _self:-m_indent, "Real const & ", %, " = ", label, "_", j-1, "; // ",
            vars_info[j][1], "\n");
        end if;
      end do;
    else
      out := cat(_self:-m_indent, "// No inputs' alias\n");
    end if;
    return out;
  end proc: # GenerateInputs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExtractElements := proc(
    _self::SandalsCodegen,
    name::string,
    func::{list, Vector, Matrix, Array},
    dims::list(nonnegint),
    {
    label::string  := "out"
    }, $)::list, string;

    description "Extract elements for a n-dimensional function <func> with name <name>, and "
      "dimensions <dims>. A label <label> is used to prefix the variables.";

    local i, j, idx, lst, out, cur, tmp, str;

    lst := [];
    out := "";
    if has(evalb~(func =~ 0), false) and (mul(dims) > 0) then
      out := cat(out, _self:-m_indent, label, "_", name, " <<\n", _self:-m_indent, _self:-m_indent);
      for i from 1 to mul(dims) do
        cur := [];
        idx := i-1;
        for j from 1 to nops(dims) do
          cur := [irem(idx, dims[j], 'idx')+1, op(cur)];
        end do;
        tmp := func[op(ListTools:-Reverse(cur))];
        map(x -> (x, "_"), cur -~ 1); str := cat(label, "_", op(%))[1..-2];
        lst := [op(lst), convert(str, symbol) = tmp];
        out := cat(out, str, ", ");
      end do;
      out := cat(StringTools:-Delete(out, -2..-1), ";\n");
    else
      out := cat(out, _self:-m_indent, label, "_", name, ".setZero();\n");
    end if;
    return lst, out;
  end proc: # ExtractElements

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateHeader := proc(
    _self::SandalsCodegen,
    name::string,
    vars::list(list(symbol)),
    vars_type::list(string),
    vars_bool::list(boolean),
    out_type::string,
    {
    info::string := "No description provided.",
    spec::string := ""
    }, $)::string;

    description "Generate a function header for a function <name> with variables <vars>, variables "
      "types <vars_type>, output type <out_type>, and variables boolean flag <vars_bool>. The "
      "optional arguments are function description <info>, and function specifier <spec>.";

    local i, out;

    if (nops(vars_bool) <> nops(vars)) then
      error("the number of variables' boolean must be equal to the number of variables.");
    elif (nops(vars_type) <> nops(vars)) then
      error("the number of variables' types must be equal to the number of variables.");
    end if;

    out := cat(
      "// ", info, "\n",
      out_type, " ", name, "("
    );
    for i from 1 to nops(vars) do
      out := cat(out,
        vars_type[i], " ", `if`(vars_bool[i], cat("in_", i-1), cat("/*in_", i-1, "*/")), ", "
      );
    end do;
    return cat(`if`(nops(vars) > 0, StringTools:-Delete(out, -2..-1), out), ") ", spec, "\n{\n");
  end proc: # GenerateHeader

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateElements := proc(
    _self::SandalsCodegen,
    func::{list, Vector, Matrix, Array},
    $)::string;

    description "Generate code for elements <func>.";

    _self:-TranslateToCpp(_self, func);
    return _self:-ApplyIndent(_self,
      `if`(% = "", "// No body\n", cat(%%, %))
    );
  end proc: # GenerateElements

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateBody := proc(
    _self::SandalsCodegen,
    name::string,
    dims::list(nonnegint),
    {
    header::string   := "No header",
    data::string     := "No data",
    inputs::string   := "No inputs",
    elements::string := "No elements",
    outputs::string  := "No outputs"
    }, $)::string;

    description "Generate code for function body for a function <name> with dimensions <dims>, "
      "header <header>, data <data>, inputs <inputs>, elements <elements>, and outputs "
      "<outputs>.";

    local tmp;

    tmp := cat(seq(cat(dims[i], ", "), i = 1..nops(dims)-1), dims[-1]);
    if (nops(dims) <= 1) then
      tmp := cat(tmp, ", 1");
    end if;

    return cat(
      header,
      _self:-m_indent, "// Create data' aliases\n", data, "\n",
      _self:-m_indent, "// Create inputs' aliases\n", inputs, "\n",
      _self:-m_indent, "// Evaluate function\n", elements, "\n",
      _self:-m_indent, "// Store outputs\n", outputs, "}\n"
    );
  end proc: # GenerateBody

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export VectorToCpp := proc(
    _self::SandalsCodegen,
    name::string,
    vars::list(list(symbol)),
    vec::{list, Vector},
    vars_type::list(string),
    vars_info::list(list(string)),
    out_type::string,
    {
    info::string := "No info",
    spec::string := ""
    }, $)::string;

    description "Translate the vector <vec> with variables <vars> into a C++ function named <name> "
      "and return it as a string. Variable types and return type are specified with <vars_type> and "
      "<out_type>. The optional arguments are function description <info>, and function specifier "
      "<spec>.";

    local vec_tmp, data_bool, vars_bool, i, j, data, header, inputs, elements, outputs, dims, lst,
      comp_sequence;

    if (nops(vars) <> nops(vars_type)) then
      error("the number of variables' types must be equal to the number of variables.");
    end if;

    if not type(vec, Vector) then
      vec_tmp := convert(vec, Vector)
    else
      vec_tmp := vec;
    end if;

    # Create boolean flags for the variables
    data_bool := [seq(has(vec_tmp, _self:-m_data[i]), i = 1..nops(_self:-m_data))];
    vars_bool := [seq([seq(has(vec_tmp, vars[i][j]), j = 1..nops(vars[i]))], i = 1..nops(vars))];

    # Extract the function data
    data := _self:-GenerateData(_self, data_bool);

    # Extract the function inputs
    inputs := _self:-GenerateInputs(_self, vars, vars_bool, vars_info);

    # Extract the function elements
    dims := [LinearAlgebra:-Dimension(vec_tmp)];
    lst, outputs := _self:-ExtractElements(_self, name, vec_tmp, dims);
    outputs := cat(_self:-m_indent, out_type, " ", label, "_", name, ";\n", outputs);
    outputs := cat(outputs, _self:-m_indent, "return ", label, "_", name, ";\n");

    # Generate the method header
    vars_bool := [seq(has(vars_bool[i], true), i = 1..nops(vars))];
    header := _self:-GenerateHeader(
      _self, name, vars, vars_type, vars_bool, out_type, parse("info") = info, parse("spec") = spec
    );

    # Generate the evaluation sequence
    # comp_sequence := _self:-GenerateElements(
    #   _self, selectremove[flatten](x -> has(vec, x), _self:-m_comp_sequence)
    # );

    # Generate the elements
    elements := _self:-GenerateElements(_self, lst);

    # Generate the generated code
    return _self:-GenerateBody(
      _self, name, dims, parse("header") = header, parse("data") = data, parse("inputs") = inputs,
      (*parse("comp_sequence") = comp_sequence,*) parse("elements") = elements, parse("outputs") = outputs
    );
  end proc: # VectorToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MatrixToCpp := proc(
    _self::SandalsCodegen,
    name::string,
    vars::list(list(symbol)),
    mat::Matrix,
    vars_type::list(string),
    vars_info::list(list(string)),
    out_type::string,
    {
    info::string := "No info",
    spec::string := ""
    }, $)::string;

    description "Translate the matrix <mat> with variables <vars> into a C++ function named <name> "
      "and return it as a string. Variable types and return type are specified with <vars_type> and "
      "<out_type>. The optional arguments are function description <info>, and function specifier "
      "<spec>.";

    local data_bool, vars_bool, i, j, data, header, inputs, elements, outputs, dims, lst,
      comp_sequence;

    if (nops(vars) <> nops(vars_type)) then
      error("the number of variables' types must be equal to the number of variables.");
    end if;

    # Create boolean flags for the variables
    data_bool := [seq(has(mat, _self:-m_data[i]), i = 1..nops(_self:-m_data))];
    vars_bool := [seq([seq(has(mat, vars[i][j]), j = 1..nops(vars[i]))], i = 1..nops(vars))];

    # Extract the function data
    data := _self:-GenerateData(_self, data_bool);

    # Extract the function inputs
    inputs := _self:-GenerateInputs(_self, vars, vars_bool, vars_info);

    # Extract the function elements
    dims := [LinearAlgebra:-Dimensions(mat)];
    lst, outputs := _self:-ExtractElements(_self, name, mat, dims);
    outputs := cat(_self:-m_indent, out_type, " ", label, "_", name, ";\n", outputs);
    outputs := cat(outputs, _self:-m_indent, "return ", label, "_", name, ";\n");

    # Generate the method header
    vars_bool := [seq(has(vars_bool[i], true), i = 1..nops(vars))];
    header := _self:-GenerateHeader(
      _self, name, vars, vars_type, vars_bool, out_type, parse("info") = info, parse("spec") = spec
    );

    # Generate the evaluation sequence
    # comp_sequence := _self:-GenerateElements(
    #   _self, selectremove[flatten](x -> has(mat, x), _self:-m_comp_sequence)
    # );

    # Generate the elements
    elements := _self:-GenerateElements(_self, lst);

    # Generate the generated code
    return _self:-GenerateBody(
      _self, name, dims, parse("header") = header, parse("data") = data, parse("inputs") = inputs,
      (*parse("comp_sequence") = comp_sequence,*) parse("elements") = elements, parse("outputs") = outputs
    );
  end proc: # MatrixToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TensorToCpp := proc(
    _self::SandalsCodegen,
    name::string,
    vars::list(list(symbol)),
    ten::Array,
    vars_type::list(string),
    vars_info::list(list(string)),
    out_type::string,
    {
    info::string := "No info",
    spec::string := ""
    }, $)::string;

    description "Translate the tensor <ten> with variables <vars> into a C++ function named <name> "
      "and return it as a string. Variable types and return type are specified with <vars_type> and "
      "<out_type>. The optional arguments are function description <info>, and function specifier "
      "<spec>.";

    local data_bool, vars_bool, i, j, data, header, inputs, elements, outputs_tmp, lst_tmp, outputs,
      mat, lst, dims, comp_sequence;

    if (nops(vars) <> nops(vars_type)) then
      error("the number of variables' types must be equal to the number of variables.");
    end if;

    # Check the tensor dimensions
    dims := op~(2, [ArrayDims(ten)]);
    if (nops(dims) < 3) then
      error("the tensor must have at least 3 dimensions.");
    elif (nops(dims) > 3) then
      error("only 3-dimensional tensors are supported.");
    end if;

    # Create boolean flags for the variables
    data_bool := [seq(has(ten, _self:-m_data[i]), i = 1..nops(_self:-m_data))];
    vars_bool := [seq([seq(has(ten, vars[i][j]), j = 1..nops(vars[i]))], i = 1..nops(vars))];

    # Extract the function data
    data := _self:-GenerateData(_self, data_bool);

    # Extract the function inputs
    inputs := _self:-GenerateInputs(_self, vars, vars_bool, vars_info);

    # Extract the function elements
    lst := []; outputs := "";
    for i from 1 to dims[-1] do
      mat := convert(ten[1..-1, 1..-1, i], Matrix);
      lst_tmp, outputs_tmp := _self:-ExtractElements(
        _self, cat(name, "[", i-1, "]"), mat, dims[1..-2]
      );
      lst     := [op(lst), op(lst_tmp)];
      outputs := cat(outputs, outputs_tmp);
    end do;
    outputs := cat(_self:-m_indent, out_type, " ", label, "_", name, "(", dims[-1], ");\n", outputs);
    outputs := cat(outputs, _self:-m_indent, "return ", label, "_", name, ";\n");

    # Generate the method header
    vars_bool := [seq(has(vars_bool[i], true), i = 1..nops(vars))];
    header := _self:-GenerateHeader(
      _self, name, vars, vars_type, vars_bool, out_type, parse("info") = info, parse("spec") = spec
    );

    # Generate the evaluation sequence
    # comp_sequence := _self:-GenerateElements(
    #   _self, selectremove[flatten](x -> has(ten, x), _self:-m_comp_sequence)
    # );

    # Generate the elements
    elements := _self:-GenerateElements(_self, lst);

    # Generate the generated code
    return _self:-GenerateBody(
      _self, name, dims, parse("header") = header, parse("data") = data, parse("inputs") = inputs,
      (*parse("comp_sequence") = comp_sequence,*) parse("elements") = elements, parse("outputs") = outputs
    );
  end proc: # TensorToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateConstructor := proc(
    _self::SandalsCodegen,
    name::string,
    type::string,
    {
    num_eqns::nonnegint := 0,
    num_invs::nonnegint := 0,
    info::string        := "Class constructor."
    }, $)::string;

    description "Generate a constructor for a system named <name> with system type <type>, number "
      "of functions <num_eqns>, number of invariants <num_invs>, and description <info>.";

    return cat(
      "// ", info, "\n",
      name, "() : ", type, "<", num_eqns, ", ", num_invs, ">(""", name, """) {}\n",
      "\n",
      `if`(nops(_self:-m_data) > 0, cat(
        "// ", info, "\n",
        name, "(", StringTools:-Delete(cat~(op(cat~("Real ", convert~(_self:-m_data, string), ", "))), -2..-1), ")\n",
        _self:-m_indent, ": ", type, "<", num_eqns, ", ", num_invs, ">(""", name, """)",
        cat~(op(cat~(", m_", convert~(lhs~(_self:-m_data), string), "(", convert~(rhs~(_self:-m_data), string), ")"))), " {}\n"),
        "")
    );
  end proc: # GenerateConstructor

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateIcs := proc(
    _self::SandalsCodegen,
    $)::string;

    description "Generate initial conditions function.";

    if (nops(_self:-m_ics) > 0) then
      return _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "ics", [], _self:-m_ics,
        [], [], "VectorF",
        parse("info") = "Evaluate the initial conditions.",
        parse("spec") = "const"
      ));
    else
      return "";
    end if;
  end proc: # GenerateIcs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateTime := proc(
    _self::SandalsCodegen,
    $)

    description "Generate integration time data and functions.";

    if (nops(_self:-m_time) > 0) then
      return _self:-ApplyIndent(_self, cat(
        "Real m_t_ini{" , _self:-m_time[1], "};\n"
        "Real m_t_end{" , _self:-m_time[2], "};\n"
      )), _self:-ApplyIndent(_self, cat(
        "// Return the initial integration time\n"
        "Real t_ini() const {return this->m_t_ini;}\n"
        "\n"
        "// Return the final integration time\n"
        "Real t_end() const {return this->m_t_end;}\n"
        "\n"
      ));
    else
      return cat(indent, "// No integration time data\n"), "";
    end if;
    return time_data, time_fun;
  end proc: # GenerateTime

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateDomain := proc(
    _self::SandalsCodegen,
    $)::string;

    description "Generate domain function.";

    if (nops(_self:-m_domain) > 0) then
      error("domain set not implemented yet.");
      return _self:-ApplyIndent(_self, _self:-VectorToCpp(
        "in_domain", [_self:-m_vars, [t]], _self:-m_domain,
        ["VectorF const &", "Real"], [_self:-m_vars_info, _self:-m_vars_info ||~ " derivative", "Time"],
        "VectorF",
        parse("spec") = "const",
        parse("info") = "Check if the system is in the defined domain."
      ));
    else
      return _self:-ApplyIndent(_self, cat(
        "// Check if the system is in the defined domain.\n",
        "bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}\n",
        "\n"
      ));
    end if;
  end proc: # GenerateDomain

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   ___                 _ _      _ _
  #  |_ _|_ __ ___  _ __ | (_) ___(_) |_
  #   | || '_ ` _ \| '_ \| | |/ __| | __|
  #   | || | | | | | |_) | | | (__| | |_
  #  |___|_| |_| |_| .__/|_|_|\___|_|\__|
  #                |_|
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ImplicitToCpp := proc(
    _self::SandalsCodegen,
    name::string,
    F::Vector(algebraic),
    h::Vector(algebraic),
    $)::string;

    description "Generate an implicit system for the first-order differential equations F(x,x',t) "
      "= 0, with function <F>, and invariants <h>.";

    local x, x_dot, mk_x_dot, rm_x_deps, rm_x_dot_deps, F_tmp, JF_x, JF_x_dot, h_tmp, Jh_x, rm_deps,
      i, data_str, data_lhs, num_eqns, num_invs, time_data, time_fun, ics_fun, domain_fun, vars_info_dot,
      vars_info;

    # Store system states
    x     := _self:-m_vars;
    x_dot := map(i -> convert(cat(op(0, i), "_dot"), symbol)(op(1..-1, i)), x);

    # Prepare veriables for substitution
    # diff(x, t) -> x_dot(t)
    mk_x_dot := convert(diff(x, t) =~ x_dot, list);
    # x(t) -> x
    rm_x_deps := convert(x =~ op~(0, x), list);
    # x_dot(t) -> x_dot
    rm_x_dot_deps := convert(x_dot =~ op~(0, x_dot), list);

    # Compose substitutions
    rm_deps := [
      op(rm_x_deps),    # x(t) -> x
      op(rm_x_dot_deps) # x_dot(t) -> x_dot
    ];

    # Calculate Jacobians
    F_tmp    := subs(op(mk_x_dot), F);
    JF_x     := SandalsUtils:-DoJacobian(F, x);
    JF_x_dot := SandalsUtils:-DoJacobian(F, x_dot);
    h_tmp    := subs(op(mk_x_dot), h);
    Jh_x     := SandalsUtils:-DoJacobian(h, x);

    # Generate expressions with proper variables dependencices
    x         := convert(subs(op(rm_deps), x), list);
    x_dot     := convert(subs(op(rm_deps), x_dot), list);
    F_tmp     := subs(op(rm_deps), F_tmp);
    JF_x      := subs(op(rm_deps), JF_x);
    JF_x_dot  := subs(op(rm_deps), JF_x_dot);
    h_tmp     := subs(op(rm_deps), h_tmp);
    Jh_x      := subs(op(rm_deps), Jh_x);

    # Function utilities strings
    i := _self:-m_indent;

    # Generate internal class data
    if (nops(_self:-m_data) > 0) then
      data_lhs := lhs~(_self:-m_data);
      data_str := cat~(op(cat~(
        i, "Real m_", convert~(data_lhs, string), "{", convert~(rhs~(_self:-m_data), string), "}; // ",
        _self:-m_data_info, "\n"
      )));
    else
      data_str := cat(i, "// No Properties\n");
      data_lhs := [];
    end if;

    # Initial conditions function code generation
    ics_fun := _self:-GenerateIcs(_self);

    # Integration time function code generation
    time_data, time_fun := _self:-GenerateTime(_self);

    # Domain function code generation
    domain_fun := _self:-GenerateDomain(_self);

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(F_tmp);
    num_invs := LinearAlgebra:-Dimension(h_tmp);

    # Variables information
    vars_info_dot := [_self:-m_vars_info, _self:-m_vars_info ||~ " derivative", ["Time"]];
    vars_info     := [_self:-m_vars_info, ["Time"]];

    # Return output string
    return cat(
      "/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\\\n"
      " * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *\n"
      " *                                                                                               *\n"
      " * The Sandals project is distributed under the BSD 2-Clause License.                            *\n"
      " *                                                                                               *\n"
      " * Davide Stocco                                                               Enrico Bertolazzi *\n"
      " * University of Trento                                                     University of Trento *\n"
      " * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *\n"
      "\\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */\n"
      "\n",
      "// C++ generated code for implicit system: ", name, "\n",
      "// This file has been automatically generated by Sandals.\n",
      "// DISCLAIMER: If you need to edit it, do it wisely!\n",
      "\n",
      "#ifndef SANDALS_MAPLE_", StringTools:-UpperCase(name), "_IMPLICIT_HH\n",
      "#define SANDALS_MAPLE_", StringTools:-UpperCase(name), "_IMPLICIT_HH\n",
      "\n",
      "using namespace Sandals;\n"
      "using namespace std;\n"
      "\n",
      "// ", _self:-m_info, "\n",
      "class ", name, " : public Implicit<", num_eqns, ", ", num_invs, ">\n",
      "{\n",
      i, "// User data\n",
      data_str,
      "\n",
      "public:\n",
      i, "using VectorF  = typename Implicit<", num_eqns, ", ", num_invs, ">::VectorF;\n",
      i, "using MatrixJF = typename Implicit<", num_eqns, ", ", num_invs, ">::MatrixJF;\n",
      i, "using VectorH  = typename Implicit<", num_eqns, ", ", num_invs, ">::VectorH;\n",
      i, "using MatrixJH = typename Implicit<", num_eqns, ", ", num_invs, ">::MatrixJH;\n",
      "\n",
      _self:-ApplyIndent(_self, _self:-GenerateConstructor(
        _self, name, "Implicit",
        parse("num_eqns") = num_eqns,
        parse("num_invs") = num_invs,
        parse("info")     = "Class constructor."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "F", [x, x_dot, [t]], F_tmp,
        ["VectorF const &", "VectorF const &", "Real"], vars_info_dot, "VectorF",
        parse("spec") = "const override",
        parse("info") = "Evaluate the function F."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "JF_x", [x, x_dot, [t]], JF_x,
        ["VectorF const &", "VectorF const &", "Real"], vars_info_dot, "MatrixJF",
        parse("spec") = "const override",
        parse("info") = "Evaluate the Jacobian of F with respect to x."
      )),
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "JF_x_dot", [x, x_dot, [t]], JF_x_dot,
        ["VectorF const &", "VectorF const &", "Real"], vars_info_dot, "MatrixJF",
        parse("spec") = "const override",
        parse("info") = "Evaluate the Jacobian of F with respect to x_dot."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "h", [x, [t]], h_tmp,
        ["VectorF const &", "Real"], vars_info,
        "VectorH",
        parse("spec") = "const override",
        parse("info") = "Calculate the vector h of the invariants."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "Jh_x", [x, [t]], Jh_x,
        ["VectorF const &", "Real"], vars_info,
        "MatrixJH",
        parse("spec") = "const override",
        parse("info") = "Calculate the Jacobian of h with respect to x."
      )),
      "\n",
      time_fun,
      ics_fun,
      domain_fun,
      "}; // class ", name, "\n",
      "\n",
      "#endif // SANDALS_MAPLE_", StringTools:-UpperCase(name), "_IMPLICIT_HH\n"
    );
  end proc: # ImplicitToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   _____            _ _      _ _
  #  | ____|_  ___ __ | (_) ___(_) |_
  #  |  _| \ \/ / '_ \| | |/ __| | __|
  #  | |___ >  <| |_) | | | (__| | |_
  #  |_____/_/\_\ .__/|_|_|\___|_|\__|
  #             |_|
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExplicitToCpp := proc(
    _self::SandalsCodegen,
    name::string,
    f::Vector(algebraic),
    h::Vector(algebraic),
    $)::string;

    description "Generate an explicit system for the first-order differential equations x' = f(x,t), "
      "withstates variables <x>, right-hand side <f>, and invariants <h>.";

    local x, f_tmp, Jf_x, h_tmp, Jh_x, rm_deps, i, data_str, data_lhs, num_eqns, num_invs,
      time_data, time_fun, ics_fun, domain_fun, vars_info;

    # Store system states
    x := _self:-m_vars;

    # Prepare veriables for substitution (x(t) -> x)
    rm_deps := convert(x =~ op~(0, x), list);

    # Generate expressions with proper variables dependencices
    x := convert(subs(op(rm_deps), x), list);
    f_tmp := subs(op(rm_deps), f);
    Jf_x  := subs(op(rm_deps), SandalsUtils:-DoJacobian(f, x));
    h_tmp := subs(op(rm_deps), h);
    Jh_x  := subs(op(rm_deps), SandalsUtils:-DoJacobian(h, x));

    # Function utilities strings
    i := _self:-m_indent;

    # Generate internal class data
    if (nops(_self:-m_data) > 0) then
      data_lhs := lhs~(_self:-m_data);
      data_str := cat~(op(cat~(
        i, "Real m_", convert~(data_lhs, string), "{", convert~(rhs~(_self:-m_data), string), "}; // ",
        _self:-m_data_info, "\n"
      )));
    else
      data_str   := cat(i, "// No Properties\n");
      data_lhs := [];
    end if;

    # Initial conditions function code generation
    ics_fun := _self:-GenerateIcs(_self);

    # Integration time function code generation
    time_data, time_fun := _self:-GenerateTime(_self);

    # Domain function code generation
    domain_fun := _self:-GenerateDomain(_self);

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(f_tmp);
    num_invs := LinearAlgebra:-Dimension(h_tmp);

    # Variables information
    vars_info := [_self:-m_vars_info, ["Time"]];

    # Return output string
    return cat(
      "/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\\\n"
      " * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *\n"
      " *                                                                                               *\n"
      " * The Sandals project is distributed under the BSD 2-Clause License.                            *\n"
      " *                                                                                               *\n"
      " * Davide Stocco                                                               Enrico Bertolazzi *\n"
      " * University of Trento                                                     University of Trento *\n"
      " * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *\n"
      "\\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */\n"
      "\n",
      "// C++ generated code for explicit system: ", name, "\n",
      "// This file has been automatically generated by Sandals.\n",
      "// DISCLAIMER: If you need to edit it, do it wisely!\n",
      "\n",
      "#ifndef SANDALS_MAPLE_", StringTools:-UpperCase(name), "_EXPLICIT_HH\n",
      "#define SANDALS_MAPLE_", StringTools:-UpperCase(name), "_EXPLICIT_HH\n",
      "\n",
      "using namespace Sandals;\n"
      "using namespace std;\n"
      "\n",
      "// ", _self:-m_info, "\n",
      "class ", name, " : public Explicit<", num_eqns, ", ", num_invs, ">\n",
      "{\n",
      i, "// User data\n",
      data_str,
      "\n",
      "public:\n",
      i, "using VectorF  = typename Explicit<", num_eqns, ", ", num_invs, ">::VectorF;\n",
      i, "using MatrixJF = typename Explicit<", num_eqns, ", ", num_invs, ">::MatrixJF;\n",
      i, "using VectorH  = typename Explicit<", num_eqns, ", ", num_invs, ">::VectorH;\n",
      i, "using MatrixJH = typename Explicit<", num_eqns, ", ", num_invs, ">::MatrixJH;\n",
      "\n",
      _self:-ApplyIndent(_self, _self:-GenerateConstructor(
        _self, name, "Explicit",
        parse("num_eqns") = num_eqns,
        parse("num_invs") = num_invs,
        parse("info")     = "Class constructor."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "f", [x, [t]], f_tmp,
        ["VectorF const &", "Real"], vars_info, "VectorF",
        parse("spec") = "const override",
        parse("info") = "Evaluate the function f."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "Jf_x", [x, [t]], Jf_x,
        ["VectorF const &", "Real"], vars_info, "MatrixJF",
        parse("spec") = "const override",
        parse("info") = "Evaluate the Jacobian of f with respect to x."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "h", [x, [t]], h_tmp,
        ["VectorF const &", "Real"], vars_info, "VectorH",
        parse("spec") = "const override",
        parse("info") = "Calculate the vector h of the invariants."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "Jh_x", [x, [t]], Jh_x,
        ["VectorF const &", "Real"], vars_info, "MatrixJH",
        parse("spec") = "const override",
        parse("info") = "Calculate the Jacobian of h with respect to x."
      )),
      "\n",
      time_fun,
      ics_fun,
      domain_fun,
      "}; // class ", name, "\n",
      "\n",
      "#endif // SANDALS_MAPLE_", StringTools:-UpperCase(name), "_EXPLICIT_HH\n"
    );
  end proc: # ExplicitToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   ____                 _ _____            _ _      _ _
  #  / ___|  ___ _ __ ___ (_) ____|_  ___ __ | (_) ___(_) |_
  #  \___ \ / _ \ '_ ` _ \| |  _| \ \/ / '_ \| | |/ __| | __|
  #   ___) |  __/ | | | | | | |___ >  <| |_) | | | (__| | |_
  #  |____/ \___|_| |_| |_|_|_____/_/\_\ .__/|_|_|\___|_|\__|
  #                                    |_|
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SemiExplicitToCpp := proc(
    _self::SandalsCodegen,
    name::string,
    A::Matrix(algebraic),
    b::Vector(algebraic),
    h::Vector(algebraic),
    $)::string;

    description "Generate a semi-explicit system for the first-order differential equations A(x,t)."
      "x' = b(x,t), with states variables <x>, mass matrix <A>, right-hand side <b>, and invariants "
      "<h>.";

    local x, A_tmp, TA_x, b_tmp, Jb_x, h_tmp, Jh_x, rm_deps, i, data_str, data_lhs, num_eqns,
      num_invs, time_data, time_fun, ics_fun, domain_fun, vars_info;

    # Store system states
    x := _self:-m_vars;

    # Prepare veriables for substitution (x(t) -> x)
    rm_deps := convert(x =~ op~(0, x), list);

    # Store system right-hand side function and calculate Jacobian
    TA_x := SandalsUtils:-DoTensor(A, x);
    Jb_x := SandalsUtils:-DoJacobian(b, x);

    # Store system invariants and calculate Jacobian
    Jh_x := SandalsUtils:-DoJacobian(h, x);

    # Generate expressions with proper variables dependencices
    x     := convert(subs(op(rm_deps), x), list);
    A_tmp := subs(op(rm_deps), A);
    TA_x  := subs(op(rm_deps), TA_x);
    b_tmp := subs(op(rm_deps), b);
    Jb_x  := subs(op(rm_deps), Jb_x);
    h_tmp := subs(op(rm_deps), h);
    Jh_x  := subs(op(rm_deps), Jh_x);

    # Function utilities strings
    i := _self:-m_indent;

    # Generate internal class data
    if (nops(_self:-m_data) > 0) then
      data_lhs := lhs~(_self:-m_data);
      data_str := cat~(op(cat~(
        i, "Real m_", convert~(data_lhs, string), "{", convert~(rhs~(_self:-m_data), string), "}; // ",
        _self:-m_data_info, "\n"
      )));
    else
      data_str   := cat(i, "// No Properties\n");
      data_lhs := [];
    end if;

    # Initial conditions function code generation
    ics_fun := _self:-GenerateIcs(_self);

    # Integration time function code generation
    time_data, time_fun := _self:-GenerateTime(_self);

    # Domain function code generation
    domain_fun := _self:-GenerateDomain(_self);

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(b);
    num_invs := LinearAlgebra:-Dimension(h);

    # Variables information
    vars_info := [_self:-m_vars_info, ["Time"]];

    # Return output string
    return cat(
      "/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\\\n"
      " * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                      *\n"
      " *                                                                                               *\n"
      " * The Sandals project is distributed under the BSD 2-Clause License.                            *\n"
      " *                                                                                               *\n"
      " * Davide Stocco                                                               Enrico Bertolazzi *\n"
      " * University of Trento                                                     University of Trento *\n"
      " * e-mail: davide.stocco@unitn.it                             e-mail: enrico.bertolazzi@unitn.it *\n"
      "\\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */\n"
      "\n",
      "// C++ generated code for semi-explicit system: ", name, "\n",
      "// This file has been automatically generated by Sandals.\n",
      "// DISCLAIMER: If you need to edit it, do it wisely!\n",
      "\n",
      "#ifndef SANDALS_MAPLE_", StringTools:-UpperCase(name), "_SEMIEXPLICIT_HH\n",
      "#define SANDALS_MAPLE_", StringTools:-UpperCase(name), "_SEMIEXPLICIT_HH\n",
      "\n",
      "using namespace Sandals;\n"
      "using namespace std;\n"
      "\n",
      "// ", _self:-m_info, "\n",
      "class ", name, " : public SemiExplicit<", num_eqns, ", ", num_invs, ">\n",
      "{\n",
      i, "// User data\n",
      data_str,
      "\n",
      i, "// Integration time\n",
      time_data,
      "\n",
      "public:\n",
      i, "using VectorF  = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::VectorF;\n",
      i, "using MatrixA  = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::MatrixA;\n",
      i, "using TensorTA = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::TensorTA;\n",
      i, "using VectorB  = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::VectorB;\n",
      i, "using MatrixJB = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::MatrixJB;\n",
      i, "using VectorH  = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::VectorH;\n",
      i, "using MatrixJH = typename SemiExplicit<", num_eqns, ", ", num_invs, ">::MatrixJH;\n",
      "\n",
      _self:-ApplyIndent(_self, _self:-GenerateConstructor(
        _self, name, "SemiExplicit",
        parse("num_eqns") = num_eqns,
        parse("num_invs") = num_invs,
        parse("info")     = "Class constructor."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "A", [x, [t]], A_tmp,
        ["VectorF const &", "Real"], vars_info, "MatrixA",
        parse("spec") = "const override",
        parse("info") = "Evaluate the mass matrix A."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-TensorToCpp(
        _self, "TA_x", [x, [t]], TA_x,
        ["VectorF const &", "Real"], vars_info, "TensorTA",
        parse("spec") = "const override",
        parse("info") = "Evaluate the tensor of A with respect to x."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "b", [x, [t]], b_tmp,
        ["VectorF const &", "Real"], vars_info, "VectorB",
        parse("spec") = "const override",
        parse("info") = "Evaluate the function b."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "Jb_x", [x, [t]], Jb_x,
        ["VectorF const &", "Real"], vars_info, "MatrixJB",
        parse("spec") = "const override",
        parse("info") = "Evaluate the Jacobian of b with respect to x."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-VectorToCpp(
        _self, "h", [x, [t]], h_tmp,
        ["VectorF const &", "Real"], vars_info, "VectorH",
        parse("spec") = "const override",
        parse("info") = "Calculate the vector h of the invariants."
      )),
      "\n",
      _self:-ApplyIndent(_self, _self:-MatrixToCpp(
        _self, "Jh_x", [x, [t]], Jh_x,
        ["VectorF const &", "Real"], vars_info, "MatrixJH",
        parse("spec") = "const override",
        parse("info") = "Calculate the Jacobian of h with respect to x."
      )),
      "\n",
      time_fun,
      ics_fun,
      domain_fun,
      "}; // class ", name, "\n",
      "\n",
      "#endif // SANDALS_MAPLE_", StringTools:-UpperCase(name), "_SEMIEXPLICIT_HH\n"
    );
  end proc: # SemiExplicitToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # SandalsCodegen
