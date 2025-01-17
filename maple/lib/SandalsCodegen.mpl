# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                                        #
#                                                                                                 #
# The Sandals project is distributed under the BSD 2-Clause License.                              #
#                                                                                                 #
# Davide Stocco                                                                 Enrico Bertolazzi #
# University of Trento                                                       University of Trento #
# e-mail: davide.stocco@unitn.it                               e-mail: enrico.bertolazzi@unitn.it #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

SandalsCodegen := module()

  description "Code Generation for the 'SandalsCodegen' module.";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  local m_WorkingDirectory;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ModuleLoad := proc()
    description "'SandalsCodegen' module load procedure.";
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ModuleUnload := proc()
    description "'SandalsCodegen' module unload procedure.";
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetWorkingDirectory := proc(
    dname::string,
    $)

    description "Set the working directory to <dname> for the code generation.";

    m_WorkingDirectory := cat(currentdir(), "/", dname, "/");
    try
      mkdir(m_WorkingDirectory);
    catch:
      printf("Working directory already exists.\n");
    end try;
    printf("Working directory: %s\n", m_WorkingDirectory);
    return NULL;
  end proc: # SetWorkingDirectory

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GenerateFile := proc(
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

  export ClearFile := proc(
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

  export TranslateToCpp := proc(
    expr_list::list,
    $)

    description "Convert a list of expressions <expr_list> into C++ code.";

    local i, user_fun, out;

    # Define language functions
    user_fun := [];

    # Define new language
    CodeGeneration:-LanguageDefinition:-Define("NewC", extend = "C",
      seq(AddFunction(user_fun[i], anything::numeric, user_fun[i]), i = 1..nops(user_fun)));

    # Generate code
    return SandalsCodegen:-ApplyIndent("Real ", CodeGeneration:-Translate(
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

  local BoundCode := proc(
    code::string,
    {
    char::nonnegint := 120
    }, $)::string;

    description "Bound the code <code> with a character limit <char>.";

    return cat(StringTools:-LengthSplit(code, char) ||~ "\n");

  end proc: # BoundCode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ApplyIndent := proc(
    ind::string,
    str::string,
    $)::string;

    description "Apply indentation <ind> to string <str>.";

    local out;

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

  local GenerateProperties := proc(
    data::list(symbol),
    {
    propbool::list(boolean) := [seq(true, i = 1..nops(data))],
    indent::string          := "  "
    }, $)::string;

    description "Generate properties code from a list of data <data>, boolean flags <propbool>, "
      "and optional indentation <indent>.";

    local i, out;

    if (nops(propbool) <> nops(data)) then
      error("the number of properties' boolean must be equal to the number of properties.");
    end if;

    if has(propbool, true) then
      out := "";
      for i from 1 to nops(data) do
        if propbool[i] then
          convert(data[i], string);
          out := cat(out, indent, "Real const & ", %, " = this->m_", %, ";\n");
        end if;
      end do;
    else
      out := cat(indent, "// No properties\n");
    end if;
    return out;
  end proc: # GenerateProperties

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateInputs := proc(
    vars::list(list(symbol)),
    {
    varsbool::list(list(boolean)) := [seq([seq(true, j = 1..nops(vars[i]))], i = 1..nops(vars))],
    indent::string := "  ",
    label::string  := "in"
    }, $)::string;

    description "Generate inputs code from a list of variables <vars>, label <label>, and "
      "indentation <indent>.";

    local i, j, out;

    if (nops(varsbool) <> nops(vars)) then
      error("the number of variables' boolean must be equal to the number of variables.");
    end if;

    if has(varsbool, true) then
      out := "";
      for j from 1 to nops(vars) do
        for i from 1 to nops(vars[j]) do
          if varsbool[j][i] then
            convert(vars[j][i], string);
            out := cat(out, indent, "Real const & ", %, " = ", label, "_", j-1, "(", i-1, ");\n");
          end if;
        end do;
      end do;
    else
      out := cat(indent, "// No inputs' alias\n");
    end if;
    return out;
  end proc: # GenerateInputs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ExtractElements := proc(
    name::string,
    func::{list, Vector, Matrix, Array},
    dims::list(nonnegint),
    {
    label::string  := "out",
    indent::string := "  "
    }, $)::list, string;

    description "Extract elements for a n-dimensional function <func> with name <name>, dimensions "
      "<dims>, label <label>, and indentation <indent>.";

    local i, j, idx, lst, out, cur, tmp, str;

    lst := [];
    out := "";
    if has(evalb~(func =~ 0), false) and (mul(dims) > 0) then
      out := cat(out, indent, label, "_", name, " <<\n", indent, indent);
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
      out := cat(out, indent, label, "_", name, ".setZero();\n");
    end if;
    return lst, out;
  end proc: # ExtractElements

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateHeader := proc(
    name::string,
    vars::list(list(symbol)),
    varstype::list(string) := [seq("undefined", i = 1..nops(vars))],
    outtype::string        := "undefined",
    {
    varsbool::list(boolean) := [seq(true, i = 1..nops(vars))],
    timebool::boolean       := true,
    info::string            := "No description provided.",
    spec::string            := "",
    indent::string          := "  "
    }, $)::string;

    description "Generate a function header for a function <name> with variables <vars>, description "
      "<info>, variables types <varstype>, output type <outtype>, function specifier <spec>, "
      "variables boolean flag <varsbool>, time boolean flag <timebool>, and indentation <indent>.";

    local i, out;

    if (nops(varsbool) <> nops(vars)) then
      error("the number of variables' boolean must be equal to the number of variables.");
    elif (nops(varstype) <> nops(vars)) then
      error("the number of variables' types must be equal to the number of variables.");
    end if;

    out := cat(
      "// ", info, "\n",
      outtype, " ", name, "("
    );
    for i from 1 to nops(vars) do
      out := cat(out,
        varstype[i], " ", `if`(varsbool[i], cat("in_", i-1), cat("/*in_", i-1, "*/")), ", "
      );
    end do;
    return cat(out, "Real ", `if`(timebool, "t", "/*t*/"), ") ", spec, "\n{\n");
  end proc: # GenerateHeader

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateElements := proc(
    func::{list, Vector, Matrix, Array},
    {
    indent::string := "  "
    }, $)::string;

    description "Generate code for elements <func> with indentation <indent>.";

    SandalsCodegen:-TranslateToCpp(func);
    return SandalsCodegen:-ApplyIndent(indent,
      `if`(% = "", "// No body\n", cat(%%, %))
    );
  end proc: # GenerateElements

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateBody := proc(
    name::string,
    dims::list(nonnegint),
    {
    header::string     := "No header",
    properties::string := "No properties",
    inputs::string     := "No inputs",
    elements::string   := "No elements",
    outputs::string    := "No outputs",
    indent::string     := "  "
    }, $)::string;

    description "Generate code for function body for a function <name> with dimensions <dims>, "
      "header <header>, properties <properties>, inputs <inputs>, elements <elements>, output "
      "type <type>, and indentation <indent>.";

    local tmp;

    tmp := cat(seq(cat(dims[i], ", "), i = 1..nops(dims)-1), dims[-1]);
    if (nops(dims) <= 1) then
      tmp := cat(tmp, ", 1");
    end if;

    return cat(
      header,
      indent, "// Create properties' aliases\n", properties, "\n",
      indent, "// Create inputs' aliases\n",     inputs,     "\n",
      indent, "// Evaluate function\n",          elements,   "\n",
      indent, "// Store outputs\n",              outputs,   "}\n"
    );
  end proc: # GenerateBody

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local VectorToCpp := proc(
    name::string,
    vars::list(list(symbol)),
    vec::{list, Vector},
    varstype::list(string) := [seq("undefined", i = 1..nops(vars))],
    outtype::string        := "undefined",
    {
    data::list(symbol) := [],
    info::string       := "No info",
    label::string      := "out",
    spec::string       := "",
    indent::string     := "  "
    }, $)::string;

    description "Translate the vector <vec> with variables <vars> into a C++ function named <name> "
      "and return it as a string. The optional arguments are class properties <data>, function "
      "description <info>, output type <outtype>, label <label>, function specifier <spec>, and "
      "indentation string <indent>.";

    local header, properties, inputs, elements, outputs, dims, lst, i, j, propbool, varsbool;

    # Extract the function properties
    propbool := [seq(has(vec, data[i]), i = 1..nops(data))];
    properties := SandalsCodegen:-GenerateProperties(
      data, parse("indent") = indent, parse("propbool") = propbool
    );

    # Extract the function inputs
    varsbool := [seq([seq(has(vec, vars[i][j]), j = 1..nops(vars[i]))], i = 1..nops(vars))];
    inputs := SandalsCodegen:-GenerateInputs(
      vars, parse("indent") = indent, parse("varsbool") = varsbool
    );

    # Extract the function elements
    dims := [LinearAlgebra:-Dimension(vec)];
    lst, outputs := SandalsCodegen:-ExtractElements(
      name, vec, dims,
      parse("indent") = indent, parse("label") = label
    );
    outputs := cat(indent, outtype, " ", label, "_", name, ";\n", outputs);
    outputs := cat(outputs, indent, "return ", label, "_", name, ";\n");

    # Generate the method header
    header := SandalsCodegen:-GenerateHeader(
      name, vars, varstype, outtype,
      parse("info") = info, parse("indent") = indent, parse("spec") = spec,
      parse("timebool") = has(vec, t), parse("varsbool") = [seq(has(varsbool[i], true), i = 1..nops(vars))]
    );

    # Generate the elements
    elements := SandalsCodegen:-GenerateElements(lst);

    # Generate the generated code
    return SandalsCodegen:-GenerateBody(
      name, dims,
      parse("header")   = header,   parse("properties") = properties, parse("inputs")  = inputs,
      parse("elements") = elements, parse("indent")     = indent,     parse("outputs") = outputs
    );
  end proc: # VectorToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local MatrixToCpp := proc(
    name::string,
    vars::list(list(symbol)),
    mat::Matrix,
    varstype::list(string) := [seq("undefined", i = 1..nops(vars))],
    outtype::string        := "undefined",
    {
    data::list(symbol) := [],
    info::string       := "No info",
    label::string      := "out",
    spec::string       := "",
    indent::string     := "  "
    }, $)::string;

    description "Translate the matrix <mat> with variables <vars> into a C++ function named <name> "
      "and return it as a string. The optional arguments are class properties <data>, function "
      "description <info>, output type <outtype>, label <label>, function specifier <spec>, and "
      "indentation string <indent>.";

    local header, properties, inputs, elements, outputs, dims, lst, i, j, propbool, varsbool;

    # Extract the function properties
    propbool := [seq(has(mat, data[i]), i = 1..nops(data))];
    properties := SandalsCodegen:-GenerateProperties(
      data, parse("indent") = indent, parse("propbool") = propbool
    );

    # Extract the function inputs
    varsbool := [seq([seq(has(mat, vars[i][j]), j = 1..nops(vars[i]))], i = 1..nops(vars))];
    inputs := SandalsCodegen:-GenerateInputs(
      vars, parse("indent") = indent, parse("varsbool") = varsbool
    );

    # Extract the function elements
    dims := [LinearAlgebra:-Dimensions(mat)];
    lst, outputs := SandalsCodegen:-ExtractElements(
      name, mat, dims,
      parse("indent") = indent, parse("label") = label
    );
    outputs := cat(indent, outtype, " ", label, "_", name, ";\n", outputs);
    outputs := cat(outputs, indent, "return ", label, "_", name, ";\n");

    # Generate the method header
    header := SandalsCodegen:-GenerateHeader(
      name, vars, varstype, outtype,
      parse("info") = info, parse("indent") = indent, parse("spec") = spec,
      parse("timebool") = has(mat, t), parse("varsbool") = [seq(has(varsbool[i], true), i = 1..nops(vars))]
    );

    # Generate the elements
    elements := SandalsCodegen:-GenerateElements(lst);

    # Store the results
    return SandalsCodegen:-GenerateBody(
      name, dims,
      parse("header")   = header,   parse("properties") = properties, parse("inputs")  = inputs,
      parse("elements") = elements, parse("indent")     = indent,     parse("outputs") = outputs
    );
  end proc: # MatrixToCpp


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local TensorToCpp := proc(
    name::string,
    vars::list(list(symbol)),
    ten::Array,
    varstype::list(string) := [seq("undefined", i = 1..nops(vars))],
    outtype::string        := "undefined",
    {
    data::list(symbol) := [],
    info::string       := "No info",
    label::string      := "out",
    spec::string       := "",
    indent::string     := "  "
    }, $)::string;

    description "Translate the tensor <ten> with variables <vars> into a C++ function named <name> "
      "and return it as a string. The optional arguments are class properties <data>, function "
      "description <info>, output type <outtype>, label <label>, function specifier <spec>, and "
      "indentation string <indent>.";

    local dims, header, properties, inputs, elements, i, mat, outputs_tmp, lst_tmp, outputs, lst,
      j, propbool, varsbool;

    # Check the tensor dimensions
    dims := op~(2, [ArrayDims(ten)]);
    if (nops(dims) < 3) then
      error("the tensor must have at least 3 dimensions.");
    elif (nops(dims) > 3) then
      error("only 3-dimensional tensors are supported.");
    end if;

    # Extract the function properties
    propbool := [seq(has(ten, data[i]), i = 1..nops(data))];
    properties := SandalsCodegen:-GenerateProperties(
      data, parse("indent") = indent, parse("propbool") = propbool
    );

    # Extract the function inputs
    varsbool := [seq([seq(has(ten, vars[i][j]), j = 1..nops(vars[i]))], i = 1..nops(vars))];
    inputs := SandalsCodegen:-GenerateInputs(
      vars, parse("indent") = indent, parse("varsbool") = varsbool
    );

    # Extract the function elements
    lst := []; outputs := "";
    for i from 1 to dims[-1] do
      mat := convert(ten[1..-1, 1..-1, i], Matrix);
      lst_tmp, outputs_tmp := SandalsCodegen:-ExtractElements(
        cat(name, "[", i-1, "]"), mat, dims[1..-2],
        parse("indent") = indent, parse("label") = label
      );
      lst     := [op(lst), op(lst_tmp)];
      outputs := cat(outputs, outputs_tmp);
    end do;
    outputs := cat(indent, outtype, " ", label, "_", name, "(", dims[-1], ");\n", outputs);
    outputs := cat(outputs, indent, "return ", label, "_", name, ";\n");

    # Generate the method header
    header := SandalsCodegen:-GenerateHeader(
      name, vars, varstype, outtype,
      parse("info") = info, parse("indent") = indent, parse("spec") = spec,
      parse("timebool") = has(ten, t), parse("varsbool") = [seq(has(varsbool[i], true), i = 1..nops(vars))]
    );

    # Generate the elements
    elements := SandalsCodegen:-GenerateElements(lst);

    # Generate the generated code
    return SandalsCodegen:-GenerateBody(
      name, dims,
      parse("header")   = header,   parse("properties") = properties, parse("inputs")  = inputs,
      parse("elements") = elements, parse("indent")     = indent,     parse("outputs") = outputs
    );
  end proc: # TensorToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateConstructor := proc(
    name::string,
    type::string,
    {
    num_eqns::nonnegint := 0,
    num_invs::nonnegint := 0,
    data::list(symbol)  := [],
    info::string        := "Class constructor.",
    indent::string      := "  "
    }, $)::string;

    description "Generate a constructor for a system named <name> with system type <type>, "
      "number of functions <num_eqns>, number of invariants <num_invs>, system data <data> and "
      "description <info>.";

    return cat(
      "// ", info, "\n",
      name, "() : ", type, "<", num_eqns, ", ", num_invs, ">(""", name, """) {}\n",
      "\n",
      `if`(nops(data) > 0, cat(
        "// ", info, "\n",
        name, "(", StringTools:-Delete(cat~(op(cat~("Real ", convert~(data, string), ", "))), -2..-1), ")\n",
        indent, ": ", type, "<", num_eqns, ", ", num_invs, ">(""", name, """)",
        cat~(op(cat~(", m_", convert~(data, string), "(", convert~(data, string), ")"))), " {}\n"),
        "")
    );
  end proc: # GenerateConstructor

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   ___                 _ _      _ _
  #  |_ _|_ __ ___  _ __ | (_) ___(_) |_
  #   | || '_ ` _ \| '_ \| | |/ __| | __|
  #   | || | | | | | |_) | | | (__| | |_
  #  |___|_| |_| |_| .__/|_|_|\___|_|\__|
  #                |_|
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ImplicitToCpp := proc(
    name::string,
    x::{Vector(algebraic), list(algebraic)},
    F::{Vector(algebraic), list(algebraic)},
    h::{Vector(algebraic), list(algebraic)},
    {
    data::list(symbol = algebraic)               := [],
    time::{list(numeric), range(numeric)}        := [],
    ics::{Vector(algebraic), list(algebraic)}    := [],
    domain::{Vector(algebraic), list(algebraic)} := [],
    info::string                                 := "No class description provided.",
    indent::string                               := "  "
    }, $)::string;

    description "Generate an implicit system for the first-order differential equations F(x,x',t) "
      "= 0, with function <F>, invariants <h>, system data <data>, domain set <domain>, and "
      "description <info>. The indentation of the output code is set to <indent>.";

    local x_tmp, x_dot_tmp, mk_x_dot, rm_x_deps, rm_x_dot_deps, F_tmp, JF_x, JF_x_dot, h_tmp, Jh_x,
      rm_deps, i, data_str, properties, num_eqns, num_invs;

    # Store system states
    if not type(x, Vector) then
      x_tmp := convert(x, Vector);
    else
      x_tmp := x;
    end if;
    x_dot_tmp := map(i -> convert(cat(op(0, i), "_dot"), symbol)(op(1..-1, i)), x);

    # Prepare veriables for substitution
    # diff(x, t) -> x_dot(t)
    mk_x_dot := convert(diff(x_tmp, t) =~ x_dot_tmp, list);
    # x(t) -> x
    rm_x_deps := convert(x_tmp =~ op~(0, x_tmp), list);
    # x_dot(t) -> x_dot
    rm_x_dot_deps := convert(x_dot_tmp =~ op~(0, x_dot_tmp), list);

    # Compose substitutions
    rm_deps := [
      op(rm_x_deps),    # x(t) -> x
      op(rm_x_dot_deps) # x_dot(t) -> x_dot
    ];

    # Store system function and calculate Jacobians
    if not type(F, Vector) then
      F_tmp := convert(F, Vector);
    else
      F_tmp := F;
    end if;
    F_tmp    := subs(op(mk_x_dot), F_tmp);
    JF_x     := SandalsUtils:-DoJacobian(F_tmp, x_tmp);
    JF_x_dot := SandalsUtils:-DoJacobian(F_tmp, x_dot_tmp);

    # Store system invariants and calculate Jacobian
    if not type(h, Vector) then
      h_tmp := convert(h, Vector);
    else
      h_tmp := h;
    end if;
    Jh_x := SandalsUtils:-DoJacobian(h_tmp, x_tmp);

    # Generate expressions with proper variables dependencices
    x_tmp     := convert(subs(op(rm_deps), x_tmp), list);
    x_dot_tmp := convert(subs(op(rm_deps), x_dot_tmp), list);
    F_tmp     := subs(op(rm_deps), F_tmp);
    JF_x      := subs(op(rm_deps), JF_x);
    JF_x_dot  := subs(op(rm_deps), JF_x_dot);
    h_tmp     := subs(op(rm_deps), h_tmp);
    Jh_x      := subs(op(rm_deps), Jh_x);

    # Function utilities strings
    i := indent;

    # Generate properties
    if (nops(data) > 0) then
      properties := lhs~(data);
      data_str   := cat~(op(cat~(
        i, "Real m_", convert~(properties, string), "{", convert~(rhs~(data), string), "};\n"
      )));
    else
      data_str   := cat(i, "// No Properties\n");
      properties := [];
    end if;

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(F_tmp);
    num_invs := LinearAlgebra:-Dimension(h_tmp);

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
      "// ", info, "\n",
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
      SandalsCodegen:-ApplyIndent(i, GenerateConstructor(
          name, "Implicit",
          parse("num_eqns") = num_eqns,
          parse("num_invs") = num_invs,
          parse("data")     = properties,
          parse("info")     = "Class constructor.",
          parse("indent")   = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "F", [x_tmp, x_dot_tmp], F_tmp,
          ["VectorF const &", "VectorF const &"], "VectorF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the function F.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "JF_x", [x_tmp, x_dot_tmp], JF_x,
          ["VectorF const &", "VectorF const &"], "MatrixJF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the Jacobian of F with respect to x.",
          parse("indent")  = i
      )),
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "JF_x_dot", [x_tmp, x_dot_tmp], JF_x_dot,
          ["VectorF const &", "VectorF const &"], "MatrixJF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the Jacobian of F with respect to x_dot.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "h", [x_tmp], h_tmp,
          ["VectorF const &"], "VectorH",
          parse("data")    = properties,
          parse("info")    = "Calculate the vector h of the invariants.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jh_x", [x_tmp], Jh_x,
          ["VectorF const &"], "MatrixJH",
          parse("data")    = properties,
          parse("info")    = "Calculate the Jacobian of h with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      i, "bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}\n",
      "\n",
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
    name::string,
    x::{Vector(algebraic), list(algebraic)},
    f::{Vector(algebraic), list(algebraic)},
    h::{Vector(algebraic), list(algebraic)},
    {
    data::list(symbol = algebraic)               := [],
    time::{list(numeric), range(numeric)}        := [],
    ics::{Vector(algebraic), list(algebraic)}    := [],
    domain::{Vector(algebraic), list(algebraic)} := [],
    info::string                                 := "No class description provided.",
    indent::string                               := "  "
    }, $)::string;

    description "Generate an explicit system for the first-order differential equations x' = f(x,t), "
      "withstates variables <x>, right-hand side <f>, invariants <h>, system data <data>, domain "
      "set <domain>, and description <info>. The indentation of the output code is set to <indent>.";

    local x_tmp, b_tmp, Jb_x, h_tmp, Jh_x, rm_deps, i, data_str, properties, num_eqns, num_invs;

    # Store system states
    if not type(x, Vector) then
      x_tmp := convert(x, Vector);
    else
      x_tmp := x;
    end if;

    # Prepare veriables for substitution (x(t) -> x)
    rm_deps := convert(x_tmp =~ op~(0, x_tmp), list);

    # Store system right-hand side function and calculate Jacobian
    if not type(f, Vector) then
      b_tmp := convert(f, Vector);
    else
      b_tmp := f;
    end if;
    Jb_x := SandalsUtils:-DoJacobian(b_tmp, x_tmp);

    # Store system invariants and calculate Jacobian
    if not type(h, Vector) then
      h_tmp := convert(h, Vector);
    else
      h_tmp := h;
    end if;
    Jh_x := SandalsUtils:-DoJacobian(h_tmp, x_tmp);

    # Generate expressions with proper variables dependencices
    x_tmp := convert(subs(op(rm_deps), x_tmp), list);
    b_tmp := subs(op(rm_deps), b_tmp);
    Jb_x  := subs(op(rm_deps), Jb_x);
    h_tmp := subs(op(rm_deps), h_tmp);
    Jh_x  := subs(op(rm_deps), Jh_x);

    # Function utilities strings
    i := indent;

    # Generate properties
    if (nops(data) > 0) then
      properties := lhs~(data);
      data_str   := cat~(op(cat~(
        i, "Real m_", convert~(properties, string), "{", convert~(rhs~(data), string), "};\n"
      )));
    else
      data_str   := cat(i, "// No Properties\n");
      properties := [];
    end if;

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(b_tmp);
    num_invs := LinearAlgebra:-Dimension(h_tmp);

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
      "// ", info, "\n",
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
      SandalsCodegen:-ApplyIndent(i, GenerateConstructor(
          name, "Explicit",
          parse("num_eqns") = num_eqns,
          parse("num_invs") = num_invs,
          parse("data")     = properties,
          parse("info")     = "Class constructor.",
          parse("indent")   = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "f", [x_tmp], b_tmp,
          ["VectorF const &"], "VectorF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the function f.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jb_x", [x_tmp], Jb_x,
          ["VectorF const &"], "MatrixJF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the Jacobian of f with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "h", [x_tmp], h_tmp,
          ["VectorF const &"], "VectorH",
          parse("data")    = properties,
          parse("info")    = "Calculate the vector h of the invariants.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jh_x", [x_tmp], Jh_x,
          ["VectorF const &"], "MatrixJH",
          parse("data")    = properties,
          parse("info")    = "Calculate the Jacobian of h with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      i, "bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}\n",
      "\n",
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
    name::string,
    x::{Vector(algebraic), list(algebraic)},
    A::{Matrix(algebraic)},
    b::{Vector(algebraic)},
    h::{Vector(algebraic), list(algebraic)},
    {
    data::list(symbol = algebraic)               := [],
    time::{list(numeric), range(numeric)}        := [],
    ics::{Vector(algebraic), list(algebraic)}    := [],
    domain::{Vector(algebraic), list(algebraic)} := [],
    info::string                                 := "No class description provided.",
    indent::string                               := "  "
    }, $)::string;

    description "Generate a semi-explicit system for the first-order differential equations A(x,t)."
      "x' = b(x,t), with states variables <x>, mass matrix <A>, right-hand side <b>, invariants <h>, "
      "system data <data>, domain set <domain>, and description <info>. The indentation of the output "
      "code is set to <indent>.";

    local x_tmp, A_tmp, TA_x, b_tmp, Jb_x, h_tmp, Jh_x, rm_deps, i, data_str, properties, num_eqns,
      num_invs, domain_tmp, domain_fun, time_tmp, time_prp, time_fun;

    # Store system states
    if not type(x, Vector) then
      x_tmp := convert(x, Vector);
    else
      x_tmp := x;
    end if;

    # Prepare veriables for substitution (x(t) -> x)
    rm_deps := convert(x_tmp =~ op~(0, x_tmp), list);

    # Store system right-hand side function and calculate Jacobian
    if not type(A, Matrix) then
      A_tmp := convert(A, Matrix);
    else
      A_tmp := A;
    end if;
    if not type(b, Vector) then
      b_tmp := convert(b, Vector);
    else
      b_tmp := b;
    end if;
    TA_x := SandalsUtils:-DoTensor(A_tmp, x_tmp);
    Jb_x := SandalsUtils:-DoJacobian(b_tmp, x_tmp);

    # Store system invariants and calculate Jacobian
    if not type(h, Vector) then
      h_tmp := convert(h, Vector);
    else
      h_tmp := h;
    end if;
    Jh_x := SandalsUtils:-DoJacobian(h_tmp, x_tmp);

    # Generate expressions with proper variables dependencices
    x_tmp := convert(subs(op(rm_deps), x_tmp), list);
    A_tmp := subs(op(rm_deps), A_tmp);
    TA_x  := subs(op(rm_deps), TA_x);
    b_tmp := subs(op(rm_deps), b_tmp);
    Jb_x  := subs(op(rm_deps), Jb_x);
    h_tmp := subs(op(rm_deps), h_tmp);
    Jh_x  := subs(op(rm_deps), Jh_x);

    # Function utilities strings
    i := indent;

    # Generate properties
    if (nops(data) > 0) then
      properties := lhs~(data);
      data_str   := cat~(op(cat~(
        i, "Real m_", convert~(properties, string), "{", convert~(rhs~(data), string), "};\n"
      )));
    else
      data_str   := cat(i, "// No Properties\n");
      properties := [];
    end if;

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(b_tmp);
    num_invs := LinearAlgebra:-Dimension(h_tmp);

    # Initial conditions function code generation

    # Integration time function code generation
    if not type(time, list) then
      time_tmp := convert(time, list)
    else
      time_tmp := time;
    end if;

    if (nops(time_tmp) > 0) then
      time_prp := SandalsCodegen:-ApplyIndent(i, cat(
        "Real m_t_ini{" , time_tmp[1], "};\n"
        "Real m_t_end{" , time_tmp[2], "};\n"
      ));
      time_fun := SandalsCodegen:-ApplyIndent(i, cat(
        "// Return the initial integration time\n"
        "Real t_ini() const {return this->m_t_ini;}\n"
        "\n"
        "// Return the final integration time\n"
        "Real t_end() const {return this->m_t_end;}\n"
        "\n"
      ));
    else
      time_prp := cat(i, "// No integration time properties\n");
      time_fun := "";
    end if;

    # Domain function code generation
    if not type(domain, list) then
      domain_tmp := convert(domain, list);
    else
      domain_tmp := domain;
    end if;
    if (nops(domain_tmp) > 0) then
      error("Domain set not supported yet.");
    #  domain_fun := SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
    #    "in_domain", [x_tmp], domain_tmp,
    #    ["VectorF const &"], "bool",
    #    parse("data")    = properties,
    #    parse("info")    = "Check if the state vector x is in the domain."
    #    parse("indent")  = i
    #  ));
    else
      domain_fun := SandalsCodegen:-ApplyIndent(i, cat(
        "// Check if the state vector x is in the domain.\n",
        "bool in_domain(VectorF const &/*x*/, Real /*t*/) const override {return true;}\n",
        "\n"
      ));
    end if;

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
      "// ", info, "\n",
      "class ", name, " : public SemiExplicit<", num_eqns, ", ", num_invs, ">\n",
      "{\n",
      i, "// User data\n",
      data_str,
      "\n",
      i, "// Integration time\n",
      time_prp,
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
      SandalsCodegen:-ApplyIndent(i, GenerateConstructor(
          name, "SemiExplicit",
          parse("num_eqns") = num_eqns,
          parse("num_invs") = num_invs,
          parse("data")     = properties,
          parse("info")     = "Class constructor.",
          parse("indent")   = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "A", [x_tmp], A_tmp,
          ["VectorF const &"], "MatrixA",
          parse("data")    = properties,
          parse("spec")    = "const override",
          parse("info")    = "Evaluate the mass matrix A.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-TensorToCpp(
          "TA_x", [x_tmp], TA_x,
          ["VectorF const &"], "TensorTA",
          parse("data")    = properties,
          parse("spec")    = "const override",
          parse("info")    = "Evaluate the tensor of A with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "b", [x_tmp], b_tmp,
          ["VectorF const &"], "VectorB",
          parse("data")    = properties,
          parse("spec")    = "const override",
          parse("info")    = "Evaluate the function b.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jb_x", [x_tmp], Jb_x,
          ["VectorF const &"], "MatrixJB",
          parse("data")    = properties,
          parse("spec")    = "const override",
          parse("info")    = "Evaluate the Jacobian of b with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "h", [x_tmp], h_tmp,
          ["VectorF const &"], "VectorH",
          parse("data")    = properties,
          parse("spec")    = "const override",
          parse("info")    = "Calculate the vector h of the invariants.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jh_x", [x_tmp], Jh_x,
          ["VectorF const &"], "MatrixJH",
          parse("data")    = properties,
          parse("spec")    = "const override",
          parse("info")    = "Calculate the Jacobian of h with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      domain_fun,
      time_fun,
      "}; // class ", name, "\n",
      "\n",
      "#endif // SANDALS_MAPLE_", StringTools:-UpperCase(name), "_SEMIEXPLICIT_HH\n"
    );
  end proc: # SemiExplicitToCpp

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # SandalsCodegen
