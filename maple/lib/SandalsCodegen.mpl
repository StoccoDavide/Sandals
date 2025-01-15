# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#          ____                  _       _       ____          _              #
#         / ___|  __ _ _ __   __| | __ _| |___ / ___|___   __| | ___          #
#         \___ \ / _` | '_ \ / _` |/ _` | / __| |   / _ \ / _` |/ _ \         #
#          ___) | (_| | | | | (_| | (_| | \__ \ |__| (_) | (_| |  __/         #
#         |____/ \__,_|_| |_|\__,_|\__,_|_|___/\____\___/ \__,_|\___|         #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors: Davide Stocco (University of Trento)
#          Enrico Bertolazzi (University of Trento)
#
# License: BSD 2-Clause License
#
# This is a module for the 'SandalsCodegen' module. It contains the functions to translate the Maple
# code into C++ code.

SandalsCodegen := module()

  description "Code Generation for the 'SandalsCodegen' module.";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  local m_CodegenOptions;
  local m_WorkingDirectory;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ModuleLoad := proc()
    description "'SandalsCodegen' module load procedure.";
    m_CodegenOptions := table([
      optimize          = true,
      digits            = 30,
      deducereturn      = false,
      coercetypes       = false,
      deducetypes       = false,
      reduceanalysis    = true,
      defaulttype       = numeric,
      functionprecision = double
    ]);
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

  export SetCodegenOptions := proc({
      optimize::boolean         := m_CodegenOptions[parse("optimize")],
      digits::posint            := m_CodegenOptions[parse("digits")],
      deducereturn::boolean     := m_CodegenOptions[parse("deducereturn")],
      coercetypes::boolean      := m_CodegenOptions[parse("coercetypes")],
      deducetypes::boolean      := m_CodegenOptions[parse("deducetypes")],
      reduceanalysis::boolean   := m_CodegenOptions[parse("reduceanalysis")],
      defaulttype::symbol       := m_CodegenOptions[parse("defaulttype")],
      functionprecision::symbol := m_CodegenOptions[parse("functionprecision")]
    }, $)

    description "Set options for code generation optimization.";

    m_CodegenOptions[parse("optimize")]          := optimize;
    m_CodegenOptions[parse("digits")]            := digits;
    m_CodegenOptions[parse("deducereturn")]      := deducereturn;
    m_CodegenOptions[parse("coercetypes")]       := coercetypes;
    m_CodegenOptions[parse("deducetypes")]       := deducetypes;
    m_CodegenOptions[parse("reduceanalysis")]    := reduceanalysis;
    m_CodegenOptions[parse("defaulttype")]       := defaulttype;
    m_CodegenOptions[parse("functionprecision")] := functionprecision;
    return NULL;
  end proc: # SetCodegenOptions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetCodegenOptions := proc(
    field::string := "all",
    $)

    description "Get options for CodeGeneration.";

    if evalb(field = "all") then
      return m_CodegenOptions;
    else
      return m_CodegenOptions[parse(field)];
    end if;
  end proc: # GetCodegenOptions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TranslateToCpp := proc(
    expr_list::list,
    $)

    description "Convert a list of expressions <expr_list> into C++ code.";

    local i, lang_fncs;

    # Define language functions
    lang_fncs := [];

    # Define new language
    CodeGeneration:-LanguageDefinition:-Define("NewC", extend = "C"
      #seq(AddFunction(lang_fncs[i], anything::numeric, lang_fncs[i]), i = 1..nops(lang_fncs)),
      #AddType(
      # 'integer' = table(['single' = "integer", 'double' = "long"]),
      # 'numeric' = table(['single' = "Real",    'double' = "Real"])
      #)
    );

    # Generate code
    return CodeGeneration:-Translate(
      expr_list,
      language = "C",
      #op(op(eval(m_CodegenOptions))),
      output = string
    );
  end proc:

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
    indent::string := "  "
    }, $)::string;

    description "Generate properties code from a list of data <data> and optional indentation <indent>.";

    local i, out;

    if evalb(nops(data) > 0) then
      out := "";
      for i from 1 to nops(data) do
        convert(data[i], string);
        out := cat(out, indent, "using ", %, " = this->m_", %, ";\n");
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
    indent::string := "  ",
    label::string  := "in"
    }, $)::string;

    description "Generate inputs code from a list of variables <vars>, label <label>, and "
      "indentation <indent>.";

    local i, j, out;

    out := "";
    for j from 1 to nops(vars) do
      for i from 1 to nops(vars[j]) do
        convert(vars[j][i], string);
        out := cat(out, indent, "using ", %, " = ", label, "_", j-1, "(", i-1, ");\n");
      end do;
    end do;
    return out;
  end proc: # GenerateInputs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local ExtractElements := proc(
    name::string,
    func::{list, Vector, Matrix, Array},
    dims::list(nonnegint),
    varstype::list(string) := [seq("undefined", i = 1..nops(vars))],
    outtype::string        := "undefined",
    {
    label::string  := "out",
    indent::string := "  "
    }, $)::list, string;

    description "Extract elements for a n-dimensional function <func> with name <name>, dimensions "
      "<dims>, output type <outtype>, label <label>, and indentation <indent>.";

    local i, j, idx, lst, out, cur, tmp, str;

    lst := [];
    out := cat(indent, outtype, " ", label, "_", name, ";\n");
    if has(false, evalb~(func = 0)) and evalb(mul(dims) > 0) then
      out := cat(out, indent, label, "_", name, " <<\n", indent, indent);
      for i from 1 to mul(dims) do
        cur := [];
        idx := i-1;
        for j from 1 to nops(dims) do
          cur := [irem(idx, dims[j], 'idx')+1, op(cur)];
        end do;
        tmp := func[op(cur)];
        map(x -> (x, "_"), cur -~ 1); str := cat(label, "_", op(%))[1..-2];
        lst := [op(lst), convert(str, symbol) = tmp];
        out := cat(out, str, ", ");
      end do;
      out := cat(StringTools:-Delete(out, -2..-1), ";\n");
    else
      out := cat(out, indent, label, "_", name, ".setZero();\n");
    end if;
    out := cat(out, indent, "return ", label, "_", name, ";\n");
    return lst, out;
  end proc: # ExtractElements

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateHeader := proc(
    name::string,
    vars::list(list(symbol)),
    varstype::list(string) := [seq("undefined", i = 1..nops(vars))],
    outtype::string        := "undefined",
    {
    info::string    := "No description provided.",
    indent::string  := "  "
    }, $)::string;

    description "Generate a function header for a function <name> with variables <vars>, "
      "description <info>, variable type <varstype>, output type <outtype>, and indentation <indent>.";

    local i, out;

    out := cat(
      "// ", info, "\n",
      outtype, " ", name, "("
    );
    for i from 1 to nops(vars) do
      out := cat(out, varstype[i], "in_", i, ", ");
    end do;
    return cat(out, "Real t) const override\n{\n");
    return out;
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
    if evalb(nops(dims) <= 1) then
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
    indent::string     := "  "
    }, $)::string;

    description "Translate the vector <vec> with variables <vars> into a C++ function named <name> "
      "and return it as a string. The optional arguments are class properties <data>, function "
      "description <info>, output type <outtype>, label <label>, and indentation string <indent>.";

    local header, properties, inputs, elements, outputs, dims, lst;

    # Extract the function properties
    properties := SandalsCodegen:-GenerateProperties(
      data, parse("indent") = indent
    );

    # Extract the function inputs
    inputs := SandalsCodegen:-GenerateInputs(
      vars, parse("indent") = indent
    );

    # Extract the function elements
    dims := [LinearAlgebra:-Dimension(vec)];
    lst, outputs := SandalsCodegen:-ExtractElements(
      name, vec, dims, varstype, outtype,
      parse("indent") = indent, parse("label") = label
      );

    # Generate the method header
    header := SandalsCodegen:-GenerateHeader(
      name, vars, varstype, outtype,
      parse("info") = info, parse("indent") = indent
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
    indent::string     := "  "
    }, $)::string;

    description "Translate the vector <vec> with variables <vars> into a C++ function named <name> "
    "and return it as a string. The optional arguments are class properties <data>, function "
    "description <info>, output type <outtype>, label <label>, and indentation string <indent>.";

    local header, properties, inputs, elements, outputs, dims, lst;

    # Extract the function properties
    properties := SandalsCodegen:-GenerateProperties(
      data, parse("indent") = indent
    );

    # Extract the function inputs
    inputs := SandalsCodegen:-GenerateInputs(
      vars, parse("indent") = indent
    );

    # Extract the function elements
    dims := [LinearAlgebra:-Dimensions(mat)];
    lst, outputs := SandalsCodegen:-ExtractElements(
      name, mat, dims, varstype, outtype,
      parse("indent") = indent, parse("label") = label
    );

    # Generate the method header
    header := SandalsCodegen:-GenerateHeader(
      name, vars, varstype, outtype,
      parse("info") = info, parse("indent") = indent
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  local GenerateConstructor := proc(
    name::string,
    sys_type::string,
    {
    num_eqns::nonnegint := 0,
    num_invs::nonnegint := 0,
    data::list(symbol)  := [],
    info::string        := "Class constructor.",
    indent::string      := "  "
    }, $)::string;

    description "Generate a constructor for a system named <name> with system type <sys_type>, "
      "number of functions <num_eqns>, number of invariants <num_invs>, system data <data> and "
      "description <info>.";

    return cat(
      "// ", info, "\n",
      name, "() : Explicit<", num_eqns, ", ", num_invs, ">(""", name, """) {}\n",
      "\n",
      `if`(nops(data) > 0, cat(
        "// ", info, "\n",
        name, "(", StringTools:-Delete(cat~(op(cat~("Real ", convert~(data, string), ", "))), -2..-1), ")\n",
        indent, ": Explicit<", num_eqns, ", ", num_invs, ">(""", name, """)",
        cat~(op(cat~(", m_", convert~(data, string), "(", convert~(data, string), ")"))), " {}\n"),
        "")
    );

  end proc: # GenerateConstructor

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #   _____            _ _      _ _
  #  | ____|_  ___ __ | (_) ___(_) |_
  #  |  _| \ \/ / '_ \| | |/ __| | __|
  #  | |___ >  <| |_) | | | (__| | |_
  #  |_____/_/\_\ .__/|_|_|\___|_|\__|
  #             |_|
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ExplicitSystemToCpp := proc(
    name::string,
    x::{Vector(algebraic), list(algebraic)},
    f::{Vector(algebraic), list(algebraic)},
    h::{Vector(algebraic), list(algebraic)},
    {
    data::list(symbol = algebraic)   := [],
    domain::list(symbol = algebraic) := [],
    info::string                     := "No class description provided.",
    indent::string                   := "  "
    }, $)::string;

    description "Generate an explicit system for the first-order differential equations <f> with "
      "states variables <x>, invariants <h>, system data <data>, domain set <domain>, and description "
      "<info>.";

    local x_tmp, f_tmp, Jf_x, h_tmp, Jh_x, rm_deps, i, bar, data_str, properties, num_eqns, num_invs;

    # Store system states
    if not type(x, Vector) then
      x_tmp := convert(x, Vector);
    else
      x_tmp := x;
    end if;

    # Prepare veriables for substitution (x(t) -> x)
    rm_deps := convert(x_tmp =~ op~(0, x_tmp), list);

    # Store system function and calculate Jacobian
    if not type(f, Vector) then
      f_tmp := convert(f, Vector);
    else
      f_tmp := f;
    end if;
    Jf_x := SandalsUtils:-DoJacobian(f_tmp, x_tmp);

    # Store system invariants and calculate Jacobian
    if not type(h, Vector) then
      h_tmp := convert(h, Vector);
    else
      h_tmp := h;
    end if;
    Jh_x := SandalsUtils:-DoJacobian(h_tmp, x_tmp);

    # Generate expressions with proper variables dependencices
    x_tmp := convert(subs(op(rm_deps), x_tmp), list);
    f_tmp := subs(op(rm_deps), f_tmp);
    Jf_x  := subs(op(rm_deps), Jf_x);
    h_tmp := subs(op(rm_deps), h_tmp);
    Jh_x  := subs(op(rm_deps), Jh_x);

    # Function utilities strings
    i   := indent;
    bar := "// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";

    # Generate properties
    if evalb(nops(data) > 0) then
      properties := lhs~(data);
      data_str   := cat~(op(cat~(i, "Real m_", convert~(properties, string), "{", convert~(rhs~(data), string), "};\n")));
    else
      data_str   := cat(i, "// No Properties\n");
      properties := [];
    end if;

    # Compute template arguments
    num_eqns := LinearAlgebra:-Dimension(f_tmp);
    num_invs := LinearAlgebra:-Dimension(h_tmp);

    # Return output string
    return cat(
      "/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\\\n"
      " * Copyright (c) 2025, Davide Stocco and Enrico Bertolazzi.                  *\n"
      " *                                                                           *\n"
      " * The Sandals project is distributed under the BSD 2-Clause License.        *\n"
      " *                                                                           *\n"
      " * Davide Stocco                                           Enrico Bertolazzi *\n"
      " * University of Trento                                 University of Trento *\n"
      " * e-mail: davide.stocco@unitn.it         e-mail: enrico.bertolazzi@unitn.it *\n"
      "\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */\n"
      "\n",
      "// C++ generated code for implicit system: ", name, "\n",
      "// This file has been automatically generated by Sandals.\n",
      "// DISCLAIMER: If you need to edit it, do it wisely!\n",
      "\n",
      "#ifndef SANDALS_CODEGEN_", StringTools:-UpperCase(name), "_HH\n"
      "#define SANDALS_CODEGEN_", StringTools:-UpperCase(name), "_HH\n"
      "\n",
      "using namespace Sandals;\n"
      "\n",
      "// ", info, "\n",
      "class ", name, " : public Explicit<", num_eqns, ", ", num_invs, ">\n",
      "{\n",
      i, "// User data\n",
      data_str,
      "\n",
      "public:\n",
      i, "using VectorF  = Eigen::Vector<Real, ", num_eqns, ">;\n",
      i, "using MatrixJF = Eigen::Matrix<Real, ", num_eqns, ", ", num_eqns, ">;\n",
      i, "using VectorH  = Eigen::Vector<Real, ", num_invs, ">;\n",
      i, "using MatrixJH = Eigen::Matrix<Real, ", num_invs, ", ", num_eqns, ">;\n",
      "\n",
      SandalsCodegen:-ApplyIndent(i, GenerateConstructor(
          name, "Explicit",
          parse("num_eqns") = LinearAlgebra:-Dimension(f_tmp),
          parse("num_invs") = LinearAlgebra:-Dimension(h_tmp),
          parse("data")     = properties,
          parse("info")     = "Class constructor.",
          parse("indent")   = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "f", [x_tmp], f_tmp,
          ["Vector2 const &"], "VectorF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the function f.",
          parse("indent")  = i
      )),
      "\n",
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jf_x", [x_tmp], Jf_x,
          ["Vector2 const &"], "MatrixJF",
          parse("data")    = properties,
          parse("info")    = "Evaluate the Jacobian of f with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-VectorToCpp(
          "h", [x_tmp], h_tmp,
          ["Vector2 const &"], "VectorH",
          parse("data")    = properties,
          parse("info")    = "Calculate the vector h of the invariants.",
          parse("indent")  = i
      )),
      "\n",
      SandalsCodegen:-ApplyIndent(i, SandalsCodegen:-MatrixToCpp(
          "Jh_x", [x_tmp], Jh_x,
          ["Vector2 const &"], "MatrixJH",
          parse("data")    = properties,
          parse("info")    = "Calculate the Jacobian of h with respect to x.",
          parse("indent")  = i
      )),
      "\n",
      i, "bool in_domain(Vector2 const &/*x*/, Real /*t*/) const override {return true;}\n",
      "\n",
      "} // class ", name, "\n",
      "\n",
      "#endif // SANDALS_CODEGEN_", StringTools:-UpperCase(name), "_HH\n"
    );
    end proc: # ExplicitSystem

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # SandalsCodegen

# That's all folks!
