# Running Dynare

In order to give instructions to Dynare, the user has to write a *model
file* whose filename extension must be `.mod`. This file
contains the description of the model and the computing tasks required
by the user. Its contents are described in `The model file`

## Dynare invocation

Once the model file is written, Dynare is invoked using the `@dynare`
Julia macro (with the filename of the `.mod`
given as argument).

In practice, the handling of the model file is done in two steps: in the
first one, the model and the processing instructions written by the user
in a *model file* are interpreted and the proper Julia 
instructions are generated; in the second step, the program actually
runs the computations. Both steps are triggered automatically by the
`@dynare` macro:

```
context = @dynare "FILENAME [.mod ]" [OPTIONS... ]";
```

This command launches Dynare and executes the instructions included in
`FILENAME.mod`. This user-supplied file contains the model and the
processing instructions, as described in `The model file`. The 
options, listed below, can be passed on the command line, following the 
name of the `.mod` file or in the first line of the `.mod` file itself 
(see below).

Dynare begins by launching the preprocessor on the `.mod file`.
### Options

- `debug`

  Instructs the preprocessor to write some debugging information about
  the scanning and parsing of the `.mod` file.

- `notmpterms`

  Instructs the preprocessor to omit temporary terms in the static and
  dynamic files; this generally decreases performance, but is used for
  debugging purposes since it makes the static and dynamic files more
  readable.

- `savemacro\[=FILENAME\]`

  Instructs Dynare to save the intermediary file which is obtained
  after macro processing (see (@ref "Macro processing language")); the 
  saved output will go in the file specified, or if no file is specified 
  in `FILENAME-macroexp.mod`. See the (@ref "note on quotes") for info on 
  passing a `FILENAME` argument containing spaces.

- `onlymacro`

  Instructs the preprocessor to only perform the macro processing step,
  and stop just after. Useful for debugging purposes or for using the
  macro processor independently of the rest of Dynare toolbox.

- `linemacro`

  Instructs the macro preprocessor include `@#line` directives
  specifying the line on which macro directives were encountered and
  expanded from. Only useful in conjunction with `savemacro
  <savemacro[=FILENAME]>`.

- `onlymodel`

  Instructs the preprocessor to print only information about the model
  in the driver file; no Dynare commands (other than the shocks
  statement and parameter initializations) are printed and hence no
  computational tasks performed. The same ancillary files are created as
  would otherwise be created (dynamic, static files, etc.).

- `nolog`

  Instructs Dynare to no create a logfile of this run in `FILENAME.log.`
  The default is to create the logfile.

- `output=second\|third`

  Instructs the preprocessor to output derivatives of the dynamic model
  at least up to the given order.

- `language=matlab\|julia`
 
  Instructs the preprocessor to write output for MATLAB or Julia.
  Default: MATLAB

- `params\_derivs\_order=0\|1\|2`

  When (@ref "identification"), (@ref "dynare_sensitivity") (with
  identification), or (@ref "estimation_cmd") are present, this option 
  is used to limit the order of the derivatives with respect to the 
  parameters that are calculated by the preprocessor. 0 means no derivatives, 
  1 means first derivatives, and 2 means second derivatives. Default: 2

- `transform\_unary\_ops`

  Transform the following operators in the model block into auxiliary
  variables: `exp`, `log`, `log10`, `cos`, `sin`, `tan`, `acos`, `asin`,
  `atan`, `cosh`, `sinh`, `tanh`, `acosh`, `asinh`, `atanh`, `sqrt`,
  `cbrt`, `abs`, `sign`, `erf`. Default: no obligatory transformation

- `json = parse\|transform\|compute`

  Causes the preprocessor to output a version of the `.mod` file in JSON
  format to `<<M_.fname>>/model/json/`. When the JSON output is created
  depends on the value passed. These values represent various steps of
  processing in the preprocessor.

  If `parse` is passed, the output will be written after the parsing of
  the `.mod` file to a file called `FILENAME.json` but before file has
  been checked (e.g. if there are unused exogenous in the model block,
  the JSON output will be created before the preprocessor exits).
 
  If `check` is passed, the output will be written to a file called
  `FILENAME.json` after the model has been checked.
 
  If `transform` is passed, the JSON output of the transformed model
  (maximum lead of 1, minimum lag of -1, expectation operators
  substituted, etc.) will be written to a file called `FILENAME.json`
  and the original, untransformed model will be written in
  `FILENAME_original.json`.
 
  And if `compute` is passed, the output is written after the computing
  pass. In this case, the transformed model is written to
  `FILENAME.json`, the original model is written to
  `FILENAME_original.json`, and the dynamic and static files are written
  to `FILENAME_dynamic.json` and `FILENAME_static.json`.

- `jsonstdout`

  Instead of writing output requested by `json` to files, write to
  standard out, i.e. to the Julia command window (and the
  log-file).

- `onlyjson`

  Quit processing once the output requested by `json` has been written.

- `jsonderivsimple`

  Print a simplified version (excluding variable name(s) and lag
  information) of the static and dynamic files in `FILENAME_static.json`
  and `FILENAME_dynamic.`.

- `warn\_uninit`
 
  Display a warning for each variable or parameter which is not
  initialized. See (@ref "Parameter initialization"), or
  (@ref "load_params_and_steady_state") for initialization of parameters. 
  See (@ref "Initial and Terminal conditions"), or (@ref "load_params_and_steady_state")
  for initialization of endogenous and exogenous variables.

- `nopreprocessoroutput`

  Prevent Dynare from printing the output of the steps leading up to the
  preprocessor as well as the preprocessor output itself.


- `-DMACRO\_VARIABLE\[=MACRO\_EXPRESSION\]`

  Defines a macro-variable from the command line (the same effect as
  using the Macro directive `@#define` in a model file, see
  (@ref "Macro processing language")). See the (@ref "note on quotes") for 
  info on passing a `MACRO_EXPRESSION` argument containing spaces. Note 
  that an expression passed on the command line can reference variables 
  defined before it. If `MACRO_EXPRESSION` is omitted, the variable is assigned
  the `true` logical value.

  *Example*

  Call dynare with command line defines

   ```julia
   julia> @dynare <<modfile.mod>> -DA=true '-DB="A string with space"' -DC=[1,2,3] '-DD=[ i in C when i > 1 ]' -DE;
   ```

- `-I\<\<path\>\>`

  Defines a path to search for files to be included by the macro
  processor (using the `@#include` command). Multiple `-I` flags can be
  passed on the command line. The paths will be searched in the order
  that the `-I` flags are passed and the first matching file will be
  used. The flags passed here take priority over those passed to
  `@#includepath`. See the (@ref  "note on quotes") for info on
  passing a `<<path>>` argument containing spaces.

- `nostrict`

  Allows Dynare to issue a warning and continue processing when
 
  1.  there are more endogenous variables than equations.
  2.  an undeclared symbol is assigned in `initval` or `endval`.
  3.  an undeclared symbol is found in the `model` block in this case,
      it is automatically declared exogenous.
  4.  exogenous variables were declared but not used in the `model`
      block.

- `stochastic`

  Tells Dynare that the model to be solved is stochastic. If no Dynare
  commands related to stochastic models (`stoch_simul`, `estimation`,
  \...) are present in the `.mod` file, Dynare understands by default
  that the model to be solved is deterministic.


- `exclude\_eqs=\<\<equation\_tags\_to\_exclude\>\>`

  Tells Dynare to exclude all equations specified by the argument. As a
  `.mod` file must have the same number of endogenous variables as
  equations, when ***exclude\_eqs*** is passed, certain rules
  are followed for excluding endogenous variables. If the `endogenous`
  tag has been set for the excluded equation, the variable it specifies
  is excluded. Otherwise, if the left hand side of the excluded equation
  is an expression that contains only one endogenous variable, that
  variable is excluded. If neither of these conditions hold, processing
  stops with an error. If an endogenous variable has been excluded by
  the ***exclude\_eqs*** option and it exists in an equation that has not 
  been excluded, it is transformed into an exogenous variable.

  To specify which equations to exclude, you must pass the argument
  `<<equation_tags_to_exclude>>`. This argument takes either a list of
  equation tags specifying the equations to be excluded or a filename
  that contains those tags.
 
  If `<<equation_tags_to_exclude>>` is a list of equation tags, it can
  take one of the following forms:
 
  1.  Given a single argument, e.g. `exclude_eqs=eq1`, the equation with
      the tag `[name='eq1']` will be excluded. Note that if there is a
      file called `eq1` in the current directory, Dynare will instead
      try to open this and read equations to exclude from it (see info
      on filename argument to `exclude_eqs` below). Further note that if
      the tag value contains a space, you must use the variant specified
      in 2 below, i.e. `exclude_eqs=[eq 1]`.
  2.  Given two or more arguments, e.g. `exclude_eqs=[eq1, eq 2]`, the
      equations with the tags `[name='eq1']` and `[name='eq 2']` will be
      excluded.
  3.  If you\'d like to exclude equations based on another tag name (as
      opposed to the default `name`), you can pass the argument as
      either e.g. `exclude_eqs=[tagname=a tag]` if a single equation
      with tag `[tagname='a tag']` is to be excluded or as e.g.
      `exclude_eqs=[tagname=(a tag, 'a tag with a, comma')]` if more
      than one equation with tags `[tagname='a tag']` and
      `[tagname='a tag with a, comma']` will be excluded (note the
      parenthesis, which are required when more than one equation is
      specified). Note that if the value of a tag contains a comma, it
      must be included inside single quotes.

  If `<<equation_tags_to_exclude>>` is a filename, the file can take one
  of the following forms:

  1.  One equation per line of the file, where every line represents the
      value passed to the `name` tag. e.g., a file such as:
        ```julia    
          eq1
          eq 2
        ```
      would exclude equations with tags `[name='eq1']` and
      `[name='eq 2']`.

  2.  One equation per line of the file, where every line after the
      first line represents the value passed to the tag specified by the
      first line. e.g., a file such as:
        ``` julia
          tagname=
          a tag
          a tag with a, comma
        ```  
      would exclude equations with tags `[tagname='a tag']` and
      `[tagname='a tag with a, comma']`. Here note that the first line
      must end in an equal sign.

- `include\_eqs=\<\<equation\_tags\_to\_include\>\>`

  Tells Dynare to run with only those equations specified by the
  argument; in other words, Dynare will exclude all equations not
  specified by the argument. The argument `<<equation_tags_to_include>>`
  is specified in the same way as the argument to `exclude_eqs
  <exclude_eqs>`{.interpreted-text role="ref"}. The functionality of
  `include_eqs` is to find which equations to exclude then take actions
  in accord with (@ref "exclude_eqs").

- `nocommutativity`
 
  This option tells the preprocessor not to use the commutativity of
  addition and multiplication when looking for common subexpressions. As
  a consequence, when using this option, equations in various outputs
  (LaTeX, JSON...) will appear as the user entered them (without terms
  or factors swapped). Note that using this option may have a
  performance impact on the preprocessing stage, though it is likely to
  be small.

  These options can be passed to the preprocessor by listing them after
  the name of the `.mod` file. They can alternatively be defined in the
  first line of the `.mod` file, this avoids typing them on the command
  line each time a `.mod` file is to be run. This line must be a Dynare
  one-line comment (i.e. must begin with `//`) and the options must be
  whitespace separated between `--+ options:` and `+--`. Note that any
  text after the `+--` will be discarded. As in the command line, if an
  option admits a value the equal symbol must not be surrounded by
  spaces. For instance `json = compute` is not correct, and should be
  written `json=compute`. The `nopathchange` option cannot be specified
  in this way, it must be passed on the command-line.


## Understanding Preprocessor Error Messages

If the preprocessor runs into an error while processing your `.mod`
file, it will issue an error. Due to the way that a parser works,
sometimes these errors can be misleading. Here, we aim to demystify
these error messages.

The preprocessor issues error messages of the form:

  1.  `ERROR: <<file.mod>>: line A, col B: <<error message>>`
  2.  `ERROR: <<file.mod>>: line A, cols B-C: <<error message>>`
  3.  `ERROR: <<file.mod>>: line A, col B - line C, col D: <<error message>>`

The first two errors occur on a single line, with error two spanning
multiple columns. Error three spans multiple rows.

Often, the line and column numbers are precise, leading you directly to
the offending syntax. Infrequently however, because of the way the
parser works, this is not the case. The most common example of
misleading line and column numbers (and error message for that matter)
is the case of a missing semicolon, as seen in the following example:

```
varexo a, b
parameters c, ...;
```

In this case, the parser doesn't know a semicolon is missing at the end
of the `varexo` command until it begins parsing the second line and
bumps into the `parameters` command. This is because we allow commands
to span multiple lines and, hence, the parser cannot know that the
second line will not have a semicolon on it until it gets there. Once
the parser begins parsing the second line, it realizes that it has
encountered a keyword, `parameters`, which it did not expect. Hence, it
throws an error of the form:
`ERROR: <<file.mod>>: line 2, cols 0-9: syntax error, unexpected PARAMETERS`.
In this case, you would simply place a semicolon at the end of line one
and the parser would continue processing.

