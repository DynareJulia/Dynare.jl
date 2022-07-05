::: {.default-domain}
dynare
:::

Running Dynare
==============

In order to give instructions to Dynare, the user has to write a *model
file* whose filename extension must be `.mod` or `.dyn`. This file
contains the description of the model and the computing tasks required
by the user. Its contents are described in
`model-file`{.interpreted-text role="ref"}.

Dynare invocation {#dyn-invoc}
-----------------

Once the model file is written, Dynare is invoked using the `dynare`
command at the MATLAB or Octave prompt (with the filename of the `.mod`
given as argument).

In practice, the handling of the model file is done in two steps: in the
first one, the model and the processing instructions written by the user
in a *model file* are interpreted and the proper MATLAB or Octave
instructions are generated; in the second step, the program actually
runs the computations. Both steps are triggered automatically by the
`dynare` command.

> ::: {.matcomm}
> context = @dynare "FILENAME\[.mod\]" \["OPTIONS"...\];
> :::
>
> This command launches Dynare and executes the instructions included in
> `FILENAME.mod`. This user-supplied file contains the model and the
> processing instructions, as described in
> `model-file`{.interpreted-text role="ref"}. The options, listed below,
> can be passed on the command line, following the name of the `.mod`
> file or in the first line of the `.mod` file itself (see below).
>
> dynare begins by launching the preprocessor on the `.mod file`. By
> default (unless the `use_dll`{.interpreted-text role="opt"} option has
> been given to `model`), the preprocessor creates three intermediary
> files:
>
> -   `+FILENAME/driver.m`
>
>     > Contains variable declarations, and computing tasks.
>
> -   `+FILENAME/dynamic.m`
>
>     > Contains the dynamic model equations. Note that Dynare might
>     > introduce auxiliary equations and variables (see
>     > `aux-variables`{.interpreted-text role="ref"}). Outputs are the
>     > residuals of the dynamic model equations in the order the
>     > equations were declared and the Jacobian of the dynamic model
>     > equations. For higher order approximations also the Hessian and
>     > the third-order derivatives are provided. When computing the
>     > Jacobian of the dynamic model, the order of the endogenous
>     > variables in the columns is stored in `M_.lead_lag_incidence`.
>     > The rows of this matrix represent time periods: the first row
>     > denotes a lagged (time t-1) variable, the second row a
>     > contemporaneous (time t) variable, and the third row a leaded
>     > (time t+1) variable. The columns of the matrix represent the
>     > endogenous variables in their order of declaration. A zero in
>     > the matrix means that this endogenous does not appear in the
>     > model in this time period. The value in the
>     > `M_.lead_lag_incidence` matrix corresponds to the column of that
>     > variable in the Jacobian of the dynamic model. Example: Let the
>     > second declared variable be `c` and the `(3,2)` entry of
>     > `M_.lead_lag_incidence` be 15. Then the 15th column of the
>     > Jacobian is the derivative with respect to `c(+1)`.
>
> -   `+FILENAME/static.m`
>
>     > Contains the long run static model equations. Note that Dynare
>     > might introduce auxiliary equations and variables (see
>     > `aux-variables`{.interpreted-text role="ref"}). Outputs are the
>     > residuals of the static model equations in the order the
>     > equations were declared and the Jacobian of the static
>     > equations. Entry `(i,j)` of the Jacobian represents the
>     > derivative of the ith static model equation with respect to the
>     > jth model variable in declaration order.
>
> These files may be looked at to understand errors reported at the
> simulation stage.
>
> `dynare` will then run the computing tasks by executing
> `+FILENAME/driver.m`. If a user needs to rerun the computing tasks
> without calling the preprocessor (or without calling the
> `dynare`{.interpreted-text role="mcomm"} command), for instance
> because he has modified the script, he just have to type the following
> on the command line:
>
> ``` {.sourceCode .matlab}
> >> FILENAME.driver
> ```
>
> A few words of warning are warranted here: under Octave the filename
> of the `.mod` file should be chosen in such a way that the generated
> `.m` files described above do not conflict with `.m` files provided by
> Octave or by Dynare. Not respecting this rule could cause crashes or
> unexpected behaviour. In particular, it means that the `.mod` file
> cannot be given the name of an Octave or Dynare command. For instance,
> under Octave, it also means that the `.mod` file cannot be named
> `test.mod` or `example.mod`.
>
> ::: {#quote-note}
> ::: {.note}
> ::: {.admonition-title}
> Note
> :::
>
> Note on Quotes
>
> When passing command line options that contains a space (or, under
> Octave, a double quote), you must surround the entire option (keyword
> and argument) with single quotes, as in the following example.
>
> *Example*
>
> Call Dynare with options containing spaces
>
> ``` {.sourceCode .matlab}
> >> dynare <<modfile.mod>> '-DA=[ i in [1,2,3] when i > 1 ]' 'conffile=C:\User\My Documents\config.txt'
> ```
> :::
> :::
>
> *Options*
>
> ::: {.option}
> noclearall
>
> By default, `dynare` will issue a `clear all` command to MATLAB
> (\<R2015b) or Octave, thereby deleting all workspace variables and
> functions; this option instructs `dynare` not to clear the workspace.
> Note that starting with MATLAB 2015b `dynare` only deletes the global
> variables and the functions using persistent variables, in order to
> benefit from the JIT (Just In Time) compilation. In this case the
> option instructs `dynare` not to clear the globals and functions.
> :::
>
> ::: {.option}
> onlyclearglobals
>
> By default, `dynare` will issue a `clear all` command to MATLAB
> versions before 2015b and to Octave, thereby deleting all workspace
> variables; this option instructs `dynare` to clear only the global
> variables (i.e. `M_, options_, oo_, estim_params_, bayestopt_`, and
> `dataset_`), leaving the other variables in the workspace.
> :::
>
> ::: {.option}
> debug
>
> Instructs the preprocessor to write some debugging information about
> the scanning and parsing of the `.mod` file.
> :::
>
> ::: {.option}
> notmpterms
>
> Instructs the preprocessor to omit temporary terms in the static and
> dynamic files; this generally decreases performance, but is used for
> debugging purposes since it makes the static and dynamic files more
> readable.
> :::
>
> ::: {.option}
> savemacro\[=FILENAME\]
>
> Instructs `dynare` to save the intermediary file which is obtained
> after macro processing (see `macro-proc-lang`{.interpreted-text
> role="ref"}); the saved output will go in the file specified, or if no
> file is specified in `FILENAME-macroexp.mod`. See the
> `note on quotes<quote-note>`{.interpreted-text role="ref"} for info on
> passing a `FILENAME` argument containing spaces.
> :::
>
> ::: {.option}
> onlymacro
>
> Instructs the preprocessor to only perform the macro processing step,
> and stop just after. Useful for debugging purposes or for using the
> macro processor independently of the rest of Dynare toolbox.
> :::
>
> ::: {.option}
> linemacro
>
> Instructs the macro preprocessor include `@#line` directives
> specifying the line on which macro directives were encountered and
> expanded from. Only useful in conjunction with `savemacro
> <savemacro[=FILENAME]>`{.interpreted-text role="opt"}.
> :::
>
> ::: {.option}
> onlymodel
>
> Instructs the preprocessor to print only information about the model
> in the driver file; no Dynare commands (other than the shocks
> statement and parameter initializations) are printed and hence no
> computational tasks performed. The same ancillary files are created as
> would otherwise be created (dynamic, static files, etc.).
> :::
>
> ::: {.option}
> nolog
>
> Instructs Dynare to no create a logfile of this run in `FILENAME.log.`
> The default is to create the logfile.
> :::
>
> ::: {.option}
> output=second\|third
>
> Instructs the preprocessor to output derivatives of the dynamic model
> at least up to the given order.
> :::
>
> ::: {.option}
> language=matlab\|julia
>
> Instructs the preprocessor to write output for MATLAB or Julia.
> Default: MATLAB
> :::
>
> ::: {.option}
> params\_derivs\_order=02
>
> When `identification`{.interpreted-text role="comm"},
> `dynare_sensitivity`{.interpreted-text role="comm"} (with
> identification), or `estimation_cmd <estim-comm>`{.interpreted-text
> role="ref"} are present, this option is used to limit the order of the
> derivatives with respect to the parameters that are calculated by the
> preprocessor. 0 means no derivatives, 1 means first derivatives, and 2
> means second derivatives. Default: 2
> :::
>
> ::: {.option}
> nowarn
>
> Suppresses all warnings.
> :::
>
> ::: {.option}
> notime
>
> Do not print the total computing time at the end of the driver, and do
> not save that total computing time to `oo_.time`.
> :::
>
> ::: {.option}
> transform\_unary\_ops
>
> Transform the following operators in the model block into auxiliary
> variables: `exp`, `log`, `log10`, `cos`, `sin`, `tan`, `acos`, `asin`,
> `atan`, `cosh`, `sinh`, `tanh`, `acosh`, `asinh`, `atanh`, `sqrt`,
> `cbrt`, `abs`, `sign`, `erf`. Default: no obligatory transformation
> :::
>
> ::: {.option}
> json = parsetransform\|compute
>
> Causes the preprocessor to output a version of the `.mod` file in JSON
> format to `<<M_.fname>>/model/json/`. When the JSON output is created
> depends on the value passed. These values represent various steps of
> processing in the preprocessor.
>
> If `parse` is passed, the output will be written after the parsing of
> the `.mod` file to a file called `FILENAME.json` but before file has
> been checked (e.g. if there are unused exogenous in the model block,
> the JSON output will be created before the preprocessor exits).
>
> If `check` is passed, the output will be written to a file called
> `FILENAME.json` after the model has been checked.
>
> If `transform` is passed, the JSON output of the transformed model
> (maximum lead of 1, minimum lag of -1, expectation operators
> substituted, etc.) will be written to a file called `FILENAME.json`
> and the original, untransformed model will be written in
> `FILENAME_original.json`.
>
> And if `compute` is passed, the output is written after the computing
> pass. In this case, the transformed model is written to
> `FILENAME.json`, the original model is written to
> `FILENAME_original.json`, and the dynamic and static files are written
> to `FILENAME_dynamic.json` and `FILENAME_static.json`.
> :::
>
> ::: {.option}
> jsonstdout
>
> Instead of writing output requested by `json` to files, write to
> standard out, i.e. to the MATLAB/Octave command window (and the
> log-file).
> :::
>
> ::: {.option}
> onlyjson
>
> Quit processing once the output requested by `json` has been written.
> :::
>
> ::: {.option}
> jsonderivsimple
>
> Print a simplified version (excluding variable name(s) and lag
> information) of the static and dynamic files in `FILENAME_static.json`
> and `FILENAME_dynamic.`.
> :::
>
> ::: {.option}
> warn\_uninit
>
> Display a warning for each variable or parameter which is not
> initialized. See `param-init`{.interpreted-text role="ref"}, or
> `load_params_and_steady_state
> <load_params_and_steady_state>`{.interpreted-text role="comm"} for
> initialization of parameters. See `init-term-cond`{.interpreted-text
> role="ref"}, or `load_params_and_steady_state
> <load_params_and_steady_state>`{.interpreted-text role="comm"} for
> initialization of endogenous and exogenous variables.
> :::
>
> ::: {.option}
> console
>
> Activate console mode. In addition to the behavior of `nodisplay`,
> Dynare will not use graphical waitbars for long computations.
> :::
>
> ::: {.option noindex=""}
> nograph
>
> Activate the `nograph` option (see `nograph`{.interpreted-text
> role="opt"}), so that Dynare will not produce any graph.
> :::
>
> ::: {.option}
> nointeractive
>
> Instructs Dynare to not request user input.
> :::
>
> ::: {.option}
> nopathchange
>
> By default Dynare will change MATLAB/Octave's path if `dynare/matlab`
> directory is not on top and if Dynare's routines are overriden by
> routines provided in other toolboxes. If one wishes to override
> Dynare's routines, the `nopathchange` options can be used.
> Alternatively, the path can be temporarly modified by the user at the
> top of the `.mod` file (using MATLAB/Octave's `addpath` command).
> :::
>
> ::: {.option}
> nopreprocessoroutput
>
> Prevent Dynare from printing the output of the steps leading up to the
> preprocessor as well as the preprocessor output itself.
> :::
>
> ::: {.option}
> mexext=mexmexw64mexa64
>
> The mex extension associated with your platform to be used when
> compiling output associated with `use_dll`{.interpreted-text
> role="opt"}. Dynare is able to set this automatically, so you should
> not need to set it yourself.
> :::
>
> ::: {.option}
> matlabroot=\<\<path\>\>
>
> The path to the MATLAB installation for use with
> `use_dll`{.interpreted-text role="opt"}. Dynare is able to set this
> automatically, so you should not need to set it yourself. See the
> `note on quotes<quote-note>`{.interpreted-text role="ref"} for info on
> passing a `<<path>>` argument containing spaces.
> :::
>
> ::: {.option}
> parallel\[=CLUSTER\_NAME\]
>
> Tells Dynare to perform computations in parallel. If CLUSTER\_NAME is
> passed, Dynare will use the specified cluster to perform parallel
> computations. Otherwise, Dynare will use the first cluster specified
> in the configuration file. See `conf-file`{.interpreted-text
> role="ref"}, for more information about the configuration file.
> :::
>
> ::: {.option}
> conffile=FILENAME
>
> Specifies the location of the configuration file if it differs from
> the default. See `conf-file`{.interpreted-text role="ref"}, for more
> information about the configuration file and its default location. See
> the `note on
> quotes<quote-note>`{.interpreted-text role="ref"} for info on passing
> a `FILENAME` argument containing spaces.
> :::
>
> ::: {.option}
> parallel\_follower\_open\_mode
>
> Instructs Dynare to leave the connection to the follower node open
> after computation is complete, closing this connection only when
> Dynare finishes processing.
> :::
>
> ::: {.option}
> parallel\_test
>
> Tests the parallel setup specified in the configuration file without
> executing the `.mod` file. See `conf-file`{.interpreted-text
> role="ref"}, for more information about the configuration file.
> :::
>
> ::: {.option}
> -DMACRO\_VARIABLE\[=MACRO\_EXPRESSION\]
>
> Defines a macro-variable from the command line (the same effect as
> using the Macro directive `@#define` in a model file, see
> `macro-proc-lang`{.interpreted-text role="ref"}). See the
> `note on quotes<quote-note>`{.interpreted-text role="ref"} for info on
> passing a `MACRO_EXPRESSION` argument containing spaces. Note that an
> expression passed on the command line can reference variables defined
> before it. If `MACRO_EXPRESSION` is omitted, the variable is assigned
> the `true` logical value.
>
> *Example*
>
> Call dynare with command line defines
>
> > ``` {.sourceCode .matlab}
> > >> dynare <<modfile.mod>> -DA=true '-DB="A string with space"' -DC=[1,2,3] '-DD=[ i in C when i > 1 ]' -DE
> > ```
> :::
>
> ::: {.option}
> -I\<\<path\>\>
>
> Defines a path to search for files to be included by the macro
> processor (using the `@#include` command). Multiple `-I` flags can be
> passed on the command line. The paths will be searched in the order
> that the `-I` flags are passed and the first matching file will be
> used. The flags passed here take priority over those passed to
> `@#includepath`. See the
> `note on quotes<quote-note>`{.interpreted-text role="ref"} for info on
> passing a `<<path>>` argument containing spaces.
> :::
>
> ::: {.option}
> nostrict
>
> Allows Dynare to issue a warning and continue processing when
>
> 1.  there are more endogenous variables than equations.
> 2.  an undeclared symbol is assigned in `initval` or `endval`.
> 3.  an undeclared symbol is found in the `model` block in this case,
>     it is automatically declared exogenous.
> 4.  exogenous variables were declared but not used in the `model`
>     block.
> :::
>
> ::: {.option}
> fast
>
> Only useful with model option `use_dll`{.interpreted-text role="opt"}.
> Don't recompile the MEX files when running again the same model file
> and the lists of variables and the equations haven't changed. We use a
> 32 bit checksum, stored in `<model filename>/checksum`. There is a
> very small probability that the preprocessor misses a change in the
> model. In case of doubt, re-run without the fast option.
> :::
>
> ::: {.option}
> minimal\_workspace
>
> Instructs Dynare not to write parameter assignments to parameter names
> in the .m file produced by the preprocessor. This is potentially
> useful when running `dynare` on a large `.mod` file that runs into
> workspace size limitations imposed by MATLAB.
> :::
>
> ::: {.option}
> compute\_xrefs
>
> Tells Dynare to compute the equation cross references, writing them to
> the output `.m` file.
> :::
>
> ::: {.option}
> stochastic
>
> Tells Dynare that the model to be solved is stochastic. If no Dynare
> commands related to stochastic models (`stoch_simul`, `estimation`,
> \...) are present in the `.mod` file, Dynare understands by default
> that the model to be solved is deterministic.
> :::
>
> ::: {#exclude_eqs}
> ::: {.option}
> exclude\_eqs=\<\<equation\_tags\_to\_exclude\>\>
>
> Tells Dynare to exclude all equations specified by the argument. As a
> `.mod` file must have the same number of endogenous variables as
> equations, when [exclude\_eqs]{.title-ref} is passed, certain rules
> are followed for excluding endogenous variables. If the `endogenous`
> tag has been set for the excluded equation, the variable it specifies
> is excluded. Otherwise, if the left hand side of the excluded equation
> is an expression that contains only one endogenous variable, that
> variable is excluded. If neither of these conditions hold, processing
> stops with an error. If an endogenous variable has been excluded by
> the [exclude\_eqs]{.title-ref} option and it exists in an equation
> that has not been excluded, it is transformed into an exogenous
> variable.
>
> To specify which equations to exclude, you must pass the argument
> `<<equation_tags_to_exclude>>`. This argument takes either a list of
> equation tags specifying the equations to be excluded or a filename
> that contains those tags.
>
> If `<<equation_tags_to_exclude>>` is a list of equation tags, it can
> take one of the following forms:
>
> 1.  Given a single argument, e.g. `exclude_eqs=eq1`, the equation with
>     the tag `[name='eq1']` will be excluded. Note that if there is a
>     file called `eq1` in the current directory, Dynare will instead
>     try to open this and read equations to exclude from it (see info
>     on filename argument to `exclude_eqs` below). Further note that if
>     the tag value contains a space, you must use the variant specified
>     in 2 below, i.e. `exclude_eqs=[eq 1]`.
> 2.  Given two or more arguments, e.g. `exclude_eqs=[eq1, eq 2]`, the
>     equations with the tags `[name='eq1']` and `[name='eq 2']` will be
>     excluded.
> 3.  If you\'d like to exclude equations based on another tag name (as
>     opposed to the default `name`), you can pass the argument as
>     either e.g. `exclude_eqs=[tagname=a tag]` if a single equation
>     with tag `[tagname='a tag']` is to be excluded or as e.g.
>     `exclude_eqs=[tagname=(a tag, 'a tag with a, comma')]` if more
>     than one equation with tags `[tagname='a tag']` and
>     `[tagname='a tag with a, comma']` will be excluded (note the
>     parenthesis, which are required when more than one equation is
>     specified). Note that if the value of a tag contains a comma, it
>     must be included inside single quotes.
>
> If `<<equation_tags_to_exclude>>` is a filename, the file can take one
> of the following forms:
>
> 1.  One equation per line of the file, where every line represents the
>     value passed to the `name` tag. e.g., a file such as:
>
>         eq1
>         eq 2
>
>     would exclude equations with tags `[name='eq1']` and
>     `[name='eq 2']`.
>
> 2.  One equation per line of the file, where every line after the
>     first line represents the value passed to the tag specified by the
>     first line. e.g., a file such as:
>
>         tagname=
>         a tag
>         a tag with a, comma
>
>     would exclude equations with tags `[tagname='a tag']` and
>     `[tagname='a tag with a, comma']`. Here note that the first line
>     must end in an equal sign.
> :::
> :::
>
> ::: {.option}
> include\_eqs=\<\<equation\_tags\_to\_include\>\>
>
> Tells Dynare to run with only those equations specified by the
> argument; in other words, Dynare will exclude all equations not
> specified by the argument. The argument `<<equation_tags_to_include>>`
> is specified in the same way as the argument to `exclude_eqs
> <exclude_eqs>`{.interpreted-text role="ref"}. The functionality of
> `include_eqs` is to find which equations to exclude then take actions
> in accord with `exclude_eqs
> <exclude_eqs>`{.interpreted-text role="ref"}.
> :::
>
> ::: {.option noindex=""}
> use\_dll
>
> Instructs the preprocessor to create dynamic loadable libraries (DLL)
> containing the model equations and derivatives, instead of writing
> those in M-files. This is equivalent to the
> `use_dll`{.interpreted-text role="opt"} option of the `model` block.
> :::
>
> ::: {.option}
> nocommutativity
>
> This option tells the preprocessor not to use the commutativity of
> addition and multiplication when looking for common subexpressions. As
> a consequence, when using this option, equations in various outputs
> (LaTeX, JSON...) will appear as the user entered them (without terms
> or factors swapped). Note that using this option may have a
> performance impact on the preprocessing stage, though it is likely to
> be small.
> :::
>
> These options can be passed to the preprocessor by listing them after
> the name of the `.mod` file. They can alternatively be defined in the
> first line of the `.mod` file, this avoids typing them on the command
> line each time a `.mod` file is to be run. This line must be a Dynare
> one-line comment (i.e. must begin with `//`) and the options must be
> whitespace separated between `--+ options:` and `+--`. Note that any
> text after the `+--` will be discarded. As in the command line, if an
> option admits a value the equal symbol must not be surrounded by
> spaces. For instance `json = compute` is not correct, and should be
> written `json=compute`. The `nopathchange` option cannot be specified
> in this way, it must be passed on the command-line.
>
> *Output*
>
> Depending on the computing tasks requested in the `.mod` file,
> executing the `dynare` command will leave variables containing results
> in the workspace available for further processing. More details are
> given under the relevant computing tasks. The `M_`,`oo_`, and
> `options_` structures are saved in a file called
> `FILENAME_results.mat` located in the `MODFILENAME/Output` folder. If
> they exist, `estim_params_`, `bayestopt_`, `dataset_`, `oo_recursive_`
> and `estimation_info` are saved in the same file. Note that Matlab by
> default only allows `.mat`-files up to 2GB. You can lift this
> restriction by enabling the `save -v7.3`-option in
> `Preferences -> General -> MAT-Files`.
>
> ::: {.matvar}
> [M]()
>
> Structure containing various information about the model.
> :::
>
> ::: {.matvar}
> [options]()
>
> Structure contains the values of the various options used by Dynare
> during the computation.
> :::
>
> ::: {.matvar}
> [oo]()
>
> Structure containing the various results of the computations.
> :::
>
> ::: {.matvar}
> [dataset]()
>
> A `dseries` object containing the data used for estimation.
> :::
>
> ::: {.matvar}
> [oo\_recursive]()
>
> Cell array containing the `oo_` structures obtained when estimating
> the model for the different samples when performing recursive
> estimation and forecasting. The `oo_` structure obtained for the
> sample ranging to the [i]{.title-ref} -th observation is saved in the
> [i]{.title-ref} -th field. The fields for non-estimated endpoints are
> empty.
> :::
>
> ::: {.matvar}
> [oo]().time
>
> Total computing time of the Dynare run, in seconds. This field is not
> set if the `notime`{.interpreted-text role="opt"} option has been
> used.
> :::
>
> *Example*
>
> Call dynare from the MATLAB or Octave prompt, without or with options:
>
> > ``` {.sourceCode .matlab}
> > >> dynare ramst
> > >> dynare ramst.mod savemacro
> > ```
>
> Alternatively the options can be passed in the first line of
> `ramst.mod`:
>
> > ``` {.sourceCode .dynare}
> > // --+ options: savemacro, json=compute +--
> > ```
>
> and then dynare called without passing options on the command line:
>
> > ``` {.sourceCode .matlab}
> > >> dynare ramst
> > ```

Dynare hooks
------------


Understanding Preprocessor Error Messages
-----------------------------------------

If the preprocessor runs into an error while processing your `.mod`
file, it will issue an error. Due to the way that a parser works,
sometimes these errors can be misleading. Here, we aim to demystify
these error messages.

The preprocessor issues error messages of the form:

> 1.  `ERROR: <<file.mod>>: line A, col B: <<error message>>`
> 2.  `ERROR: <<file.mod>>: line A, cols B-C: <<error message>>`
> 3.  `ERROR: <<file.mod>>: line A, col B - line C, col D: <<error message>>`

The first two errors occur on a single line, with error two spanning
multiple columns. Error three spans multiple rows.

Often, the line and column numbers are precise, leading you directly to
the offending syntax. Infrequently however, because of the way the
parser works, this is not the case. The most common example of
misleading line and column numbers (and error message for that matter)
is the case of a missing semicolon, as seen in the following example:

    varexo a, b
    parameters c, ...;

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

It is also helpful to keep in mind that any piece of code that does not
violate Dynare syntax, but at the same time is not recognized by the
parser, is interpreted as native MATLAB code. This code will be directly
passed to the `driver` script. Investigating `driver.m` file then helps
with debugging. Such problems most often occur when defined variable or
parameter names have been misspelled so that Dynare\'s parser is unable
to recognize them.
