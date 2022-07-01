::: {.default-domain}
dynare
:::

The model file {#model-file}
==============

Conventions {#conv}
-----------

A model file contains a list of commands and of blocks. Each command and
each element of a block is terminated by a semicolon (;). Blocks are
terminated by `end;`.

If Dynare encounters an unknown expression at the beginning of a line or
after a semicolon, it will parse the rest of that line as native MATLAB
code, even if there are more statements separated by semicolons present.
To prevent cryptic error messages, it is strongly recommended to always
only put one statement/command into each line and start a new line after
each semicolon.[^1]

Lines of codes can be commented out line by line or as a block.
Single-line comments begin with `//` and stop at the end of the line.
Multiline comments are introduced by `/*` and terminated by `*/`.

*Examples*

>     // This is a single line comment
>
>     var x; // This is a comment about x
>
>     /* This is another inline comment about alpha */  alpha = 0.3;
>
>     /*
>      This comment is spanning
>      two lines.
>     */

Note that these comment marks should not be used in native MATLAB code
regions where the [%]{.title-ref} should be preferred instead to
introduce a comment. In a `verbatim` block, see
`verbatim`{.interpreted-text role="ref"}, this would result in a crash
since `//` is not a valid MATLAB statement).

Most Dynare commands have arguments and several accept options,
indicated in parentheses after the command keyword. Several options are
separated by commas.

In the description of Dynare commands, the following conventions are
observed:

-   Optional arguments or options are indicated between square brackets:
    '\[\]';
-   Repeated arguments are indicated by ellipses: "\...";
-   Mutually exclusive arguments are separated by vertical bars: '\|';
-   INTEGER indicates an integer number;
-   INTEGER\_VECTOR indicates a vector of integer numbers separated by
    spaces, enclosed by square brackets;
-   DOUBLE indicates a double precision number. The following syntaxes
    are valid: `1.1e3`, `1.1E3`, `1.1d3`, `1.1D3`. In some places,
    infinite Values `Inf` and `-Inf` are also allowed;
-   NUMERICAL\_VECTOR indicates a vector of numbers separated by spaces,
    enclosed by square brackets;
-   EXPRESSION indicates a mathematical expression valid outside the
    model description (see `expr`{.interpreted-text role="ref"});
-   MODEL\_EXPRESSION (sometimes MODEL\_EXP) indicates a mathematical
    expression valid in the model description (see
    `expr`{.interpreted-text role="ref"} and
    `model-decl`{.interpreted-text role="ref"});
-   MACRO\_EXPRESSION designates an expression of the macro processor
    (see `macro-exp`{.interpreted-text role="ref"});
-   VARIABLE\_NAME (sometimes VAR\_NAME) indicates a variable name
    starting with an alphabetical character and can't contain:
    '()+-\*/\^=!;:@\#.' or accentuated characters;
-   PARAMETER\_NAME (sometimes PARAM\_NAME) indicates a parameter name
    starting with an alphabetical character and can't contain:
    '()+-\*/\^=!;:@\#.' or accentuated characters;
-   LATEX\_NAME (sometimes TEX\_NAME) indicates a valid LaTeX expression
    in math mode (not including the dollar signs);
-   FUNCTION\_NAME indicates a valid MATLAB function name;
-   FILENAME indicates a filename valid in the underlying operating
    system; it is necessary to put it between quotes when specifying the
    extension or if the filename contains a non-alphanumeric character;
-   QUOTED\_STRING indicates an arbitrary string enclosed between
    (single) quotes.

Variable declarations {#var-decl}
---------------------

While Dynare allows the user to choose their own variable names, there
are some restrictions to be kept in mind. First, variables and
parameters must not have the same name as Dynare commands or built-in
functions. In this respect, Dynare is not case-sensitive. For example,
do not use `Ln` or `Sigma_e` to name your variable. Not conforming to
this rule might yield hard-to-debug error messages or crashes. Second,
when employing user-defined steady state files it is recommended to
avoid using the name of MATLAB functions as this may cause conflicts. In
particular, when working with user-defined steady state files, do not
use correctly-spelled greek names like [alpha]{.title-ref}, because
there are MATLAB functions of the same name. Rather go for `alppha` or
`alph`. Lastly, please do not name a variable or parameter `i`. This may
interfere with the imaginary number i and the index in many loops.
Rather, name investment `invest`. Using `inv` is also not recommended as
it already denotes the inverse operator. Commands for declaring
variables and parameters are described below.

::: {.command}
var VAR\_NAME \[\$TEX\_NAME\$\]
\[(long\_name=QUOTED\_STRINGNAME=QUOTED\_STRING)\]\...;
var(deflator=MODEL\_EXPR) VAR\_NAME (\... same options apply) var(log,
deflator=MODEL\_EXPR) VAR\_NAME (\... same options apply)
var(log\_deflator=MODEL\_EXPR) VAR\_NAME (\... same options apply)

This required command declares the endogenous variables in the model.
See `conv`{.interpreted-text role="ref"} for the syntax of *VAR\_NAME*
and *MODEL\_EXPR*. Optionally it is possible to give a LaTeX name to the
variable or, if it is nonstationary, provide information regarding its
deflator. The variables in the list can be separated by spaces or by
commas. `var` commands can appear several times in the file and Dynare
will concatenate them. Dynare stores the list of declared parameters, in
the order of declaration, in a column cell array `M_.endo_names`.

If the model is nonstationary and is to be written as such in the
`model` block, Dynare will need the trend deflator for the appropriate
endogenous variables in order to stationarize the model. The trend
deflator must be provided alongside the variables that follow this
trend.

*Options*

::: {.option}
log

In addition to the endogenous variable(s) thus declared, this option
also triggers the creation of auxiliary variable(s) equal to the log of
the corresponding endogenous variable(s). For example, given a
`var(log) y` statement, two endogenous will be created (`y` and
`LOG_y`), and an auxiliary equation linking the two will also be added
(equal to `LOG_y = log(y)`). Moreover, every occurence of `y` in the
model will be replaced by `exp(LOG_y)`. This option is for example
useful when one wants to perform a loglinear approximation of some
variable(s) in the context of a first-order stochastic approximation; or
when one wants to ensure the variable(s) stay(s) in the definition
domain of the function defining the steady state or the dynamic
residuals when the nonlinear solver is used.
:::

::: {.option}
deflator = MODEL\_EXPR

The expression used to detrend an endogenous variable. All trend
variables, endogenous variables and parameters referenced in MODEL\_EXPR
must already have been declared by the `trend_var, log_trend_var, var`
and `parameters` commands. The deflator is assumed to be multiplicative;
for an additive deflator, use `log_deflator`. This option can be used
together with the `log` option (the latter must come first).
:::

::: {.option}
log\_deflator = MODEL\_EXPR

Same as `deflator`, except that the deflator is assumed to be additive
instead of multiplicative (or, to put it otherwise, the declared
variable is equal to the log of a variable with a multiplicative trend).
This option cannot be used together with the `log` option, because it
would not make much sense from an economic point of view (the
corresponding auxiliary variable would correspond to the log taken two
times on a variable with a multiplicative trend).
:::

::: {#long-name}
::: {.option}
long\_name = QUOTED\_STRING

This is the long version of the variable name. Its value is stored in
`M_.endo_names_long` (a column cell array, in the same order as
`M_.endo_names`). In case multiple `long_name` options are provided, the
last one will be used. Default: `VAR_NAME`.
:::
:::

::: {#partitioning}
::: {.option}
NAME = QUOTED\_STRING

This is used to create a partitioning of variables. It results in the
direct output in the `.m` file analogous to:
`M_.endo_partitions.NAME = QUOTED_STRING`;.
:::
:::

*Example (variable partitioning)*

>     var c gnp cva (country=`US', state=`VA')
>               cca (country=`US', state=`CA', long_name=`Consumption CA');
>     var(deflator=A) i b;
>     var c $C$ (long_name=`Consumption');
:::

::: {.command}
varexo\_det VAR\_NAME \[\$TEX\_NAME\$\]
\[(long\_name=QUOTED\_STRING\|NAME=QUOTED\_STRING)\...\];

This optional command declares exogenous deterministic variables in a
stochastic model. See `conv`{.interpreted-text role="ref"} for the
syntax of VARIABLE\_NAME. Optionally it is possible to give a LaTeX name
to the variable. The variables in the list can be separated by spaces or
by commas. `varexo_det` commands can appear several times in the file
and Dynare will concatenate them.

It is possible to mix deterministic and stochastic shocks to build
models where agents know from the start of the simulation about future
exogenous changes. In that case `stoch_simul` will compute the rational
expectation solution adding future information to the state space
(nothing is shown in the output of `stoch_simul`) and forecast will
compute a simulation conditional on initial conditions and future
information.

Note that exogenous deterministic variables cannot appear with a lead or
a lag in the model.

*Options*

::: {.option}
long\_name = QUOTED\_STRING

Like `long_name <long-name>`{.interpreted-text role="ref"} but value
stored in `M_.exo_det_names_long`.
:::

::: {.option}
NAME = QUOTED\_STRING

Like `partitioning <partitioning>`{.interpreted-text role="ref"} but
QUOTED\_STRING stored in `M_.exo_det_partitions.NAME`.
:::

*Example*

>     varexo m gov;
>     varexo_det tau;
:::

::: {.command}
var\_remove VAR\_NAME \| PARAM\_NAME\...;

Removes the listed variables (or parameters) from the model. Removing a
variable that has already been used in a model equation or elsewhere
will lead to an error.
:::

::: {.command}
predetermined\_variables VAR\_NAME\...;

In Dynare, the default convention is that the timing of a variable
reflects when this variable is decided. The typical example is for
capital stock: since the capital stock used at current period is
actually decided at the previous period, then the capital stock entering
the production function is `k(-1)`, and the law of motion of capital
must be written:

    k = i + (1-delta)*k(-1)

Put another way, for stock variables, the default in Dynare is to use a
"stock at the end of the period" concept, instead of a "stock at the
beginning of the period" convention.

The `predetermined_variables` is used to change that convention. The
endogenous variables declared as predetermined variables are supposed to
be decided one period ahead of all other endogenous variables. For stock
variables, they are supposed to follow a "stock at the beginning of the
period" convention.

Note that Dynare internally always uses the "stock at the end of the
period" concept, even when the model has been entered using the
`predetermined_variables` command. Thus, when plotting, computing or
simulating variables, Dynare will follow the convention to use variables
that are decided in the current period. For example, when generating
impulse response functions for capital, Dynare will plot `k`, which is
the capital stock decided upon by investment today (and which will be
used in tomorrow's production function). This is the reason that capital
is shown to be moving on impact, because it is `k` and not the
predetermined `k(-1)` that is displayed. It is important to remember
that this also affects simulated time series and output from smoother
routines for predetermined variables. Compared to non-predetermined
variables they might otherwise appear to be falsely shifted to the
future by one period.

*Example*

> The following two program snippets are strictly equivalent.
>
> Using default Dynare timing convention:
>
>     var y, k, i;
>     ...
>     model;
>     y = k(-1)^alpha;
>     k = i + (1-delta)*k(-1);
>     ...
>     end;
>
> Using the alternative timing convention:
>
>     var y, k, i;
>     predetermined_variables k;
>     ...
>     model;
>     y = k^alpha;
>     k(+1) = i + (1-delta)*k;
>     ...
>     end;
:::

::: {.command}
trend\_var (growth\_factor = MODEL\_EXPR) VAR\_NAME
\[\$LATEX\_NAME\$\]\...;

This optional command declares the trend variables in the model. See
`conv`{.interpreted-text role="ref"} for the syntax of MODEL\_EXPR and
VAR\_NAME. Optionally it is possible to give a LaTeX name to the
variable.

The variable is assumed to have a multiplicative growth trend. For an
additive growth trend, use `log_trend_var` instead.

Trend variables are required if the user wants to be able to write a
nonstationary model in the `model` block. The `trend_var` command must
appear before the var command that references the trend variable.

`trend_var` commands can appear several times in the file and Dynare
will concatenate them.

If the model is nonstationary and is to be written as such in the
`model` block, Dynare will need the growth factor of every trend
variable in order to stationarize the model. The growth factor must be
provided within the declaration of the trend variable, using the
`growth_factor` keyword. All endogenous variables and parameters
referenced in MODEL\_EXPR must already have been declared by the var and
parameters commands.

*Example*

>     trend_var (growth_factor=gA) A;
:::

::: {.command}
model\_local\_variable VARIABLE\_NAME \[LATEX\_NAME\]\... ;

This optional command declares a model local variable. See
`conv`{.interpreted-text role="ref"} for the syntax of VARIABLE\_NAME.
As you can create model local variables on the fly in the model block
(see `model-decl`{.interpreted-text role="ref"}), the interest of this
command is primarily to assign a LATEX\_NAME to the model local
variable.

*Example*

>     model_local_variable GDP_US $GDPUS$;
:::

### On-the-fly Model Variable Declaration {#on-the-fly-declaration}

Endogenous variables, exogenous variables, and parameters can also be
declared inside the model block. You can do this in two different ways:
either via the equation tag or directly in an equation.

To declare a variable on-the-fly in an equation tag, simply state the
type of variable to be declared (`endogenous`, `exogenous`, or
`parameter` followed by an equal sign and the variable name in single
quotes. Hence, to declare a variable `c` as endogenous in an equation
tag, you can type `[endogenous='c']`.

To perform on-the-fly variable declaration in an equation, simply follow
the symbol name with a vertical line (`|`, pipe character) and either an
`e`, an `x`, or a `p`. For example, to declare a parameter named
`alphaa` in the model block, you could write `alphaa|p` directly in an
equation where it appears. Similarly, to declare an endogenous variable
`c` in the model block you could write `c|e`. Note that in-equation
on-the-fly variable declarations must be made on contemporaneous
variables.

On-the-fly variable declarations do not have to appear in the first
place where this variable is encountered.

*Example*

> The following two snippets are equivalent:
>
> >     model;
> >       [endogenous='k',name='law of motion of capital']
> >       k(+1) = i|e + (1-delta|p)*k;
> >       y|e = k^alpha|p;
> >       ...
> >     end;
> >     delta = 0.025;
> >     alpha = 0.36;
> >
> >     var k, i, y;
> >     parameters delta, alpha;
> >     delta = 0.025;
> >     alpha = 0.36;
> >     ...
> >     model;
> >       [name='law of motion of capital']
> >       k(1) = i|e + (1-delta|p)*k;
> >       y|e = k|e^alpha|p;
> >       ...
> >     end;

Expressions {#expr}
-----------

Dynare distinguishes between two types of mathematical expressions:
those that are used to describe the model, and those that are used
outside the model block (e.g. for initializing parameters or variables,
or as command options). In this manual, those two types of expressions
are respectively denoted by MODEL\_EXPRESSION and EXPRESSION.

Unlike MATLAB or Octave expressions, Dynare expressions are necessarily
scalar ones: they cannot contain matrices or evaluate to matrices.[^2]

Expressions can be constructed using integers (INTEGER), floating point
numbers (DOUBLE), parameter names (PARAMETER\_NAME), variable names
(VARIABLE\_NAME), operators and functions.

The following special constants are also accepted in some contexts:

::: {.constant}
inf

Represents infinity.
:::

::: {.constant}
nan

"Not a number": represents an undefined or unrepresentable value.
:::

### Parameters and variables

Parameters and variables can be introduced in expressions by simply
typing their names. The semantics of parameters and variables is quite
different whether they are used inside or outside the model block.

#### Inside the model

Parameters used inside the model refer to the value given through
parameter initialization (see `param-init`{.interpreted-text
role="ref"}) or `homotopy_setup` when doing a simulation, or are the
estimated variables when doing an estimation.

Variables used in a MODEL\_EXPRESSION denote current period values when
neither a lead or a lag is given. A lead or a lag can be given by
enclosing an integer between parenthesis just after the variable name: a
positive integer means a lead, a negative one means a lag. Leads or lags
of more than one period are allowed. For example, if `c` is an
endogenous variable, then `c(+1)` is the variable one period ahead, and
`c(-2)` is the variable two periods before.

When specifying the leads and lags of endogenous variables, it is
important to respect the following convention: in Dynare, the timing of
a variable reflects when that variable is decided. A control variable
--- which by definition is decided in the current period --- must have
no lead. A predetermined variable --- which by definition has been
decided in a previous period --- must have a lag. A consequence of this
is that all stock variables must use the "stock at the end of the
period" convention.

Leads and lags are primarily used for endogenous variables, but can be
used for exogenous variables. They have no effect on parameters and are
forbidden for local model variables (see Model declaration).

#### Outside the model

When used in an expression outside the model block, a parameter or a
variable simply refers to the last value given to that variable. More
precisely, for a parameter it refers to the value given in the
corresponding parameter initialization (see
`param-init`{.interpreted-text role="ref"}); for an endogenous or
exogenous variable, it refers to the value given in the most recent
`initval` or `endval` block.

### Operators

The following operators are allowed in both MODEL\_EXPRESSION and
EXPRESSION:

-   Binary arithmetic operators: `+`, `-`, `*`, `/`, `^`
-   Unary arithmetic operators: `+`, `-`
-   Binary comparison operators (which evaluate to either 0 or 1): `<`,
    `>`, `<=`, `>=`, `==`, `!=`

Note the binary comparison operators are differentiable everywhere
except on a line of the 2-dimensional real plane. However for
facilitating convergence of Newton-type methods, Dynare assumes that, at
the points of non-differentiability, the partial derivatives of these
operators with respect to both arguments is equal to 0 (since this is
the value of the partial derivatives everywhere else).

The following special operators are accepted in MODEL\_EXPRESSION (but
not in EXPRESSION):

::: {.operator}
STEADY\_STATE (MODEL\_EXPRESSION)

This operator is used to take the value of the enclosed expression at
the steady state. A typical usage is in the Taylor rule, where you may
want to use the value of GDP at steady state to compute the output gap.

Exogenous and exogenous deterministic variables may not appear in
MODEL\_EXPRESSION.

::: {.warning}
::: {.admonition-title}
Warning
:::

The concept of a steady state is ambiguous in a perfect foresight
context with permament and potentially anticipated shocks occuring.
Dynare will use the contents of `oo_.steady_state` as its reference for
calls to the `STEADY_STATE()`-operator. In the presence of `endval`,
this implies that the terminal state provided by the user is used. This
may be a steady state computed by Dynare (if `endval` is followed by
`steady`) or simply the terminal state provided by the user (if `endval`
is not followed by `steady`). Put differently, Dynare will not
automatically compute the steady state conditional on the specificed
value of the exogenous variables in the respective periods.
:::
:::

::: {.operator}
EXPECTATION (INTEGER) (MODEL\_EXPRESSION)

This operator is used to take the expectation of some expression using a
different information set than the information available at current
period. For example, `EXPECTATION(-1)(x(+1))` is equal to the expected
value of variable x at next period, using the information set available
at the previous period. See `aux-variables`{.interpreted-text
role="ref"} for an explanation of how this operator is handled
internally and how this affects the output.
:::

### Functions

#### Built-in functions

The following standard functions are supported internally for both
MODEL\_EXPRESSION and EXPRESSION:

::: {.function}
exp(x)

Natural exponential.
:::

::: {.function}
log(x)
:::

::: {.function}
ln(x)

Natural logarithm.
:::

::: {.function}
log10(x)

Base 10 logarithm.
:::

::: {.function}
sqrt(x)

Square root.
:::

::: {.function}
cbrt(x)

Cube root.
:::

::: {.function}
sign(x)

Signum function, defined as:

> $$\begin{aligned}
> \textrm{sign}(x) =
>        \begin{cases}
>        -1 &\quad\text{if }x<0\\
>        0 &\quad\text{if }x=0\\
>        1 &\quad\text{if }x>0
>        \end{cases}
> \end{aligned}$$

Note that this function is not continuous, hence not differentiable, at
$x=0$. However, for facilitating convergence of Newton-type methods,
Dynare assumes that the derivative at $x=0$ is equal to $0$. This
assumption comes from the observation that both the right- and
left-derivatives at this point exist and are equal to $0$, so we can
remove the singularity by postulating that the derivative at $x=0$ is
$0$.
:::

::: {.function}
abs(x)

Absolute value.

Note that this continuous function is not differentiable at $x=0$.
However, for facilitating convergence of Newton-type methods, Dynare
assumes that the derivative at $x=0$ is equal to $0$ (even if the
derivative does not exist). The rational for this mathematically
unfounded definition, rely on the observation that the derivative of
$\mathrm{abs}(x)$ is equal to $\mathrm{sign}(x)$ for any $x\neq 0$ in
$\mathbb R$ and from the convention for the value of $\mathrm{sign}(x)$
at $x=0$).
:::

::: {.function}
sin(x)
:::

::: {.function}
cos(x)
:::

::: {.function}
tan(x)
:::

::: {.function}
asin(x)
:::

::: {.function}
acos(x)
:::

::: {.function}
atan(x)

Trigonometric functions.
:::

::: {.function}
sinh(x)
:::

::: {.function}
cosh(x)
:::

::: {.function}
tanh(x)
:::

::: {.function}
asinh(x)
:::

::: {.function}
acosh(x)
:::

::: {.function}
atanh(x)

Hyperbolic functions.
:::

::: {.function}
max(a, b)
:::

::: {.function}
min(a, b)

Maximum and minimum of two reals.

Note that these functions are differentiable everywhere except on a line
of the 2-dimensional real plane defined by $a=b$. However for
facilitating convergence of Newton-type methods, Dynare assumes that, at
the points of non-differentiability, the partial derivative of these
functions with respect to the first (resp. the second) argument is equal
to $1$ (resp. to $0$) (i.e. the derivatives at the kink are equal to the
derivatives observed on the half-plane where the function is equal to
its first argument).
:::

::: {.function}
normcdf(x) normcdf(x, mu, sigma)

Gaussian cumulative density function, with mean *mu* and standard
deviation *sigma*. Note that `normcdf(x)` is equivalent to
`normcdf(x,0,1)`.
:::

::: {.function}
normpdf(x) normpdf(x, mu, sigma)

Gaussian probability density function, with mean *mu* and standard
deviation *sigma*. Note that `normpdf(x)` is equivalent to
`normpdf(x,0,1)`.
:::

::: {.function}
erf(x)

Gauss error function.
:::

::: {.function}
erfc(x)

Complementary error function, *i.e.*
$\mathrm{erfc}(x) = 1-\mathrm{erf}(x)$.
:::

#### External functions

Any other user-defined (or built-in) MATLAB or Octave function may be
used in both a MODEL\_EXPRESSION and an EXPRESSION, provided that this
function has a scalar argument as a return value.

To use an external function in a MODEL\_EXPRESSION, one must declare the
function using the `external_function` statement. This is not required
for external functions used in an EXPRESSION outside of a `model` block
or `steady_state_model` block.

::: {.command}
external\_function (OPTIONS\...);

This command declares the external functions used in the model block. It
is required for every unique function used in the model block.

`external_function` commands can appear several times in the file and
must come before the model block.

*Options*

::: {.option}
name = NAME

The name of the function, which must also be the name of the M-/MEX file
implementing it. This option is mandatory.
:::

::: {.option}
nargs = INTEGER

The number of arguments of the function. If this option is not provided,
Dynare assumes `nargs = 1`.
:::

::: {.option}
first\_deriv\_provided \[= NAME\]

If NAME is provided, this tells Dynare that the Jacobian is provided as
the only output of the M-/MEX file given as the option argument. If NAME
is not provided, this tells Dynare that the M-/MEX file specified by the
argument passed to NAME returns the Jacobian as its second output
argument. When this option is not provided, Dynare will use finite
difference approximations for computing the derivatives of the function,
whenever needed.
:::

::: {.option}
second\_deriv\_provided \[= NAME\]

If NAME is provided, this tells Dynare that the Hessian is provided as
the only output of the M-/MEX file given as the option argument. If NAME
is not provided, this tells Dynare that the M-/MEX file specified by the
argument passed to NAME returns the Hessian as its third output
argument. NB: This option can only be used if the `first_deriv_provided`
option is used in the same `external_function` command. When this option
is not provided, Dynare will use finite difference approximations for
computing the Hessian derivatives of the function, whenever needed.
:::

*Example*

>     external_function(name = funcname);
>     external_function(name = otherfuncname, nargs = 2, first_deriv_provided, second_deriv_provided);
>     external_function(name = yetotherfuncname, nargs = 3, first_deriv_provided = funcname_deriv);
:::

### A few words of warning in stochastic context

The use of the following functions and operators is strongly discouraged
in a stochastic context: `max`, `min`, `abs`, `sign`, `<`, `>`, `<=`,
`>=`, `==`, `!=`.

The reason is that the local approximation used by `stoch_simul` or
`estimation` will by nature ignore the non-linearities introduced by
these functions if the steady state is away from the kink. And, if the
steady state is exactly at the kink, then the approximation will be
bogus because the derivative of these functions at the kink is bogus (as
explained in the respective documentations of these functions and
operators).

Note that `extended_path` is not affected by this problem, because it
does not rely on a local approximation of the mode.

Parameter initialization {#param-init}
------------------------

When using Dynare for computing simulations, it is necessary to
calibrate the parameters of the model. This is done through parameter
initialization.

The syntax is the following:

    PARAMETER_NAME = EXPRESSION;

Here is an example of calibration:

    parameters alpha, beta;

    beta = 0.99;
    alpha = 0.36;
    A = 1-alpha*beta;

Internally, the parameter values are stored in `M_.params`:

::: {.matvar}
[M]().params

Contains the values of model parameters. The parameters are in the order
that was used in the `parameters` command, hence ordered as in
`M_.param_names`.
:::

The parameter names are stored in `M_.param_names`:

::: {.matvar}
[M]().param\_names

Cell array containing the names of the model parameters.
:::

::: {.matcomm}
get\_param\_by\_name (\'PARAMETER\_NAME\');

Given the name of a parameter, returns its calibrated value as it is
stored in `M_.params`.
:::

::: {.matcomm}
set\_param\_value (\'PARAMETER\_NAME\', MATLAB\_EXPRESSION);

Sets the calibrated value of a parameter to the provided expression.
This does essentially the same as the parameter initialization syntax
described above, except that it accepts arbitrary MATLAB/Octave
expressions, and that it works from MATLAB/Octave scripts.
:::

Model declaration {#model-decl}
-----------------

The model is declared inside a `model` block:

::: {.block}
model ; model (OPTIONS\...);

The equations of the model are written in a block delimited by `model`
and `end` keywords.

There must be as many equations as there are endogenous variables in the
model, except when computing the unconstrained optimal policy with
`ramsey_model`, `ramsey_policy` or `discretionary_policy`.

The syntax of equations must follow the conventions for
MODEL\_EXPRESSION as described in `expr`{.interpreted-text role="ref"}.
Each equation must be terminated by a semicolon (';'). A normal equation
looks like:

> MODEL\_EXPRESSION = MODEL\_EXPRESSION;

When the equations are written in homogenous form, it is possible to
omit the '=0' part and write only the left hand side of the equation. A
homogenous equation looks like:

> MODEL\_EXPRESSION;

Inside the model block, Dynare allows the creation of *model-local
variables*, which constitute a simple way to share a common expression
between several equations. The syntax consists of a pound sign (\#)
followed by the name of the new model local variable (which must **not**
be declared as in `var-decl`{.interpreted-text role="ref"}, but may have
been declared by `model_local_variable`{.interpreted-text role="comm"}),
an equal sign, and the expression for which this new variable will
stand. Later on, every time this variable appears in the model, Dynare
will substitute it by the expression assigned to the variable. Note that
the scope of this variable is restricted to the model block; it cannot
be used outside. To assign a LaTeX name to the model local variable, use
the declaration syntax outlined by
`model_local_variable`{.interpreted-text role="comm"}. A model local
variable declaration looks like:

> \#VARIABLE\_NAME = MODEL\_EXPRESSION;

It is possible to tag equations written in the model block. A tag can
serve different purposes by allowing the user to attach arbitrary
informations to each equation and to recover them at runtime. For
instance, it is possible to name the equations with a `name`-tag, using
a syntax like:

    model;

    [name = 'Budget constraint'];
    c + k = k^theta*A;

    end;

Here, `name` is the keyword indicating that the tag names the equation.
If an equation of the model is tagged with a name, the `resid` command
will display the name of the equations (which may be more informative
than the equation numbers) in addition to the equation number. Several
tags for one equation can be separated using a comma:

    model;

    [name='Taylor rule',mcp = 'r > -1.94478']
    r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;

    end;

More information on tags is available at
<https://git.dynare.org/Dynare/dynare/-/wikis/Equations-Tags>.

There can be several `model` blocks, in which case they are simply
concatenated. The set of effective options is also the concatenation of
the options declared in all the blocks, but in that case you may rather
want to use the `model_options`{.interpreted-text role="comm"} command.

*Options*

::: {.option}
linear

Declares the model as being linear. It spares oneself from having to
declare initial values for computing the steady state of a stationary
linear model. This option can't be used with non-linear models, it will
NOT trigger linearization of the model.
:::

::: {.option}
use\_dll

Instructs the preprocessor to create dynamic loadable libraries (DLL)
containing the model equations and derivatives, instead of writing those
in M-files. You need a working compilation environment, i.e. a working
`mex` command (see `compil-install`{.interpreted-text role="ref"} for
more details). Using this option can result in faster simulations or
estimations, at the expense of some initial compilation time.
Alternatively, this option can be given to the `dynare` command (see
`dyn-invoc`{.interpreted-text role="ref"}).[^3]
:::

::: {.option}
block

Perform the block decomposition of the model, and exploit it in
computations (steady-state, deterministic simulation, stochastic
simulation with first order approximation and estimation). See
<https://archives.dynare.org/DynareWiki/FastDeterministicSimulationAndSteadyStateComputation>
for details on the algorithms used in deterministic simulation and
steady-state computation.
:::

::: {.option}
bytecode

Instead of M-files, use a bytecode representation of the model, i.e. a
binary file containing a compact representation of all the equations.
:::

::: {.option}
cutoff = DOUBLE

Threshold under which a jacobian element is considered as null during
the model normalization. Only available with option `block`. Default:
`1e-15`
:::

::: {.option}
mfs = INTEGER

Controls the handling of minimum feedback set of endogenous variables.
Only available with option `block`. Possible values:

`0`

> All the endogenous variables are considered as feedback variables
> (Default).

`1`

> The endogenous variables assigned to equation naturally normalized
> (i.e. of the form $x=f(Y)$ where $x$ does not appear in $Y$) are
> potentially recursive variables. All the other variables are forced to
> belong to the set of feedback variables.

`2`

> In addition of variables with `mfs = 1` the endogenous variables
> related to linear equations which could be normalized are potential
> recursive variables. All the other variables are forced to belong to
> the set of feedback variables.

`3`

> In addition of variables with `mfs = 2` the endogenous variables
> related to non-linear equations which could be normalized are
> potential recursive variables. All the other variables are forced to
> belong to the set of feedback variables.
:::

::: {.option}
no\_static

Don't create the static model file. This can be useful for models which
don't have a steady state.
:::

::: {.option}
differentiate\_forward\_vars differentiate\_forward\_vars = (
VARIABLE\_NAME \[VARIABLE\_NAME \...\] )

Tells Dynare to create a new auxiliary variable for each endogenous
variable that appears with a lead, such that the new variable is the
time differentiate of the original one. More precisely, if the model
contains `x(+1)`, then a variable `AUX_DIFF_VAR` will be created such
that `AUX_DIFF_VAR=x-x(-1)`, and `x(+1)` will be replaced with
`x+AUX_DIFF_VAR(+1)`.

The transformation is applied to all endogenous variables with a lead if
the option is given without a list of variables. If there is a list, the
transformation is restricted to endogenous with a lead that also appear
in the list.

This option can useful for some deterministic simulations where
convergence is hard to obtain. Bad values for terminal conditions in the
case of very persistent dynamics or permanent shocks can hinder correct
solutions or any convergence. The new differentiated variables have
obvious zero terminal conditions (if the terminal condition is a steady
state) and this in many cases helps convergence of simulations.
:::

::: {.option}
parallel\_local\_files = ( FILENAME \[, FILENAME\]\... )

Declares a list of extra files that should be transferred to follower
nodes when doing a parallel computation (see
`paral-conf`{.interpreted-text role="ref"}).
:::

::: {.option}
balanced\_growth\_test\_tol = DOUBLE

Tolerance used for determining whether cross-derivatives are zero in the
test for balanced growth path (the latter is documented on
<https://archives.dynare.org/DynareWiki/RemovingTrends>). Default:
`1e-6`
:::

*Example* (Elementary RBC model)

>     var c k;
>     varexo x;
>     parameters aa alph bet delt gam;
>
>     model;
>     c =  - k + aa*x*k(-1)^alph + (1-delt)*k(-1);
>     c^(-gam) = (aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam)/(1+bet);
>     end;

*Example* (Use of model local variables)

> The following program:
>
>     model;
>     # gamma = 1 - 1/sigma;
>     u1 = c1^gamma/gamma;
>     u2 = c2^gamma/gamma;
>     end;
>
> \...is formally equivalent to:
>
>     model;
>     u1 = c1^(1-1/sigma)/(1-1/sigma);
>     u2 = c2^(1-1/sigma)/(1-1/sigma);
>     end;

*Example* (A linear model)

>     model(linear);
>     x = a*x(-1)+b*y(+1)+e_x;
>     y = d*y(-1)+e_y;
>     end;
:::

::: {.command}
model\_options (OPTIONS\...);

This command accepts the same options as the `model`{.interpreted-text
role="bck"} block.

The purpose of this statement is to specify the options that apply to
the whole model, when there are several `model` blocks, so as to restore
the symmetry between those blocks (since otherwise one `model` block
would typically bear the options, while the other ones would typically
have no option).
:::

::: {.command}
model\_remove (TAGS\...);

This command removes equations that appeared in a previous
`model`{.interpreted-text role="bck"} block.

The equations must be specified by a list of tag values, separated by
commas. Each element of the list is either a simple quoted string, in
which case it designates an equation by its `name` tag; or a tag name
(without quotes), followed by an equal sign, then by the tag value
(within quotes).

Each removed equation must either have an `endogenous` tag, or have a
left hand side containing a single endogenous variable. The
corresponding endogenous variable will be either turned into an
exogenous (if it is still used in somewhere in the model at that point),
otherwise it will be removed from the model.

*Example*

>     var c k dummy1 dummy2;
>
>     model;
>       c + k - aa*x*k(-1)^alph - (1-delt)*k(-1) + dummy1;
>       c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
>       [ name = 'eq:dummy1', endogenous = 'dummy1' ]
>       c*k = dummy1;
>       [ foo = 'eq:dummy2' ]
>       log(dummy2) = k + 2;
>     end;
>
>     model_remove('eq:dummy1', foo = 'eq:dummy2');
>
> In the above example, the last two equations will be removed, `dummy1`
> will be turned into an exogenous, and `dummy2` will be removed.
:::

::: {.block}
model\_replace (TAGS\...);

This block replaces several equations in the model. It removes the
equations given by the tags list (with the same syntax as in
`model_remove`{.interpreted-text role="comm"}), and it adds equations
given within the block (with the same syntax as
`model`{.interpreted-text role="bck"}).

No variable is removed or has its type changed in the process.

*Example*

>     var c k;
>
>     model;
>       c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
>       [ name = 'dummy' ]
>       c*k = 1;
>     end;
>
>     model_replace('dummy');
>       c^(-gam) = (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
>     end;
>
> In the above example, the dummy equation is replaced by a proper Euler
> equation.
:::

Dynare has the ability to output the original list of model equations to
a LaTeX file, using the `write_latex_original_model` command, the list
of transformed model equations using the
`write_latex_dynamic_model command`, and the list of static model
equations using the `write_latex_static_model` command.

::: {.command}
write\_latex\_original\_model (OPTIONS);

This command creates two LaTeX files: one containing the model as
defined in the model block and one containing the LaTeX document header
information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/original.tex`, which includes a file called
`FILENAME/latex/original_content.tex` (also created by Dynare)
containing the list of all the original model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Time subscripts (`t`, `t+1`, `t-1`, \...) will be appended to the
variable names, as LaTeX subscripts.

Compiling the TeX file requires the following LaTeX packages:
`geometry, fullpage, breqn`.

*Options*

::: {.option}
write\_equation\_tags

Write the equation tags in the LaTeX output. The equation tags will be
interpreted with LaTeX markups.
:::
:::

::: {.command}
write\_latex\_dynamic\_model ; write\_latex\_dynamic\_model (OPTIONS);

This command creates two LaTeX files: one containing the dynamic model
and one containing the LaTeX document header information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/dynamic.tex`, which includes a file called
`FILENAME/latex/dynamic_content.tex` (also created by Dynare) containing
the list of all the dynamic model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Time subscripts (`t`, `t+1`, `t-1`, \...) will be appended to the
variable names, as LaTeX subscripts.

Note that the model written in the TeX file will differ from the model
declared by the user in the following dimensions:

> -   The timing convention of predetermined variables (see
>     `predetermined_variables`{.interpreted-text role="comm"}) will
>     have been changed to the default Dynare timing convention; in
>     other words, variables declared as predetermined will be lagged on
>     period back,
> -   The `EXPECTATION` operators will have been removed, replaced by
>     auxiliary variables and new equations (as explained in the
>     documentation of
>     `EXPECTATION <EXPECTATION (INTEGER) (MODEL_EXPRESSION)>`{.interpreted-text
>     role="op"}),
> -   Endogenous variables with leads or lags greater or equal than two
>     will have been removed, replaced by new auxiliary variables and
>     equations,
> -   Exogenous variables with leads or lags will also have been
>     replaced by new auxiliary variables and equations.

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.

*Options*

::: {.option}
write\_equation\_tags

See `write_equation_tags`{.interpreted-text role="opt"}
:::
:::

::: {.command}
write\_latex\_static\_model (OPTIONS);

This command creates two LaTeX files: one containing the static model
and one containing the LaTeX document header information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/static.tex`, which includes a file called
`FILENAME/latex/static_content.tex` (also created by Dynare) containing
the list of all the steady state model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Note that the model written in the TeX file will differ from the model
declared by the user in the some dimensions (see
`write_latex_dynamic_model`{.interpreted-text role="comm"} for details).

Also note that this command will not output the contents of the optional
`steady_state_model` block (see `steady_state_model`{.interpreted-text
role="bck"}); it will rather output a static version (i.e. without leads
and lags) of the dynamic `model` declared in the model block. To write
the LaTeX contents of the `steady_state_model` see
`write_latex_steady_state_model`{.interpreted-text role="comm"}.

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.

*Options*

::: {.option}
write\_equation\_tags

See `write_equation_tags`{.interpreted-text role="opt"}.
:::
:::

::: {.command}
write\_latex\_steady\_state\_model

This command creates two LaTeX files: one containing the steady state
model and one containing the LaTeX document header information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/steady_state.tex`, which includes a file called
`FILENAME/latex/steady_state_content.tex` (also created by Dynare)
containing the list of all the steady state model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Note that the model written in the `.tex` file will differ from the
model declared by the user in some dimensions (see
`write_latex_dynamic_model`{.interpreted-text role="comm"} for details).

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.
:::

Auxiliary variables {#aux-variables}
-------------------

The model which is solved internally by Dynare is not exactly the model
declared by the user. In some cases, Dynare will introduce auxiliary
endogenous variables---along with corresponding auxiliary
equations---which will appear in the final output.

The main transformation concerns leads and lags. Dynare will perform a
transformation of the model so that there is only one lead and one lag
on endogenous variables and no leads/lags on exogenous variables.

This transformation is achieved by the creation of auxiliary variables
and corresponding equations. For example, if `x(+2)` exists in the
model, Dynare will create one auxiliary variable
`AUX_ENDO_LEAD = x(+1)`, and replace `x(+2)` by `AUX_ENDO_LEAD(+1)`.

A similar transformation is done for lags greater than 2 on endogenous
(auxiliary variables will have a name beginning with `AUX_ENDO_LAG`),
and for exogenous with leads and lags (auxiliary variables will have a
name beginning with `AUX_EXO_LEAD` or `AUX_EXO_LAG` respectively).

Another transformation is done for the `EXPECTATION` operator. For each
occurrence of this operator, Dynare creates an auxiliary variable
defined by a new equation, and replaces the expectation operator by a
reference to the new auxiliary variable. For example, the expression
`EXPECTATION(-1)(x(+1))` is replaced by `AUX_EXPECT_LAG_1(-1)`, and the
new auxiliary variable is declared as `AUX_EXPECT_LAG_1 = x(+2)`.

Auxiliary variables are also introduced by the preprocessor for the
`ramsey_model` and `ramsey_policy` commands. In this case, they are used
to represent the Lagrange multipliers when first order conditions of the
Ramsey problem are computed. The new variables take the form `MULT_i`,
where *i* represents the constraint with which the multiplier is
associated (counted from the order of declaration in the model block).

Auxiliary variables are also introduced by the
`differentiate_forward_vars` option of the model block. The new
variables take the form `AUX_DIFF_FWRD_i`, and are equal to `x-x(-1)`
for some endogenous variable `x`.

Finally, auxiliary variables will arise in the context of employing the
`diff`-operator.

Once created, all auxiliary variables are included in the set of
endogenous variables. The output of decision rules (see below) is such
that auxiliary variable names are replaced by the original variables
they refer to.

The number of endogenous variables before the creation of auxiliary
variables is stored in `M_.orig_endo_nbr`, and the number of endogenous
variables after the creation of auxiliary variables is stored in
`M_.endo_nbr`.

See <https://git.dynare.org/Dynare/dynare/-/wikis/Auxiliary-variables>
for more technical details on auxiliary variables.

Initial and terminal conditions {#init-term-cond}
-------------------------------

For most simulation exercises, it is necessary to provide initial (and
possibly terminal) conditions. It is also necessary to provide initial
guess values for non-linear solvers. This section describes the
statements used for those purposes.

In many contexts (deterministic or stochastic), it is necessary to
compute the steady state of a non-linear model: `initval` then specifies
numerical initial values for the non-linear solver. The command `resid`
can be used to compute the equation residuals for the given initial
values.

Used in perfect foresight mode, the types of forward-looking models for
which Dynare was designed require both initial and terminal conditions.
Most often these initial and terminal conditions are static equilibria,
but not necessarily.

One typical application is to consider an economy at the equilibrium at
time 0, trigger a shock in first period, and study the trajectory of
return to the initial equilibrium. To do that, one needs `initval` and
`shocks` (see `shocks-exo`{.interpreted-text role="ref"}).

Another one is to study how an economy, starting from arbitrary initial
conditions at time 0 converges towards equilibrium. In this case models,
the command `histval` permits to specify different historical initial
values for variables with lags for the periods before the beginning of
the simulation. Due to the design of Dynare, in this case `initval` is
used to specify the terminal conditions.

::: {.block}
initval ; initval(OPTIONS\...);

The `initval` block has two main purposes: providing guess values for
non-linear solvers in the context of perfect foresight simulations and
providing guess values for steady state computations in both perfect
foresight and stochastic simulations. Depending on the presence of
`histval` and `endval` blocks it is also used for declaring the initial
and terminal conditions in a perfect foresight simulation exercise.
Because of this interaction of the meaning of an `initval` block with
the presence of `histval` and `endval` blocks in perfect foresight
simulations, it is strongly recommended to check that the constructed
`oo_.endo_simul` and `oo_.exo_simul` variables contain the desired
values after running `perfect_foresight_setup` and before running
`perfect_foresight_solver`. In the presence of leads and lags, these
subfields of the results structure will store the historical values for
the lags in the first column/row and the terminal values for the leads
in the last column/row.

The `initval` block is terminated by `end;` and contains lines of the
form:

> VARIABLE\_NAME = EXPRESSION;

*In a deterministic (i.e. perfect foresight) model*

First, both the `oo_.endo_simul` and `oo_.exo_simul` variables storing
the endogenous and exogenous variables will be filled with the values
provided by this block. If there are no other blocks present, it will
therefore provide the initial and terminal conditions for all the
endogenous and exogenous variables, because it will also fill the last
column/row of these matrices. For the intermediate simulation periods it
thereby provides the starting values for the solver. In the presence of
a `histval` block (and therefore absence of an `endval` block), this
`histval` block will provide/overwrite the historical values for the
state variables (lags) by setting the first column/row of
`oo_.endo_simul` and `oo_.exo_simul`. This implies that the `initval`
block in the presence of `histval` only sets the terminal values for the
variables with leads and provides initial values for the perfect
foresight solver.

Because of these various functions of `initval` it is often necessary to
provide values for all the endogenous variables in an `initval` block.
Initial and terminal conditions are strictly necessary for lagged/leaded
variables, while feasible starting values are required for the solver.
It is important to be aware that if some variables, endogenous or
exogenous, are not mentioned in the `initval` block, a zero value is
assumed. It is particularly important to keep this in mind when
specifying exogenous variables using `varexo` that are not allowed to
take on the value of zero, like e.g. TFP.

Note that if the `initval` block is immediately followed by a `steady`
command, its semantics are slightly changed. The `steady` command will
compute the steady state of the model for all the endogenous variables,
assuming that exogenous variables are kept constant at the value
declared in the `initval` block. These steady state values conditional
on the declared exogenous variables are then written into
`oo_.endo_simul` and take up the potential roles as historical and
terminal conditions as well as starting values for the solver. An
`initval` block followed by `steady` is therefore formally equivalent to
an `initval` block with the specified values for the exogenous
variables, and the endogenous variables set to the associated steady
state values conditional on the exogenous variables.

*In a stochastic model*

The main purpose of `initval` is to provide initial guess values for the
non-linear solver in the steady state computation. Note that if the
`initval` block is not followed by `steady`, the steady state
computation will still be triggered by subsequent commands
(`stoch_simul`, `estimation`\...).

As such, `initval` allows specifying the initial instrument value for
steady state finding when providing an analytical conditional steady
state file for `ramsey_model`-computations.

It is not necessary to declare 0 as initial value for exogenous
stochastic variables, since it is the only possible value.

The subsequently computed steady state (not the initial values, use
histval for this) will be used as the initial condition at all the
periods preceeding the first simulation period for the three possible
types of simulations in stochastic mode:

> -   `stoch_simul`{.interpreted-text role="comm"}, if the `periods`
>     option is specified.
> -   `forecast`{.interpreted-text role="comm"} as the initial point at
>     which the forecasts are computed.
> -   `conditional_forecast`{.interpreted-text role="comm"} as the
>     initial point at which the conditional forecasts are computed.

To start simulations at a particular set of starting values that are not
a computed steady state, use `histval`{.interpreted-text role="bck"}.

*Options*

::: {.option}
all\_values\_required

Issues an error and stops processing the .mod file if there is at least
one endogenous or exogenous variable that has not been set in the
initval block.
:::

*Example*

:   initval;
        c = 1.2;
        k = 12;
        x = 1;
        end;

        steady;
:::

::: {.block}
endval ; endval (OPTIONS\...);

This block is terminated by `end;` and contains lines of the form:

> VARIABLE\_NAME = EXPRESSION;

The `endval` block makes only sense in a deterministic model and cannot
be used together with `histval`. Similar to the `initval` command, it
will fill both the `oo_.endo_simul` and `oo_.exo_simul` variables
storing the endogenous and exogenous variables with the values provided
by this block. If no `initval` block is present, it will fill the whole
matrices, therefore providing the initial and terminal conditions for
all the endogenous and exogenous variables, because it will also fill
the first and last column/row of these matrices. Due to also filling the
intermediate simulation periods it will provide the starting values for
the solver as well.

If an `initval` block is present, `initval` will provide the historical
values for the variables (if there are states/lags), while `endval` will
fill the remainder of the matrices, thereby still providing *i*) the
terminal conditions for variables entering the model with a lead and
*ii*) the initial guess values for all endogenous variables at all the
simulation dates for the perfect foresight solver.

Note that if some variables, endogenous or exogenous, are NOT mentioned
in the `endval` block, the value assumed is that of the last `initval`
block or `steady` command (if present). Therefore, in contrast to
`initval`, omitted variables are not automatically assumed to be 0 in
this case. Again, it is strongly recommended to check the constructed
`oo_.endo_simul` and `oo_.exo_simul` variables after running
`perfect_foresight_setup` and before running `perfect_foresight_solver`
to see whether the desired outcome has been achieved.

Like `initval`, if the `endval` block is immediately followed by a
`steady` command, its semantics are slightly changed. The `steady`
command will compute the steady state of the model for all the
endogenous variables, assuming that exogenous variables are kept
constant to the value declared in the `endval` block. These steady state
values conditional on the declared exogenous variables are then written
into `oo_.endo_simul` and therefore take up the potential roles as
historical and terminal conditions as well as starting values for the
solver. An `endval` block followed by `steady` is therefore formally
equivalent to an `endval` block with the specified values for the
exogenous variables, and the endogenous variables set to the associated
steady state values.

*Options*

::: {.option}
all\_values\_required

See `all_values_required`{.interpreted-text role="opt"}.
:::

*Example*

>     var c k;
>     varexo x;
>
>     model;
>     c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
>     c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
>     end;
>
>     initval;
>     c = 1.2;
>     k = 12;
>     x = 1;
>     end;
>
>     steady;
>
>     endval;
>     c = 2;
>     k = 20;
>     x = 2;
>     end;
>
>     steady;
>
>     perfect_foresight_setup(periods=200);
>     perfect_foresight_solver;
>
> In this example, the problem is finding the optimal path for
> consumption and capital for the periods $t=1$ to $T=200$, given the
> path of the exogenous technology level `x`. `c` is a forward-looking
> variable and the exogenous variable `x` appears with a lead in the
> expected return of physical capital, while `k` is a purely
> backward-looking (state) variable.
>
> The initial equilibrium is computed by `steady` conditional on `x=1`,
> and the terminal one conditional on `x=2`. The `initval` block sets
> the initial condition for `k` (since it is the only backward-looking
> variable), while the `endval` block sets the terminal condition for
> `c` (since it is the only forward-looking endogenous variable). The
> starting values for the perfect foresight solver are given by the
> `endval` block. See below for more details.

*Example*

>     var c k;
>     varexo x;
>
>     model;
>     c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
>     c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
>     end;
>
>     initval;
>     k = 12;
>     end;
>
>     endval;
>     c = 2;
>     x = 1.1;
>     end;
>
>     perfect_foresight_setup(periods=200);
>     perfect_foresight_solver;
>
> In this example, there is no [steady]{.title-ref} command, hence the
> conditions are exactly those specified in the [initval]{.title-ref}
> and [endval]{.title-ref} blocks. We need terminal conditions for `c`
> and `x`, since both appear with a lead, and an initial condition for
> `k`, since it appears with a lag.
>
> Setting `x=1.1` in the `endval` block without a `shocks` block implies
> that technology is at $1.1$ in $t=1$ and stays there forever, because
> `endval` is filling all entries of `oo_.endo_simul` and
> `oo_.exo_simul` except for the very first one, which stores the
> initial conditions and was set to $0$ by the `initval` block when not
> explicitly specifying a value for it.
>
> Because the law of motion for capital is backward-looking, we need an
> initial condition for `k` at time $0$. Due to the presence of
> `endval`, this cannot be done via a `histval` block, but rather must
> be specified in the `initval` block. Similarly, because the Euler
> equation is forward-looking, we need a terminal condition for `c` at
> $t=201$, which is specified in the `endval` block.
>
> As can be seen, it is not necessary to specify `c` and `x` in the
> `initval` block and `k` in the `endval` block, because they have no
> impact on the results. Due to the optimization problem in the first
> period being to choose `c,k` at $t=1$ given the predetermined capital
> stock `k` inherited from $t=0$ as well as the current and future
> values for technology `x`, the values for `c` and `x` at time $t=0$
> play no role. The same applies to the choice of `c,k` at time $t=200$,
> which does not depend on `k` at $t=201$. As the Euler equation shows,
> that choice only depends on current capital as well as future
> consumption `c` and technology `x`, but not on future capital `k`. The
> intuitive reason is that those variables are the consequence of
> optimization problems taking place in at periods $t=0$ and $t=201$,
> respectively, which are not modeled here.

*Example*

>     initval;
>     c = 1.2;
>     k = 12;
>     x = 1;
>     end;
>
>     endval;
>     c = 2;
>     k = 20;
>     x = 1.1;
>     end;
>
> In this example, initial conditions for the forward-looking variables
> `x` and `c` are provided, together with a terminal condition for the
> backward-looking variable `k`. As shown in the previous example, these
> values will not affect the simulation results. Dynare simply takes
> them as given and basically assumes that there were realizations of
> exogenous variables and states that make those choices equilibrium
> values (basically initial/terminal conditions at the unspecified time
> periods $t<0$ and $t>201$).
>
> The above example suggests another way of looking at the use of
> `steady` after `initval` and `endval`. Instead of saying that the
> implicit unspecified conditions before and after the simulation range
> have to fit the initial/terminal conditions of the endogenous
> variables in those blocks, steady specifies that those conditions at
> $t<0$ and $t>201$ are equal to being at the steady state given the
> exogenous variables in the `initval` and `endval` blocks. The
> endogenous variables at $t=0$ and $t=201$ are then set to the
> corresponding steady state equilibrium values.
>
> The fact that `c` at $t=0$ and `k` at $t=201$ specified in `initval`
> and `endval` are taken as given has an important implication for
> plotting the simulated vector for the endogenous variables, i.e. the
> rows of `oo_.endo_simul`: this vector will also contain the initial
> and terminal conditions and thus is 202 periods long in the example.
> When you specify arbitrary values for the initial and terminal
> conditions for forward- and backward-looking variables, respectively,
> these values can be very far away from the endogenously determined
> values at $t=1$ and $t=200$. While the values at $t=0$ and $t=201$ are
> unrelated to the dynamics for $0<t<201$, they may result in
> strange-looking large jumps. In the example above, consumption will
> display a large jump from $t=0$ to $t=1$ and capital will jump from
> $t=200$ to $t=201$ when using `rplot`{.interpreted-text role="comm"}
> or manually plotting `oo_.endo_val`.
:::

::: {.block}
histval ; histval (OPTIONS\...);

*In a deterministic perfect foresight context*

In models with lags on more than one period, the `histval` block permits
to specify different historical initial values for different periods of
the state variables. In this case, the `initval` block takes over the
role of specifying terminal conditions and starting values for the
solver. Note that the `histval` block does not take non-state variables.

This block is terminated by `end;` and contains lines of the form:

> VARIABLE\_NAME(INTEGER) = EXPRESSION;

EXPRESSION is any valid expression returning a numerical value and can
contain already initialized variable names.

By convention in Dynare, period 1 is the first period of the simulation.
Going backward in time, the first period before the start of the
simulation is period 0, then period -1, and so on.

State variables not initialized in the `histval` block are assumed to
have a value of zero at period 0 and before. Note that `histval` cannot
be followed by `steady`.

*Example*

>     model;
>     x=1.5*x(-1)-0.6*x(-2)+epsilon;
>     log(c)=0.5*x+0.5*log(c(+1));
>     end;
>
>     histval;
>     x(0)=-1;
>     x(-1)=0.2;
>     end;
>
>     initval;
>     c=1;
>     x=1;
>     end;
>
> In this example, `histval` is used to set the historical conditions
> for the two lags of the endogenous variable `x`, stored in the first
> column of `oo_.endo_simul`. The `initval` block is used to set the
> terminal condition for the forward looking variable `c`, stored in the
> last column of `oo_.endo_simul`. Moreover, the `initval` block defines
> the starting values for the perfect foresight solver for both
> endogenous variables `c` and `x`.

*In a stochastic simulation context*

In the context of stochastic simulations, `histval` allows setting the
starting point of those simulations in the state space. As for the case
of perfect foresight simulations, all not explicitly specified variables
are set to 0. Moreover, as only states enter the recursive policy
functions, all values specified for control variables will be ignored.
This can be used

> -   In `stoch_simul`{.interpreted-text role="comm"}, if the `periods`
>     option is specified. Note that this only affects the starting
>     point for the simulation, but not for the impulse response
>     functions. When using the `loglinear <logl>`{.interpreted-text
>     role="ref"} option, the `histval` block nevertheless takes the
>     unlogged starting values.
> -   In `forecast`{.interpreted-text role="comm"} as the initial point
>     at which the forecasts are computed. When using the `loglinear
>     <logl>`{.interpreted-text role="ref"} option, the `histval` block
>     nevertheless takes the unlogged starting values.
> -   In `conditional_forecast`{.interpreted-text role="comm"} for a
>     calibrated model as the initial point at which the conditional
>     forecasts are computed. When using the
>     `loglinear <logl>`{.interpreted-text role="ref"} option, the
>     histval-block nevertheless takes the unlogged starting values.
> -   In `Ramsey policy <ramsey_model>`{.interpreted-text role="comm"},
>     where it also specifies the values of the endogenous states
>     (including lagged exogenous) at which the objective function of
>     the planner is computed. Note that the initial values of the
>     Lagrange multipliers associated with the planner's problem cannot
>     be set (see `evaluate_planner_objective`{.interpreted-text
>     role="comm"}).

*Options*

::: {.option}
all\_values\_required

See `all_values_required`{.interpreted-text role="opt"}.
:::

*Example*

>     var x y;
>     varexo e;
>
>     model;
>     x = y(-1)^alpha*y(-2)^(1-alpha)+e;
>
>     end;
>
>     initval;
>     x = 1;
>     y = 1;
>     e = 0.5;
>     end;
>
>     steady;
>
>     histval;
>     y(0) = 1.1;
>     y(-1) = 0.9;
>     end;
>
>     stoch_simul(periods=100);
:::

::: {.command}
resid ;

This command will display the residuals of the static equations of the
model, using the values given for the endogenous in the last `initval`
or `endval` block (or the steady state file if you provided one, see
`st-st`{.interpreted-text role="ref"}).

*Options*

::: {.option}
non\_zero

Only display non-zero residuals.
:::
:::

::: {.command}
initval\_file (OPTIONS\...);

In a deterministic setup, this command is used to specify a path for all
endogenous and exogenous variables. The length of these paths must be
equal to the number of simulation periods, plus the number of leads and
the number of lags of the model (for example, with 50 simulation
periods, in a model with 2 lags and 1 lead, the paths must have a length
of 53). Note that these paths cover two different things:

> -   The constraints of the problem, which are given by the path for
>     exogenous and the initial and terminal values for endogenous
> -   The initial guess for the non-linear solver, which is given by the
>     path for endogenous variables for the simulation periods
>     (excluding initial and terminal conditions)

In perfect foresight and stochastic contexts, `steady` uses the first
observation loaded by `initval_file` as guess value to solve for the
steady state of the model. This first observation is determined by the
`first_obs` option when it is used.

Don't mix `initval_file` with `initval` statements. However, after
`initval_file`, you can modify the historical initial values with
`histval` or `histval_file` statement.

There can be several `initval_file` statements in a model file. Each
statement resets `oo_.initval_series`.

*Options*

::: {.option}
datafile = FILENAME filename = FILENAME (deprecated)

The name of the file containing the data. It must be included in quotes
if the filename contains a path or an extension. The command accepts the
following file formats:

-   M-file (extension `.m`): for each endogenous and exogenous variable,
    the file must contain a row or column vector of the same name.
-   MAT-file (extension `.mat`): same as for M-files.
-   Excel file (extension `.xls` or `.xlsx`): for each endogenous and
    exogenous variable, the file must contain a column of the same name.
    NB: Octave only supports the `.xlsx` file extension and must have
    the [io](https://octave.sourceforge.io/io/) package installed
    (easily done via octave by typing '`pkg install -forge io`'). The
    first column may contain the date
:::

> of each observation.
>
> :   -   CSV files (extension `.csv`): for each endogenous and
>         exogenous variable, the file must contain a column of the same
>         name. The first column may contain the date of each
>         observation.
>
::: {.option}
first\_obs = {INTEGER \| DATE}

The observation number or the date (see
:::

`dates-members`{.interpreted-text role="ref"}) of the first observation
to be used in the file

::: {.option}
first\_simulation\_period = {INTEGER \| DATE}
:::

The observation number in the file or the date (see
`dates <dates-members>`{.interpreted-text role="ref"}) at which the
simulation (or the forecast) is starting. This option avoids to have to
compute the maximum number of lags in the model. The observation
corresponding to the first period of simulation doesn't need to exist in
the file as the only dates necessary for initialization are before that
date.

::: {.option}
last\_simulation\_period = {INTEGER \| DATE}
:::

The observation number in the file or the date (see
`dates <dates-members>`{.interpreted-text role="ref"}) at which the
simulation (or the forecast) is ending. This option avoids to have to
compute the maximum number of leads in the model.

::: {.option}
last\_obs = {INTEGER \| DATE}

The observaton number or the date (see
:::

`dates-members`{.interpreted-text role="ref"}) of the last observation
to be used in the file.

::: {.option}
nobs = INTEGER
:::

The number of observations to be used in the file (starting with first
of `first_obs` observation).

::: {.option}
series = DSERIES NAME

The name of a DSERIES containing the data (see
`dseries-members`{.interpreted-text role="ref"})
:::

*Example 1*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

The initial and terminal values are taken from file `mydata.csv`
(nothing guarantees that these vales are the steady state of the model).
The guess value for the trajectories are also taken from the file. The
file must contain at least 203 observations of variables `c`, `x` and
`e`. If there are more than 203 observations available in the file, the
first 203 are used by `perfect_foresight_setup(periods=200)`. Note that
the values for the auxiliary variable corresponding to `x(-2)` are
automatically computed by `initval_file`.

*Example 2*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   first\_obs=10);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

The initial and terminal values are taken from file `mydata.csv`
starting with the 10th observation in the file. There must be at least
212 observations in the file.

*Example 3*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> ds = dseries(mydata.csv); lds = log(ds);
>
> initval\_file(series=lds,
>
> :   first\_obs=2010Q1);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

The initial and terminal values are taken from dseries `lds`. All
observations are loaded starting with the 1st quarter of 2010 until the
end of the file. There must be data available at least until 2050Q3.

*Example 4*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   first\_simulation\_period=2010Q1);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

The initial and terminal values are taken from file `mydata.csv`. The
observations in the file must have dates. All observations are loaded
from the 3rd quarter of 2009 until the end of the file. There must be
data available in the file at least until 2050Q1.

*Example 5*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   last\_obs = 212);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

The initial and terminal values are taken from file `mydata.csv`. The
first 212 observations are loaded and the first 203 observations will be
used by `perfect_foresight_setup(periods=200)`.

*Example 6*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   first\_obs = 10,
>
> > nobs = 203);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

The initial and terminal values are taken from file `mydata.csv`.
Observations 10 to 212 are loaded.

*Example 7*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   first\_obs = 10);
>
> steady;
>
> The values of the 10th observation of `mydata.csv` are used

as guess value to compute the steady state. The exogenous variables are
set to values found in the file or zero if these variables aren\'t
present.
:::

::: {.command}
histval\_file (OPTIONS\...);

This command is equivalent to `histval`, except that it reads its input
from a file, and is typically used in conjunction with
`smoother2histval`.

> *Options*

::: {.option}
datafile = FILENAME filename = FILENAME (deprecated)

The name of the file containing the data. The command accepts
:::

the following file formats:

> -   M-file (extension `.m`): for each endogenous and exogenous
>     variable, the file must contain a row or column vector of the same
>     name.
> -   MAT-file (extension `.mat`): same as for M-files.
> -   Excel file (extension `.xls` or `.xlsx`): for each endogenous and
>     exogenous variable, the file must contain a column of the same
>     name. NB: Octave only supports the `.xlsx` file extension and must
>     have the [io](https://octave.sourceforge.io/io/) package installed
>     (easily done via octave by typing '`pkg install -forge io`'). The
>     first column may contain the date of each observation.
> -   CSV files (extension `.csv`): for each endogenous and exogenous
>     variable, the file must contain a column of the same name. The
>     first column may contain the date of each observation.

::: {.option}
first\_obs = {INTEGER \| DATE}

The observation number or the date (see
`dates-members`{.interpreted-text role="ref"}) of
:::

the first observation to be used in the file

::: {.option}
first\_simulation\_period = {INTEGER \| DATE}
:::

The observation number in the file or the date (see
`dates-members`{.interpreted-text role="ref"}) at which the simulation
(or the forecast) is starting. This option avoids to have to compute the
maximum number of lags in the model. The observation corresponding to
the first period of simulation doesn't need to exist in the file as the
only dates necessary for initialization are before that date.

::: {.option}
last\_simulation\_period = {INTEGER \| DATE}
:::

The observation number in the file or the date (see
`dates <dates-members>`{.interpreted-text role="ref"}) at which the
simulation (or the forecast) is ending. This option avoids to have to
compute the maximum number of leads in the model.

::: {.option}
last\_obs = {INTEGER \| DATE}

The observation number or the date (see
`dates-members`{.interpreted-text role="ref"}) of the
:::

last observation to be used in the file.

::: {.option}
nobs = INTEGER
:::

The number of observations to be used in the file (starting with first
of `first_obs` observation).

::: {.option}
series = DSERIES NAME

The name of a DSERIES containing the data (see
`dseries-members`{.interpreted-text role="ref"})
:::

*Example 1*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> steady\_state\_model; x = 0; c = exp(c\*x/(1 - d)); end;
>
> histval\_file(datafile=mydata.csv);
>
> stoch\_simul(order=1,periods=100);

The initial values for the stochastic simulation are taken from the two
first rows of file `mydata.csv`.

*Example 2*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> histval\_file(datafile=mydata.csv,
>
> :   first\_obs=10);
>
> stoch\_simul(order=1,periods=100);

The initial values for the stochastic simulation are taken from rows 10
and 11 of file `mydata.csv`.

*Example 3*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> histval\_file(datafile=mydata.csv,
>
> :   first\_obs=2010Q1);
>
> stoch\_simul(order=1,periods=100);

The initial values for the stochastic simulation are taken from
observations 2010Q1 and 2010Q2 of file `mydata.csv`.

*Example 4*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> histval\_file(datafile=mydata.csv,
>
> :   first\_simulation\_period=2010Q1)
>
> stoch\_simul(order=1,periods=100);

The initial values for the stochastic simulation are taken from
observations 2009Q3 and 2009Q4 of file `mydata.csv`.

*Example 5*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> histval\_file(datafile=mydata.csv,
>
> :   last\_obs = 4);
>
> stoch\_simul(order=1,periods=100);

The initial values for the stochastic simulation are taken from the two
first rows of file `mydata.csv`.

*Example 6*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   first\_obs = 10,
>
> > nobs = 4);
>
> stoch\_simul(order=1,periods=100);

The initial values for the stochastic simulation are taken from rows 10
and 11 of file `mydata.csv`.

*Example 7*

>     var c x;
>     varexo e;
>
> parameters a b c d;
>
> a = 1.5; b = -0,6; c = 0.5; d = 0.5;
>
> > model; x = a\*x(-1) + b\*x(-2) + e; log(c) = c\*x + d\*log(c(+1));
> > end;
>
> initval\_file(datafile=mydata.csv,
>
> :   first\_obs=10);
>
> > histval\_file(datafile=myotherdata.csv);
>
> perfect\_foresight\_setup(periods=200); perfect\_foresight\_solver;

Historical initial values for the simulation are taken from the two
first rows of file `myotherdata.csv`.

> Terminal values and guess values for the simulation are taken

from file `mydata.csv` starting with the 12th observation in the file.
There must be at least 212 observations in the file.
:::

Shocks on exogenous variables {#shocks-exo}
-----------------------------

In a deterministic context, when one wants to study the transition of
one equilibrium position to another, it is equivalent to analyze the
consequences of a permanent shock and this in done in Dynare through the
proper use of `initval` and `endval`.

Another typical experiment is to study the effects of a temporary shock
after which the system goes back to the original equilibrium (if the
model is stable\...). A temporary shock is a temporary change of value
of one or several exogenous variables in the model. Temporary shocks are
specified with the command `shocks`.

In a stochastic framework, the exogenous variables take random values in
each period. In Dynare, these random values follow a normal distribution
with zero mean, but it belongs to the user to specify the variability of
these shocks. The non-zero elements of the matrix of variance-covariance
of the shocks can be entered with the `shocks` command. Or, the entire
matrix can be directly entered with `Sigma_e` (this use is however
deprecated).

If the variance of an exogenous variable is set to zero, this variable
will appear in the report on policy and transition functions, but isn't
used in the computation of moments and of Impulse Response Functions.
Setting a variance to zero is an easy way of removing an exogenous
shock.

Note that, by default, if there are several `shocks` or `mshocks` blocks
in the same `.mod` file, then they are cumulative: all the shocks
declared in all the blocks are considered; however, if a `shocks` or
`mshocks` block is declared with the `overwrite` option, then it
replaces all the previous `shocks` and `mshocks` blocks.

::: {.block}
shocks ; shocks(overwrite);

See above for the meaning of the `overwrite` option.

*In deterministic context*

For deterministic simulations, the `shocks` block specifies temporary
changes in the value of exogenous variables. For permanent shocks, use
an `endval` block.

The block should contain one or more occurrences of the following group
of three lines:

    var VARIABLE_NAME;
    periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
    values DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;

It is possible to specify shocks which last several periods and which
can vary over time. The `periods` keyword accepts a list of several
dates or date ranges, which must be matched by as many shock values in
the `values` keyword. Note that a range in the `periods` keyword can be
matched by only one value in the `values` keyword. If `values`
represents a scalar, the same value applies to the whole range. If
`values` represents a vector, it must have as many elements as there are
periods in the range.

Note that shock values are not restricted to numerical constants:
arbitrary expressions are also allowed, but you have to enclose them
inside parentheses.

The feasible range of `periods` is from 0 to the number of `periods`
specified in `perfect_foresight_setup`.

::: {.warning}
::: {.admonition-title}
Warning
:::

Note that the first endogenous simulation period is period 1. Thus, a
shock value specified for the initial period 0 may conflict with (i.e.
may overwrite or be overwritten by) values for the initial period
specified with `initval` or `endval` (depending on the exact context).
Users should always verify the correct setting of `oo_.exo_simul` after
`perfect_foresight_setup`.
:::

*Example* (with scalar values)

    shocks;

    var e;
    periods 1;
    values 0.5;
    var u;
    periods 4:5;
    values 0;
    var v;
    periods 4:5 6 7:9;
    values 1 1.1 0.9;
    var w;
    periods 1 2;
    values (1+p) (exp(z));

    end;

*Example* (with vector values)

    xx = [1.2; 1.3; 1];

    shocks;
    var e;
    periods 1:3;
    values (xx);
    end;

*In stochastic context*

For stochastic simulations, the `shocks` block specifies the non zero
elements of the covariance matrix of the shocks of exogenous variables.

You can use the following types of entries in the block:

-   Specification of the standard error of an exogenous variable.

        var VARIABLE_NAME; stderr EXPRESSION;

-   Specification of the variance of an exogenous variable.

        var VARIABLE_NAME = EXPRESSION;

-   Specification the covariance of two exogenous variables.

        var VARIABLE_NAME, VARIABLE_NAME = EXPRESSION;

-   Specification of the correlation of two exogenous variables.

        corr VARIABLE_NAME, VARIABLE_NAME = EXPRESSION;

In an estimation context, it is also possible to specify variances and
covariances on endogenous variables: in that case, these values are
interpreted as the calibration of the measurement errors on these
variables. This requires the `varobs` command to be specified before the
`shocks` block.

*Example*

    shocks;
    var e = 0.000081;
    var u; stderr 0.009;
    corr e, u = 0.8;
    var v, w = 2;
    end;

*In stochastic optimal policy context*

When computing conditional welfare in a `ramsey_model` or
`discretionary_policy` context, welfare is conditional on the state
values inherited by planner when making choices in the first period. The
information set of the first period includes the respective exogenous
shock realizations. Thus, their known value can be specified using the
perfect foresight syntax. Note that i) all other values specified for
periods than period 1 will be ignored and ii) the value of lagged shocks
(e.g. in the case of news shocks) is specified with `histval`.

*Example*

    shocks;
     var u; stderr 0.008;
     var u;
     periods 1;
     values 1;
    end;

*Mixing deterministic and stochastic shocks*

It is possible to mix deterministic and stochastic shocks to build
models where agents know from the start of the simulation about future
exogenous changes. In that case `stoch_simul` will compute the rational
expectation solution adding future information to the state space
(nothing is shown in the output of `stoch_simul`) and `forecast` will
compute a simulation conditional on initial conditions and future
information.

*Example*

    varexo_det tau;
    varexo e;
    ...
    shocks;
    var e; stderr 0.01;
    var tau;
    periods 1:9;
    values -0.15;
    end;

    stoch_simul(irf=0);

    forecast;
:::

::: {.block}
mshocks ; mshocks(overwrite);

The purpose of this block is similar to that of the `shocks` block for
deterministic shocks, except that the numeric values given will be
interpreted in a multiplicative way. For example, if a value of `1.05`
is given as shock value for some exogenous at some date, it means 5%
above its steady state value (as given by the last `initval` or `endval`
block).

The syntax is the same as `shocks` in a deterministic context.

This command is only meaningful in two situations:

-   on exogenous variables with a non-zero steady state, in a
    deterministic setup,
-   on deterministic exogenous variables with a non-zero steady state,
    in a stochastic setup.

See above for the meaning of the `overwrite` option.
:::

::: {.block}
heteroskedastic\_shocks ; heteroskedastic\_shocks(overwrite);

In *estimation context*, it implements heteroskedastic filters, where
the standard error of shocks may unexpectedly change in every period.
The standard deviation of shocks may be either provided directly or
set/modified in each observed period by a scale factor. If `std0` is the
usual standard error for `shock1`, then:

-   using a scale factor in period `t` implies:
    `std(shock1|t)=std0(shock1)*scale(t)`
-   using a provided value in period `t` implies:
    `std(shock1|t)=value(t)`.

The block has a similar syntax as the `shocks` block in a perfect
foresight context. It should contain one or more occurrences of the
following group of three lines (for setting values):

    var VARIABLE_NAME;
    periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
    values DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;

OR (for setting scale factors):

    var VARIABLE_NAME;
    periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
    scales DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;

NOTE: `scales` and `values` cannot be simultaneously set for the same
shock in the same period, but it is possible to set `values` for some
periods and `scales` for other periods for the same shock. There can be
only one `scales` and `values` directive each for a given shock, so all
affected periods must be set in one statement.

*Example*

    heteroskedastic_shocks;

    var e1;
    periods 86:87, 89:97;
    scales 0.5, 0;

    var e1;
    periods 88;
    values 0.1;

    var e2;
    periods 86:87 88:97;
    values 0.04 0.01;

    end;
:::

::: {.specvar}
Sigma\_e

This special variable specifies directly the covariance matrix of the
stochastic shocks, as an upper (or lower) triangular matrix. Dynare
builds the corresponding symmetric matrix. Each row of the triangular
matrix, except the last one, must be terminated by a semi-colon ;. For a
given element, an arbitrary *EXPRESSION* is allowed (instead of a simple
constant), but in that case you need to enclose the expression in
parentheses. The order of the covariances in the matrix is the same as
the one used in the `varexo` declaration.

*Example*

    varexo u, e;

    Sigma_e = [ 0.81 (phi*0.9*0.009);
                0.000081];

This sets the variance of `u` to 0.81, the variance of `e` to 0.000081,
and the correlation between `e` and `u` to `phi`.

::: {.warning}
::: {.admonition-title}
Warning
:::

**The use of this special variable is deprecated and is strongly
discouraged**. You should use a `shocks` block instead.
:::
:::

::: {.matcomm}
get\_shock\_stderr\_by\_name (\'EXOGENOUS\_NAME\');

Given the name of an exogenous variable, returns its standard deviation,
as set by a previous `shocks` block.
:::

::: {.matcomm}
set\_shock\_stderr\_value (\'EXOGENOUS\_NAME\', MATLAB\_EXPRESSION);

Sets the standard deviation of an exgonous variable. This does
essentially the same as setting the standard error via a `shocks` block,
except that it accepts arbitrary MATLAB/Octave expressions, and that it
works from MATLAB/Octave scripts.
:::

Other general declarations
--------------------------

::: {.command}
dsample INTEGER \[INTEGER\];

Reduces the number of periods considered in subsequent output commands.
:::

::: {.command}
periods INTEGER

This command is now deprecated (but will still work for older model
files). It is not necessary when no simulation is performed and is
replaced by an option `periods` in `perfect_foresight_setup`, `simul`
and `stoch_simul`.

This command sets the number of periods in the simulation. The periods
are numbered from 1 to INTEGER. In perfect foresight simulations, it is
assumed that all future events are perfectly known at the beginning of
period 1.

*Example*

    periods 100;
:::

Steady state {#st-st}
------------

There are two ways of computing the steady state (i.e. the static
equilibrium) of a model. The first way is to let Dynare compute the
steady state using a nonlinear Newton-type solver; this should work for
most models, and is relatively simple to use. The second way is to give
more guidance to Dynare, using your knowledge of the model, by providing
it with a method to compute the steady state, either using a
[steady\_state\_model]{.title-ref} block or writing matlab routine.

### Finding the steady state with Dynare nonlinear solver

::: {.command}
steady ; steady (OPTIONS\...);

This command computes the steady state of a model using a nonlinear
Newton-type solver and displays it. When a steady state file is used
`steady` displays the steady state and checks that it is a solution of
the static model.

More precisely, it computes the equilibrium value of the endogenous
variables for the value of the exogenous variables specified in the
previous `initval` or `endval` block.

`steady` uses an iterative procedure and takes as initial guess the
value of the endogenous variables set in the previous `initval` or
`endval` block.

For complicated models, finding good numerical initial values for the
endogenous variables is the trickiest part of finding the equilibrium of
that model. Often, it is better to start with a smaller model and add
new variables one by one.

*Options*

::: {#steady_maxit}
::: {.option}
maxit = INTEGER

Determines the maximum number of iterations used in the non-linear
solver. The default value of `maxit` is 50.
:::
:::

::: {#steady_tolf}
::: {.option}
tolf = DOUBLE

Convergence criterion for termination based on the function value.
Iteration will cease when the residuals are smaller than `tolf`.
Default: `eps^(1/3)`
:::
:::

::: {.option}
tolx = DOUBLE

Convergence criterion for termination based on the step tolerance along.
Iteration will cease when the attempted step size is smaller than
`tolx`. Default: `eps^(2/3)`
:::

::: {#solvalg}
::: {.option}
solve\_algo = INTEGER

Determines the non-linear solver to use. Possible values for the option
are:

> `0`
>
> > Use `fsolve` (under MATLAB, only available if you have the
> > Optimization Toolbox; always available under Octave).
>
> `1`
>
> > Use Dynare's own nonlinear equation solver (a Newton-like algorithm
> > with line-search).
>
> `2`
>
> > Splits the model into recursive blocks and solves each block in turn
> > using the same solver as value 1.
>
> `3`
>
> > Use Chris Sims' solver.
>
> `4`
>
> > Splits the model into recursive blocks and solves each block in turn
> > using a trust-region solver with autoscaling.
>
> `5`
>
> > Newton algorithm with a sparse Gaussian elimination (SPE) (requires
> > `bytecode` option, see `model-decl`{.interpreted-text role="ref"}).
>
> `6`
>
> > Newton algorithm with a sparse LU solver at each iteration (requires
> > `bytecode` and/or `block` option, see `model-decl`{.interpreted-text
> > role="ref"}).
>
> `7`
>
> > Newton algorithm with a Generalized Minimal Residual (GMRES) solver
> > at each iteration (requires `bytecode` and/or `block` option, see
> > `model-decl`{.interpreted-text role="ref"}).
>
> `8`
>
> > Newton algorithm with a Stabilized Bi-Conjugate Gradient (BICGSTAB)
> > solver at each iteration (requires bytecode and/or block option, see
> > `model-decl`{.interpreted-text role="ref"}).
>
> `9`
>
> > Trust-region algorithm on the entire model.
>
> `10`
>
> > Levenberg-Marquardt mixed complementarity problem (LMMCP) solver
> > (*Kanzow and Petra (2004)*).
>
> `11`
>
> > PATH mixed complementarity problem solver of *Ferris and Munson
> > (1999)*. The complementarity conditions are specified with an `mcp`
> > equation tag, see `lmmcp`{.interpreted-text role="opt"}. Dynare only
> > provides the interface for using the solver. Due to licence
> > restrictions, you have to download the solver's most current version
> > yourself from <http://pages.cs.wisc.edu/~ferris/path.html> and place
> > it in MATLAB's search path.
>
> `12`
>
> > Specialized version of `2` for models where all the equations have
> > one endogenous variable on the left hand side and where each
> > equation determines a different endogenous variable. Only
> > expressions allowed on the left hand side are the natural logarithm
> > of an endogenous variable, the first difference of an endogenous
> > variable (with the `diff` operator), or the first difference of the
> > logarithm of an endogenous variable. Univariate blocks are solved by
> > evaluating the expression on the right hand side.
>
> `14`
>
> > Specialized version of `4` for models where all the equations have
> > one endogenous variable on the left hand side and where each
> > equation determines a different endogenous variable. Only
> > expressions allowed on the left hand side are the natural logarithm
> > of an endogenous variable, the first difference of an endogenous
> > variable (with the `diff` operator), or the first difference of the
> > logarithm of an endogenous variable.. Univariate blocks are solved
> > by evaluating the expression on the right hand side.

Default value is `4`.
:::
:::

::: {.option}
homotopy\_mode = INTEGER

Use a homotopy (or divide-and-conquer) technique to solve for the steady
state. If you use this option, you must specify a `homotopy_setup`
block. This option can take three possible values:

> `1`
>
> > In this mode, all the parameters are changed simultaneously, and the
> > distance between the boundaries for each parameter is divided in as
> > many intervals as there are steps (as defined by the
> > `homotopy_steps` option); the problem is solved as many times as
> > there are steps.
>
> `2`
>
> > Same as mode `1`, except that only one parameter is changed at a
> > time; the problem is solved as many times as steps times number of
> > parameters.
>
> `3`
>
> > Dynare tries first the most extreme values. If it fails to compute
> > the steady state, the interval between initial and desired values is
> > divided by two for all parameters. Every time that it is impossible
> > to find a steady state, the previous interval is divided by two.
> > When it succeeds to find a steady state, the previous interval is
> > multiplied by two. In that last case `homotopy_steps` contains the
> > maximum number of computations attempted before giving up.
:::

::: {.option}
homotopy\_steps = INTEGER

Defines the number of steps when performing a homotopy. See
`homotopy_mode` option for more details.
:::

::: {.option}
homotopy\_force\_continue = INTEGER

This option controls what happens when homotopy fails.

> `0`
>
> > `steady` fails with an error message
>
> `1`
>
> > `steady` keeps the values of the last homotopy step that was
> > successful and continues. **BE CAREFUL**: parameters and/or
> > exogenous variables are NOT at the value expected by the user

Default is `0`.
:::

::: {.option}
nocheck

Don't check the steady state values when they are provided explicitly
either by a steady state file or a `steady_state_model` block. This is
useful for models with unit roots as, in this case, the steady state is
not unique or doesn't exist.
:::

::: {#steady_markowitz}
::: {.option}
markowitz = DOUBLE

Value of the Markowitz criterion (:math:(0,infty)[) used to select the
pivot with sparse Gaussian elimination (]{.title-ref}[solve\_algo =
5]{.title-ref}[). This criterion governs the tradeoff between selecting
the pivot resulting in the most accurate solution (low
]{.title-ref}[markowitz]{.title-ref}[ values) and the one that preserves
maximum sparsity (high ]{.title-ref}[markowitz]{.title-ref}\` values).
Default: 0.5.
:::
:::

*Example*

See `init-term-cond`{.interpreted-text role="ref"}.
:::

After computation, the steady state is available in the following
variable:

::: {.matvar}
[oo]().steady\_state

Contains the computed steady state. Endogenous variables are ordered in
the order of declaration used in the `var` command (which is also the
order used in `M_.endo_names`).
:::

::: {.matcomm}
get\_mean (\'ENDOGENOUS\_NAME\' \[, \'ENDOGENOUS\_NAME\'\]\... );

Returns the steady of state of the given endogenous variable(s), as it
is stored in `oo_.steady_state`. Note that, if the steady state has not
yet been computed with `steady`, it will first try to compute it.
:::

::: {.block}
homotopy\_setup ;

This block is used to declare initial and final values when using a
homotopy method. It is used in conjunction with the option
`homotopy_mode` of the steady command.

The idea of homotopy (also called divide-and-conquer by some authors) is
to subdivide the problem of finding the steady state into smaller
problems. It assumes that you know how to compute the steady state for a
given set of parameters, and it helps you finding the steady state for
another set of parameters, by incrementally moving from one to another
set of parameters.

The purpose of the `homotopy_setup` block is to declare the final (and
possibly also the initial) values for the parameters or exogenous that
will be changed during the homotopy. It should contain lines of the
form:

    VARIABLE_NAME, EXPRESSION, EXPRESSION;

This syntax specifies the initial and final values of a given
parameter/exogenous.

There is an alternative syntax:

    VARIABLE_NAME, EXPRESSION;

Here only the final value is specified for a given parameter/exogenous;
the initial value is taken from the preceeding `initval` block.

A necessary condition for a successful homotopy is that Dynare must be
able to solve the steady state for the initial parameters/exogenous
without additional help (using the guess values given in the `initval`
block).

If the homotopy fails, a possible solution is to increase the number of
steps (given in `homotopy_steps` option of `steady`).

*Example*

In the following example, Dynare will first compute the steady state for
the initial values (`gam=0.5` and `x=1`), and then subdivide the problem
into 50 smaller problems to find the steady state for the final values
(`gam=2` and `x=2`):

    var c k;
    varexo x;

    parameters alph gam delt bet aa;
    alph=0.5;
    delt=0.02;
    aa=0.5;
    bet=0.05;

    model;
    c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
    c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
    end;

    initval;
    x = 1;
    k = ((delt+bet)/(aa*x*alph))^(1/(alph-1));
    c = aa*x*k^alph-delt*k;
    end;

    homotopy_setup;
    gam, 0.5, 2;
    x, 2;
    end;

    steady(homotopy_mode = 1, homotopy_steps = 50);
:::

### Providing the steady state to Dynare

If you know how to compute the steady state for your model, you can
provide a MATLAB/Octave function doing the computation instead of using
`steady`. Again, there are two options for doing that:

> -   The easiest way is to write a `steady_state_model` block, which is
>     described below in more details. See also `fs2000.mod` in the
>     `examples` directory for an example. The steady state file
>     generated by Dynare will be called `+FILENAME/steadystate.m.`
> -   You can write the corresponding MATLAB function by hand. If your
>     MOD-file is called `FILENAME.mod`, the steady state file must be
>     called `FILENAME_steadystate.m`. See `NK_baseline_steadystate.m`
>     in the examples directory for an example. This option gives a bit
>     more flexibility (loops and conditional structures can be used),
>     at the expense of a heavier programming burden and a lesser
>     efficiency.

Note that both files allow to update parameters in each call of the
function. This allows for example to calibrate a model to a labor supply
of 0.2 in steady state by setting the labor disutility parameter to a
corresponding value (see `NK_baseline_steadystate.m` in the `examples`
directory). They can also be used in estimation where some parameter may
be a function of an estimated parameter and needs to be updated for
every parameter draw. For example, one might want to set the capital
utilization cost parameter as a function of the discount rate to ensure
that capacity utilization is 1 in steady state. Treating both parameters
as independent or not updating one as a function of the other would lead
to wrong results. But this also means that care is required. Do not
accidentally overwrite your parameters with new values as it will lead
to wrong results.

::: {.block}
steady\_state\_model ;

When the analytical solution of the model is known, this command can be
used to help Dynare find the steady state in a more efficient and
reliable way, especially during estimation where the steady state has to
be recomputed for every point in the parameter space.

Each line of this block consists of a variable (either an endogenous, a
temporary variable or a parameter) which is assigned an expression
(which can contain parameters, exogenous at the steady state, or any
endogenous or temporary variable already declared above). Each line
therefore looks like:

    VARIABLE_NAME = EXPRESSION;

Note that it is also possible to assign several variables at the same
time, if the main function in the right hand side is a MATLAB/Octave
function returning several arguments:

    [ VARIABLE_NAME, VARIABLE_NAME... ] = EXPRESSION;

Dynare will automatically generate a steady state file (of the form
`+FILENAME/steadystate.m`) using the information provided in this block.

*Steady state file for deterministic models*

The `steady_state_model` block also works with deterministic models. An
`initval` block and, when necessary, an `endval` block, is used to set
the value of the exogenous variables. Each `initval` or `endval` block
must be followed by `steady` to execute the function created by
`steady_state_model` and set the initial, respectively terminal, steady
state.

*Example*

>     var m P c e W R k d n l gy_obs gp_obs y dA;
>     varexo e_a e_m;
>
>     parameters alp bet gam mst rho psi del;
>
>     ...
>     // parameter calibration, (dynamic) model declaration, shock calibration...
>     ...
>
>     steady_state_model;
>       dA = exp(gam);
>       gst = 1/dA; // A temporary variable
>       m = mst;
>
>       // Three other temporary variables
>       khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
>       xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
>       nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
>
>       n  = xist/(nust+xist);
>       P  = xist + nust;
>       k  = khst*n;
>
>       l  = psi*mst*n/( (1-psi)*(1-n) );
>       c  = mst/P;
>       d  = l - mst + 1;
>       y  = k^alp*n^(1-alp)*gst^alp;
>       R  = mst/bet;
>
>       // You can use MATLAB functions which return several arguments
>       [W, e] = my_function(l, n);
>
>       gp_obs = m/dA;
>       gy_obs = dA;
>     end;
>
>     steady;
:::

### Replace some equations during steady state computations {#eq-tag-ss}

When there is no steady state file, Dynare computes the steady state by
solving the static model, i.e. the model from the `.mod` file from which
leads and lags have been removed.

In some specific cases, one may want to have more control over the way
this static model is created. Dynare therefore offers the possibility to
explicitly give the form of equations that should be in the static
model.

More precisely, if an equation is prepended by a `[static]` tag, then it
will appear in the static model used for steady state computation, but
that equation will not be used for other computations. For every
equation tagged in this way, you must tag another equation with
`[dynamic]`: that equation will not be used for steady state
computation, but will be used for other computations.

This functionality can be useful on models with a unit root, where there
is an infinity of steady states. An equation (tagged `[dynamic]`) would
give the law of motion of the nonstationary variable (like a random
walk). To pin down one specific steady state, an equation tagged
`[static]` would affect a constant value to the nonstationary variable.
Another situation where the `[static]` tag can be useful is when one has
only a partial closed form solution for the steady state.

*Example*

This is a trivial example with two endogenous variables. The second
equation takes a different form in the static model:

    var c k;
    varexo x;
    ...
    model;
    c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
    [dynamic] c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
    [static] k = ((delt+bet)/(x*aa*alph))^(1/(alph-1));
    end;

Getting information about the model
-----------------------------------

::: {.command}
check ; check (OPTIONS\...);

Computes the eigenvalues of the model linearized around the values
specified by the last `initval`, `endval` or `steady` statement.
Generally, the eigenvalues are only meaningful if the linearization is
done around a steady state of the model. It is a device for local
analysis in the neighborhood of this steady state.

A necessary condition for the uniqueness of a stable equilibrium in the
neighborhood of the steady state is that there are as many eigenvalues
larger than one in modulus as there are forward looking variables in the
system. An additional rank condition requires that the square submatrix
of the right Schur vectors corresponding to the forward looking
variables (jumpers) and to the explosive eigenvalues must have full
rank.

Note that the outcome may be different from what would be suggested by
`sum(abs(oo_.dr.eigval))` when eigenvalues are very close to
`qz_criterium <qz_criterium = DOUBLE>`{.interpreted-text role="opt"}.

*Options*

::: {.option}
solve\_algo = INTEGER

See `solve_algo <solvalg>`{.interpreted-text role="ref"}, for the
possible values and their meaning.
:::

::: {.option}
qz\_zero\_threshold = DOUBLE

Value used to test if a generalized eigenvalue is $0/0$ in the
generalized Schur decomposition (in which case the model does not admit
a unique solution). Default: `1e-6`.
:::

*Output*

`check` returns the eigenvalues in the global variable `oo_.dr.eigval`.
:::

::: {.matvar}
[oo]().dr.eigval

Contains the eigenvalues of the model, as computed by the `check`
command.
:::

::: {.command}
model\_diagnostics ;

This command performs various sanity checks on the model, and prints a
message if a problem is detected (missing variables at current period,
invalid steady state, singular Jacobian of static model).
:::

::: {.command}
model\_info ; model\_info (OPTIONS\...);

This command provides information about the model.

When used outside the context of the `block` option of the `model`
block, it will provide a list of predetermined state variables,
forward-looking variables, and purely static variables.

When used in conjunction with the `block` option of the `model` block,
it displays:

-   The normalization of the model: an endogenous variable is attributed
    to each equation of the model;
-   The block structure of the model: for each block `model_info`
    indicates its type, the equations number and endogenous variables
    belonging to this block.

There are five different types of blocks depending on the simulation
method used:

-   `EVALUATE FORWARD`

    In this case the block contains only equations where the endogenous
    variable attributed to the equation appears at current period on the
    left hand side and where no forward looking endogenous variables
    appear. The block has the form:
    $y_{j,t} = f_j(y_t, y_{t-1}, \ldots, y_{t-k})$.

-   `EVALUATE BACKWARD`

    The block contains only equations where the endogenous variable
    attributed to the equation appears at current period on the left
    hand side and where no backward looking endogenous variables appear.
    The block has the form: $y_{j,t} = f_j(y_t,
    y_{t+1}, \ldots, y_{t+k})$.

-   `SOLVE BACKWARD x`

    The block contains only equations where the endogenous variable
    attributed to the equation does not appear at current period on the
    left hand side and where no forward looking endogenous variables
    appear. The block has the form: $g_j(y_{j,t},
    y_t, y_{t-1}, \ldots, y_{t-k})=0$. `x` is equal to `SIMPLE` if the
    block has only one equation. If several equations appear in the
    block, `x` is equal to `COMPLETE`.

-   `SOLVE FORWARD x`

    The block contains only equations where the endogenous variable
    attributed to the equation does not appear at current period on the
    left hand side and where no backward looking endogenous variables
    appear. The block has the form: $g_j(y_{j,t},
    y_t, y_{t+1}, \ldots, y_{t+k})=0$. `x` is equal to `SIMPLE` if the
    block has only one equation. If several equations appear in the
    block, `x` is equal to `COMPLETE`.

-   `SOLVE TWO BOUNDARIES x`

    The block contains equations depending on both forward and backward
    variables. The block looks like: $g_j(y_{j,t},
    y_t, y_{t-1}, \ldots, y_{t-k} ,y_t, y_{t+1}, \ldots,
    y_{t+k})=0$. `x` is equal to `SIMPLE` if the block has only one
    equation. If several equations appear in the block, `x` is equal to
    `COMPLETE`.

*Options*

::: {.option}
static

Prints out the block decomposition of the static model. Without the
`static` option, `model_info` displays the block decomposition of the
dynamic model.
:::

::: {.option}
incidence

Displays the gross incidence matrix and the reordered incidence matrix
of the block decomposed model.
:::
:::

::: {.command}
print\_bytecode\_dynamic\_model ;

Prints the equations and the Jacobian matrix of the dynamic model stored
in the bytecode binary format file. Can only be used in conjunction with
the `bytecode` option of the `model` block.
:::

::: {.command}
print\_bytecode\_static\_model ;

Prints the equations and the Jacobian matrix of the static model stored
in the bytecode binary format file. Can only be used in conjunction with
the `bytecode` option of the `model` block.
:::

Deterministic simulation {#det-simul}
------------------------

### Perfect foresight

When the framework is deterministic, Dynare can be used for models with
the assumption of perfect foresight. Typically, the system is supposed
to be in a state of equilibrium before a period `1` when the news of a
contemporaneous or of a future shock is learned by the agents in the
model. The purpose of the simulation is to describe the reaction in
anticipation of, then in reaction to the shock, until the system returns
to the old or to a new state of equilibrium. In most models, this return
to equilibrium is only an asymptotic phenomenon, which one must
approximate by an horizon of simulation far enough in the future.
Another exercise for which Dynare is well suited is to study the
transition path to a new equilibrium following a permanent shock. For
deterministic simulations, the numerical problem consists of solving a
nonlinear system of simultaneous equations in `n` endogenous variables
in `T` periods. Dynare offers several algorithms for solving this
problem, which can be chosen via the `stack_solve_algo` option. By
default (`stack_solve_algo=0`), Dynare uses a Newton-type method to
solve the simultaneous equation system. Because the resulting Jacobian
is in the order of `n` by `T` and hence will be very large for long
simulations with many variables, Dynare makes use of the sparse matrix
capacities of MATLAB/Octave. A slower but potentially less memory
consuming alternative (`stack_solve_algo=1`) is based on a Newton-type
algorithm first proposed by *Laffargue (1990)* and *Boucekkine (1995)*,
which avoids ever storing the full Jacobian. The details of the
algorithm can be found in *Juillard (1996)*. The third type of
algorithms makes use of block decomposition techniques
(divide-and-conquer methods) that exploit the structure of the model.
The principle is to identify recursive and simultaneous blocks in the
model structure and use this information to aid the solution process.
These solution algorithms can provide a significant speed-up on large
models.

::: {.warning}
::: {.admonition-title}
Warning
:::

Be careful when employing auxiliary variables in the context of perfect
foresight computations. The same model may work for stochastic
simulations, but fail for perfect foresight simulations. The issue
arises when an equation suddenly only contains variables dated `t+1` (or
`t-1` for that matter). In this case, the derivative in the last (first)
period with respect to all variables will be 0, rendering the stacked
Jacobian singular.

*Example*

> Consider the following specification of an Euler equation with log
> utility:
>
>     Lambda = beta*C(-1)/C;
>     Lambda(+1)*R(+1)= 1;
>
> Clearly, the derivative of the second equation with respect to all
> endogenous variables at time `t` is zero, causing
> `perfect_foresight_solver` to generally fail. This is due to the use
> of the Lagrange multiplier `Lambda` as an auxiliary variable. Instead,
> employing the identical
>
>     beta*C/C(+1)*R(+1)= 1;
>
> will work.
:::

::: {.command}
perfect\_foresight\_setup ; perfect\_foresight\_setup (OPTIONS\...);

Prepares a perfect foresight simulation, by extracting the information
in the `initval`, `endval` and `shocks` blocks and converting them into
simulation paths for exogenous and endogenous variables.

This command must always be called before running the simulation with
`perfect_foresight_solver`.

*Options*

::: {.option}
periods = INTEGER

Number of periods of the simulation.
:::

::: {.option}
datafile = FILENAME

Used to specify path for all endogenous and exogenous variables.
Strictly equivalent to `initval_file`{.interpreted-text role="comm"}.
:::

*Output*

The paths for the exogenous variables are stored into `oo_.exo_simul`.

The initial and terminal conditions for the endogenous variables and the
initial guess for the path of endogenous variables are stored into
`oo_.endo_simul`.
:::

::: {.command}
perfect\_foresight\_solver ; perfect\_foresight\_solver (OPTIONS\...);

Computes the perfect foresight (or deterministic) simulation of the
model.

Note that `perfect_foresight_setup` must be called before this command,
in order to setup the environment for the simulation.

*Options*

::: {.option}
maxit = INTEGER

Determines the maximum number of iterations used in the non-linear
solver. The default value of `maxit` is `50`.
:::

::: {.option}
tolf = DOUBLE

Convergence criterion for termination based on the function value.
Iteration will cease when it proves impossible to improve the function
value by more than `tolf`. Default: `1e-5`
:::

::: {.option}
tolx = DOUBLE

Convergence criterion for termination based on the change in the
function argument. Iteration will cease when the solver attempts to take
a step that is smaller than `tolx`. Default: `1e-5`
:::

::: {.option}
noprint

Don't print anything. Useful for loops.
:::

::: {.option}
print

Print results (opposite of `noprint`).
:::

::: {.option}
stack\_solve\_algo = INTEGER

Algorithm used for computing the solution. Possible values are:

> `0`
>
> > Newton method to solve simultaneously all the equations for every
> > period, using sparse matrices (Default).
>
> `1`
>
> > Use the Laffargue-Boucekkine-Juillard (LBJ) algorithm proposed in
> > *Juillard (1996)*. It is slower than `stack_solve_algo=0`, but may
> > be less memory consuming on big models. Note that if the `block`
> > option is used (see `model-decl`{.interpreted-text role="ref"}), a
> > simple Newton algorithm with sparse matrices is used for blocks
> > which are purely backward or forward (of type `SOLVE BACKWARD` or
> > `SOLVE FORWARD`, see `model_info`{.interpreted-text role="comm"}),
> > since LBJ only makes sense on blocks with both leads and lags (of
> > type `SOLVE TWO BOUNDARIES`).
>
> `2`
>
> > Use a Newton algorithm with a Generalized Minimal Residual (GMRES)
> > solver at each iteration (requires `bytecode` and/or `block` option,
> > see `model-decl`{.interpreted-text role="ref"})
>
> `3`
>
> > Use a Newton algorithm with a Stabilized Bi-Conjugate Gradient
> > (BICGSTAB) solver at each iteration (requires `bytecode` and/or
> > `block` option, see `model-decl`{.interpreted-text role="ref"}).
>
> `4`
>
> > Use a Newton algorithm with a optimal path length at each iteration
> > (requires `bytecode` and/or `block` option, see
> > `model-decl`{.interpreted-text role="ref"}).
>
> `5`
>
> > Use a Newton algorithm with a sparse Gaussian elimination (SPE)
> > solver at each iteration (requires `bytecode` option, see
> > `model-decl`{.interpreted-text role="ref"}).
>
> `6`
>
> > Synonymous for `stack_solve_algo=1`. Kept for historical reasons.
>
> `7`
>
> > Allows the user to solve the perfect foresight model with the
> > solvers available through option `solve_algo` (See
> > `solve_algo <solvalg>`{.interpreted-text role="ref"} for a list of
> > possible values, note that values 5, 6, 7 and 8, which require
> > `bytecode` and/or `block` options, are not allowed). For instance,
> > the following commands:
> >
> >     perfect_foresight_setup(periods=400);
> >     perfect_foresight_solver(stack_solve_algo=7, solve_algo=9)
> >
> > trigger the computation of the solution with a trust region
> > algorithm.
:::

::: {.option}
robust\_lin\_solve

Triggers the use of a robust linear solver for the default
`stack_solve_algo=0`.
:::

::: {.option}
solve\_algo

See `solve_algo <solvalg>`{.interpreted-text role="ref"}. Allows
selecting the solver used with `stack_solve_algo=7`.
:::

::: {.option}
no\_homotopy

By default, the perfect foresight solver uses a homotopy technique if it
cannot solve the problem. Concretely, it divides the problem into
smaller steps by diminishing the size of shocks and increasing them
progressively until the problem converges. This option tells Dynare to
disable that behavior. Note that the homotopy is not implemented for
purely forward or backward models.
:::

::: {.option}
markowitz = DOUBLE

Value of the Markowitz criterion, used to select the pivot. Only used
when `stack_solve_algo = 5`. Default: `0.5`.
:::

::: {.option}
minimal\_solving\_periods = INTEGER

Specify the minimal number of periods where the model has to be solved,
before using a constant set of operations for the remaining periods.
Only used when `stack_solve_algo = 5`. Default: `1`.
:::

::: {.option}
lmmcp

Solves the perfect foresight model with a Levenberg-Marquardt mixed
complementarity problem (LMMCP) solver (*Kanzow and Petra (2004)*),
which allows to consider inequality constraints on the endogenous
variables (such as a ZLB on the nominal interest rate or a model with
irreversible investment). This option is equivalent to
`stack_solve_algo=7` **and** `solve_algo=10`. Using the LMMCP solver
requires a particular model setup as the goal is to get rid of any
min/max operators and complementary slackness conditions that might
introduce a singularity into the Jacobian. This is done by attaching an
equation tag (see `model-decl`{.interpreted-text role="ref"}) with the
`mcp` keyword to affected equations. This tag states that the equation
to which the tag is attached has to hold unless the expression within
the tag is binding. For instance, a ZLB on the nominal interest rate
would be specified as follows in the model block:

    model;
       ...
       [mcp = 'r > -1.94478']
       r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;
       ...
    end;

where `1.94478` is the steady state level of the nominal interest rate
and `r` is the nominal interest rate in deviation from the steady state.
This construct implies that the Taylor rule is operative, unless the
implied interest rate `r<=-1.94478`, in which case the `r` is fixed at
`-1.94478` (thereby being equivalent to a complementary slackness
condition). By restricting the value of `r` coming out of this equation,
the `mcp` tag also avoids using `max(r,-1.94478)` for other occurrences
of `r` in the rest of the model. It is important to keep in mind that,
because the `mcp` tag effectively replaces a complementary slackness
condition, it cannot be simply attached to any equation. Rather, it must
be attached to the correct affected equation as otherwise the solver
will solve a different problem than originally intended. Also, since the
problem to be solved is nonlinear, the sign of the residuals of the
dynamic equation matters. In the previous example, for the nominal
interest rate rule, if the LHS and RHS are reversed the sign of the
residuals (the difference between the LHS and the RHS) will change and
it may happen that solver fails to identify the solution path. More
generally, convergence of the nonlinear solver is not guaranteed when
using mathematically equivalent representations of the same equation.

Note that in the current implementation, the content of the `mcp`
equation tag is not parsed by the preprocessor. The inequalities must
therefore be as simple as possible: an endogenous variable, followed by
a relational operator, followed by a number (not a variable, parameter
or expression).
:::

::: {.option}
endogenous\_terminal\_period

The number of periods is not constant across Newton iterations when
solving the perfect foresight model. The size of the nonlinear system of
equations is reduced by removing the portion of the paths (and
associated equations) for which the solution has already been identified
(up to the tolerance parameter). This strategy can be interpreted as a
mix of the shooting and relaxation approaches. Note that round off
errors are more important with this mixed strategy (user should check
the reported value of the maximum absolute error). Only available with
option `stack_solve_algo==0`.
:::

::: {.option}
linear\_approximation

Solves the linearized version of the perfect foresight model. The model
must be stationary. Only available with option `stack_solve_algo==0` or
`stack_solve_algo==7`.
:::

*Output*

The simulated endogenous variables are available in global matrix
`oo_.endo_simul`.
:::

::: {.command}
simul ; simul (OPTIONS\...);

Short-form command for triggering the computation of a deterministic
simulation of the model. It is strictly equivalent to a call to
`perfect_foresight_setup` followed by a call to
`perfect_foresight_solver`.

*Options*

Accepts all the options of `perfect_foresight_setup` and
`perfect_foresight_solver`.
:::

::: {.matvar}
[oo]().endo\_simul

This variable stores the result of a deterministic simulation (computed
by `perfect_foresight_solver` or `simul`) or of a stochastic simulation
(computed by `stoch_simul` with the periods option or by
`extended_path`). The variables are arranged row by row, in order of
declaration (as in `M_.endo_names`). Note that this variable also
contains initial and terminal conditions, so it has more columns than
the value of `periods` option.
:::

::: {.matvar}
[oo]().exo\_simul

This variable stores the path of exogenous variables during a simulation
(computed by `perfect_foresight_solver`, `simul`, `stoch_simul` or
`extended_path`). The variables are arranged in columns, in order of
declaration (as in `M_.exo_names`). Periods are in rows. Note that this
convention regarding columns and rows is the opposite of the convention
for `oo_.endo_simul`!
:::

### Perfect foresight with expectation errors

The solution under perfect foresight that was presented in the previous
section makes the assumption that agents learn the complete path of
future shocks in period 1, without making any expectation errors.

One may however want to study a scenario where it turns out that agents
make expectation errors, in the sense that the path they had anticipated
in period 1 does not realize exactly. More precisely, in some simulation
periods, they may receive new information that makes them revise their
anticipation for the path of future shocks. Also, under this scenario,
it is assumed that agents behave as under perfect foresight, *i.e.* they
take their decisions as if there was no uncertainty and they knew
exactly the path of future shocks; the new information that they may
receive comes as a total surprise to them.

Such a scenario can be solved by Dynare using the
`perfect_foresight_with_expectation_errors_setup` and
`perfect_foresight_with_expectation_errors_solver` commands, alongside
`shocks` and `endval` blocks which are given a special `learnt_in`
option.

::: {.block}
shocks(learnt\_in=INTEGER) ; shocks(learnt\_in=INTEGER,overwrite) ;

The `shocks(learnt_in=INTEGER)` can be used to specify temporary shocks
that are learnt in a specific period. It should contain one or more
occurences of the following group of three lines, with the same
semantics as a regular `shocks`{.interpreted-text role="bck"} block:

    var VARIABLE_NAME;
    periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
    values DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;

If the period in which information is learnt is greater or equal than 2,
then it is possible to specify the shock values in deviation with
respect to the values that were expected from the perspective of the
previous period. If the new information consists of an addition to the
previously-anticipated value, the `values` keyword can be replaced by
the `add` keyword; similarly, if the new information consists of an
addition to the previously-anticipated value, the `values` keyword can
be replaced by the `multiply` keyword.

The `overwrite` option says that this block cancels and replaces
previous `shocks` blocks that have the same `learnt_in` option.

Note that a `shocks(learnt_in=1)` block is equivalent to a regular
`shocks`{.interpreted-text role="bck"} block.

*Example*

    shocks(learnt_in=1);
      var x;
      periods 1:2 3:4 5;
      values 1 1.2 1.4;
    end;

    shocks(learnt_in=2);
      var x;
      periods 3:4;
      add 0.1;
    end;

    shocks(learnt_in=4);
      var x;
      periods 5;
      multiply 2;
    end;

This syntax means that:

> -   from the perspective of period 1, `x` is expected to be equal to 1
>     in periods 1 and 2, to 1.2 in periods 3 and 4, and to 1.4 in
>     period 5;
> -   from the perspective of periods 2 (and 3), `x` is expected to be
>     equal to 1 in period 2, to 1.3 in periods 3 and 4, and to 1.4 in
>     period 5;
> -   from the perspective of periods 4 (and following), `x` is expected
>     to be equal to 1.3 in period 4, and to 2.8 in period 5.
:::

::: {.block}
endval(learnt\_in=INTEGER) ;

The `endval(learnt_in=INTEGER)` can be used to specify terminal
conditions that are learnt in a specific period.

Note that an `endval(learnt_in=1)` block is equivalent to a regular
`endval`{.interpreted-text role="bck"} block.

It is possible to express the terminal condition by specifying the level
of the exogenous variable (using an equal symbol, as in a regular
`endval`{.interpreted-text role="bck"} blocks without the `learnt_in`
option). But it is also possible to express the terminal condition as an
addition to the value expected from the perspective of the previous
previous period (using the `+=` operator), or as a multiplicative factor
over that previously expected value (using the `*=` operator).

*Example*

    endval(learnt_in = 3);
      x = 1.1;
      y += 0.1;
      z *= 2;
    end;

This syntax means that, in period 3, the agents learn that:

> -   the terminal condition for `x` will be 1.1;
> -   the terminal condition for `y` will be 0.1 above the terminal
>     condition for `y` that was expected from the perspective of period
>     2;
> -   the terminal condition for `z` will be 2 times the terminal
>     condition for `z` that was expected from the perspective of
>     period 2.

Those values will be the realized ones, unless there is another
`endval(learnt_in=p)` block with `p>3`.
:::

::: {.command}
perfect\_foresight\_with\_expectation\_errors\_setup ;
perfect\_foresight\_with\_expectation\_errors\_setup (OPTIONS\...);

Prepares a perfect foresight simulation with expectation errors, by
extracting the contents of the `initval`, `endval` and `shocks` blocks
(the latter two types of blocks typically used with the `learnt_in`
option); alternatively, the information about future shocks can be given
in a CSV file using the `datafile` option.

This command must always be called before running the simulation with
`perfect_foresight_with_expectation_errors_solver`.

Note that this command makes the assumption that the terminal condition
is always a steady state. Hence, it will recompute the terminal steady
state as many times as the anticipation about the terminal condition
changes. In particular, the information about endogenous variables that
may be given in the `endval` block is ignored.

*Options*

::: {.option}
periods = INTEGER

Number of periods of the simulation.
:::

::: {.option}
datafile = FILENAME

Used to specify the information about future shocks and their
anticipation, as an alternative to `shocks` and `endval` blocks.

The file has the following format:

-   the first column is ignored (can be used to add descriptive labels)
-   the first line contains names of exogenous variables
-   the second line contains, in columns, indices of periods *at which*
    expectations are formed; the information set used in a given period
    is described by all the columns for which that line is equal to the
    period index
-   the subsequent lines correspond to the periods *for which*
    expectations are formed, one period per line; each line gives the
    values of present and future exogenous variables, as seen from the
    period given in the second line
-   the last line corresponds to the terminal condition for exogenous
    variables, as anticipated in the various informational periods

If `p` is the value of the `periods` option and `k` is the number of
exogenous variables, then the CSV file has `p+3` lines and `k×p+1`
columns.

Concretely, the value of a given exogenous in period `t`, as anticipated
from period `s`, is given in line `t+2`, and in the column which has the
name of the variable on the first line and `s` on the second line. Of
course, values in cells corresponding to `t<s` are ignored.
:::

::: {.option}
solve\_algo = INTEGER

See `solve_algo <solvalg>`{.interpreted-text role="ref"}. Used when
computing the terminal steady state.
:::

::: {.option}
tolf = DOUBLE

See `tolf <steady_tolf>`{.interpreted-text role="ref"}. Used when
computing the terminal steady state.
:::

::: {.option}
tolx = DOUBLE

See `tolx <steady_tolx>`{.interpreted-text role="ref"}. Used when
computing the terminal steady state.
:::

::: {.option}
maxit = INTEGER

See `maxit <steady_maxit>`{.interpreted-text role="ref"}. Used when
computing the terminal steady state.
:::

::: {.option}
markowitz = DOUBLE

See `markowitz <steady_markowitz>`{.interpreted-text role="ref"}. Used
when computing the terminal steady state.
:::

*Output*

`oo_.exo_simul` and `oo_.endo_simul` are initialized before the
simulation. Temporary shocks are stored in `oo_.pfwee.shocks_info`,
terminal conditions for exogenous variables are stored in
`oo_.pfwee.terminal_info`, and terminal steady states are stored in
`oo_.pfwee.terminal_steady_state`.

*Example*

Here is a CSV file example that could be given to the `datafile` option
(adding some extra padding space for clarity):

>     Exogenous      ,   x,   x,   x,   x,   x,   x,   x
>     Period   (info),   1,   2,   3,   4,   5,   6,   7
>     Period 1 (real), 1.2,    ,    ,    ,    ,    ,
>     Period 2 (real),   1, 1.3,    ,    ,    ,    ,
>     Period 3 (real),   1,   1, 1.4,    ,    ,    ,
>     Period 4 (real),   1,   1,   1,   1,    ,    ,
>     Period 5 (real),   1,   1,   1,   1,   1,    ,
>     Period 6 (real),   1,   1,   1,   1,   1, 1.1,
>     Period 7 (real),   1,   1,   1,   1,   1, 1.1, 1.1
>     Terminal (real),   1, 1.1, 1.2, 1.2, 1.2, 1.1, 1.1

In this example, there is only one exogenous variable (`x`), and 7
simulation periods. In the first period, agents learn a contemporary
shock (1.2), but anticipate no further shock. In period 2, they learn an
unexpected contemporary shock (1.3), and also a change in the terminal
condition (1.1). In period 3 again there is an unexpected contemporary
shock and a change in the terminal condition. No new information comes
in period 4 and 5. In period 6, an unexpected permanent shock is learnt.
No new information comes in period 7.

Alternatively, instead of using a CSV file, the same sequence of
information sets could be described using the following blocks:

>     initval;
>       x = 1;
>     end;
>
>     steady;
>
>     shocks(learnt_in = 1);
>       var x;
>       periods 1;
>       values 1.2;
>     end;
>
>     shocks(learnt_in = 2);
>       var x;
>       periods 2;
>       values 1.3;
>     end;
>
>     endval(learnt_in = 2);
>       x = 1.1;
>     end;
>
>     shocks(learnt_in = 3);
>       var x;
>       periods 3;
>       values 1.4;
>     end;
>
>     endval(learnt_in = 3);
>       x = 1.2;
>     end;
>
>     shocks(learnt_in = 6);
>       var x;
>       periods 6:7;
>       values 1.1;
>     end;
>
>     endval(learnt_in = 6);
>       x = 1.1;
>     end;
:::

::: {.command}
perfect\_foresight\_with\_expectation\_solver ;
perfect\_foresight\_with\_expectation\_solver (OPTIONS\...);

Computes the perfect foresight simulation with expectation errors of the
model.

Note that `perfect_foresight_with_expectation_errors_setup` must be
called before this command, in order to setup the environment for the
simulation.

*Options*

This command accepts all the options of
`perfect_foresight_solver`{.interpreted-text role="comm"}, with the same
semantics, plus the following ones:

::: {.option}
terminal\_steady\_state\_as\_guess\_value

By default, the initial guess for the computation of the path of
endogenous is the initial steady state (when using the information set
from period 1) or the previously simulated path (when using an
information set that is different from that of period 1). When this
option is given, the initial guess is instead the terminal steady state.
:::

::: {.option}
constant\_simulation\_length

By default, every time the information set changes, the simulation with
the new information set is shorter than the previous one (because the
terminal date is getting closer). When this option is set, every new
simulation has the same length (as specified by the
[periods]{.title-ref}[ option of
:comm:\`perfect\_foresight\_with\_expectation\_errors\_setup]{.title-ref};
as a consequence, the simulated paths as stored in `oo_.endo_simul` will
be longer when this option is set (if [s]{.title-ref} is the last period
in which the information set is modified, then they will contain
[s+periods-1]{.title-ref} periods, excluding initial and terminal
conditions).
:::

*Output*

The simulated paths of endogenous variables are available in
`oo_.endo_simul`.
:::

::: {.matvar}
[oo]().pfwee.shocks\_info

This variable stores the temporary shocks used during perfect foresight
simulations with expectation errors, after
`perfect_foresight_with_expectation_errors_setup`{.interpreted-text
role="comm"} has been run. It is a three-dimensional matrix: first
dimension correspond to exogenous variables (in declaration order);
second dimension corresponds to real time; third dimension corresponds
to informational time. In other words, the value of exogenous indexed
`k` in period `t`, as anticipated from period `s`, is stored in
`oo_.pfwee.shocks_info(k,t,s)`.
:::

::: {.matvar}
[oo]().pfwee.terminal\_info

This variable stores the terminal conditions for exogenous variables
used during perfect foresight simulations with expectation errors, after
`perfect_foresight_with_expectation_errors_setup`{.interpreted-text
role="comm"} has been run. It is a matrix, whose lines correspond to
exogenous variables (in declaration order), and whose columns correspond
to informational time. In other words, the terminal condition for
exogenous indexed `k`, as anticipated from period `s`, is stored in
`oo_.pfwee.terminal_info(k,s)`.
:::

::: {.matvar}
[oo]().pfwee.terminal\_steady\_state

This variable stores the terminal steady states for endogenous variables
used during perfect foresight simulations with expectation errors, after
`perfect_foresight_with_expectation_errors_setup`{.interpreted-text
role="comm"} has been run. It is a matrix, whose lines correspond to
endogenous variables (in declaration order), and whose columns
correspond to informational time. In other words, the terminal steady
state for endogenous indexed `k`, as anticipated from period `s`, is
stored in `oo_.pfwee.terminal_steady_state(k,s)`.
:::

Stochastic solution and simulation {#stoch-sol}
----------------------------------

In a stochastic context, Dynare computes one or several simulations
corresponding to a random draw of the shocks.

The main algorithm for solving stochastic models relies on a Taylor
approximation, up to third order, of the expectation functions (see
*Judd (1996)*, *Collard and Juillard (2001a, 2001b)*, and *Schmitt-Grohé
and Uríbe (2004)*). The details of the Dynare implementation of the
first order solution are given in *Villemot (2011)*. Such a solution is
computed using the `stoch_simul` command.

As an alternative, it is possible to compute a simulation to a
stochastic model using the *extended path* method presented by *Fair and
Taylor (1983)*. This method is especially useful when there are strong
nonlinearities or binding constraints. Such a solution is computed using
the `extended_path` command.

### Computing the stochastic solution {#stoch-sol-simul}

::: {.command}
stoch\_simul \[VARIABLE\_NAME\...\]; stoch\_simul (OPTIONS\...)
\[VARIABLE\_NAME\...\];

Solves a stochastic (i.e. rational expectations) model, using
perturbation techniques.

More precisely, `stoch_simul` computes a Taylor approximation of the
model around the deterministic steady state and solves of the the
decision and transition functions for the approximated model. Using
this, it computes impulse response functions and various descriptive
statistics (moments, variance decomposition, correlation and
autocorrelation coefficients). For correlated shocks, the variance
decomposition is computed as in the VAR literature through a Cholesky
decomposition of the covariance matrix of the exogenous variables. When
the shocks are correlated, the variance decomposition depends upon the
order of the variables in the `varexo` command.

The Taylor approximation is computed around the steady state (see
`st-st`{.interpreted-text role="ref"}).

The IRFs are computed as the difference between the trajectory of a
variable following a shock at the beginning of period `1` and its steady
state value. More details on the computation of IRFs can be found at
<https://archives.dynare.org/DynareWiki/IrFs>.

Variance decomposition, correlation, autocorrelation are only displayed
for variables with strictly positive variance. Impulse response
functions are only plotted for variables with response larger than
$10^{-10}$.

Variance decomposition is computed relative to the sum of the
contribution of each shock. Normally, this is of course equal to
aggregate variance, but if a model generates very large variances, it
may happen that, due to numerical error, the two differ by a significant
amount. Dynare issues a warning if the maximum relative difference
between the sum of the contribution of each shock and aggregate variance
is larger than `0.01%`.

The covariance matrix of the shocks is specified with the `shocks`
command (see `shocks-exo`{.interpreted-text role="ref"}).

When a list of `VARIABLE_NAME` is specified, results are displayed only
for these variables.

The `stoch_simul` command with a first order approximation can benefit
from the block decomposition of the model (see `block`{.interpreted-text
role="opt"}).

*Options*

::: {.option}
ar = INTEGER

Order of autocorrelation coefficients to compute and to print. Default:
`5`.
:::

::: {.option}
drop = INTEGER

Number of points (burnin) dropped at the beginning of simulation before
computing the summary statistics. Note that this option does not affect
the simulated series stored in `oo_.endo_simul` and the workspace. Here,
no periods are dropped. Default: `100`.
:::

::: {.option}
hp\_filter = DOUBLE

Uses HP filter with $\lambda =$ `DOUBLE` before computing moments. If
theoretical moments are requested, the spectrum of the model solution is
filtered following the approach outlined in Uhlig (2001). Default: no
filter.
:::

::: {.option}
one\_sided\_hp\_filter = DOUBLE

Uses the one-sided HP filter with $\lambda =$ `DOUBLE` described in
*Stock and Watson (1999)* before computing moments. This option is only
available with simulated moments. Default: no filter.
:::

::: {.option}
bandpass\_filter

Uses a bandpass filter with the default passband before computing
moments. If theoretical moments are requested, the spectrum of the model
solution is filtered using an ideal bandpass filter. If empirical
moments are requested, the *Baxter and King (1999)* filter is used.
Default: no filter.
:::

::: {.option}
bandpass\_filter = \[HIGHEST\_PERIODICITY LOWEST\_PERIODICITY\]

Uses a bandpass filter before computing moments. The passband is set to
a periodicity of to LOWEST\_PERIODICITY, e.g. $6$ to $32$ quarters if
the model frequency is quarterly. Default: `[6,32]`.
:::

::: {.option}
filtered\_theoretical\_moments\_grid = INTEGER

When computing filtered theoretical moments (with either option
`hp_filter` or option `bandpass_filter`), this option governs the number
of points in the grid for the discrete Inverse Fast Fourier Transform.
It may be necessary to increase it for highly autocorrelated processes.
Default: `512`.
:::

::: {.option}
irf = INTEGER

Number of periods on which to compute the IRFs. Setting `irf=0`
suppresses the plotting of IRFs. Default: `40`.
:::

::: {.option}
irf\_shocks = ( VARIABLE\_NAME \[\[,\] VARIABLE\_NAME \...\] )

The exogenous variables for which to compute IRFs. Default: all.
:::

::: {.option}
relative\_irf

Requests the computation of normalized IRFs. At first order, the normal
shock vector of size one standard deviation is divided by the standard
deviation of the current shock and multiplied by 100. The impulse
responses are hence the responses to a unit shock of size 1 (as opposed
to the regular shock size of one standard deviation), multiplied by 100.
Thus, for a loglinearized model where the variables are measured in
percent, the IRFs have the interpretation of the percent responses to a
100 percent shock. For example, a response of 400 of output to a TFP
shock shows that output increases by 400 percent after a 100 percent TFP
shock (you will see that TFP increases by 100 on impact). Given
linearity at `order=1`, it is straightforward to rescale the IRFs stored
in `oo_.irfs` to any desired size. At higher order, the interpretation
is different. The `relative_irf` option then triggers the generation of
IRFs as the response to a 0.01 unit shock (corresponding to 1 percent
for shocks measured in percent) and no multiplication with 100 is
performed. That is, the normal shock vector of size one standard
deviation is divided by the standard deviation of the current shock and
divided by 100. For example, a response of 0.04 of log output (thus
measured in percent of the steady state output level) to a TFP shock
also measured in percent then shows that output increases by 4 percent
after a 1 percent TFP shock (you will see that TFP increases by 0.01 on
impact).
:::

::: {.option}
irf\_plot\_threshold = DOUBLE

Threshold size for plotting IRFs. All IRFs for a particular variable
with a maximum absolute deviation from the steady state smaller than
this value are not displayed. Default: `1e-10`.
:::

::: {.option}
nocorr

Don't print the correlation matrix (printing them is the default).
:::

::: {.option}
nodecomposition

Don't compute (and don't print) unconditional variance decomposition.
:::

::: {.option}
nofunctions

Don't print the coefficients of the approximated solution (printing them
is the default).
:::

::: {.option}
nomoments

Don't print moments of the endogenous variables (printing them is the
default).
:::

::: {.option}
nograph

Do not create graphs (which implies that they are not saved to the disk
nor displayed). If this option is not used, graphs will be saved to disk
(to the format specified by `graph_format` option, except if
`graph_format=none`) and displayed to screen (unless `nodisplay` option
is used).
:::

::: {.option}
graph

Re-enables the generation of graphs previously shut off with `nograph`.
:::

::: {.option}
nodisplay

Do not display the graphs, but still save them to disk (unless `nograph`
is used).
:::

::: {.option}
graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )

Specify the file format(s) for graphs saved to disk. Possible values are
`eps` (the default), `pdf`, `fig` and `none` (under Octave, `fig` is
unavailable). If the file format is set equal to `none`, the graphs are
displayed but not saved to the disk.
:::

::: {.option}
noprint

See `noprint`{.interpreted-text role="opt"}.
:::

::: {.option}
print

See `print`{.interpreted-text role="opt"}.
:::

::: {.option}
order = INTEGER

Order of Taylor approximation. Note that for third order and above, the
`k_order_solver` option is implied and only empirical moments are
available (you must provide a value for `periods` option). Default: `2`
(except after an `estimation` command, in which case the default is the
value used for the estimation).
:::

::: {.option}
k\_order\_solver

Use a k-order solver (implemented in C++) instead of the default Dynare
solver. This option is not yet compatible with the `bytecode` option
(see `model-decl`{.interpreted-text role="ref"}). Default: disabled for
order 1 and 2, enabled for order 3 and above.
:::

::: {.option}
periods = INTEGER

If different from zero, empirical moments will be computed instead of
theoretical moments. The value of the option specifies the number of
periods to use in the simulations. Values of the initval block, possibly
recomputed by `steady`, will be used as starting point for the
simulation. The simulated endogenous variables are made available to the
user in a vector for each variable and in the global matrix
`oo_.endo_simul` (see `oo_.endo_simul`{.interpreted-text role="mvar"}).
The simulated exogenous variables are made available in `oo_.exo_simul`
(see `oo_.exo_simul`{.interpreted-text role="mvar"}). Default: `0`.
:::

::: {.option}
qz\_criterium = DOUBLE

Value used to split stable from unstable eigenvalues in reordering the
Generalized Schur decomposition used for solving first order problems.
Default: `1.000001` (except when estimating with `lik_init` option equal
to `1`: the default is `0.999999` in that case; see
`estim`{.interpreted-text role="ref"}).
:::

::: {.option}
qz\_zero\_threshold = DOUBLE

See `qz_zero_threshold <qz_zero_threshold = DOUBLE>`{.interpreted-text
role="opt"}.
:::

::: {.option}
replic = INTEGER

Number of simulated series used to compute the IRFs. Default: `1` if
`order=1`, and `50` otherwise.
:::

::: {.option}
simul\_replic = INTEGER

Number of series to simulate when empirical moments are requested (i.e.
`periods` $>$ 0). Note that if this option is greater than 1, the
additional series will not be used for computing the empirical moments
but will simply be saved in binary form to the file `FILENAME_simul` in
the `FILENAME/Output`-folder. Default: `1`.
:::

::: {.option}
solve\_algo = INTEGER

See `solve_algo <solvalg>`{.interpreted-text role="ref"}, for the
possible values and their meaning.
:::

::: {.option}
aim\_solver

Use the Anderson-Moore Algorithm (AIM) to compute the decision rules,
instead of using Dynare's default method based on a generalized Schur
decomposition. This option is only valid for first order approximation.
See [AIM website](https://www.federalreserve.gov/econres/ama-index.htm)
for more details on the algorithm.
:::

::: {.option}
conditional\_variance\_decomposition = INTEGER
conditional\_variance\_decomposition = \[INTEGER1:INTEGER2\]
conditional\_variance\_decomposition = \[INTEGER1 INTEGER2 \...\]

Computes a conditional variance decomposition for the specified
period(s). The periods must be strictly positive. Conditional variances
are given by $var(y_{t+k}\vert t)$. For period 1, the conditional
variance decomposition provides the decomposition of the effects of
shocks upon impact.

The results are stored in `oo_.conditional_variance_decomposition` (see
`oo_.conditional_variance_decomposition`{.interpreted-text
role="mvar"}). In the presence of measurement error, the
`oo_.conditional_variance_decomposition` field will contain the variance
contribution after measurement error has been taken out, i.e. the
decomposition will be conducted of the actual as opposed to the measured
variables. The variance decomposition of the measured variables will be
stored in `oo_.conditional_variance_decomposition_ME` (see
`oo_.conditional_variance_decomposition_ME`{.interpreted-text
role="mvar"}). The variance decomposition is only conducted, if
theoretical moments are requested, *i.e.* using the `periods=0`-option.
Only available at `order<3` and without
`pruning''. In case of`order=2`, Dynare provides a second-order accurate approximation to the true second moments based on the linear terms of the second-order solution (see *Kim, Kim, Schaumburg and Sims (2008)*). Note that the unconditional variance decomposition *i.e.* at horizon infinity) is automatically conducted if theoretical moments are requested and if`nodecomposition\`[
is not set (see :mvar:\`oo\_.variance\_decomposition]{.title-ref}).
:::

::: {.option}
pruning

Discard higher order terms when iteratively computing simulations of the
solution. At second order, Dynare uses the algorithm of *Kim, Kim,
Schaumburg and Sims (2008)*, while at third order its generalization by
*Andreasen, Fernández-Villaverde and Rubio-Ramírez (2018)* is used. Not
available above third order. When specified, theoretical moments are
based on the pruned state space, i.e. the computation of second moments
uses all terms as in *Andreasen, Fernández-Villaverde and Rubio-Ramírez
(2018), page 10* as opposed to simply providing a second-order accurate
result based on the linear solution as in *Kim, Kim, Schaumburg and Sims
(2008)*.
:::

::: {.option}
partial\_information

Computes the solution of the model under partial information, along the
lines of *Pearlman, Currie and Levine (1986)*. Agents are supposed to
observe only some variables of the economy. The set of observed
variables is declared using the `varobs` command. Note that if `varobs`
is not present or contains all endogenous variables, then this is the
full information case and this option has no effect. More references can
be found
[here](https://archives.dynare.org/DynareWiki/PartialInformation) .
:::

::: {.option}
sylvester = OPTION

Determines the algorithm used to solve the Sylvester equation for block
decomposed model. Possible values for OPTION are:

> `default`
>
> > Uses the default solver for Sylvester equations (`gensylv`) based on
> > Ondra Kamenik's algorithm (see
> > [here](https://www.dynare.org/assets/dynare++/sylvester.pdf) for
> > more information).
>
> `fixed_point`
>
> > Uses a fixed point algorithm to solve the Sylvester equation
> > (`gensylv_fp`). This method is faster than the default one for large
> > scale models.

Default value is `default`.
:::

::: {.option}
sylvester\_fixed\_point\_tol = DOUBLE

The convergence criterion used in the fixed point Sylvester solver. Its
default value is `1e-12`.
:::

::: {.option}
dr = OPTION

Determines the method used to compute the decision rule. Possible values
for OPTION are:

> `default`
>
> > Uses the default method to compute the decision rule based on the
> > generalized Schur decomposition (see *Villemot (2011)* for more
> > information).
>
> `cycle_reduction`
>
> > Uses the cycle reduction algorithm to solve the polynomial equation
> > for retrieving the coefficients associated to the endogenous
> > variables in the decision rule. This method is faster than the
> > default one for large scale models.
>
> `logarithmic_reduction`
>
> > Uses the logarithmic reduction algorithm to solve the polynomial
> > equation for retrieving the coefficients associated to the
> > endogenous variables in the decision rule. This method is in general
> > slower than the `cycle_reduction`.

Default value is `default`.
:::

::: {.option}
dr\_cycle\_reduction\_tol = DOUBLE

The convergence criterion used in the cycle reduction algorithm. Its
default value is `1e-7`.
:::

::: {.option}
dr\_logarithmic\_reduction\_tol = DOUBLE

The convergence criterion used in the logarithmic reduction algorithm.
Its default value is `1e-12`.
:::

::: {.option}
dr\_logarithmic\_reduction\_maxiter = INTEGER

The maximum number of iterations used in the logarithmic reduction
algorithm. Its default value is `100`.
:::

::: {.option}
loglinear

See `loglinear <logl>`{.interpreted-text role="ref"}. Note that ALL
variables are log-transformed by using the Jacobian transformation, not
only selected ones. Thus, you have to make sure that your variables have
strictly positive steady states. `stoch_simul` will display the moments,
decision rules, and impulse responses for the log-linearized variables.
The decision rules saved in `oo_.dr` and the simulated variables will
also be the ones for the log-linear variables.
:::

::: {.option}
tex

Requests the printing of results and graphs in TeX tables and graphics
that can be later directly included in LaTeX files.
:::

::: {.option}
dr\_display\_tol = DOUBLE

Tolerance for the suppression of small terms in the display of decision
rules. Rows where all terms are smaller than `dr_display_tol` are not
displayed. Default value: `1e-6`.
:::

::: {.option}
contemporaneous\_correlation

Saves the contemporaneous correlation between the endogenous variables
in `oo_.contemporaneous_correlation`. Requires the `nocorr` option not
to be set.
:::

::: {.option}
spectral\_density

Triggers the computation and display of the theoretical spectral density
of the (filtered) model variables. Results are stored in
`oo_.SpectralDensity`, defined below. Default: do not request spectral
density estimates.
:::

::: {.option}
hp\_ngrid = INTEGER

Deprecated option. It has the same effect as
`filtered_theoretical_moments_grid <filtered_theoretical_moments_grid = INTEGER>`{.interpreted-text
role="opt"}.
:::

*Output*

This command sets `oo_.dr`, `oo_.mean`, `oo_.var`, `oo_.var_list`, and
`oo_.autocorr`, which are described below.

If the `periods` option is present, sets `oo_.skewness`, `oo_.kurtosis`,
and `oo_.endo_simul` (see `oo_.endo_simul`{.interpreted-text
role="mvar"}), and also saves the simulated variables in MATLAB/Octave
vectors of the global workspace with the same name as the endogenous
variables.

If option `irf` is different from zero, sets `oo_.irfs` (see below) and
also saves the IRFs in MATLAB/Octave vectors of the global workspace
(this latter way of accessing the IRFs is deprecated and will disappear
in a future version).

If the option `contemporaneous_correlation` is different from `0`, sets
`oo_.contemporaneous_correlation`, which is described below.

*Example*

>     shocks;
>     var e;
>     stderr 0.0348;
>     end;
>
>     stoch_simul;
>
> Performs the simulation of the 2nd-order approximation of a model with
> a single stochastic shock `e`, with a standard error of `0.0348`.

*Example*

>     stoch_simul(irf=60) y k;
>
> Performs the simulation of a model and displays impulse response
> functions on 60 periods for variables `y` and `k`.
:::

::: {.matvar}
[oo]().mean

After a run of `stoch_simul`, contains the mean of the endogenous
variables. Contains theoretical mean if the `periods` option is not
present, and simulated mean otherwise. The variables are arranged in
declaration order.
:::

::: {.matvar}
[oo]().var

After a run of `stoch_simul`, contains the variance-covariance of the
endogenous variables. Contains theoretical variance if the `periods`
option is not present and simulated variance otherwise. Only available
for `order<4`. At `order=2` it will be be a second-order accurate
approximation (i.e. ignoring terms of order 3 and 4 that would arise
when using the full second-order policy function). At `order=3`,
theoretical moments are only available with `pruning`. The variables are
arranged in declaration order.
:::

::: {.matvar}
[oo]().var\_list

The list of variables for which results are displayed.
:::

::: {.matvar}
[oo]().skewness

After a run of `stoch_simul` contains the skewness (standardized third
moment) of the simulated variables if the `periods` option is present.
The variables are arranged in declaration order.
:::

::: {.matvar}
[oo]().kurtosis

After a run of `stoch_simul` contains the excess kurtosis (standardized
fourth moment) of the simulated variables if the `periods` option is
present. The variables are arranged in declaration order.
:::

::: {.matvar}
[oo]().autocorr

After a run of `stoch_simul`, contains a cell array of the
autocorrelation matrices of the endogenous variables. The element number
of the matrix in the cell array corresponds to the order of
autocorrelation. The option ar specifies the number of autocorrelation
matrices available. Contains theoretical autocorrelations if the
`periods` option is not present and simulated autocorrelations
otherwise. Only available for `order<4`. At `order=2` it will be be a
second-order accurate approximation. At `order=3`, theoretical moments
are only available with `pruning`. The field is only created if
stationary variables are present.

The element `oo_.autocorr{i}(k,l)` is equal to the correlation between
$y^k_t$ and $y^l_{t-i}$, where $y^k$ (resp. $y^l$) is the $k$-th (resp.
$l$-th) endogenous variable in the declaration order.

Note that if theoretical moments have been requested, `oo_.autocorr{i}`
is the same than `oo_.gamma_y{i+1}`.
:::

::: {.matvar}
[oo]().gamma\_y

After a run of `stoch_simul`, if theoretical moments have been requested
(i.e. if the `periods` option is not present), this variable contains a
cell array with the following values (where `ar` is the value of the
option of the same name):

> `oo_.gamma{1}`
>
> > Variance/covariance matrix.
>
> `oo_.gamma{i+1}` (for i=1:ar)
>
> > Autocorrelation function. See `oo_.autocorr`{.interpreted-text
> > role="mvar"} for more details. **Beware**, this is the
> > autocorrelation function, not the autocovariance function.
>
> `oo_.gamma{ar+2}`
>
> > Unconditional variance decomposition, see
> > `oo_.variance_decomposition`{.interpreted-text role="mvar"}.
>
> `oo_.gamma{ar+3}`
>
> > If a second order approximation has been requested, contains the
> > vector of the mean correction terms.
> >
> > Only available at `order<4`. In case `order=2`, the theoretical
> > second moments are a second order accurate approximation of the true
> > second moments. See conditional\_variance\_decomposition. At
> > `order=3`, theoretical moments are only available with `pruning`.
:::

::: {.matvar}
[oo]().variance\_decomposition

After a run of `stoch_simul` when requesting theoretical moments
(`periods=0`), contains a matrix with the result of the unconditional
variance decomposition (i.e. at horizon infinity). The first dimension
corresponds to the endogenous variables (in the order of declaration
after the command or in `M_.endo_names`) and the second dimension
corresponds to exogenous variables (in the order of declaration).
Numbers are in percent and sum up to 100 across columns. In the presence
of measurement error, the field will contain the variance contribution
after measurement error has been taken out, *i.e.* the decomposition
will be conducted of the actual as opposed to the measured variables.
:::

::: {.matvar}
[oo]().variance\_decomposition\_ME

Field set after a run of `stoch_simul` when requesting theoretical
moments (`periods=0`) if measurement error is present. It is similar to
`oo_.variance_decomposition`{.interpreted-text role="mvar"}, but the
decomposition will be conducted of the measured variables. The field
contains a matrix with the result of the unconditional variance
decomposition (*i.e.* at horizon infinity). The first dimension
corresponds to the observed endoogenous variables (in the order of
declaration after the command) and the second dimension corresponds to
exogenous variables (in the order of declaration), with the last column
corresponding to the contribution of measurement error. Numbers are in
percent and sum up to 100 across columns.
:::

::: {.matvar}
[oo]().conditional\_variance\_decomposition

After a run of `stoch_simul` with the
`conditional_variance_decomposition` option, contains a
three-dimensional array with the result of the decomposition. The first
dimension corresponds to the endogenous variables (in the order of
declaration after the command or in `M_.endo_names` if not specified),
the second dimension corresponds to the forecast horizons (as declared
with the option), and the third dimension corresponds to the exogenous
variables (in the order of declaration). In the presence of measurement
error, the field will contain the variance contribution after
measurement error has been taken out, *i.e.* the decomposition will be
conductedof the actual as opposed to the measured variables.
:::

::: {.matvar}
[oo]().conditional\_variance\_decomposition\_ME

Field set after a run of `stoch_simul` with the
`conditional_variance_decomposition` option if measurement error is
present. It is similar to
`oo_.conditional_variance_decomposition`{.interpreted-text role="mvar"},
but the decomposition will be conducted of the measured variables. It
contains a three-dimensional array with the result of the decomposition.
The first dimension corresponds to the endogenous variables (in the
order of declaration after the command or in `M_.endo_names` if not
specified), the second dimension corresponds to the forecast horizons
(as declared with the option), and the third dimension corresponds to
the exogenous variables (in the order of declaration), with the last
column corresponding to the contribution of the measurement error.
:::

::: {.matvar}
[oo]().contemporaneous\_correlation

After a run of `stoch_simul` with the
`contemporaneous_correlation option`, contains theoretical
contemporaneous correlations if the `periods` option is not present, and
simulated contemporaneous correlations otherwise. Only available for
`order<4`. At `order=2` it will be be a second-order accurate
approximation. At `order=3`, theoretical moments are only available with
`pruning`. The variables are arranged in declaration order.
:::

::: {.matvar}
[oo]().SpectralDensity

After a run of `stoch_simul` with option `spectral_density`, contains
the spectral density of the model variables. There will be a `nvars` by
`nfrequencies` subfield `freqs` storing the respective frequency grid
points ranging from $0$ to $2\pi$ and a same sized subfield `density`
storing the corresponding density.
:::

::: {.matvar}
[oo]().irfs

After a run of `stoch_simul` with option `irf` different from zero,
contains the impulse responses, with the following naming convention:
[VARIABLE\_NAME\_SHOCK\_NAME]{.title-ref}.

> For example, `oo_.irfs.gnp_ea` contains the effect on `gnp` of a
> one-standard deviation shock on `ea`.
:::

::: {.matcomm}
get\_irf (\'EXOGENOUS\_NAME\' \[, \'ENDOGENOUS\_NAME\'\]\... );

Given the name of an exogenous variables, returns the IRFs for the
requested endogenous variable(s), as they are stored in `oo_.irfs`.
:::

The approximated solution of a model takes the form of a set of decision
rules or transition equations expressing the current value of the
endogenous variables of the model as function of the previous state of
the model and shocks observed at the beginning of the period. The
decision rules are stored in the structure `oo_.dr` which is described
below.

::: {.matvar}
[oo]().dr

Structure storing the decision rules. The subfields for different orders
of approximation are explained below.
:::

::: {.command}
extended\_path ; extended\_path (OPTIONS\...);

Simulates a stochastic (i.e. rational expectations) model, using the
extended path method presented by *Fair and Taylor (1983)*. Time series
for the endogenous variables are generated by assuming that the agents
believe that there will no more shocks in the following periods.

This function first computes a random path for the exogenous variables
(stored in `oo_.exo_simul`, see `oo_.exo_simul`{.interpreted-text
role="mvar"}) and then computes the corresponding path for endogenous
variables, taking the steady state as starting point. The result of the
simulation is stored in `oo_.endo_simul` (see
`oo_.endo_simul`{.interpreted-text role="mvar"}). Note that this
simulation approach does not solve for the policy and transition
equations but for paths for the endogenous variables.

*Options*

::: {.option}
periods = INTEGER

The number of periods for which the simulation is to be computed. No
default value, mandatory option.
:::

::: {.option}
solver\_periods = INTEGER

The number of periods used to compute the solution of the perfect
foresight at every iteration of the algorithm. Default: `200`.
:::

::: {.option}
order = INTEGER

If order is greater than `0` Dynare uses a gaussian quadrature to take
into account the effects of future uncertainty. If `order` $=S$ then the
time series for the endogenous variables are generated by assuming that
the agents believe that there will no more shocks after period $t+S$.
This is an experimental feature and can be quite slow. A non-zero value
is not compatible with either the `bytecode` or the `block` option of
the `model` block. Default: `0`.
:::

::: {.option}
hybrid

Use the constant of the second order perturbation reduced form to
correct the paths generated by the (stochastic) extended path algorithm.
:::

::: {.option}
lmmcp

Solves the perfect foresight model with a Levenberg-Marquardt mixed
complementarity problem (LMMCP) solver (*Kanzow and Petra (2004)*),
which allows to consider inequality constraints on the endogenous
variables (such as a ZLB on the nominal interest rate or a model with
irreversible investment). For specifying the necessary `mcp`-tag, see
`lmmcp`{.interpreted-text role="opt"}.
:::
:::

### Typology and ordering of variables

Dynare distinguishes four types of endogenous variables:

*Purely backward (or purely predetermined) variables*

> Those that appear only at current and past period in the model, but
> not at future period (i.e. at $t$ and $t-1$ but not $t+1$). The number
> of such variables is equal to `M_.npred`.

*Purely forward variables*

> Those that appear only at current and future period in the model, but
> not at past period (i.e. at $t$ and $t+1$ but not $t-1$). The number
> of such variables is stored in `M_.nfwrd`.

*Mixed variables*

> Those that appear at current, past and future period in the model
> (i.e. at $t$, $t+1$ and $t-1$). The number of such variables is stored
> in `M_.nboth`.

*Static variables*

> Those that appear only at current, not past and future period in the
> model (i.e. only at $t$, not at $t+1$ or $t-1$). The number of such
> variables is stored in `M_.nstatic`.

Note that all endogenous variables fall into one of these four
categories, since after the creation of auxiliary variables (see
`aux-variables`{.interpreted-text role="ref"}), all endogenous have at
most one lead and one lag. We therefore have the following identity:

> ``` {.sourceCode .matlab}
> M_.npred + M_.both + M_.nfwrd + M_.nstatic = M_.endo_nbr
> ```

::: {.matvar}
[M]().state\_var

Vector of numerical indices identifying the state variables in the
vector of declared variables. `M_.endo_names(M_.state_var)` therefore
yields the name of all variables that are states in the model
declaration, i.e. that show up with a lag.
:::

Internally, Dynare uses two orderings of the endogenous variables: the
order of declaration (which is reflected in `M_.endo_names`), and an
order based on the four types described above, which we will call the
DR-order ("DR" stands for decision rules). Most of the time, the
declaration order is used, but for elements of the decision rules, the
DR-order is used.

The DR-order is the following: static variables appear first, then
purely backward variables, then mixed variables, and finally purely
forward variables. Inside each category, variables are arranged
according to the declaration order.

::: {.matvar}
[oo]().dr.order\_var

This variables maps DR-order to declaration order.
:::

::: {.matvar}
[oo]().dr.inv\_order\_var

This variable contains the inverse map.
:::

In other words, the k-th variable in the DR-order corresponds to the
endogenous variable numbered `oo_.dr.order_var(k)` in declaration order.
Conversely, k-th declared variable is numbered `oo_.dr.inv_order_var(k)`
in DR-order.

Finally, the state variables of the model are the purely backward
variables and the mixed variables. They are ordered in DR-order when
they appear in decision rules elements. There are
`M_.nspred = M_.npred + M_.nboth` such variables. Similarly, one has
`M_.nsfwrd = M_.nfwrd + M_.nboth`, and
`M_.ndynamic = M_.nfwrd + M_.nboth + M_.npred`.

### First-order approximation

The approximation has the stylized form:

$$y_t = y^s + A y^h_{t-1} + B u_t$$

where $y^s$ is the steady state value of $y$ and $y^h_t=y_t-y^s$.

::: {.matvar}
oo.dr.state\_var

Vector of numerical indices identifying the state variables in the
vector of declared variables, *given the current parameter values* for
which the decision rules have been computed. It may differ from
`M_.state_var` in case a state variable drops from the model given the
current parameterization, because it only gets 0 coefficients in the
decision rules. See `M_.state_var`{.interpreted-text role="mvar"}.
:::

The coefficients of the decision rules are stored as follows:

-   $y^s$ is stored in `oo_.dr.ys`. The vector rows correspond to all
    endogenous in the declaration order.
-   $A$ is stored in `oo_.dr.ghx`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to state
    variables in DR-order, as given by `oo_.dr.state_var`. (N.B.: if the
    `block` option to the `model` block has been specified, then rows
    are in declaration order, and columns are ordered according to
    `oo_.dr.state_var` which may differ from DR-order.)
-   $B$ is stored `oo_.dr.ghu`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to exogenous
    variables in declaration order. (N.B.: if the `block` option to the
    `model` block has been specified, then rows are in declaration
    order.)

Of course, the shown form of the approximation is only stylized, because
it neglects the required different ordering in $y^s$ and $y^h_t$. The
precise form of the approximation that shows the way Dynare deals with
differences between declaration and DR-order, is

$$y_t(\mathrm{oo\_.dr.order\_var}) =
y^s(\mathrm{oo\_.dr.order\_var}) + A \cdot
y_{t-1}(\mathrm{oo\_.dr.order\_var(k2)}) -
y^s(\mathrm{oo\_.dr.order\_var(k2)}) + B\cdot u_t$$

where $\mathrm{k2}$ selects the state variables, $y_t$ and $y^s$ are in
declaration order and the coefficient matrices are in DR-order.
Effectively, all variables on the right hand side are brought into DR
order for computations and then assigned to $y_t$ in declaration order.

### Second-order approximation

The approximation has the form:

$$y_t = y^s + 0.5 \Delta^2 + A y^h_{t-1} + B u_t + 0.5 C (y^h_{t-1}\otimes y^h_{t-1}) + 0.5 D (u_t \otimes u_t) + E (y^h_{t-1} \otimes u_t)$$

where $y^s$ is the steady state value of $y$, $y^h_t=y_t-y^s$, and
$\Delta^2$ is the shift effect of the variance of future shocks. For the
reordering required due to differences in declaration and DR order, see
the first order approximation.

The coefficients of the decision rules are stored in the variables
described for first order approximation, plus the following variables:

-   $\Delta^2$ is stored in `oo_.dr.ghs2`. The vector rows correspond to
    all endogenous in DR-order.
-   $C$ is stored in `oo_.dr.ghxx`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of the vector of state variables in DR-order.
-   $D$ is stored in `oo_.dr.ghuu`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of exogenous variables in declaration order.
-   $E$ is stored in `oo_.dr.ghxu`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of the vector of state variables (in DR-order) by
    the vector of exogenous variables (in declaration order).

### Third-order approximation

The approximation has the form:

$$y_t = y^s + G_0 + G_1 z_t + G_2 (z_t \otimes z_t) + G_3 (z_t \otimes z_t \otimes z_t)$$

where $y^s$ is the steady state value of $y$, and $z_t$ is a vector
consisting of the deviation from the steady state of the state variables
(in DR-order) at date $t-1$ followed by the exogenous variables at date
$t$ (in declaration order). The vector $z_t$ is therefore of size $n_z$
= `M_.nspred` + `M_.exo_nbr`.

The coefficients of the decision rules are stored as follows:

-   $y^s$ is stored in `oo_.dr.ys`. The vector rows correspond to all
    endogenous in the declaration order.
-   $G_0$ is stored in `oo_.dr.g_0`. The vector rows correspond to all
    endogenous in DR-order.
-   $G_1$ is stored in `oo_.dr.g_1`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to state
    variables in DR-order, followed by exogenous in declaration order.
-   $G_2$ is stored in `oo_.dr.g_2`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of state variables (in DR-order), followed by
    exogenous (in declaration order). Note that the Kronecker product is
    stored in a folded way, i.e. symmetric elements are stored only
    once, which implies that the matrix has $n_z(n_z+1)/2$ columns. More
    precisely, each column of this matrix corresponds to a pair
    $(i_1, i_2)$ where each index represents an element of $z_t$ and is
    therefore between $1$ and $n_z$. Only non-decreasing pairs are
    stored, i.e. those for which $i_1
    \leq i_2$. The columns are arranged in the lexicographical order of
    non-decreasing pairs. Also note that for those pairs where
    $i_1 \neq i_2$, since the element is stored only once but appears
    two times in the unfolded $G_2$ matrix, it must be multiplied by 2
    when computing the decision rules.
-   $G_3$ is stored in `oo_.dr.g_3`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the third
    Kronecker power of state variables (in DR-order), followed by
    exogenous (in declaration order). Note that the third Kronecker
    power is stored in a folded way, i.e. symmetric elements are stored
    only once, which implies that the matrix has $n_z(n_z+1)(n_z+2)/6$
    columns. More precisely, each column of this matrix corresponds to a
    tuple $(i_1, i_2, i_3)$ where each index represents an element of
    $z_t$ and is therefore between $1$ and $n_z$. Only non-decreasing
    tuples are stored, i.e. those for which $i_1 \leq i_2 \leq i_3$. The
    columns are arranged in the lexicographical order of non-decreasing
    tuples. Also note that for tuples that have three distinct indices
    (i.e. $i_1 \neq i_2$ and $i_1 \neq i_3$ and $i_2
    \neq i_3$), since these elements are stored only once but appears
    six times in the unfolded $G_3$ matrix, they must be multiplied by 6
    when computing the decision rules. Similarly, for those tuples that
    have two equal indices (i.e. of the form $(a,a,b)$ or $(a,b,a)$ or
    $(b,a,a)$), since these elements are stored only once but appears
    three times in the unfolded $G_3$ matrix, they must be multiplied by
    3 when computing the decision rules.

### Higher-order approximation

Higher-order approximations are simply a generalization of what is done
at order 3.

The steady state is stored in `oo_.dr.ys` and the constant correction is
stored in `oo_.dr.g_0`. The coefficient for orders 1, 2, 3, 4... are
respectively stored in `oo_.dr.g_0`, `oo_.dr.g_1`, `oo_.dr.g_2`,
`oo_.dr.g_3`, `oo_.dr.g_4`... The columns of those matrices correspond
to multidimensional indices of state variables, in such a way that
symmetric elements are never repeated (for more details, see the
description of `oo_.dr.g_3` in the third-order case).

Occasionally binding constraints (OCCBIN)
-----------------------------------------

Dynare allows simulating models with up to two occasionally-binding
constraints by relying on a piecewise linear solution as in *Guerrieri
and Iacoviello (2015)*. It also allows estimating such models employing
either the inversion filter of *Cuba-Borda, Guerrieri, Iacoviello, and
Zhong (2019)* or the piecewise Kalman filter of *Giovannini, Pfeiffer,
and Ratto (2021)*. To trigger computations involving
occasionally-binding constraints requires

1.  defining and naming the occasionally-binding constraints using an
    `occbin_constraints`-block
2.  specifying the model equations for the respective regimes in the
    `model`-block using appropriate equation tags.
3.  potentially specifying a sequence of surprise shocks using a
    `shocks(surprise)`-block
4.  setting up Occbin simulations or estimation with `occbin_setup`
5.  triggering a simulation with `occbin_solver` or running `estimation`
    or `calib_smoother`.

All of these elements are discussed in the following.

::: {.block}
occbin\_constraints ;

The `occbin_constraints`-block specifies the occasionally-binding
constraints. It contains one or two of the following lines:

> name \'STRING\'; bind EXPRESSION; \[relax EXPRESSION;\] \[error\_bind
> EXPRESSION;\] \[error\_relax EXPRESSION;\]

`STRING` is the name of constraint that is used to reference the
constraint in `relax/bind` equation-tags to identify the respective
regime (see below). The `bind`-expresssion is mandatory and defines a
logical condition that is evaluated in the baseline/steady state regime
to check whether the specified constraint becomes binding. In contrast,
the `relax`-expression is optional and specifies a logical condition
that is evaluated in the binding regime to check whether the regime
returns to the baseline/steady state regime. If not specified, Dynare
will simply check in the binding regime whether the `bind`-expression
evaluates to false. However, there are cases where the `bind`-expression
cannot be evaluated in the binding regime(s), because the variables
involved are constant by definition so that e.g. the value of the
Lagrange multiplier on the complementary slackness condition needs to be
checked. In these cases, it is necessary to provide an explicit
condition that can be evaluated in the binding regime that allows to
check whether it should be left.

Note that the baseline regime denotes the steady state of the model
where the economy will settle in the long-run without shocks. For that
matter, it may be one where e.g. a borrowing constraint is binding. In
that type of setup, the `bind`-condition is used to specify the
condition when this borrowing constraint becomes non-binding so that the
alternative regime is entered.

Three things are important to keep in mind when specifying the
expressions. First, feasible expressions may only contain
contemporaneous endogenous variables. If you want to include leads/lags
or exogenous variables, you need to define an auxiliary variable.
Second, Dynare will at the current stage not linearly approximate the
entered expressions. Because Occbin will work with a linearized model,
consistency will often require the user to enter a linearized
constraint. Otherwise, the condition employed for checking constraint
violations may differ from the one employed within model simulations
based on the piecewise-linear model solution. Third, in contrast to the
original Occbin replication codes, the variables used in expressions are
not automatically demeaned, i.e. they refer to the levels, not
deviations from the steady state. To access the steady state level of a
variable, the `STEADY_STATE()`-operator can be used.

Finally, it\'s worth keeping in mind that for each simulation period,
Occbin will check the respective conditions for whether the current
regime should be left. Small numerical differences from the cutoff point
for a regime can sometimes lead to oscillations between regimes and
cause a spurious periodic solution. Such cases may be prevented by
introducing a small buffer between the two regimes, e.g.

    occbin_constraints;
    name 'ELB'; bind inom <= iss-1e8; relax inom > iss+1e-8;
    end;

The `error_bind` and `error_relax`-options are optional and allow
specifying numerical criteria for the size of the respective constraint
violations employed in numerical routines. By default, Dynare will
simply use the absolute value of the `bind` and `relax` inequalities.
But occasionnally, user-specified expressions perform better.

*Example*

>     occbin_constraints;
>         name 'IRR'; bind log_Invest-log(steady_state(Invest))<log(phi); relax Lambda<0;
>         name 'INEG'; bind log_Invest-log(steady_state(Invest))<0;
>     end;
>
> IRR is a constraint for irreversible investment that becomes binding
> if investment drops below its steady state by more than 0.025 percent
> in the non-binding regime. The constraint will be relaxed whenever the
> associated Lagrange multiplier `Lambda` in the binding regime becomes
> negative. Note that the constraint here takes on a linear form to be
> consistent with a piecewise linear model solution

The specification of the model equations belonging to the respective
regimes is done in the `model`-block, with equation tags indicating to
which regime a particular equation belongs. All equations that differ
across regimes must have a `name`-tag attached to them that allows
uniquely identifying different versions of the same equation. The name
of the constraints specified is then used in conjunction with a `bind`
or `relax` tag to indicate to which regime a particular equation
belongs. In case of more than one occasionally-binding constraint, if an
equation belongs to several regimes (e.g. both constraints binding), the
constraint name tags must be separated by a comma. If only one name tag
is present, the respective equation is assumed to hold for both states
of the other constraint.

*Example*

>     [name='investment',bind='IRR,INEG']
>     (log_Invest - log(phi*steady_state(Invest))) = 0;
>     [name='investment',relax='IRR']
>     Lambda=0;
>     [name='investment',bind='IRR',relax='INEG']
>     (log_Invest - log(phi*steady_state(Invest))) = 0;
>
> The three entered equations for the investment condition define the
> model equation for all four possible combinations of the two
> constraints. The first equation defines the model equation in the
> regime where both the IRR and INEG constraint are binding. The second
> equation defines the model equation for the regimes where the IRR
> constraint is non-binding, regardless of whether the INEG constraint
> is binding or not. Finally, the last equation defines the model
> equation for the final regime where the IRR constraint is binding, but
> the INEG one is not.
:::

::: {.block}
shocks(surprise) ; shocks(surprise,overwrite);

The `shocks(surprise)`-block allows specifying a sequence of temporary
changes in the value of exogenous variables that in each period come as
a surprise to agents, i.e. are not anticipated. Note that to actually
use the specified shocks in subsequent commands like `occbin_solver`,
the block needs to be followed by a call to `occbin_setup`.

The block mirrors the perfect foresight syntax in that it should contain
one or more occurrences of the following group of three lines:

    var VARIABLE_NAME;
    periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
    values DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;

*Example* (with vector values and overwrite option)

    shockssequence = randn(100,1)*0.02;

    shocks(surprise,overwrite);
    var epsilon;
    periods 1:100;
    values (shockssequence);
    end;
:::

::: {.command}
occbin\_setup ; occbin\_setup (OPTIONS\...);

Prepares a simulation with occasionally binding constraints. This
command will also translate the contents of a `shocks(surprise)`-block
for use in subsequent commands.

In order to conduct `estimation` with occasionally binding constraints,
it needs to be prefaced by a call to `occbin_setup` to trigger the use
of either the inversion filter or the piecewise Kalman filter (default).
An issue that can arise in the context of estimation is a structural
shock dropping out of the model in a particular regime. For example, at
the zero lower bound on interest rates, the monetary policy shock in the
Taylor rule will not appear anymore. This may create a problem of
stochastic singularity if there are then more observables than shocks.
To avoid this issue, the data points for the zero interest rate should
be set to NaN and the standard deviation of the associated shock set to
0 for the corresponding periods using the
`heteroskedastic_shocks`-block.

Note that models with unit roots will require the user to specify the
`diffuse_filter`-option as otherwise Blanchard-Kahn errors will be
triggered. For the piecewise Kalman filter, the initialization steps in
the diffuse filter will always rely on the model solved for the baseline
regime, without checking whether this is the actual regime in the first
period(s).

*Example*

>     occbin_setup(likelihood_inversion_filter,smoother_inversion_filter);
>     estimation(smoother,heteroskedastic_filter,...);

The above piece of code sets up an estimation employing the inversion
filter for both the likelihood evaluation and the smoother, while also
accounting for `heteroskedastic_shocks` using the
`heteroskedastic_filter`-option.

Be aware that Occbin has largely command-specific options, i.e. there
are separate options to control the behavior of Occbin when called by
the smoother or when computing the likelihood. These latter commands
will not inherit the options potentially previously set for simulations.

*Options*

::: {.option}
simul\_periods = INTEGER

Number of periods of the simulation. Default: 100.
:::

::: {.option}
simul\_maxit = INTEGER

Maximum number of iterations when trying to find the regimes of the
piecewise solution. Default: 30.
:::

::: {.option}
simul\_check\_ahead\_periods = INTEGER

Number of periods for which to check ahead for return to the baseline
regime. This number should be chosen large enough, because Occbin
requires the simulation to return to the baseline regime at the end of
time. Default: 200.
:::

::: {.option}
simul\_curb\_retrench

Instead of basing the initial regime guess for the current iteration on
the last iteration, update the guess only one period at a time. This
will slow down the iterations, but may lead to more robust convergence
behavior. Default: not enabled.
:::

::: {.option}
simul\_periodic\_solution

Accept a periodic solution where the solution alternates between two
sets of results across iterations, i.e. is not found to be unique. This
is sometimes caused by spurious numerical errors that lead to
oscillations between regiems and may be prevented by allowing for a
small buffer in regime transitions. Default: not enabled.
:::

::: {.option}
simul\_debug

Provide additional debugging information during solving. Default: not
enabled.
:::

::: {.option}
smoother\_periods = INTEGER

Number of periods employed during the simulation when called by the
smoother (equivalent of `simul_periods`). Default: 100.
:::

::: {.option}
smoother\_maxit = INTEGER

Maximum number of iterations employed during the simulation when called
by the smoother (equivalent of `simul_maxit`). Default: 30.
:::

::: {.option}
smoother\_check\_ahead\_periods = INTEGER

Number of periods for which to check ahead for return to the baseline
regime during the simulation when called by the smoother (equivalent of
`simul_check_ahead_periods`). Default: 200.
:::

::: {.option}
smoother\_curb\_retrench

Have the smoother invoke the `simul_curb_retrench`-option during
simulations. Default: not enabled.
:::

::: {.option}
smoother\_periodic\_solution

Accept periodic solution where solution alternates between two sets of
results (equivalent of `simul_periodic_solution`). Default: not enabled.
:::

::: {.option}
likelihood\_periods = INTEGER

Number of periods employed during the simulation when computing the
likelihood (equivalent of `simul_periods`). Default: 100.
:::

::: {.option}
likelihood\_maxit = INTEGER

Maximum number of iterations employed during the simulation when
computing the likelihood (equivalent of `simul_maxit`). Default: 30.
:::

::: {.option}
likelihood\_check\_ahead\_periods = INTEGER

Number of periods for which to check ahead for return to the baseline
regime during the simulation when computing the likelihood (equivalent
of `simul_check_ahead_periods`). Default: 200.
:::

::: {.option}
likelihood\_curb\_retrench

Have the likelihood computation invoke the `simul_curb_retrench`-option
during simulations. Default: not enabled.
:::

::: {.option}
likelihood\_periodic\_solution

Accept periodic solution where solution alternates between two sets of
results (equivalent of `simul_periodic_solution`). Default: not enabled.
:::

::: {.option}
likelihood\_inversion\_filter

Employ the inversion filter of *Cuba-Borda, Guerrieri, Iacoviello, and
Zhong (2019)* when estimating the model. Default: not enabled.
:::

::: {.option}
likelihood\_piecewise\_kalman\_filter

Employ the piecewise Kalman filter of *Giovannini, Pfeiffer, and Ratto
(2021)* when estimating the model. Note that this filter is incompatible
with univariate Kalman filters, i.e. `kalman_algo=2,4`. Default:
enabled.
:::

::: {.option}
likelihood\_max\_kalman\_iterations

Maximum number of iterations of the outer loop for the piecewise Kalman
filter. Default: 10.
:::

::: {.option}
smoother\_inversion\_filter

Employ the inversion filter of *Cuba-Borda, Guerrieri, Iacoviello, and
Zhong (2019)* when running the smoother. The underlying assumption is
that the system starts at the steady state. In this case, the inversion
filter will provide the required smoother output. Default: not enabled.
:::

::: {.option}
smoother\_piecewise\_kalman\_filter

Employ the piecewise Kalman filter of *Giovannini, Pfeiffer, and Ratto
(2021)* when running the smoother. Default: enabled.
:::

::: {.option}
filter\_use\_relaxation

Triggers relaxation within the guess and verify algorithm used in the
update step of the piecewise Kalman filter. When old and new guess
regime differ to much, use a new guess closer to the previous guess. In
case of multiple solutions, tends to provide an occasionally binding
regime with a shorter duration (typically preferable). Specifying this
option may slow down convergence. Default: not enabled.
:::
:::

> *Output*
>
> > The paths for the exogenous variables are stored into
> > `options_.occbin.simul.SHOCKS`.

::: {.command}
occbin\_solver ; occbin\_solver (OPTIONS\...);

Computes a simulation with occasionally-binding constraints based on a
piecewise-linear solution.

Note that `occbin_setup` must be called before this command in order for
the simulation to take into account previous
`shocks(surprise)`-commands.

*Options*

::: {.option}
simul\_periods = INTEGER

See `simul_periods <simul_periods = INTEGER>`{.interpreted-text
role="opt"}.
:::

::: {.option}
simul\_maxit = INTEGER

See `simul_maxit <simul_maxit = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
simul\_check\_ahead\_periods = INTEGER

See
`simul_check_ahead_periods <simul_check_ahead_periods = INTEGER>`{.interpreted-text
role="opt"}.
:::

::: {.option}
simul\_curb\_retrench

See `simul_curb_retrench`{.interpreted-text role="opt"}.
:::

::: {.option}
simul\_debug

See `simul_debug`{.interpreted-text role="opt"}.
:::

*Output*

The command outputs various objects into `oo_.occbin`.
:::

::: {.matvar}
[oo]().occbin.simul.piecewise

Matrix storing the simulations based on the piecewise-linear solution.
The variables are arranged column by column, in order of declaration (as
in `M_.endo_names`), while the the rows correspond to the
`simul_periods`.
:::

::: {.matvar}
[oo]().occbin.simul.linear

Matrix storing the simulations based on the linear solution, i.e.
ignoring the occasionally binding constraint(s). The variables are
arranged column by column, in order of declaration (as in
`M_.endo_names`), while the the rows correspond to the `simul_periods`.
:::

::: {.matvar}
[oo]().occbin.simul.shocks\_sequence

Matrix storing the shock sequence employed during the simulation. The
shocks are arranged column by column, with their order in `M_.exo_names`
stored in `oo_.occbin.exo_pos`. The the rows correspond to the number of
shock periods specified in a [surprise(shocks)]{.title-ref}-block, which
may be smaller than `simul_periods`.
:::

::: {.matvar}
[oo]().occbin.simul.regime\_history

Structure storing information on the regime history, conditional on the
shock that happened in the respective period (stored along the rows).
`type` is equal to either `smoother` or `simul`, depending on whether
the output comes from a run of simulations or the smoother. The subfield
`regime` contains a vector storing the regime state, while the the
subfield `regimestart` indicates the expected start of the respective
regime state. For example, if row 40 contains `[1,0]` for `regime2` and
`[1,6]` for `regimestart2`, it indicates that - after the shock in
period 40 has occurred - the second constraint became binding (1) and is
expected to revert to non-binding (0) after six periods including the
current one, i.e. period 45.
:::

::: {.matvar}
[oo]().occbin.simul.ys

Vector of steady state values
:::

::: {.command}
occbin\_graph \[VARIABLE\_NAME\...\]; occbin\_graph (OPTIONS\...)
\[VARIABLE\_NAME\...\];

Plots a graph comparing the simulation results of the piecewise-linear
solution with the occasionally binding contraints to the linear solution
ignoring the constraint.

*Options*

::: {.option}
noconstant

Omit the steady state in the graphs.
:::
:::

::: {.command}
occbin\_write\_regimes ; occbin\_write\_regimes (OPTIONS\...);

Write the information on the regime history stored in
`oo_.occbin.simul.regime_history` or
``[oo]().occbin.smoother.regime\_history`into an Excel file stored in the`FILENAME/Output`-folder.  *Options*  .. option:: periods = INTEGER     Number of periods for which to write the expected regime durations. Default: write all    available periods.  .. option:: filename = FILENAME     Name of the Excel-file to write. Default:`FILENAME\_occbin\_regimes`.  .. option:: simul     Selects the regime history from the last run of simulations. Default: enabled.  .. option:: smoother     Selects the regime history from the last run of the smoother. Default: use`simul\`\`.
:::

Estimation based on likelihood {#estim}
------------------------------

Provided that you have observations on some endogenous variables, it is
possible to use Dynare to estimate some or all parameters. Both maximum
likelihood (as in *Ireland (2004)*) and Bayesian techniques (as in
*Fernández-Villaverde and Rubio-Ramírez (2004)*, *Rabanal and
Rubio-Ramirez (2003)*, *Schorfheide (2000)* or *Smets and Wouters
(2003)*) are available. Using Bayesian methods, it is possible to
estimate DSGE models, VAR models, or a combination of the two techniques
called DSGE-VAR.

Note that in order to avoid stochastic singularity, you must have at
least as many shocks or measurement errors in your model as you have
observed variables.

The estimation using a first order approximation can benefit from the
block decomposition of the model (see `block`{.interpreted-text
role="opt"}).

::: {#varobs}
::: {.command}
varobs VARIABLE\_NAME\...;

This command lists the name of observed endogenous variables for the
estimation procedure. These variables must be available in the data file
(see `estimation_cmd <estim-comm>`{.interpreted-text role="ref"}).

Alternatively, this command is also used in conjunction with the
`partial_information` option of `stoch_simul`, for declaring the set of
observed variables when solving the model under partial information.

Only one instance of `varobs` is allowed in a model file. If one needs
to declare observed variables in a loop, the macro processor can be used
as shown in the second example below.

*Example*

>     varobs C y rr;
>
> Declares endogenous variables `C`, `y` and `rr` as observed variables.

*Example* (with a macro processor loop)

>     varobs
>     @#for co in countries
>     GDP_@{co}
>     @#endfor
>     ;
:::
:::

::: {.block}
observation\_trends ;

This block specifies linear trends for observed variables as functions
of model parameters. In case the `loglinear` option is used, this
corresponds to a linear trend in the logged observables, i.e. an
exponential trend in the level of the observables.

Each line inside of the block should be of the form:

    VARIABLE_NAME(EXPRESSION);

In most cases, variables shouldn't be centered when `observation_trends`
is used.

*Example*

>     observation_trends;
>     Y (eta);
>     P (mu/eta);
>     end;
:::

::: {.block}
estimated\_params ; estimated\_params (overwrite) ;

This block lists all parameters to be estimated and specifies bounds and
priors as necessary.

Each line corresponds to an estimated parameter.

In a maximum likelihood or a method of moments estimation, each line
follows this syntax:

    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME
    , INITIAL_VALUE [, LOWER_BOUND, UPPER_BOUND ];

In a Bayesian MCMC or a penalized method of moments estimation, each
line follows this syntax:

    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME | DSGE_PRIOR_WEIGHT
    [, INITIAL_VALUE [, LOWER_BOUND, UPPER_BOUND]], PRIOR_SHAPE,
    PRIOR_MEAN, PRIOR_STANDARD_ERROR [, PRIOR_3RD_PARAMETER [,
    PRIOR_4TH_PARAMETER [, SCALE_PARAMETER ] ] ];

The first part of the line consists of one of the four following
alternatives:

-   `stderr VARIABLE_NAME`

    Indicates that the standard error of the exogenous variable
    VARIABLE\_NAME, or of the observation error/measurement errors
    associated with endogenous observed variable VARIABLE\_NAME, is to
    be estimated.

-   `corr VARIABLE_NAME1, VARIABLE_NAME2`

    Indicates that the correlation between the exogenous variables
    VARIABLE\_NAME1 and VARIABLE\_NAME2, or the correlation of the
    observation errors/measurement errors associated with endogenous
    observed variables VARIABLE\_NAME1 and VARIABLE\_NAME2, is to be
    estimated. Note that correlations set by previous shocks-blocks or
    estimation-commands are kept at their value set prior to estimation
    if they are not estimated again subsequently. Thus, the treatment is
    the same as in the case of deep parameters set during model
    calibration and not estimated.

-   `PARAMETER_NAME`

    The name of a model parameter to be estimated

-   `DSGE_PRIOR_WEIGHT`

    Special name for the weigh of the DSGE model in DSGE-VAR model.

The rest of the line consists of the following fields, some of them
being optional:

::: {.option}
INITIAL\_VALUE

Specifies a starting value for the posterior mode optimizer or the
maximum likelihood estimation. If unset, defaults to the prior mean.
:::

::: {.option}
LOWER\_BOUND

Specifies a lower bound for the parameter value in maximum likelihood
estimation. In a Bayesian estimation context, sets a lower bound only
effective while maximizing the posterior kernel. This lower bound does
not modify the shape of the prior density, and is only aimed at helping
the optimizer in identifying the posterior mode (no consequences for the
MCMC). For some prior densities (namely inverse gamma, gamma, uniform,
beta or Weibull) it is possible to shift the support of the prior
distributions to the left or the right using
`prior_3rd_parameter <PRIOR_3RD_PARAMETER>`{.interpreted-text
role="opt"}. In this case the prior density is effectively modified
(note that the truncated Gaussian density is not implemented in Dynare).
If unset, defaults to minus infinity (ML) or the natural lower bound of
the prior (Bayesian estimation).
:::

::: {.option}
UPPER\_BOUND

Same as `lower_bound`, but specifying an upper bound instead.
:::

::: {.option}
PRIOR\_SHAPE

A keyword specifying the shape of the prior density. The possible values
are: `beta_pdf`, `gamma_pdf`, `normal_pdf`, `uniform_pdf`,
`inv_gamma_pdf`, `inv_gamma1_pdf`, `inv_gamma2_pdf` and `weibull_pdf`.
Note that `inv_gamma_pdf` is equivalent to `inv_gamma1_pdf`.
:::

::: {.option}
PRIOR\_MEAN

The mean of the prior distribution.
:::

::: {.option}
PRIOR\_STANDARD\_ERROR

The standard error of the prior distribution.
:::

::: {.option}
PRIOR\_3RD\_PARAMETER

A third parameter of the prior used for generalized beta distribution,
generalized gamma, generalized Weibull and for the uniform distribution.
Default: `0`.
:::

::: {.option}
PRIOR\_4TH\_PARAMETER

A fourth parameter of the prior used for generalized beta distribution
and for the uniform distribution. Default: `1`.
:::

::: {.option}
SCALE\_PARAMETER

A parameter specific scale parameter for the jumping distribution's
covariance matrix of the Metropolis-Hasting algorithm.
:::

Note that INITIAL\_VALUE, LOWER\_BOUND, UPPER\_BOUND, PRIOR\_MEAN,
PRIOR\_STANDARD\_ERROR, PRIOR\_3RD\_PARAMETER, PRIOR\_4TH\_PARAMETER and
SCALE\_PARAMETER can be any valid EXPRESSION. Some of them can be empty,
in which Dynare will select a default value depending on the context and
the prior shape.

In case of the uniform distribution, it can be specified either by
providing an upper and a lower bound using
`PRIOR_3RD_PARAMETER`{.interpreted-text role="opt"} and
`PRIOR_4TH_PARAMETER`{.interpreted-text role="opt"} or via mean and
standard deviation using `PRIOR_MEAN`{.interpreted-text role="opt"},
`PRIOR_STANDARD_ERROR`{.interpreted-text role="opt"}. The other two will
automatically be filled out. Note that providing both sets of
hyperparameters will yield an error message.

As one uses options more towards the end of the list, all previous
options must be filled: for example, if you want to specify
SCALE\_PARAMETER, you must specify `PRIOR_3RD_PARAMETER` and
`PRIOR_4TH_PARAMETER`. Use empty values, if these parameters don't
apply.

*Example*

>     corr eps_1, eps_2, 0.5,  ,  , beta_pdf, 0, 0.3, -1, 1;
>
> Sets a generalized beta prior for the correlation between `eps_1` and
> `eps_2` with mean `0` and variance `0.3`. By setting
> `PRIOR_3RD_PARAMETER` to `-1` and `PRIOR_4TH_PARAMETER` to `1` the
> standard beta distribution with support `[0,1]` is changed to a
> generalized beta with support `[-1,1]`. Note that LOWER\_BOUND and
> UPPER\_BOUND are left empty and thus default to `-1` and `1`,
> respectively. The initial value is set to `0.5`.

*Example*

>     corr eps_1, eps_2, 0.5,  -0.5,  1, beta_pdf, 0, 0.3, -1, 1;
>
> Sets the same generalized beta distribution as before, but now
> truncates this distribution to `[-0.5,1]` through the use of
> LOWER\_BOUND and UPPER\_BOUND.

*Parameter transformation*

Sometimes, it is desirable to estimate a transformation of a parameter
appearing in the model, rather than the parameter itself. It is of
course possible to replace the original parameter by a function of the
estimated parameter everywhere is the model, but it is often
unpractical.

In such a case, it is possible to declare the parameter to be estimated
in the parameters statement and to define the transformation, using a
pound sign (\#) expression (see `model-decl`{.interpreted-text
role="ref"}).

*Example*

>     parameters bet;
>
>     model;
>     # sig = 1/bet;
>     c = sig*c(+1)*mpk;
>     end;
>
>     estimated_params;
>     bet, normal_pdf, 1, 0.05;
>     end;

It is possible to have several `estimated_params` blocks. By default,
subsequent blocks are concatenated with the previous ones; this can be
useful when building models in a modular fashion (see also
`estimated_params_remove`{.interpreted-text role="bck"} for that use
case). However, if an `estimated_params` block has the `overwrite`
option, its contents becomes the new list of estimated parameters,
cancelling previous blocks; this can be useful when doing several
estimations in a single `.mod` file.
:::

::: {.block}
estimated\_params\_init ; estimated\_params\_init (OPTIONS\...);

This block declares numerical initial values for the optimizer when
these ones are different from the prior mean. It should be specified
after the `estimated_params` block as otherwise the specified starting
values are overwritten by the latter.

Each line has the following syntax:

    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME, INITIAL_VALUE;

*Options*

::: {.option}
use\_calibration

For not specifically initialized parameters, use the deep parameters and
the elements of the covariance matrix specified in the `shocks` block
from calibration as starting values for estimation. For components of
the `shocks` block that were not explicitly specified during calibration
or which violate the prior, the prior mean is used.
:::

See `estimated_params`{.interpreted-text role="bck"}, for the meaning
and syntax of the various components.
:::

::: {.block}
estimated\_params\_bounds ;

This block declares lower and upper bounds for parameters in maximum
likelihood estimation.

Each line has the following syntax:

    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME, LOWER_BOUND, UPPER_BOUND;

See `estimated_params`{.interpreted-text role="bck"}, for the meaning
and syntax of the various components.
:::

::: {.block}
estimated\_params\_remove ;

This block partially undoes the effect of a previous
`estimated_params`{.interpreted-text role="bck"} block, by removing some
parameters from the estimation.

Each line has the following syntax:

    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME;
:::

::: {#estim-comm}
::: {.command}
estimation \[VARIABLE\_NAME\...\]; estimation (OPTIONS\...)
\[VARIABLE\_NAME\...\];

This command runs Bayesian or maximum likelihood estimation.

The following information will be displayed by the command:

-   Results from posterior optimization (also for maximum likelihood)
-   Marginal log data density
-   Posterior mean and highest posterior density interval (shortest
    credible set) from posterior simulation
-   Convergence diagnostic table when only one MCM chain is used or
    Metropolis-Hastings convergence graphs documented in
    *Pfeifer (2014)* in case of multiple MCM chains
-   Table with numerical inefficiency factors of the MCMC
-   Graphs with prior, posterior, and mode
-   Graphs of smoothed shocks, smoothed observation errors, smoothed and
    historical variables

Note that the posterior moments, smoothed variables, k-step ahead
filtered variables and forecasts (when requested) will only be computed
on the variables listed after the `estimation` command. Alternatively,
one can choose to compute these quantities on all endogenous or on all
observed variables (see `consider_all_endogenous`,
`consider_all_endogenous_and_auxiliary`, and `consider_only_observed`
options below). If no variable is listed after the estimation command,
then Dynare will interactively ask which variable set to use.

Also, during the MCMC (Bayesian estimation with `mh_replic` $>0$) a
(graphical or text) waiting bar is displayed showing the progress of the
Monte-Carlo and the current value of the acceptance ratio. Note that if
the `load_mh_file` option is used (see below) the reported acceptance
ratio does not take into account the draws from the previous MCMC. In
the literature there is a general agreement for saying that the
acceptance ratio should be close to one third or one quarter. If this
not the case, you can stop the MCMC (`Ctrl-C`) and change the value of
option `mh_jscale` (see below).

Note that by default Dynare generates random numbers using the algorithm
`mt199937ar` (i.e. Mersenne Twister method) with a seed set equal to
`0`. Consequently the MCMCs in Dynare are deterministic: one will get
exactly the same results across different Dynare runs (*ceteris
paribus*). For instance, the posterior moments or posterior densities
will be exactly the same. This behaviour allows to easily identify the
consequences of a change on the model, the priors or the estimation
options. But one may also want to check that across multiple runs, with
different sequences of proposals, the returned results are almost
identical. This should be true if the number of iterations (i.e. the
value of `mh_replic`) is important enough to ensure the convergence of
the MCMC to its ergodic distribution. In this case the default behaviour
of the random number generators in not wanted, and the user should set
the seed according to the system clock before the estimation command
using the following command:

    set_dynare_seed('clock');

so that the sequence of proposals will be different across different
runs.

Finally, Dynare does not always properly distinguish between maximum
likelihood and Bayesian estimation in its field names. While there is an
important conceptual distinction between frequentist confidence
intervals and Bayesian highest posterior density intervals (HPDI) as
well as between posterior density and likelilhood, Dynare sometimes uses
the Bayesian terms as a stand-in in its display of maximum likelihood
results. An example is the storage of the output of the
`forecast`-option of `estimation` with ML, which will use
`HPDinf/HPDsup` to denote the confidence interval.

*Algorithms*

The Monte Carlo Markov Chain (MCMC) diagnostics are generated by the
estimation command if
`mh_replic <mh_replic = INTEGER>`{.interpreted-text role="opt"} is
larger than 2000 and if option `nodiagnostic`{.interpreted-text
role="opt"} is not used. If
`mh_nblocks <mh_nblocks = INTEGER>`{.interpreted-text role="opt"} is
equal to one, the convergence diagnostics of *Geweke (1992,1999)* is
computed. It uses a chi-square test to compare the means of the first
and last draws specified by `geweke_interval
<geweke_interval = [DOUBLE DOUBLE]>`{.interpreted-text role="opt"} after
discarding the burn-in of `mh_drop <mh_drop = DOUBLE>`{.interpreted-text
role="opt"}. The test is computed using variance estimates under the
assumption of no serial correlation as well as using tapering windows
specified in `taper_steps
<taper_steps = [INTEGER1 INTEGER2 ...]>`{.interpreted-text role="opt"}.
If `mh_nblocks
<mh_nblocks = INTEGER>`{.interpreted-text role="opt"} is larger than 1,
the convergence diagnostics of *Brooks and Gelman (1998)* are used
instead. As described in section 3 of *Brooks and Gelman (1998)* the
univariate convergence diagnostics are based on comparing pooled and
within MCMC moments (Dynare displays the second and third order moments,
and the length of the Highest Probability Density interval covering 80%
of the posterior distribution). Due to computational reasons, the
multivariate convergence diagnostic does not follow *Brooks and Gelman
(1998)* strictly, but rather applies their idea for univariate
convergence diagnostics to the range of the posterior likelihood
function instead of the individual parameters. The posterior kernel is
used to aggregate the parameters into a scalar statistic whose
convergence is then checked using the *Brooks and Gelman (1998)*
univariate convergence diagnostic.

The inefficiency factors are computed as in *Giordano et al.(2011)*
based on Parzen windows as in e.g. *Andrews (1991)*.

*Options*

::: {#dataf}
::: {.option}
datafile = FILENAME

The datafile: a `.m` file, a `.mat` file, a `.csv` file, or a
`.xls/.xlsx` file (under Octave, the
[io](https://octave.sourceforge.io/io/) package from Octave-Forge is
required for the `.csv` and `.xlsx` formats and the `.xls` file
extension is not supported). Note that the base name (i.e. without
extension) of the datafile has to be different from the base name of the
model file. If there are several files named FILENAME, but with
different file endings, the file name must be included in quoted strings
and provide the file ending like:

    estimation(datafile='../fsdat_simul.mat',...);
:::
:::

::: {.option}
dirname = FILENAME

Directory in which to store `estimation` output. To pass a subdirectory
of a directory, you must quote the argument. Default: `<mod_file>`.
:::

::: {.option}
xls\_sheet = QUOTED\_STRING

The name of the sheet with the data in an Excel file.
:::

::: {.option}
xls\_range = RANGE

The range with the data in an Excel file. For example,
`xls_range=B2:D200`.
:::

::: {.option}
nobs = INTEGER

The number of observations following `first_obs <first_obs
= [INTEGER1:INTEGER2]>`{.interpreted-text role="opt"} to be used.
Default: all observations in the file after `first_obs`.
:::

::: {.option}
nobs = \[INTEGER1:INTEGER2\]

Runs a recursive estimation and forecast for samples of size ranging of
`INTEGER1` to `INTEGER2`. Option `forecast` must also be specified. The
forecasts are stored in the `RecursiveForecast` field of the results
structure (see
`RecursiveForecast <oo_.RecursiveForecast>`{.interpreted-text
role="mvar"}). The respective results structures `oo_` are saved in
`oo_recursive_` (see `oo_recursive_`{.interpreted-text role="mvar"}) and
are indexed with the respective sample length.
:::

::: {.option}
first\_obs = INTEGER

The number of the first observation to be used. In case of estimating a
DSGE-VAR, `first_obs` needs to be larger than the number of lags.
Default: `1`.
:::

::: {.option}
first\_obs = \[INTEGER1:INTEGER2\]

Runs a rolling window estimation and forecast for samples of fixed size
`nobs` starting with the first observation ranging from `INTEGER1` to
`INTEGER2`. Option `forecast` must also be specified. This option is
incompatible with requesting recursive forecasts using an expanding
window (see `nobs
<nobs = [INTEGER1:INTEGER2]>`{.interpreted-text role="opt"}). The
respective results structures `oo_` are saved in `oo_recursive_` (see
`oo_recursive_`{.interpreted-text role="mvar"}) and are indexed with the
respective first observation of the rolling window.
:::

::: {.option}
prefilter = INTEGER

A value of 1 means that the estimation procedure will demean each data
series by its empirical mean. If the `loglinear
<logl>`{.interpreted-text role="ref"} option without the
`logdata`{.interpreted-text role="opt"} option is requested, the data
will first be logged and then demeaned. Default: `0`, i.e. no
prefiltering.
:::

::: {.option}
presample = INTEGER

The number of observations after `first_obs <first_obs =
[INTEGER1:INTEGER2]>`{.interpreted-text role="opt"} to be skipped before
evaluating the likelihood. These presample observations do not enter the
likelihood, but are used as a training sample for starting the Kalman
filter iterations. This option is incompatible with estimating a
DSGE-VAR. Default: `0`.
:::

::: {#logl}
::: {.option}
loglinear

Computes a log-linear approximation of the model instead of a linear
approximation. As always in the context of estimation, the data must
correspond to the definition of the variables used in the model (see
*Pfeifer (2013)* for more details on how to correctly specify
observation equations linking model variables and the data). If you
specify the loglinear option, Dynare will take the logarithm of both
your model variables and of your data as it assumes the data to
correspond to the original non-logged model variables. The displayed
posterior results like impulse responses, smoothed variables, and
moments will be for the logged variables, not the original un-logged
ones. Default: computes a linear approximation.
:::
:::

::: {.option}
logdata

Dynare applies the $log$ transformation to the provided data if a
log-linearization of the model is requested
(`loglinear`{.interpreted-text role="opt"}) unless `logdata` option is
used. This option is necessary if the user provides data already in
logs, otherwise the $log$ transformation will be applied twice (this may
result in complex data).
:::

::: {.option}
plot\_priors = INTEGER

Control the plotting of priors.

> `0`
>
> > No prior plot.
>
> `1`
>
> > Prior density for each estimated parameter is plotted. It is
> > important to check that the actual shape of prior densities matches
> > what you have in mind. Ill-chosen values for the prior standard
> > density can result in absurd prior densities.

Default value is `1`.
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}.
:::

::: {.option}
posterior\_nograph

Suppresses the generation of graphs associated with Bayesian IRFs
(`bayesian_irf`{.interpreted-text role="opt"}), posterior smoothed
objects (`smoother`{.interpreted-text role="opt"}), and posterior
forecasts (`forecast`{.interpreted-text role="opt"}).
:::

::: {.option}
posterior\_graph

Re-enables the generation of graphs previously shut off with
`posterior_nograph`{.interpreted-text role="opt"}.
:::

::: {.option}
nodisplay

See `nodisplay`{.interpreted-text role="opt"}.
:::

::: {.option}
graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )

See
`graph_format <graph_format = ( FORMAT, FORMAT... )>`{.interpreted-text
role="opt"}.
:::

::: {.option}
no\_init\_estimation\_check\_first\_obs

Do not check for stochastic singularity in first period. If used,
[ESTIMATION CHECKS]{.title-ref} does not return an error if the check
fails only in first observation. This should only be used when observing
stock variables (e.g. capital) in first period, on top of their
associated flow (e.g. investment). Using this option may lead to a crash
or provide undesired/wrong results for badly specified problems (e.g.
the additional variable observed in first period is not predetermined).

For advanced use only.
:::

::: {.option}
lik\_init = INTEGER

Type of initialization of Kalman filter:

> `1`
>
> > For stationary models, the initial matrix of variance of the error
> > of forecast is set equal to the unconditional variance of the state
> > variables.
>
> `2`
>
> > For nonstationary models: a wide prior is used with an initial
> > matrix of variance of the error of forecast diagonal with 10 on the
> > diagonal (follows the suggestion of *Harvey and Phillips(1979)*).
>
> `3`
>
> > For nonstationary models: use a diffuse filter (use rather the
> > `diffuse_filter` option).
>
> `4`
>
> > The filter is initialized with the fixed point of the Riccati
> > equation.
>
> `5`
>
> > Use i) option 2 for the non-stationary elements by setting their
> > initial variance in the forecast error matrix to 10 on the diagonal
> > and all covariances to 0 and ii) option 1 for the stationary
> > elements.

Default value is 1. For advanced use only.
:::

::: {.option}
lik\_algo = INTEGER

For internal use and testing only.
:::

::: {.option}
conf\_sig = DOUBLE

Level of significance of the confidence interval used for classical
forecasting after estimation. Default: 0.9.
:::

::: {.option}
mh\_conf\_sig = DOUBLE

Confidence/HPD interval used for the computation of prior and posterior
statistics like: parameter distributions, prior/posterior moments,
conditional variance decomposition, impulse response functions, Bayesian
forecasting. Default: `0.9`.
:::

::: {.option}
mh\_replic = INTEGER

Number of replications for each chain of the Metropolis-Hastings
algorithm. The number of draws should be sufficient to achieve
convergence of the MCMC and to meaningfully compute posterior objects.
Default: `20000`.
:::

::: {.option}
sub\_draws = INTEGER

Number of draws from the MCMC that are used to compute posterior
distribution of various objects (smoothed variable, smoothed shocks,
forecast, moments, IRF). The draws used to compute these posterior
moments are sampled uniformly in the estimated empirical posterior
distribution (i.e. draws of the MCMC). `sub_draws` should be smaller
than the total number of MCMC draws available. Default:
`min(posterior_max_subsample_draws, (Total number of draws)*(number of chains) )`.
:::

::: {.option}
posterior\_max\_subsample\_draws = INTEGER

Maximum number of draws from the MCMC used to compute posterior
distribution of various objects (smoothed variable, smoothed shocks,
forecast, moments, IRF), if not overriden by option `sub_draws`.
Default: `1200`.
:::

::: {.option}
mh\_nblocks = INTEGER

Number of parallel chains for Metropolis-Hastings algorithm. Default:
`2`.
:::

::: {.option}
mh\_drop = DOUBLE

The fraction of initially generated parameter vectors to be dropped as a
burn-in before using posterior simulations. Default: `0.5`.
:::

::: {.option}
mh\_jscale = DOUBLE

The scale parameter of the jumping distribution's covariance matrix
(Metropolis-Hastings or TaRB-algorithm). The default value is rarely
satisfactory. This option must be tuned to obtain, ideally, an
acceptance ratio of 25%-33%. Basically, the idea is to increase the
variance of the jumping distribution if the acceptance ratio is too
high, and decrease the same variance if the acceptance ratio is too low.
In some situations it may help to consider parameter-specific values for
this scale parameter. This can be done in the
`estimated_params`{.interpreted-text role="bck"} block.

Note that `mode_compute=6` will tune the scale parameter to achieve an
acceptance rate of `AcceptanceRateTarget<art>`{.interpreted-text
role="ref"}. The resulting scale parameter will be saved into a file
named `MODEL_FILENAME_mh_scale.mat` in the `FILENAME/Output`-folder.
This file can be loaded in subsequent runs via the
`posterior_sampler_options` option
`scale_file <scale-file>`{.interpreted-text role="ref"}. Both
`mode_compute=6` and `scale_file` will overwrite any value specified in
`estimated_params` with the tuned value. Default: `0.2`.

Note also that for the Random Walk Metropolis Hastings algorithm, it is
possible to use option `mh_tune_jscale
<mh_tune_jscale [= DOUBLE]>`{.interpreted-text role="opt"}, to
automatically tune the value of `mh_jscale`. In this case, the
`mh_jscale` option must not be used.
:::

::: {.option}
mh\_init\_scale = DOUBLE

The scale to be used for drawing the initial value of the
Metropolis-Hastings chain. Generally, the starting points should be
overdispersed for the *Brooks and Gelman (1998)* convergence diagnostics
to be meaningful. Default: `2*mh_jscale.`

It is important to keep in mind that `mh_init_scale` is set at the
beginning of Dynare execution, i.e. the default will not take into
account potential changes in `mh_jscale` introduced by either
`mode_compute=6` or the `posterior_sampler_options` option `scale_file
<scale-file>`{.interpreted-text role="ref"}. If `mh_init_scale` is too
wide during initalization of the posterior sampler so that 100 tested
draws are inadmissible (e.g. Blanchard-Kahn conditions are always
violated), Dynare will request user input of a new `mh_init_scale` value
with which the next 100 draws will be drawn and tested. If the
`nointeractive`{.interpreted-text role="opt"} option has been invoked,
the program will instead automatically decrease `mh_init_scale` by 10
percent after 100 futile draws and try another 100 draws. This iterative
procedure will take place at most 10 times, at which point Dynare will
abort with an error message.
:::

::: {.option}
mh\_tune\_jscale \[= DOUBLE\]

Automatically tunes the scale parameter of the jumping distribution\'s
covariance matrix (Metropolis-Hastings), so that the overall acceptance
ratio is close to the desired level. Default value is `0.33`. It is not
possible to match exactly the desired acceptance ratio because of the
stochastic nature of the algorithm (the proposals and the initial
conditions of the markov chains if `mh_nblocks>1`). This option is only
available for the Random Walk Metropolis Hastings algorithm. Must not be
used in conjunction with `mh_jscale = DOUBLE`{.interpreted-text
role="opt"}.
:::

::: {.option}
mh\_tune\_guess = DOUBLE

Specifies the initial value for the `mh_tune_jscale
<mh_tune_jscale [= DOUBLE]>`{.interpreted-text role="opt"} option.
Default: `0.2`. Must not be set if
`mh_tune_jscale <mh_tune_jscale [= DOUBLE]>`{.interpreted-text
role="opt"} is not used.
:::

::: {.option}
mh\_recover

Attempts to recover a Metropolis-Hastings simulation that crashed
prematurely, starting with the last available saved `mh`-file. Shouldn't
be used together with `load_mh_file` or a different `mh_replic` than in
the crashed run. Since Dynare 4.5 the proposal density from the previous
run will automatically be loaded. In older versions, to assure a neat
continuation of the chain with the same proposal density, you should
provide the `mode_file` used in the previous run or the same
user-defined `mcmc_jumping_covariance` when using this option. Note that
under Octave, a neat continuation of the crashed chain with the
respective last random number generator state is currently not
supported.
:::

::: {.option}
mh\_posterior\_mode\_estimation

Skip optimizer-based mode-finding and instead compute the mode based on
a run of a MCMC. The MCMC will start at the prior mode and use the prior
variances to compute the inverse Hessian.
:::

::: {.option}
mode\_file = FILENAME

Name of the file containing previous value for the mode. When computing
the mode, Dynare stores the mode (`xparam1`) and the hessian (`hh`, only
if `cova_compute=1`) in a file called `MODEL_FILENAME_mode.mat` in the
`FILENAME/Output`-folder. After a successful run of the estimation
command, the `mode_file` will be disabled to prevent other function
calls from implicitly using an updated mode-file. Thus, if the `.mod`
file contains subsequent `estimation` commands, the `mode_file` option,
if desired, needs to be specified again.
:::

::: {.option}
mode\_compute = INTEGER \| FUNCTION\_NAME

Specifies the optimizer for the mode computation:

> `0`
>
> > The mode isn't computed. When the `mode_file` option is specified,
> > the mode is simply read from that file.
> >
> > When `mode_file` option is not specified, Dynare reports the value
> > of the log posterior (log likelihood) evaluated at the initial value
> > of the parameters.
> >
> > When `mode_file` is not specified and there is no `estimated_params`
> > block, but the `smoother` option is used, it is a roundabout way to
> > compute the smoothed value of the variables of a model with
> > calibrated parameters.
>
> `1`
>
> > Uses `fmincon` optimization routine (available under MATLAB if the
> > Optimization Toolbox is installed; available under Octave if the
> > [optim](https://octave.sourceforge.io/optim/) package from
> > Octave-Forge, version 1.6 or above, is installed).
>
> `2`
>
> > Uses the continuous simulated annealing global optimization
> > algorithm described in *Corana et al.(1987)* and *Goffe et
> > al.(1994)*.
>
> `3`
>
> > Uses `fminunc` optimization routine (available under MATLAB if the
> > Optimization Toolbox is installed; available under Octave if the
> > [optim](https://octave.sourceforge.io/optim/) package from
> > Octave-Forge is installed).
>
> `4`
>
> > Uses Chris Sims's `csminwel`.
>
> `5`
>
> > Uses Marco Ratto's `newrat`. This value is not compatible with non
> > linear filters or DSGE-VAR models. This is a slice optimizer: most
> > iterations are a sequence of univariate optimization step, one for
> > each estimated parameter or shock. Uses `csminwel` for line search
> > in each step.
>
> `6`
>
> > Uses a Monte-Carlo based optimization routine (see
> > <https://archives.dynare.org/DynareWiki/MonteCarloOptimization> for
> > more details).
>
> `7`
>
> > Uses `fminsearch`, a simplex-based optimization routine (available
> > under MATLAB if the Optimization Toolbox is installed; available
> > under Octave if the optim package from Octave-Forge is installed).
>
> `8`
>
> > Uses Dynare implementation of the Nelder-Mead simplex-based
> > optimization routine (generally more efficient than the MATLAB or
> > Octave implementation available with `mode_compute=7`).
>
> `9`
>
> > Uses the CMA-ES (Covariance Matrix Adaptation Evolution Strategy)
> > algorithm of *Hansen and Kern (2004)*, an evolutionary algorithm for
> > difficult non-linear non-convex optimization.
>
> `10`
>
> > Uses the `simpsa` algorithm, based on the combination of the
> > non-linear simplex and simulated annealing algorithms as proposed by
> > *Cardoso, Salcedo and Feyo de Azevedo (1996)*.
>
> `11`
>
> > This is not strictly speaking an optimization algorithm. The
> > (estimated) parameters are treated as state variables and estimated
> > jointly with the original state variables of the model using a
> > nonlinear filter. The algorithm implemented in Dynare is described
> > in *Liu and West (2001)*, and works with `k` order local
> > approximations of the model.
>
> `12`
>
> > Uses the `particleswarm` optimization routine (available under
> > MATLAB if the Global Optimization Toolbox is installed; not
> > available under Octave).
>
> `13`
>
> > Uses the `lsqnonlin` non-linear least squares optimization routine
> > (available under MATLAB if the Optimization Toolbox is installed;
> > available under Octave if the
> > [optim](https://octave.sourceforge.io/optim/) package from
> > Octave-Forge is installed). Only supported for `method_of_moments`.
>
> `101`
>
> > Uses the SolveOpt algorithm for local nonlinear optimization
> > problems proposed by *Kuntsevich and Kappel (1997)*.
>
> `102`
>
> > Uses `simulannealbnd` optimization routine (available under MATLAB
> > if the Global Optimization Toolbox is installed; not available under
> > Octave)
>
> `FUNCTION_NAME`
>
> > It is also possible to give a FUNCTION\_NAME to this option, instead
> > of an INTEGER. In that case, Dynare takes the return value of that
> > function as the posterior mode.

Default value is `4`.
:::

::: {.option}
silent\_optimizer

Instructs Dynare to run mode computing/optimization silently without
displaying results or saving files in between. Useful when running
loops.
:::

::: {.option}
mcmc\_jumping\_covariance = OPTION

Tells Dynare which covariance to use for the proposal density of the
MCMC sampler. OPTION can be one of the following:

> `hessian`
>
> > Uses the Hessian matrix computed at the mode.
>
> `prior_variance`
>
> > Uses the prior variances. No infinite prior variances are allowed in
> > this case.
>
> `identity_matrix`
>
> > Uses an identity matrix.
>
> `FILENAME`
>
> > Loads an arbitrary user-specified covariance matrix from
> > `FILENAME.mat`. The covariance matrix must be saved in a variable
> > named `jumping_covariance`, must be square, positive definite, and
> > have the same dimension as the number of estimated parameters.

Note that the covariance matrices are still scaled with
`mh_jscale <mh_jscale = DOUBLE>`{.interpreted-text role="opt"}. Default
value is `hessian`.
:::

::: {.option}
mode\_check

Tells Dynare to plot the posterior density for values around the
computed mode for each estimated parameter in turn. This is helpful to
diagnose problems with the optimizer. Note that for `order>1` the
likelihood function resulting from the particle filter is not
differentiable anymore due to the resampling step. For this reason, the
`mode_check` plot may look wiggly.
:::

::: {.option}
mode\_check\_neighbourhood\_size = DOUBLE

Used in conjunction with option `mode_check`, gives the width of the
window around the posterior mode to be displayed on the diagnostic
plots. This width is expressed in percentage deviation. The `Inf` value
is allowed, and will trigger a plot over the entire domain (see also
`mode_check_symmetric_plots`). Default:`0.5`.
:::

::: {.option}
mode\_check\_symmetric\_plots = INTEGER

Used in conjunction with option `mode_check`, if set to `1`, tells
Dynare to ensure that the check plots are symmetric around the posterior
mode. A value of `0` allows to have asymmetric plots, which can be
useful if the posterior mode is close to a domain boundary, or in
conjunction with `mode_check_neighbourhood_size = Inf` when the domain
in not the entire real line. Default: `1`.
:::

::: {.option}
mode\_check\_number\_of\_points = INTEGER

Number of points around the posterior mode where the posterior kernel is
evaluated (for each parameter). Default is `20`.
:::

::: {.option}
prior\_trunc = DOUBLE

Probability of extreme values of the prior density that is ignored when
computing bounds for the parameters. Default: `1e-32`.
:::

::: {.option}
huge\_number = DOUBLE

Value for replacing infinite values in the definition of (prior) bounds
when finite values are required for computational reasons. Default:
`1e7`.
:::

::: {.option}
load\_mh\_file

Tells Dynare to add to previous Metropolis-Hastings simulations instead
of starting from scratch. Since Dynare 4.5 the proposal density from the
previous run will automatically be loaded. In older versions, to assure
a neat continuation of the chain with the same proposal density, you
should provide the `mode_file` used in the previous run or the same
user-defined `mcmc_jumping_covariance` when using this option. Shouldn't
be used together with `mh_recover`. Note that under Octave, a neat
continuation of the chain with the last random number generator state of
the already present draws is currently not supported.
:::

::: {.option}
load\_results\_after\_load\_mh

This option is available when loading a previous MCMC run without adding
additional draws, i.e. when `load_mh_file` is specified with
`mh_replic=0`. It tells Dynare to load the previously computed
convergence diagnostics, marginal data density, and posterior statistics
from an existing `_results` file instead of recomputing them.
:::

::: {.option}
mh\_initialize\_from\_previous\_mcmc

This option allows to pick initial values for new MCMC from a previous
one, where the model specification, the number of estimated parameters,
(some) prior might have changed (so a situation where `load_mh_file`
would not work). If an additional parameter is estimated, it is
automatically initialized from prior\_draw. Note that, if this option is
used to skip the optimization step, you should use a sampling method
which does not require a proposal density, like slice. Otherwise,
optimization should always be done beforehand or a mode file with an
appropriate posterior covariance matrix should be used.
:::

::: {.option}
mh\_initialize\_from\_previous\_mcmc\_directory = FILENAME

If `mh_initialize_from_previous_mcmc` is set, users must provide here
the path to the standard FNAME folder from where to load prior
definitions and last MCMC values to be used to initialize the new MCMC.

Example: if previous project directory is `/my_previous_dir` and FNAME
is `mymodel`, users should set the option as

`mh_initialize_from_previous_mcmc_directory = '/my_previous_dir/mymodel'`

Dynare will then look for the last record file into

`/my_previous_dir/mymodel/metropolis/mymodel_mh_history_<LAST>.mat`

and for the prior definition file into

`/my_previous_dir/mymodel/prior/definition.mat`
:::

::: {.option}
mh\_initialize\_from\_previous\_mcmc\_record = FILENAME

If `mh_initialize_from_previous_mcmc` is set, and whenever the standard
file or directory tree is not applicable to load initial values, users
may directly provide here the path to the record file from which to load
values to be used to initialize the new MCMC.
:::

::: {.option}
mh\_initialize\_from\_previous\_mcmc\_prior = FILENAME

If `mh_initialize_from_previous_mcmc` is set, and whenever the standard
file or directory tree is not applicable to load initial values, users
may directly provide here the path to the prior definition file, to get
info in the priors used in previous MCMC.
:::

::: {.option}
optim = (NAME, VALUE, \...)

A list of NAME and VALUE pairs. Can be used to set options for the
optimization routines. The set of available options depends on the
selected optimization routine (i.e. on the value of option
`mode_compute <mode_compute = INTEGER |
FUNCTION_NAME>`{.interpreted-text role="opt"}):

> `1, 3, 7, 12, 13`
>
> > Available options are given in the documentation of the MATLAB
> > Optimization Toolbox or in Octave's documentation.
>
> `2`
>
> > Available options are:
> >
> > > `'initial_step_length'`
> > >
> > > > Initial step length. Default: `1`.
> > >
> > > `'initial_temperature'`
> > >
> > > > Initial temperature. Default: `15`.
> > >
> > > `'MaxIter'`
> > >
> > > > Maximum number of function evaluations. Default: `100000`.
> > >
> > > `'neps'`
> > >
> > > > Number of final function values used to decide upon termination.
> > > > Default: `10`.
> > >
> > > `'ns'`
> > >
> > > > Number of cycles. Default: `10`.
> > >
> > > `'nt'`
> > >
> > > > Number of iterations before temperature reduction. Default:
> > > > `10`.
> > >
> > > `'step_length_c'`
> > >
> > > > Step length adjustment. Default: `0.1`.
> > >
> > > `'TolFun'`
> > >
> > > > Stopping criteria. Default: `1e-8`.
> > >
> > > `'rt'`
> > >
> > > > Temperature reduction factor. Default: `0.1`.
> > >
> > > `'verbosity'`
> > >
> > > > Controls verbosity of display during optimization, ranging from
> > > > `0` (silent) to `3` (each function evaluation). Default: `1`
>
> `4`
>
> > Available options are:
> >
> > > `'InitialInverseHessian'`
> > >
> > > > Initial approximation for the inverse of the Hessian matrix of
> > > > the posterior kernel (or likelihood). Obviously this
> > > > approximation has to be a square, positive definite and
> > > > symmetric matrix. Default: `'1e-4*eye(nx)'`, where nx is the
> > > > number of parameters to be estimated.
> > >
> > > `'MaxIter'`
> > >
> > > > Maximum number of iterations. Default: `1000`.
> > >
> > > `'NumgradAlgorithm'`
> > >
> > > > Possible values are `2`, `3` and `5`, respectively,
> > > > corresponding to the two, three and five points formula used to
> > > > compute the gradient of the objective function (see *Abramowitz
> > > > and Stegun (1964)*). Values `13` and `15` are more experimental.
> > > > If perturbations on the right and the left increase the value of
> > > > the objective function (we minimize this function) then we force
> > > > the corresponding element of the gradient to be zero. The idea
> > > > is to temporarily reduce the size of the optimization problem.
> > > > Default: `2`.
> > >
> > > `'NumgradEpsilon'`
> > >
> > > > Size of the perturbation used to compute numerically the
> > > > gradient of the objective function. Default: `1e-6`.
> > >
> > > `'TolFun'`
> > >
> > > > Stopping criteria. Default: `1e-7`.
> > >
> > > `'verbosity'`
> > >
> > > > Controls verbosity of display during optimization. Set to `0` to
> > > > set to silent. Default: `1`.
> > >
> > > `'SaveFiles'`
> > >
> > > > Controls saving of intermediate results during optimization. Set
> > > > to `0` to shut off saving. Default: `1`.
>
> `5`
>
> > Available options are:
> >
> > `'Hessian'`
> >
> > > Triggers three types of Hessian computations. `0`: outer product
> > > gradient; `1`: default Dynare Hessian routine; `2`: 'mixed' outer
> > > product gradient, where diagonal elements are obtained using
> > > second order derivation formula and outer product is used for
> > > correlation structure. Both {0} and {2} options require univariate
> > > filters, to ensure using maximum number of individual densities
> > > and a positive definite Hessian. Both {0} and {2} are quicker than
> > > default Dynare numeric Hessian, but provide decent starting values
> > > for Metropolis for large models (option {2} being more accurate
> > > than {0}). Default: `1`.
> >
> > `'MaxIter'`
> >
> > > Maximum number of iterations. Default: `1000`.
> >
> > `'TolFun'`
> >
> > > Stopping criteria. Default: `1e-5` for numerical derivatives,
> > > `1e-7` for analytic derivatives.
> >
> > `'verbosity'`
> >
> > > Controls verbosity of display during optimization. Set to `0` to
> > > set to silent. Default: `1`.
> >
> > `'SaveFiles'`
> >
> > > Controls saving of intermediate results during optimization. Set
> > > to `0` to shut off saving. Default: `1`.
>
> `6`
>
> > Available options are:
> >
> > > ::: {#art}
> > > `'AcceptanceRateTarget'`
> > > :::
> > >
> > > > A real number between zero and one. The scale parameter of the
> > > > jumping distribution is adjusted so that the effective
> > > > acceptance rate matches the value of option
> > > > `'AcceptanceRateTarget'`. Default: `1.0/3.0`.
> > >
> > > `'InitialCovarianceMatrix'`
> > >
> > > > Initial covariance matrix of the jumping distribution. Default
> > > > is `'previous'` if option `mode_file` is used, `'prior'`
> > > > otherwise.
> > >
> > > `'nclimb-mh'`
> > >
> > > > Number of iterations in the last MCMC (climbing mode). Default:
> > > > `200000`.
> > >
> > > `'ncov-mh'`
> > >
> > > > Number of iterations used for updating the covariance matrix of
> > > > the jumping distribution. Default: `20000`.
> > >
> > > `'nscale-mh'`
> > >
> > > > Maximum number of iterations used for adjusting the scale
> > > > parameter of the jumping distribution. Default: `200000`.
> > >
> > > `'NumberOfMh'`
> > >
> > > > Number of MCMC run sequentially. Default: `3`.
>
> `8`
>
> > Available options are:
> >
> > > `'InitialSimplexSize'`
> > >
> > > > Initial size of the simplex, expressed as percentage deviation
> > > > from the provided initial guess in each direction. Default:
> > > > `.05`.
> > >
> > > `'MaxIter'`
> > >
> > > > Maximum number of iterations. Default: `5000`.
> > >
> > > `'MaxFunEvals'`
> > >
> > > > Maximum number of objective function evaluations. No default.
> > >
> > > `'MaxFunvEvalFactor'`
> > >
> > > > Set `MaxFunvEvals` equal to `MaxFunvEvalFactor` times the number
> > > > of estimated parameters. Default: `500`.
> > >
> > > `'TolFun'`
> > >
> > > > Tolerance parameter (w.r.t the objective function). Default:
> > > > `1e-4`.
> > >
> > > `'TolX'`
> > >
> > > > Tolerance parameter (w.r.t the instruments). Default: `1e-4`.
> > >
> > > `'verbosity'`
> > >
> > > > Controls verbosity of display during optimization. Set to `0` to
> > > > set to silent. Default: `1`.
>
> `9`
>
> > Available options are:
> >
> > > `'CMAESResume'`
> > >
> > > > Resume previous run. Requires the `variablescmaes.mat` from the
> > > > last run. Set to `1` to enable. Default: `0`.
> > >
> > > `'MaxIter'`
> > >
> > > > Maximum number of iterations.
> > >
> > > `'MaxFunEvals'`
> > >
> > > > Maximum number of objective function evaluations. Default:
> > > > `Inf`.
> > >
> > > `'TolFun'`
> > >
> > > > Tolerance parameter (w.r.t the objective function). Default:
> > > > `1e-7`.
> > >
> > > `'TolX'`
> > >
> > > > Tolerance parameter (w.r.t the instruments). Default: `1e-7`.
> > >
> > > `'verbosity'`
> > >
> > > > Controls verbosity of display during optimization. Set to `0` to
> > > > set to silent. Default: `1`.
> > >
> > > `'SaveFiles'`
> > >
> > > > Controls saving of intermediate results during optimization. Set
> > > > to `0` to shut off saving. Default: `1`.
>
> `10`
>
> > Available options are:
> >
> > > `'EndTemperature'`
> > >
> > > > Terminal condition w.r.t the temperature. When the temperature
> > > > reaches `EndTemperature`, the temperature is set to zero and the
> > > > algorithm falls back into a standard simplex algorithm. Default:
> > > > `0.1`.
> > >
> > > `'MaxIter'`
> > >
> > > > Maximum number of iterations. Default: `5000`.
> > >
> > > `'MaxFunvEvals'`
> > >
> > > > Maximum number of objective function evaluations. No default.
> > >
> > > `'TolFun'`
> > >
> > > > Tolerance parameter (w.r.t the objective function). Default:
> > > > `1e-4`.
> > >
> > > `'TolX'`
> > >
> > > > Tolerance parameter (w.r.t the instruments). Default: `1e-4`.
> > >
> > > `'verbosity'`
> > >
> > > > Controls verbosity of display during optimization. Set to `0` to
> > > > set to silent. Default: `1`.
>
> `101`
>
> > Available options are:
> >
> > > `'LBGradientStep'`
> > >
> > > > Lower bound for the stepsize used for the difference
> > > > approximation of gradients. Default: `1e-11`.
> > >
> > > `'MaxIter'`
> > >
> > > > Maximum number of iterations. Default: `15000`
> > >
> > > `'SpaceDilation'`
> > >
> > > > Coefficient of space dilation. Default: `2.5`.
> > >
> > > `'TolFun'`
> > >
> > > > Tolerance parameter (w.r.t the objective function). Default:
> > > > `1e-6`.
> > >
> > > `'TolX'`
> > >
> > > > Tolerance parameter (w.r.t the instruments). Default: `1e-6`.
> > >
> > > `'verbosity'`
> > >
> > > > Controls verbosity of display during optimization. Set to `0` to
> > > > set to silent. Default: `1`.
>
> `102`
>
> > Available options are given in the documentation of the MATLAB
> > Global Optimization Toolbox.

*Example*

> To change the defaults of `csminwel` (`mode_compute=4`):
>
>     estimation(..., mode_compute=4,optim=('NumgradAlgorithm',3,'TolFun',1e-5),...);
:::

::: {.option}
nodiagnostic

Does not compute the convergence diagnostics for Metropolis-Hastings.
Default: diagnostics are computed and displayed.
:::

::: {.option}
bayesian\_irf

Triggers the computation of the posterior distribution of IRFs. The
length of the IRFs are controlled by the `irf` option. Results are
stored in `oo_.PosteriorIRF.dsge` (see below for a description of this
variable).
:::

::: {.option}
relative\_irf

See `relative_irf`{.interpreted-text role="opt"}.
:::

::: {.option}
dsge\_var = DOUBLE

Triggers the estimation of a DSGE-VAR model, where the weight of the
DSGE prior of the VAR model is calibrated to the value passed (see *Del
Negro and Schorfheide (2004)*). It represents the ratio of dummy over
actual observations. To assure that the prior is proper, the value must
be bigger than $(k+n)/T$, where $k$ is the number of estimated
parameters, $n$ is the number of observables, and $T$ is the number of
observations.

> NB: The previous method of declaring `dsge_prior_weight` as a
> parameter and then calibrating it is now deprecated and will be
> removed in a future release of Dynare. Some of objects arising during
> estimation are stored with their values at the mode in
> `oo_.dsge_var.posterior_mode`.
:::

::: {.option}
dsge\_var

Triggers the estimation of a DSGE-VAR model, where the weight of the
DSGE prior of the VAR model will be estimated (as in *Adjemian et
al.(2008)*). The prior on the weight of the DSGE prior,
`dsge_prior_weight`, must be defined in the `estimated_params` section.

NB: The previous method of declaring `dsge_prior_weight` as a parameter
and then placing it in `estimated_params` is now deprecated and will be
removed in a future release of Dynare.
:::

::: {.option}
dsge\_varlag = INTEGER

The number of lags used to estimate a DSGE-VAR model. Default: `4`.
:::

::: {.option}
posterior\_sampling\_method = NAME

Selects the sampler used to sample from the posterior distribution
during Bayesian estimation. Default:`’random_walk_metropolis_hastings’`.

> `'random_walk_metropolis_hastings'`
>
> > Instructs Dynare to use the Random-Walk Metropolis-Hastings. In this
> > algorithm, the proposal density is recentered to the previous draw
> > in every step.
>
> `'tailored_random_block_metropolis_hastings'`
>
> > Instructs Dynare to use the Tailored randomized block (TaRB)
> > Metropolis-Hastings algorithm proposed by *Chib and Ramamurthy
> > (2010)* instead of the standard Random-Walk Metropolis-Hastings. In
> > this algorithm, at each iteration the estimated parameters are
> > randomly assigned to different blocks. For each of these blocks a
> > mode-finding step is conducted. The inverse Hessian at this mode is
> > then used as the covariance of the proposal density for a
> > Random-Walk Metropolis-Hastings step. If the numerical Hessian is
> > not positive definite, the generalized Cholesky decomposition of
> > *Schnabel and Eskow (1990)* is used, but without pivoting. The
> > TaRB-MH algorithm massively reduces the autocorrelation in the MH
> > draws and thus reduces the number of draws required to
> > representatively sample from the posterior. However, this comes at a
> > computational cost as the algorithm takes more time to run.
>
> `'independent_metropolis_hastings'`
>
> > Use the Independent Metropolis-Hastings algorithm where the proposal
> > distribution - in contrast to the Random Walk Metropolis-Hastings
> > algorithm - does not depend on the state of the chain.
>
> `'slice'`
>
> > Instructs Dynare to use the Slice sampler of *Planas, Ratto, and
> > Rossi (2015)*. Note that `'slice'` is incompatible with
> > `prior_trunc=0`.
:::

::: {.option}
posterior\_sampler\_options = (NAME, VALUE, \...)

A list of NAME and VALUE pairs. Can be used to set options for the
posterior sampling methods. The set of available options depends on the
selected posterior sampling routine (i.e. on the value of option
`posterior_sampling_method
<posterior_sampling_method = NAME>`{.interpreted-text role="opt"}):

> `'random_walk_metropolis_hastings'`
>
> > Available options are:
>
> ::: {#prop_distrib}
> `'proposal_distribution'`
> :::
>
> > Specifies the statistical distribution used for the proposal
> > density.
>
> `'rand_multivariate_normal'`
>
> > Use a multivariate normal distribution. This is the default.
>
> `'rand_multivariate_student'`
>
> > Use a multivariate student distribution.
>
> `'student_degrees_of_freedom'`
>
> > Specifies the degrees of freedom to be used with the multivariate
> > student distribution. Default: `3`.
>
> ::: {#usemhcov}
> `'use_mh_covariance_matrix'`
> :::
>
> > Indicates to use the covariance matrix of the draws from a previous
> > MCMC run to define the covariance of the proposal distribution.
> > Requires the `load_mh_file`{.interpreted-text role="opt"} option to
> > be specified. Default: `0`.
>
> ::: {#scale-file}
> `'scale_file'`
> :::
>
> > Provides the name of a `_mh_scale.mat` file storing the tuned scale
> > factor from a previous run of `mode_compute=6`.
>
> ::: {#savetmp}
> `'save_tmp_file'`
> :::
>
> > Save the MCMC draws into a `_mh_tmp_blck` file at the refresh rate
> > of the status bar instead of just saving the draws when the current
> > `_mh*_blck` file is full. Default: `0`
>
> `'independent_metropolis_hastings'`
>
> > Takes the same options as in the case of
> > `random_walk_metropolis_hastings`.
>
> `'slice'`
>
> `'rotated'`
>
> > Triggers rotated slice iterations using a covariance matrix from
> > initial burn-in iterations. Requires either
> > `use_mh_covariance_matrix` or `slice_initialize_with_mode`. Default:
> > `0`.
>
> `'mode_files'`
>
> > For multimodal posteriors, provide the name of a file containing a
> > `nparam` by `nmodes` variable called `xparams` storing the different
> > modes. This array must have one column vector per mode and the
> > estimated parameters along the row dimension. With this info, the
> > code will automatically trigger the `rotated` and `mode` options.
> > Default: `[]`.
>
> `'slice_initialize_with_mode'`
>
> > The default for slice is to set `mode_compute=0` and start the
> > chain(s) from a random location in the prior space. This option
> > first runs the mode-finder and then starts the chain from the mode.
> > Together with `rotated`, it will use the inverse Hessian from the
> > mode to perform rotated slice iterations. Default: `0`.
>
> `'initial_step_size'`
>
> > Sets the initial size of the interval in the stepping-out procedure
> > as fraction of the prior support, i.e. the size will be
> > `initial_step_size * (UB-LB)`. `initial_step_size` must be a real
> > number in the interval `[0,1]`. Default: `0.8`.
>
> `'use_mh_covariance_matrix'`
>
> > See `use_mh_covariance_matrix <usemhcov>`{.interpreted-text
> > role="ref"}. Must be used with `'rotated'`. Default: `0`.
>
> `'save_tmp_file'`
>
> > See `save_tmp_file <savetmp>`{.interpreted-text role="ref"}.
> > Default: `1`.
>
> `'tailored_random_block_metropolis_hastings'`
>
> `'proposal_distribution'`
>
> > Specifies the statistical distribution used for the proposal
> > density. See
> > `proposal_distribution <prop_distrib>`{.interpreted-text
> > role="ref"}.
>
> `new_block_probability = DOUBLE`
>
> > Specifies the probability of the next parameter belonging to a new
> > block when the random blocking in the TaRB Metropolis-Hastings
> > algorithm is conducted. The higher this number, the smaller is the
> > average block size and the more random blocks are formed during each
> > parameter sweep. Default: `0.25`.
>
> `mode_compute = INTEGER`
>
> > Specifies the mode-finder run in every iteration for every block of
> > the TaRB Metropolis-Hastings algorithm. See
> > `mode_compute <mode_compute =
> > INTEGER | FUNCTION_NAME>`{.interpreted-text role="opt"}. Default:
> > `4`.
>
> `optim = (NAME, VALUE,...)`
>
> > Specifies the options for the mode-finder used in the TaRB
> > Metropolis-Hastings algorithm. See `optim
> > <optim = (NAME, VALUE, ...)>`{.interpreted-text role="opt"}.
>
> `'scale_file'`
>
> > See `scale_file <scale-file>`{.interpreted-text role="ref"}..
>
> `'save_tmp_file'`
>
> > See `save_tmp_file <savetmp>`{.interpreted-text role="ref"}.
> > Default: `1`.
:::

::: {.option}
moments\_varendo

Triggers the computation of the posterior distribution of the
theoretical moments of the endogenous variables. Results are stored in
`oo_.PosteriorTheoreticalMoments` (see
`oo_.PosteriorTheoreticalMoments`{.interpreted-text role="mvar"}). The
number of lags in the autocorrelation function is controlled by the `ar`
option.
:::

::: {.option}
contemporaneous\_correlation

See `contemporaneous_correlation`{.interpreted-text role="opt"}. Results
are stored in `oo_.PosteriorTheoreticalMoments`. Note that the `nocorr`
option has no effect.
:::

::: {.option}
no\_posterior\_kernel\_density

Shuts off the computation of the kernel density estimator for the
posterior objects (see `density <dens>`{.interpreted-text role="ref"}
field).
:::

::: {.option}
conditional\_variance\_decomposition = INTEGER
conditional\_variance\_decomposition = \[INTEGER1:INTEGER2\]
conditional\_variance\_decomposition = \[INTEGER1 INTEGER2 \...\]

Computes the posterior distribution of the conditional variance
decomposition for the specified period(s). The periods must be strictly
positive. Conditional variances are given by $var(y_{t+k}\vert t)$. For
period 1, the conditional variance decomposition provides the
decomposition of the effects of shocks upon impact. The results are
stored in
`oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition`..
Note that this option requires the option `moments_varendo` to be
specified. In the presence of measurement error, the field will contain
the variance contribution after measurement error has been taken out,
*i.e.* the decomposition will be conducted of the actual as opposed to
the measured variables. The variance decomposition of the measured
variables will be stored in
`oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecompositionME`.
:::

::: {.option}
filtered\_vars

Triggers the computation of the posterior distribution of filtered
endogenous variables/one-step ahead forecasts, i.e. $E_{t}{y_{t+1}}$.
Results are stored in `oo_.FilteredVariables` (see below for a
description of this variable)
:::

::: {.option}
smoother

Triggers the computation of the posterior distribution of smoothed
endogenous variables and shocks, i.e. the expected value of variables
and shocks given the information available in all observations up to the
final date ($E_{T}{y_t}$). Results are stored in
`oo_.SmoothedVariables`, `oo_.SmoothedShocks` and
`oo_.SmoothedMeasurementErrors`. Also triggers the computation of
`oo_.UpdatedVariables`, which contains the estimation of the expected
value of variables given the information available at the current date
($E_{t}{y_t}$). See below for a description of all these variables.
:::

::: {.option}
smoother\_redux

Triggers a faster computation of the smoothed endogenous variables and
shocks for large models. It runs the smoother only for the state
variables (i.e. with the same representation used for likelihood
computations) and computes the remaining variables ex-post. Static
unobserved objects (filtered, smoothed, updated, k-step ahead) are
recovered, but there are exceptions to a full recovery, depending on how
static unobserved variables depend on the restricted state space
adopted. For example, lagged shocks which are ONLY used to recover
NON-observed static variables will not be recovered). For such
exceptions, only the following output is provided:

> `FilteredVariablesKStepAhead`: will be fully recovered
>
> `SmoothedVariables`, `FilteredVariables`, `UpdatedVariables`: recovered for all periods beyond period `d+1`,
>
> :   where `d` denotes the number of diffuse filtering steps.
>
> `FilteredVariablesKStepAheadVariances`, `Variance`, and
> `State_uncertainty` cannot be recovered, and ZERO is provided as
> output.

If you need variances for those variables, either do not set the option,
or declare the variable as observed, using NaNs as data points.
:::

::: {.option}
forecast = INTEGER

Computes the posterior distribution of a forecast on INTEGER periods
after the end of the sample used in estimation. If no
Metropolis-Hastings is computed, the result is stored in variable
`oo_.forecast` and corresponds to the forecast at the posterior mode. If
a Metropolis-Hastings is computed, the distribution of forecasts is
stored in variables `oo_.PointForecast` and `oo_.MeanForecast`. See
`fore`{.interpreted-text role="ref"}, for a description of these
variables.
:::

::: {.option}
tex

See `tex`{.interpreted-text role="opt"}.
:::

::: {.option}
kalman\_algo = INTEGER

`0`

> Automatically use the Multivariate Kalman Filter for stationary models
> and the Multivariate Diffuse Kalman Filter for non-stationary models.

`1`

> Use the Multivariate Kalman Filter.

`2`

> Use the Univariate Kalman Filter.

`3`

> Use the Multivariate Diffuse Kalman Filter.

`4`

> Use the Univariate Diffuse Kalman Filter.
:::

> Default value is `0`. In case of missing observations of single or all
> series, Dynare treats those missing values as unobserved states and
> uses the Kalman filter to infer their value (see e.g. *Durbin and
> Koopman (2012)*, Ch. 4.10) This procedure has the advantage of being
> capable of dealing with observations where the forecast error variance
> matrix becomes singular for some variable(s). If this happens, the
> respective observation enters with a weight of zero in the
> log-likelihood, i.e. this observation for the respective variable(s)
> is dropped from the likelihood computations (for details see *Durbin
> and Koopman (2012)*, Ch. 6.4 and 7.2.5 and *Koopman and Durbin
> (2000)*). If the use of a multivariate Kalman filter is specified and
> a singularity is encountered, Dynare by default automatically switches
> to the univariate Kalman filter for this parameter draw. This behavior
> can be changed via the
> `use_univariate_filters_if_singularity_is_detected
> <use_univariate_filters_if_singularity_is_detected = INTEGER>`{.interpreted-text
> role="opt"} option.

::: {.option}
fast\_kalman\_filter

Select the fast Kalman filter using Chandrasekhar recursions as
described by `Herbst (2015)`. This setting is only used with
`kalman_algo=1` or `kalman_algo=3`. In case of using the diffuse Kalman
filter (`kalman_algo=3/lik_init=3`), the observables must be stationary.
This option is not yet compatible with
`analytic_derivation`{.interpreted-text role="opt"}.
:::

::: {.option}
kalman\_tol = DOUBLE

Numerical tolerance for determining the singularity of the covariance
matrix of the prediction errors during the Kalman filter (minimum
allowed reciprocal of the matrix condition number). Default value is
`1e-10`.
:::

::: {.option}
diffuse\_kalman\_tol = DOUBLE

Numerical tolerance for determining the singularity of the covariance
matrix of the prediction errors ($F_{\infty}$) and the rank of the
covariance matrix of the non-stationary state variables ($P_{\infty}$)
during the Diffuse Kalman filter. Default value is `1e-6`.
:::

::: {.option}
filter\_covariance

Saves the series of one step ahead error of forecast covariance
matrices. With Metropolis, they are saved in
`oo_.FilterCovariance`{.interpreted-text role="mvar"}, otherwise in
`oo_.Smoother.Variance`{.interpreted-text role="mvar"}. Saves also
k-step ahead error of forecast covariance matrices if
`filter_step_ahead` is set.
:::

::: {.option}
filter\_step\_ahead = \[INTEGER1:INTEGER2\] filter\_step\_ahead =
\[INTEGER1 INTEGER2 \...\]

Triggers the computation k-step ahead filtered values, i.e.
$E_{t}{y_{t+k}}$. Stores results in `oo_.FilteredVariablesKStepAhead`.
Also stores 1-step ahead values in `oo_.FilteredVariables`.
`oo_.FilteredVariablesKStepAheadVariances` is stored if
`filter_covariance`.
:::

::: {.option}
filter\_decomposition

Triggers the computation of the shock decomposition of the above k-step
ahead filtered values. Stores results in
`oo_.FilteredVariablesShockDecomposition`.
:::

::: {.option}
smoothed\_state\_uncertainty

Triggers the computation of the variance of smoothed estimates, i.e.
$var_T(y_t)$. Stores results in `oo_.Smoother.State_uncertainty`.
:::

::: {.option}
diffuse\_filter

Uses the diffuse Kalman filter (as described in *Durbin and Koopman
(2012)* and *Koopman and Durbin (2003)* for the multivariate and
*Koopman and Durbin (2000)* for the univariate filter) to estimate
models with non-stationary observed variables. This option will also
reset the `qz_criterium` to count unit root variables towards the stable
variables. Trying to estimate a model with unit roots will otherwise
result in a Blanchard-Kahn error.

When `diffuse_filter` is used the `lik_init` option of `estimation` has
no effect.

When there are nonstationary exogenous variables in a model, there is no
unique deterministic steady state. For instance, if productivity is a
pure random walk:

> $$a_t = a_{t-1} + e_t$$

any value of $\bar a$ of $a$ is a deterministic steady state for
productivity. Consequently, the model admits an infinity of steady
states. In this situation, the user must help Dynare in selecting one
steady state, except if zero is a trivial model's steady state, which
happens when the `linear` option is used in the model declaration. The
user can either provide the steady state to Dynare using a
`steady_state_model` block (or writing a steady state file) if a closed
form solution is available, see `steady_state_model`{.interpreted-text
role="bck"}, or specify some constraints on the steady state, see
`equation_tag_for_conditional_steady_state <eq-tag-ss>`{.interpreted-text
role="ref"}, so that Dynare computes the steady state conditionally on
some predefined levels for the non stationary variables. In both cases,
the idea is to use dummy values for the steady state level of the
exogenous non stationary variables.

Note that the nonstationary variables in the model must be integrated
processes (their first difference or k-difference must be stationary).
:::

::: {.option}
heteroskedastic\_filter

Runs filter, likelihood, and smoother using heteroskedastic definitions
provided in a `heteroskedastic_shocks` block.
:::

::: {.option}
selected\_variables\_only

Only run the classical smoother on the variables listed just after the
`estimation` command. This option is incompatible with requesting
classical frequentist forecasts and will be overridden in this case.
When using Bayesian estimation, the smoother is by default only run on
the declared endogenous variables. Default: run the smoother on all the
declared endogenous variables.
:::

::: {.option}
cova\_compute = INTEGER

When `0`, the covariance matrix of estimated parameters is not computed
after the computation of posterior mode (or maximum likelihood). This
increases speed of computation in large models during development, when
this information is not always necessary. Of course, it will break all
successive computations that would require this covariance matrix.
Otherwise, if this option is equal to `1`, the covariance matrix is
computed and stored in variable `hh` of `MODEL_FILENAME_mode.mat`.
Default is `1`.
:::

::: {.option}
solve\_algo = INTEGER

See `solve_algo <solvalg>`{.interpreted-text role="ref"}.
:::

::: {.option}
order = INTEGER

Order of approximation around the deterministic steady state. When
greater than 1, the likelihood is evaluated with a particle or nonlinear
filter (see *Fernández-Villaverde and Rubio-Ramírez (2005)*). Default is
`1`, i.e. the likelihood of the linearized model is evaluated using a
standard Kalman filter.
:::

::: {.option}
irf = INTEGER

See `irf <irf = INTEGER>`{.interpreted-text role="opt"}. Only used if
`bayesian_irf`{.interpreted-text role="opt"} is passed.
:::

::: {.option}
irf\_shocks = ( VARIABLE\_NAME \[\[,\] VARIABLE\_NAME \...\] )

See `irf_shocks <irf_shocks = ( VARIABLE_NAME [[,]
VARIABLE_NAME ...] )>`{.interpreted-text role="opt"}. Only used if
`bayesian_irf`{.interpreted-text role="opt"} is passed.
:::

::: {.option}
irf\_plot\_threshold = DOUBLE

See `irf_plot_threshold <irf_plot_threshold =
DOUBLE>`{.interpreted-text role="opt"}. Only used if
`bayesian_irf`{.interpreted-text role="opt"} is passed.
:::

::: {.option}
aim\_solver

See `aim_solver`{.interpreted-text role="opt"}.
:::

::: {.option}
sylvester = OPTION

See `sylvester <sylvester = OPTION>`{.interpreted-text role="opt"}.
:::

::: {.option}
sylvester\_fixed\_point\_tol = DOUBLE

See `sylvester_fixed_point_tol <sylvester_fixed_point_tol
= DOUBLE>`{.interpreted-text role="opt"} .
:::

::: {.option}
lyapunov = OPTION

Determines the algorithm used to solve the Lyapunov equation to
initialized the variance-covariance matrix of the Kalman filter using
the steady-state value of state variables. Possible values for OPTION
are:

> `default`
>
> > Uses the default solver for Lyapunov equations based on
> > Bartels-Stewart algorithm.
>
> `fixed_point`
>
> > Uses a fixed point algorithm to solve the Lyapunov equation. This
> > method is faster than the `default` one for large scale models, but
> > it could require a large amount of iterations.
>
> `doubling`
>
> > Uses a doubling algorithm to solve the Lyapunov equation
> > (`disclyap_fast`). This method is faster than the two previous one
> > for large scale models.
>
> `square_root_solver`
>
> > Uses a square-root solver for Lyapunov equations (`dlyapchol`). This
> > method is fast for large scale models (available under MATLAB if the
> > Control System Toolbox is installed; available under Octave if the
> > [control](https://octave.sourceforge.io/control/) package from
> > Octave-Forge is installed)

Default value is `default`.
:::

::: {.option}
lyapunov\_fixed\_point\_tol = DOUBLE

This is the convergence criterion used in the fixed point Lyapunov
solver. Its default value is `1e-10`.
:::

::: {.option}
lyapunov\_doubling\_tol = DOUBLE

This is the convergence criterion used in the doubling algorithm to
solve the Lyapunov equation. Its default value is `1e-16`.
:::

::: {.option}
use\_penalized\_objective\_for\_hessian

Use the penalized objective instead of the objective function to compute
numerically the hessian matrix at the mode. The penalties decrease the
value of the posterior density (or likelihood) when, for some
perturbations, Dynare is not able to solve the model (issues with steady
state existence, Blanchard and Kahn conditions, \...). In pratice, the
penalized and original objectives will only differ if the posterior mode
is found to be near a region where the model is ill-behaved. By default
the original objective function is used.
:::

::: {.option}
analytic\_derivation

Triggers estimation with analytic gradient at `order=1`. The final
hessian at the mode is also computed analytically. Only works for
stationary models without missing observations, i.e. for
`kalman_algo<3`. Optimizers that rely on analytic gradients are
`mode_compute=1,3,4,5,101`.
:::

::: {.option}
ar = INTEGER

See `ar <ar = INTEGER>`{.interpreted-text role="opt"}. Only useful in
conjunction with option `moments_varendo`.
:::

::: {.option}
endogenous\_prior

Use endogenous priors as in *Christiano, Trabandt and Walentin (2011)*.
The procedure is motivated by sequential Bayesian learning. Starting
from independent initial priors on the parameters, specified in the
`estimated_params` block, the standard deviations observed in a
\"pre-sample\", taken to be the actual sample, are used to update the
initial priors. Thus, the product of the initial priors and the
pre-sample likelihood of the standard deviations of the observables is
used as the new prior (for more information, see the technical appendix
of *Christiano, Trabandt and Walentin (2011)*). This procedure helps in
cases where the regular posterior estimates, which minimize in-sample
forecast errors, result in a large overprediction of model variable
variances (a statistic that is not explicitly targeted, but often of
particular interest to researchers).
:::

::: {.option}
use\_univariate\_filters\_if\_singularity\_is\_detected = INTEGER

Decide whether Dynare should automatically switch to univariate filter
if a singularity is encountered in the likelihood computation (this is
the behaviour if the option is equal to `1`). Alternatively, if the
option is equal to `0`, Dynare will not automatically change the filter,
but rather use a penalty value for the likelihood when such a
singularity is encountered. Default: `1`.
:::

::: {.option}
keep\_kalman\_algo\_if\_singularity\_is\_detected

With the default `use_univariate_filters_if_singularity_is_detected=1
<use_univariate_filters_if_singularity_is_detected = INTEGER>`{.interpreted-text
role="opt"}, Dynare will switch to the univariate Kalman filter when it
encounters a singular forecast error variance matrix during Kalman
filtering. Upon encountering such a singularity for the first time, all
subsequent parameter draws and computations will automatically rely on
univariate filter, i.e. Dynare will never try the multivariate filter
again. Use the `keep_kalman_algo_if_singularity_is_detected` option to
have the `use_univariate_filters_if_singularity_is_detected` only affect
the behavior for the current draw/computation.
:::

::: {.option}
rescale\_prediction\_error\_covariance

Rescales the prediction error covariance in the Kalman filter to avoid
badly scaled matrix and reduce the probability of a switch to univariate
Kalman filters (which are slower). By default no rescaling is done.
:::

::: {.option}
qz\_zero\_threshold = DOUBLE

See `qz_zero_threshold <qz_zero_threshold = DOUBLE>`{.interpreted-text
role="opt"}.
:::

::: {.option}
taper\_steps = \[INTEGER1 INTEGER2 \...\]

Percent tapering used for the spectral window in the *Geweke
(1992,1999)* convergence diagnostics (requires
`mh_nblocks=1 <mh_nblocks = INTEGER>`{.interpreted-text role="opt"}).
The tapering is used to take the serial correlation of the posterior
draws into account. Default: `[4 8 15]`.
:::

::: {.option}
geweke\_interval = \[DOUBLE DOUBLE\]

Percentage of MCMC draws at the beginning and end of the MCMC chain
taken to compute the *Geweke (1992,1999)* convergence diagnostics
(requires `mh_nblocks=1 <mh_nblocks =
INTEGER>`{.interpreted-text role="opt"}) after discarding the first
`mh_drop = DOUBLE
<mh_drop>`{.interpreted-text role="opt"} percent of draws as a burnin.
Default: \[0.2 0.5\].
:::

::: {.option}
raftery\_lewis\_diagnostics

Triggers the computation of the *Raftery and Lewis (1992)* convergence
diagnostics. The goal is deliver the number of draws required to
estimate a particular quantile of the CDF `q` with precision `r` with a
probability `s`. Typically, one wants to estimate the `q=0.025`
percentile (corresponding to a 95 percent HPDI) with a precision of 0.5
percent (`r=0.005`) with 95 percent certainty (`s=0.95`). The defaults
can be changed via `raftery_lewis_qrs
<raftery_lewis_qrs = [DOUBLE DOUBLE DOUBLE]>`{.interpreted-text
role="opt"}. Based on the theory of first order Markov Chains, the
diagnostics will provide a required burn-in (`M`), the number of draws
after the burnin (`N`) as well as a thinning factor that would deliver a
first order chain (`k`). The last line of the table will also deliver
the maximum over all parameters for the respective values.
:::

::: {.option}
raftery\_lewis\_qrs = \[DOUBLE DOUBLE DOUBLE\]

Sets the quantile of the CDF `q` that is estimated with precision `r`
with a probability `s` in the *Raftery and Lewis (1992)* convergence
diagnostics. Default: `[0.025 0.005 0.95]`.
:::

::: {.option}
consider\_all\_endogenous

Compute the posterior moments, smoothed variables, k-step ahead filtered
variables and forecasts (when requested) on all the endogenous
variables. This is equivalent to manually listing all the endogenous
variables after the `estimation` command.
:::

::: {.option}
consider\_all\_endogenous\_and\_auxiliary

Compute the posterior moments, smoothed variables, k-step ahead filtered
variables and forecasts (when requested) on all the endogenous variables
and the auxiliary variables introduced by the preprocessor. This option
is useful when e.g. running `smoother2histval` on the results of the
Kalman smoother.
:::

::: {.option}
consider\_only\_observed

Compute the posterior moments, smoothed variables, k-step ahead filtered
variables and forecasts (when requested) on all the observed variables.
This is equivalent to manually listing all the observed variables after
the `estimation` command.
:::

::: {.option}
number\_of\_particles = INTEGER

Number of particles used when evaluating the likelihood of a non linear
state space model. Default: `1000`.
:::

::: {.option}
resampling = OPTION

Determines if resampling of the particles is done. Possible values for
OPTION are:

> `none`
>
> > No resampling.
>
> `systematic`
>
> > Resampling at each iteration, this is the default value.
>
> `generic`
>
> > Resampling if and only if the effective sample size is below a
> > certain level defined by
> > `resampling_threshold <resampling_threshold =
> > DOUBLE>`{.interpreted-text role="opt"} \* `number_of_particles
> > <number_of_particles = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
resampling\_threshold = DOUBLE

A real number between zero and one. The resampling step is triggered as
soon as the effective number of particles is less than this number times
the total number of particles (as set by
`number_of_particles <number_of_particles =
INTEGER>`{.interpreted-text role="opt"}). This option is effective if
and only if option `resampling <resampling = OPTION>`{.interpreted-text
role="opt"} has value `generic`.
:::

::: {.option}
resampling\_method = OPTION

Sets the resampling method. Possible values for OPTION are: `kitagawa`,
`stratified` and `smooth`.
:::

::: {.option}
filter\_algorithm = OPTION

Sets the particle filter algorithm. Possible values for OPTION are:

> `sis`
>
> > Sequential importance sampling algorithm, this is the default value.
>
> `apf`
>
> > Auxiliary particle filter.
>
> `gf`
>
> > Gaussian filter.
>
> `gmf`
>
> > Gaussian mixture filter.
>
> `cpf`
>
> > Conditional particle filter.
>
> `nlkf`
>
> > Use a standard (linear) Kalman filter algorithm with the nonlinear
> > measurement and state equations.
:::

::: {.option}
proposal\_approximation = OPTION

Sets the method for approximating the proposal distribution. Possible
values for OPTION are: `cubature`, `montecarlo` and `unscented`. Default
value is `unscented`.
:::

::: {.option}
distribution\_approximation = OPTION

Sets the method for approximating the particle distribution. Possible
values for OPTION are: `cubature`, `montecarlo` and `unscented`. Default
value is `unscented`.
:::

::: {.option}
cpf\_weights = OPTION

Controls the method used to update the weights in conditional particle
filter, possible values are `amisanotristani` (*Amisano et al. (2010)*)
or `murrayjonesparslow` (*Murray et al. (2013)*). Default value is
`amisanotristani`.
:::

::: {.option}
nonlinear\_filter\_initialization = INTEGER

Sets the initial condition of the nonlinear filters. By default the
nonlinear filters are initialized with the unconditional covariance
matrix of the state variables, computed with the reduced form solution
of the first order approximation of the model. If
`nonlinear_filter_initialization=2`, the nonlinear filter is instead
initialized with a covariance matrix estimated with a stochastic
simulation of the reduced form solution of the second order
approximation of the model. Both these initializations assume that the
model is stationary, and cannot be used if the model has unit roots
(which can be seen with the `check`{.interpreted-text role="comm"}
command prior to estimation). If the model has stochastic trends, user
must use `nonlinear_filter_initialization=3`, the filters are then
initialized with an identity matrix for the covariance matrix of the
state variables. Default value is `nonlinear_filter_initialization=1`
(initialization based on the first order approximation of the model).
:::

::: {.option}
particle\_filter\_options = (NAME, VALUE, \...)

A list of NAME and VALUE pairs. Can be used to set some fine-grained
options for the particle filter routines. The set of available options
depends on the selected filter routine.

More information on particle filter options is available at
<https://git.dynare.org/Dynare/dynare/-/wikis/Particle-filters>.

Available options are:

> `'pruning'`
>
> > Enable pruning for particle filter-related simulations. Default:
> > `false`.
>
> `'liu_west_delta'`
>
> > Set the value for delta for the Liu/West online filter. Default:
> > `0.99`.
>
> `'unscented_alpha'`
>
> > Set the value for alpha for unscented transforms. Default: `1`.
>
> `'unscented_beta'`
>
> > Set the value for beta for unscented transforms. Default: `2`.
>
> `'unscented_kappa'`
>
> > Set the value for kappa for unscented transforms. Default: `1`.
>
> `'initial_state_prior_std'`
>
> > Value of the diagonal elements for the initial covariance of the
> > state variables when employing `nonlinear_filter_initialization=3`.
> > Default: `1`.
>
> `'mixture_state_variables'`
>
> > Number of mixture components in the Gaussian-mixture filter (gmf)
> > for the state variables. Default: `5`.
>
> `'mixture_structural_shocks'`
>
> > Number of mixture components in the Gaussian-mixture filter (gmf)
> > for the structural shocks. Default: `1`.
>
> `'mixture_measurement_shocks'`
>
> > Number of mixture components in the Gaussian-mixture filter (gmf)
> > for the measurement errors. Default: `1`.
:::

*Note*

If no `mh_jscale` parameter is used for a parameter in
`estimated_params`, the procedure uses `mh_jscale` for all parameters.
If `mh_jscale` option isn't set, the procedure uses `0.2` for all
parameters. Note that if `mode_compute=6` is used or the
`posterior_sampler_option` called `scale_file` is specified, the values
set in `estimated_params` will be overwritten.

*"Endogenous" prior restrictions*

It is also possible to impose implicit "endogenous" priors about IRFs
and moments on the model during estimation. For example, one can specify
that all valid parameter draws for the model must generate fiscal
multipliers that are bigger than 1 by specifying how the IRF to a
government spending shock must look like. The prior restrictions can be
imposed via `irf_calibration` and `moment_calibration` blocks (see
`irf-momcal`{.interpreted-text role="ref"}). The way it works internally
is that any parameter draw that is inconsistent with the "calibration"
provided in these blocks is discarded, i.e. assigned a prior density of
0. When specifying these blocks, it is important to keep in mind that
one won't be able to easily do `model_comparison` in this case, because
the prior density will not integrate to 1.

*Output*

After running estimation, the parameters `M_.params` and the variance
matrix `M_.Sigma_e` of the shocks are set to the mode for maximum
likelihood estimation or posterior mode computation without Metropolis
iterations. After estimation with Metropolis iterations (option
`mh_replic > 0` or option `load_mh_file` set) the parameters `M_.params`
and the variance matrix `M_.Sigma_e` of the shocks are set to the
posterior mean.

Depending on the options, `estimation` stores results in various fields
of the `oo_` structure, described below. In the following variables, we
will adopt the following shortcuts for specific field names:

> `MOMENT_NAME`
>
> > This field can take the following values:
> >
> > `HPDinf`
> >
> > > Lower bound of a 90% HPD interval.[^4]
> >
> > `HPDsup`
> >
> > > Upper bound of a 90% HPD interval.
> >
> > `HPDinf_ME`
> >
> > > Lower bound of a 90% HPD interval[^5] for observables when taking
> > > measurement error into account (see e.g. *Christoffel et al.
> > > (2010*), p.17).
> >
> > `HPDsup_ME`
> >
> > > Upper bound of a 90% HPD interval for observables when taking
> > > measurement error into account.
> >
> > `Mean`
> >
> > > Mean of the posterior distribution.
> >
> > `Median`
> >
> > > Median of the posterior distribution.
> >
> > `Std`
> >
> > > Standard deviation of the posterior distribution.
> >
> > `Variance`
> >
> > > Variance of the posterior distribution.
> >
> > `deciles`
> >
> > > Deciles of the distribution.
> >
> > ::: {#dens}
> > `density`
> > :::
> >
> > > Non parametric estimate of the posterior density following the
> > > approach outlined in *Skoeld and Roberts (2003)*. First and second
> > > columns are respectively abscissa and ordinate coordinates.
>
> `ESTIMATED_OBJECT`
>
> > This field can take the following values:
> >
> > `measurement_errors_corr`
> >
> > > Correlation between two measurement errors.
> >
> > `measurement_errors_std`
> >
> > > Standard deviation of measurement errors.
> >
> > `parameters`
> >
> > > Parameters.
> >
> > `shocks_corr`
> >
> > > Correlation between two structural shocks.
> >
> > `shocks_std`
> >
> > > Standard deviation of structural shocks.

::: {.matvar}
[oo]().MarginalDensity.LaplaceApproximation

Variable set by the `estimation` command. Stores the marginal data
density based on the Laplace Approximation.
:::

::: {.matvar}
[oo]().MarginalDensity.ModifiedHarmonicMean

Variable set by the `estimation command`, if it is used with
`mh_replic > 0` or `load_mh_file` option. Stores the marginal data
density based on *Geweke (1999)* Modified Harmonic Mean estimator.
:::

::: {.matvar}
[oo]().posterior.optimization

Variable set by the `estimation` command if mode-finding is used. Stores
the results at the mode. Fields are of the form:

    oo_.posterior.optimization.OBJECT

where OBJECT is one of the following:

> `mode`
>
> > Parameter vector at the mode.
>
> `Variance`
>
> > Inverse Hessian matrix at the mode or MCMC jumping covariance matrix
> > when used with the `MCMC_jumping_covariance <mcmc_jumping_covariance
> > = OPTION>`{.interpreted-text role="opt"} option.
>
> `log_density`
>
> > Log likelihood (ML)/log posterior density (Bayesian) at the mode
> > when used with `mode_compute>0`.
:::

::: {.matvar}
[oo]().posterior.metropolis

Variable set by the `estimation` command if `mh_replic>0` is used.
Fields are of the form:

    oo_.posterior.metropolis.OBJECT

where OBJECT is one of the following:

> `mean`
>
> > Mean parameter vector from the MCMC.
>
> `Variance`
>
> > Covariance matrix of the parameter draws in the MCMC.
:::

::: {.matvar}
[oo]().FilteredVariables

Variable set by the `estimation` command, if it is used with the
`filtered_vars` option.

After an estimation without Metropolis, fields are of the form:

    oo_.FilteredVariables.VARIABLE_NAME

After an estimation with Metropolis, fields are of the form:

    oo_.FilteredVariables.MOMENT_NAME.VARIABLE_NAME
:::

::: {.matvar}
[oo]().FilteredVariablesKStepAhead

Variable set by the `estimation` command, if it is used with the
`filter_step_ahead` option. The k-steps are stored along the rows while
the columns indicate the respective variables. The third dimension of
the array provides the observation for which the forecast has been made.
For example, if `filter_step_ahead=[1 2 4]` and `nobs=200`, the element
(3,5,204) stores the four period ahead filtered value of variable 5
computed at time t=200 for time t=204. The periods at the beginning and
end of the sample for which no forecasts can be made, e.g. entries
(1,5,1) and (1,5,204) in the example, are set to zero. Note that in case
of Bayesian estimation the variables will be ordered in the order of
declaration after the estimation command (or in general declaration
order if no variables are specified here). In case of running the
classical smoother, the variables will always be ordered in general
declaration order. If the `selected_variables_only`{.interpreted-text
role="opt"} option is specified with the classical smoother,
non-requested variables will be simply left out in this order.
:::

::: {.matvar}
[oo]().FilteredVariablesKStepAheadVariances

Variable set by the `estimation` command, if it is used with the
`filter_step_ahead option`. It is a 4 dimensional array where the
k-steps are stored along the first dimension, while the fourth dimension
of the array provides the observation for which the forecast has been
made. The second and third dimension provide the respective variables.
For example, if `filter_step_ahead=[1 2 4]` and `nobs=200`, the element
(3,4,5,204) stores the four period ahead forecast error covariance
between variable 4 and variable 5, computed at time t=200 for time
t=204. Padding with zeros and variable ordering is analogous to
`oo_.FilteredVariablesKStepAhead`.
:::

::: {.matvar}
[oo]().Filtered\_Variables\_X\_step\_ahead

Variable set by the `estimation` command, if it is used with the
`filter_step_ahead option` in the context of Bayesian estimation. Fields
are of the form:

    oo_.Filtered_Variables_X_step_ahead.VARIABLE_NAME

The n-th entry stores the k-step ahead filtered variable computed at
time n for time n+k.
:::

::: {.matvar}
[oo]().FilteredVariablesShockDecomposition

Variable set by the `estimation` command, if it is used with the
`filter_step_ahead` option. The k-steps are stored along the rows while
the columns indicate the respective variables. The third dimension
corresponds to the shocks in declaration order. The fourth dimension of
the array provides the observation for which the forecast has been made.
For example, if `filter_step_ahead=[1 2 4]` and `nobs=200`, the element
(3,5,2,204) stores the contribution of the second shock to the four
period ahead filtered value of variable 5 (in deviations from the mean)
computed at time t=200 for time t=204. The periods at the beginning and
end of the sample for which no forecasts can be made, e.g. entries
(1,5,1) and (1,5,204) in the example, are set to zero. Padding with
zeros and variable ordering is analogous to
`oo_.FilteredVariablesKStepAhead`.
:::

::: {.matvar}
[oo]().PosteriorIRF.dsge

Variable set by the `estimation` command, if it is used with the
`bayesian_irf` option. Fields are of the form:

    oo_.PosteriorIRF.dsge.MOMENT_NAME.VARIABLE_NAME_SHOCK_NAME
:::

::: {.matvar}
[oo]().SmoothedMeasurementErrors

Variable set by the `estimation` command, if it is used with the
`smoother` option. Fields are of the form:

    oo_.SmoothedMeasurementErrors.VARIABLE_NAME
:::

::: {.matvar}
[oo]().SmoothedShocks

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command.

After an estimation without Metropolis, or if computed by
`calib_smoother`, fields are of the form:

    oo_.SmoothedShocks.VARIABLE_NAME

After an estimation with Metropolis, fields are of the form:

    oo_.SmoothedShocks.MOMENT_NAME.VARIABLE_NAME
:::

::: {.matvar}
[oo]().SmoothedVariables

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command.

After an estimation without Metropolis, or if computed by
`calib_smoother`, fields are of the form:

    oo_.SmoothedVariables.VARIABLE_NAME

After an estimation with Metropolis, fields are of the form:

    oo_.SmoothedVariables.MOMENT_NAME.VARIABLE_NAME
:::

::: {.matcomm}
get\_smooth (\'VARIABLE\_NAME\' \[, \'VARIABLE\_NAME\'\]\...);

Returns the smoothed values of the given endogenous or exogenous
variable(s), as they are stored in the `oo_.SmoothedVariables` and
`oo_.SmoothedShocks` variables.
:::

::: {.matvar}
[oo]().UpdatedVariables

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command. Contains the estimation of
the expected value of variables given the information available at the
current date.

After an estimation without Metropolis, or if computed by
`calib_smoother`, fields are of the form:

    oo_.UpdatedVariables.VARIABLE_NAME

After an estimation with Metropolis, fields are of the form:

    oo_.UpdatedVariables.MOMENT_NAME.VARIABLE_NAME
:::

::: {.matcomm}
get\_update (\'VARIABLE\_NAME\' \[, \'VARIABLE\_NAME\'\]\...);

Returns the updated values of the given variable(s), as they are stored
in the `oo_.UpdatedVariables` variable.
:::

::: {.matvar}
[oo]().FilterCovariance

Three-dimensional array set by the `estimation` command if used with the
`smoother` and Metropolis, if the `filter_covariance` option has been
requested. Contains the series of one-step ahead forecast error
covariance matrices from the Kalman smoother. The `M_.endo_nbr` times
`M_.endo_nbr` times `T+1` array contains the variables in declaration
order along the first two dimensions. The third dimension of the array
provides the observation for which the forecast has been made. Fields
are of the form:

    oo_.FilterCovariance.MOMENT_NAME

Note that density estimation is not supported.
:::

::: {.matvar}
[oo]().Smoother.Variance

Three-dimensional array set by the `estimation` command (if used with
the `smoother`) without Metropolis, or by the `calib_smoother` command,
if the `filter_covariance` option has been requested. Contains the
series of one-step ahead forecast error covariance matrices from the
Kalman smoother. The `M_.endo_nbr` times `M_.endo_nbr` times `T+1` array
contains the variables in declaration order along the first two
dimensions. The third dimension of the array provides the observation
for which the forecast has been made.
:::

::: {.matvar}
[oo]().Smoother.State\_uncertainty

Three-dimensional array set by the `estimation` command (if used with
the `smoother` option) without Metropolis, or by the `calib_smoother`
command, if the `smoothed_state_uncertainty` option has been requested.
Contains the series of covariance matrices for the state estimate given
the full data from the Kalman smoother. The `M_.endo_nbr` times
`M_.endo_nbr` times `T` array contains the variables in declaration
order along the first two dimensions. The third dimension of the array
provides the observation for which the smoothed estimate has been made.
:::

::: {.matvar}
[oo]().Smoother.SteadyState

Variable set by the `estimation` command (if used with the `smoother`)
without Metropolis, or by the `calib_smoother` command. Contains the
steady state component of the endogenous variables used in the smoother
in order of variable declaration.
:::

::: {.matvar}
[oo]().Smoother.TrendCoeffs

Variable set by the `estimation` command (if used with the `smoother`)
without Metropolis, or by the `calib_smoother` command. Contains the
trend coefficients of the observed variables used in the smoother in
order of declaration of the observed variables.
:::

::: {.matvar}
[oo]().Smoother.Trend

Variable set by the `estimation command` (if used with the `smoother`
option), or by the `calib_smoother` command. Contains the trend
component of the variables used in the smoother.

Fields are of the form:

    oo_.Smoother.Trend.VARIABLE_NAME
:::

::: {.matvar}
[oo]().Smoother.Constant

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command. Contains the constant part
of the endogenous variables used in the smoother, accounting e.g. for
the data mean when using the prefilter option.

Fields are of the form:

    oo_.Smoother.Constant.VARIABLE_NAME
:::

::: {.matvar}
[oo]().Smoother.loglinear

Indicator keeping track of whether the smoother was run with the
`loglinear <logl>`{.interpreted-text role="ref"} option and thus whether
stored smoothed objects are in logs.
:::

::: {.matvar}
[oo]().PosteriorTheoreticalMoments

Variable set by the `estimation` command, if it is used with the
`moments_varendo` option. Fields are of the form:

    oo_.PosteriorTheoreticalMoments.dsge.THEORETICAL_MOMENT.ESTIMATED_OBJECT.MOMENT_NAME.VARIABLE_NAME

where *THEORETICAL\_MOMENT* is one of the following:

> `covariance`
>
> > Variance-covariance of endogenous variables.
>
> `contemporaneous_correlation`
>
> > Contemporaneous correlation of endogenous variables when the
> > `contemporaneous_correlation`{.interpreted-text role="opt"} option
> > is specified.
>
> `correlation`
>
> > Auto- and cross-correlation of endogenous variables. Fields are
> > vectors with correlations from 1 up to order `options_.ar`.
>
> ::: {#VarianceDecomposition}
> `VarianceDecomposition`
> :::
>
> > Decomposition of variance (unconditional variance, i.e. at horizon
> > infinity).[^6]
>
> `VarianceDecompositionME`
>
> > Same as [VarianceDecomposition](), but contains the decomposition of
> > the measured as opposed to the actual variable. The joint
> > contribution of the measurement error will be saved in a field named
> > `ME`.
>
> ::: {#ConditionalVarianceDecomposition}
> `ConditionalVarianceDecomposition`
> :::
>
> > Only if the `conditional_variance_decomposition` option has been
> > specified. In the presence of measurement error, the field will
> > contain the variance contribution after measurement error has been
> > taken out, i.e. the decomposition will be conducted of the actual as
> > opposed to the measured variables.
>
> `ConditionalVarianceDecompositionME`
>
> > Only if the `conditional_variance_decomposition` option has been
> > specified. Same as [ConditionalVarianceDecomposition](), but
> > contains the decomposition of the measured as opposed to the actual
> > variable. The joint contribution of the measurement error will be
> > saved in a field names `ME`.
:::

::: {.matvar}
[oo]().posterior\_density

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_density.PARAMETER_NAME
:::

::: {.matvar}
[oo]().posterior\_hpdinf

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_hpdinf.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_hpdsup

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_hpdsup.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_mean

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_mean.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_mode

Variable set by the `estimation` command during mode-finding. Fields are
of the form:

    oo_.posterior_mode.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_std\_at\_mode

Variable set by the `estimation` command during mode-finding. It is
based on the inverse Hessian at `oo_.posterior_mode`. Fields are of the
form:

    oo_.posterior_std_at_mode.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_std

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_std.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_var

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_var.ESTIMATED_OBJECT.VARIABLE_NAME
:::

::: {.matvar}
[oo]().posterior\_median

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

    oo_.posterior_median.ESTIMATED_OBJECT.VARIABLE_NAME
:::

*Example*

Here are some examples of generated variables:

    oo_.posterior_mode.parameters.alp
    oo_.posterior_mean.shocks_std.ex
    oo_.posterior_hpdsup.measurement_errors_corr.gdp_conso

::: {.matvar}
[oo]().dsge\_var.posterior\_mode

Structure set by the `dsge_var` option of the `estimation` command after
mode\_compute.

The following fields are saved:

> `PHI_tilde`
>
> > Stacked posterior DSGE-BVAR autoregressive matrices at the mode
> > (equation (28) of *Del Negro and Schorfheide (2004)*).
>
> `SIGMA_u_tilde`
>
> > Posterior covariance matrix of the DSGE-BVAR at the mode (equation
> > (29) of *Del Negro and Schorfheide (2004)*).
>
> `iXX`
>
> > Posterior population moments in the DSGE-BVAR at the mode (
> > $inv(\lambda T \Gamma_{XX}^*+ X'X)$).
>
> `prior`
>
> > Structure storing the DSGE-BVAR prior.
>
> `PHI_star`
>
> > Stacked prior DSGE-BVAR autoregressive matrices at the mode
> > (equation (22) of *Del Negro and Schorfheide (2004)*).
>
> `SIGMA_star`
>
> > Prior covariance matrix of the DSGE-BVAR at the mode (equation (23)
> > of *Del Negro and Schorfheide (2004)*).
>
> `ArtificialSampleSize`
>
> > Size of the artifical prior sample ( $inv(\lambda T)$).
>
> `DF`
>
> > Prior degrees of freedom ( $inv(\lambda T-k-n)$).
>
> `iGXX_star`
>
> > Inverse of the theoretical prior "covariance" between X and X
> > ($\Gamma_{xx}^*$ in *Del Negro and Schorfheide (2004)*).
:::

::: {.matvar}
[oo]().RecursiveForecast

Variable set by the `forecast` option of the `estimation` command when
used with the nobs = \[INTEGER1:INTEGER2\] option (see
`nobs <nobs = [INTEGER1:INTEGER2]>`{.interpreted-text role="opt"}).

Fields are of the form:

    oo_.RecursiveForecast.FORECAST_OBJECT.VARIABLE_NAME

where `FORECAST_OBJECT` is one of the following[^7] :

`Mean`

> Mean of the posterior forecast distribution.

`HPDinf/HPDsup`

> Upper/lower bound of the 90% HPD interval taking into account only
> parameter uncertainty (corresponding to
> `oo_.MeanForecast`{.interpreted-text role="mvar"}).

`HPDTotalinf/HPDTotalsup`.

> Upper/lower bound of the 90% HPD interval taking into account both
> parameter and future shock uncertainty (corresponding to
> `oo_.PointForecast`{.interpreted-text role="mvar"})

`VARIABLE_NAME` contains a matrix of the following size: number of time
periods for which forecasts are requested using the
`nobs = [INTEGER1:INTEGER2]` option times the number of forecast
horizons requested by the forecast option. i.e., the row indicates the
period at which the forecast is performed and the column the respective
k-step ahead forecast. The starting periods are sorted in ascending
order, not in declaration order.
:::

::: {.matvar}
[oo]().convergence.geweke

Variable set by the convergence diagnostics of the `estimation` command
when used with `mh_nblocks=1` option (see
`mh_nblocks <mh_nblocks = INTEGER>`{.interpreted-text role="opt"}).

Fields are of the form:

    oo_.convergence.geweke.VARIABLE_NAME.DIAGNOSTIC_OBJECT

where *DIAGNOSTIC\_OBJECT* is one of the following:

`posteriormean`

> Mean of the posterior parameter distribution.

`posteriorstd`

> Standard deviation of the posterior parameter distribution.

`nse_iid`

> Numerical standard error (NSE) under the assumption of iid draws.

`rne_iid`

> Relative numerical efficiency (RNE) under the assumption of iid draws.

`nse_x`

> Numerical standard error (NSE) when using an x% taper.

`rne_x`

> Relative numerical efficiency (RNE) when using an x% taper.

`pooled_mean`

> Mean of the parameter when pooling the beginning and end parts of the
> chain specified in `geweke_interval
> <geweke_interval = [DOUBLE DOUBLE]>`{.interpreted-text role="opt"} and
> weighting them with their relative precision. It is a vector
> containing the results under the iid assumption followed by the ones
> using the `taper_steps` option (see `taper_steps <taper_steps
> = [INTEGER1 INTEGER2 ...]>`{.interpreted-text role="opt"}).

`pooled_nse`

> NSE of the parameter when pooling the beginning and end parts of the
> chain and weighting them with their relative precision. See
> `pooled_mean`.

`prob_chi2_test`

> p-value of a chi-squared test for equality of means in the beginning
> and the end of the MCMC chain. See `pooled_mean`. A value above 0.05
> indicates that the null hypothesis of equal means and thus convergence
> cannot be rejected at the 5 percent level. Differing values along the
> `taper_steps` signal the presence of significant autocorrelation in
> draws. In this case, the estimates using a higher tapering are usually
> more reliable.
:::
:::
:::

::: {.command}
unit\_root\_vars VARIABLE\_NAME\...;

This command is deprecated. Use `estimation` option `diffuse_filter`
instead for estimating a model with non-stationary observed variables or
`steady` option `nocheck` to prevent `steady` to check the steady state
returned by your steady state file.
:::

Dynare also has the ability to estimate Bayesian VARs:

::: {.command}
bvar\_density ;

Computes the marginal density of an estimated BVAR model, using
Minnesota priors.

See `bvar-a-la-sims.pdf`, which comes with Dynare distribution, for more
information on this command.
:::

Estimation based on moments
---------------------------

Provided that you have observations on some endogenous variables, it is
possible to use Dynare to estimate some or all parameters using a method
of moments approach. Both the Simulated Method of Moments (SMM) and the
Generalized Method of Moments (GMM) are available. The general idea is
to minimize the distance between unconditional model moments and
corresponding data moments (so called orthogonality or moment
conditions). For SMM, Dynare computes model moments via stochastic
simulations based on the perturbation approximation up to any order,
whereas for GMM model moments are computed in closed-form based on the
pruned state-space representation of the perturbation solution up to
third order. The implementation of SMM is inspired by *Born and Pfeifer
(2014)* and *Ruge-Murcia (2012)*, whereas the one for GMM is adapted
from *Andreasen, Fernández-Villaverde and Rubio-Ramírez (2018)* and
*Mutschler (2018)*. Successful estimation heavily relies on the accuracy
and efficiency of the perturbation approximation, so it is advised to
tune this as much as possible (see `stoch-sol-simul`{.interpreted-text
role="ref"}). The method of moments estimator is consistent and
asymptotically normally distributed given certain regularity conditions
(see *Duffie and Singleton (1993)* for SMM and *Hansen (1982)* for GMM).
For instance, it is required to have at least as many moment conditions
as estimated parameters (over-identified or just identified). Moreover,
the Jacobian of the moments with respect to the estimated parameters
needs to have full rank. `identification-analysis`{.interpreted-text
role="ref"} helps to check this regularity condition.

In the over-identified case of declaring more moment conditions than
estimated parameters, the choice of
`weighting_matrix <weighting_matrix = ['WM1','WM2',...,'WMn']>`{.interpreted-text
role="opt"} matters for the efficiency of the estimation, because the
estimated orthogonality conditions are random variables with unequal
variances and usually non-zero cross-moment covariances. A weighting
matrix allows to re-weight moments to put more emphasis on moment
conditions that are more informative or better measured (in the sense of
having a smaller variance). To achieve asymptotic efficiency, the
weighting matrix needs to be chosen such that, after appropriate
scaling, it has a probability limit proportional to the inverse of the
covariance matrix of the limiting distribution of the vector of
orthogonality conditions. Dynare uses a Newey-West-type estimator with a
Bartlett kernel to compute an estimate of this so-called optimal
weighting matrix. Note that in this over-identified case, it is advised
to perform the estimation in at least two stages by setting e.g.
`weighting_matrix=['DIAGONAL','DIAGONAL'] <weighting_matrix = ['WM1','WM2',...,'WMn']>`{.interpreted-text
role="opt"} so that the computation of the optimal weighting matrix
benefits from the consistent estimation of the previous stages. The
optimal weighting matrix is used to compute standard errors and the
J-test of overidentifying restrictions, which tests whether the model
and selection of moment conditions fits the data sufficiently well. If
the null hypothesis of a \"valid\" model is rejected, then something is
(most likely) wrong with either your model or selection of orthogonality
conditions.

In case the (presumed) global minimum of the moment distance function is
located in a region of the parameter space that is typically considered
unlikely ([dilemma of absurd parameters]{.title-ref}), you may opt to
choose the `penalized_estimator <penalized_estimator>`{.interpreted-text
role="opt"} option. Similar to adding priors to the likelihood, this
option incorporates prior knowledge (i.e. the prior mean) as additional
moment restrictions and weights them by their prior precision to guide
the minimization algorithm to more plausible regions of the parameter
space. Ideally, these regions are characterized by only slightly worse
values of the objective function. Note that adding prior information
comes at the cost of a loss in efficiency of the estimator.

::: {.command}
varobs VARIABLE\_NAME\...;

Required. All variables used in the `matched_moments`{.interpreted-text
role="bck"} block need to be observable. See
`varobs <varobs>`{.interpreted-text role="ref"} for more details.
:::

::: {.block}
matched\_moments ;

This block specifies the product moments which are used in estimation.
Currently, only linear product moments (e.g.
$E[y_t], E[y_t^2], E[x_t y_t], E[y_t y_{t-1}], E[y_t^3 x^2_{t-4}]$) are
supported. For other functions like $E[\log(y_t)e^{x_t}]$ you need to
declare auxiliary endogenous variables.

Each line inside of the block should be of the form:

    VARIABLE_NAME(LEAD/LAG)^POWER*VARIABLE_NAME(LEAD/LAG)^POWER*...*VARIABLE_NAME(LEAD/LAG)^POWER;

where [VARIABLE\_NAME]{.title-ref} is the name of a declared observable
variable, [LEAD/LAG]{.title-ref} is either a negative integer for lags
or a positive one for leads, and [POWER]{.title-ref} is a positive
integer indicating the exponent on the variable. You can omit
[LEAD/LAG]{.title-ref} equal to [0]{.title-ref} or [POWER]{.title-ref}
equal to [1]{.title-ref}.

*Example*

For
$E[c_t], E[y_t], E[c_t^2], E[c_t y_t], E[y_t^2], E[c_t c_{t+3}], E[y_{t+1}^2 c^3_{t-4}], E[c^3_{t-5} y_{t}^2]$
use the following block:

>     matched_moments;
>     c;
>     y;
>     c*c;
>     c*y;
>     y^2;
>     c*c(3);
>     y(1)^2*c(-4)^3;
>     c(-5)^3*y(0)^2;
>     end;

*Limitations*

1\. For GMM, Dynare can only compute the theoretical mean, covariance,
and autocovariances (i.e. first and second moments). Higher-order
moments are only supported for SMM.

2\. By default, the product moments are not demeaned, unless the
`prefilter <prefilter = INTEGER>`{.interpreted-text role="opt"} option
is set to 1. That is, by default, [c\*c]{.title-ref} corresponds to
$E[c_t^2]$ and not to $Var[c_t]=E[c_t^2]-E[c_t]^2$.

*Output*

Dynare translates the `matched_moments`{.interpreted-text role="bck"}
block into a cell array `M_.matched_moments` where:

-   the first column contains a vector of indices for the chosen
    variables in declaration order
-   the second column contains the corresponding vector of leads and
    lags
-   the third column contains the corresponding vector of powers

During the estimation phase, Dynare will eliminate all redundant or
duplicate orthogonality conditions in `M_.matched_moments` and display
which conditions were removed. In the example above, this would be the
case for the last row, which is the same as the second-to-last one. The
original block is saved in `M_.matched_moments_orig`.
:::

::: {.block}
estimated\_params ;

Required. See `estimated_params`{.interpreted-text role="bck"} for the
meaning and syntax.
:::

::: {.block}
estimated\_params\_init ;

See `estimated_params_init`{.interpreted-text role="bck"} for the
meaning and syntax.
:::

::: {.block}
estimated\_params\_bounds ;

See `estimated_params_bounds`{.interpreted-text role="bck"} for the
meaning and syntax.
:::

::: {.command}
method\_of\_moments (OPTIONS\...);

This command runs the method of moments estimation. The following
information will be displayed in the command window:

-   Overview of options chosen by the user
-   Estimation results for each stage and iteration
-   Value of minimized moment distance objective function
-   Result of the J-test
-   Table of data moments and estimated model moments

*Necessary options*

::: {.option}
mom\_method = SMM\|GMM

\"Simulated Method of Moments\" is triggered by [SMM]{.title-ref} and
\"Generalized Method of Moments\" by [GMM]{.title-ref}.
:::

::: {.option}
datafile = FILENAME

The name of the file containing the data. See
`datafile <datafile = FILENAME>`{.interpreted-text role="opt"} for the
meaning and syntax.
:::

*Options common for SMM and GMM*

::: {.option}
order = INTEGER

Order of perturbation approximation. For GMM only orders 13 are
supported. For SMM, you can choose an arbitrary order. Note that the
order set in other functions will not overwrite the default. Default:
`1`.
:::

::: {.option}
pruning

Discard higher order terms when iteratively computing simulations of the
solution. See `pruning <pruning>`{.interpreted-text role="opt"} for more
details. Default: not set for SMM, always set for GMM.
:::

::: {.option}
penalized\_estimator

This option includes deviations of the estimated parameters from the
prior mean as additional moment restrictions and weights them by their
prior precision. Default: not set.
:::

::: {.option}
weighting\_matrix = \[\'WM1\',\'WM2\',\...,\'WMn\'\]

Determines the weighting matrix used at each estimation stage. The
number of elements will define the number of stages, i.e.
`weighting_matrix = ['DIAGONAL','DIAGONAL','OPTIMAL']` performs a
three-stage estimation. Possible values for `WM` are:

> `IDENTITY_MATRIX`
>
> > Sets the weighting matrix equal to the identity matrix.
>
> `OPTIMAL`
>
> > Uses the optimal weighting matrix computed by a Newey-West-type
> > estimate with a Bartlett kernel. At the first stage, the
> > data-moments are used as initial estimate of the model moments,
> > whereas at subsequent stages the previous estimate of model moments
> > will be used when computing the optimal weighting matrix.
>
> `DIAGONAL`
>
> > Uses the diagonal of the `OPTIMAL` weighting matrix. This choice
> > puts weights on the specified moments instead of on their linear
> > combinations.
>
> `FILENAME`
>
> > The name of the mat-file (extension `.mat`) containing a
> > user-specified weighting matrix. The file must include a positive
> > definite square matrix called [weighting\_matrix]{.title-ref} with
> > both dimensions equal to the number of orthogonality conditions.

> Default value is `['DIAGONAL','OPTIMAL']`.
:::

::: {.option}
weighting\_matrix\_scaling\_factor = DOUBLE

Scaling of weighting matrix in objective function. This value should be
chosen to obtain values of the objective function in a reasonable
numerical range to prevent over- and underflows. Default: `1`.
:::

::: {.option}
bartlett\_kernel\_lag = INTEGER

Bandwidth of kernel for computing the optimal weighting matrix. Default:
`20`.
:::

::: {.option}
se\_tolx = DOUBLE

Step size for numerical differentiation when computing standard errors
with a two-sided finite difference method. Default: `1e-5`.
:::

::: {.option}
verbose

Display and store intermediate estimation results in `oo_.mom`. Default:
not set.
:::

*SMM-specific options*

::: {.option}
burnin = INTEGER

Number of periods dropped at the beginning of simulation. Default:
`500`.
:::

::: {.option}
bounded\_shock\_support

Trim shocks in simulations to $\pm 2$ standard deviations. Default: not
set.
:::

::: {.option}
seed = INTEGER

Common seed used in simulations. Default: `24051986`.
:::

::: {.option}
simulation\_multiple = INTEGER

Multiple of data length used for simulation. Default: `7`.
:::

*GMM-specific options*

::: {.option}
analytic\_standard\_errors

Compute standard errors using analytical derivatives of moments with
respect to estimated parameters. Default: not set, i.e. standard errors
are computed using a two-sided finite difference method, see
`se_tolx <se_tolx = DOUBLE>`{.interpreted-text role="opt"}.
:::

*General options*

::: {.option}
dirname = FILENAME

Directory in which to store `estimation` output. See
`dirname <dirname = FILENAME>`{.interpreted-text role="opt"} for more
details. Default: `<mod_file>`.
:::

::: {.option}
graph\_format = FORMAT

Specify the file format(s) for graphs saved to disk. See
`graph_format <graph_format = FORMAT>`{.interpreted-text role="opt"} for
more details. Default: `eps`.
:::

::: {.option}
nodisplay

See `nodisplay`{.interpreted-text role="opt"}. Default: not set.
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}. Default: not set.
:::

::: {.option}
noprint

See `noprint`{.interpreted-text role="opt"}. Default: not set.
:::

::: {.option}
plot\_priors = INTEGER

Control the plotting of priors. See
`plot_priors <plot_priors = INTEGER>`{.interpreted-text role="opt"} for
more details. Default: `1`, i.e. plot priors.
:::

::: {.option}
prior\_trunc = DOUBLE

See `prior_trunc <prior_trunc = DOUBLE>`{.interpreted-text role="opt"}
for more details. Default: `1e-10`.
:::

::: {.option}
tex

See `tex`{.interpreted-text role="opt"}. Default: not set.
:::

*Data options*

::: {.option}
first\_obs = INTEGER

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.
Default: `1`.
:::

::: {.option}
nobs = INTEGER

See `nobs <nobs = INTEGER>`{.interpreted-text role="opt"}. Default: all
observations are considered.
:::

::: {.option}
prefilter = INTEGER

A value of 1 means that the estimation procedure will demean each data
series by its empirical mean and each model moment by its theoretical
mean. See `prefilter <prefilter = INTEGER>`{.interpreted-text
role="opt"} for more details. Default: [0]{.title-ref}, i.e. no
prefiltering.
:::

::: {.option}
logdata

See `logdata <logdata>`{.interpreted-text role="opt"}. Default: not set.
:::

::: {.option}
xls\_sheet = QUOTED\_STRING

See `xls_sheet <xls_sheet = QUOTED_STRING>`{.interpreted-text
role="opt"}.
:::

::: {.option}
xls\_range = RANGE

See `xls_range <xls_range = RANGE>`{.interpreted-text role="opt"}.
:::

*Optimization options*

::: {.option}
huge\_number = DOUBLE

See `huge_number <huge_number = DOUBLE>`{.interpreted-text role="opt"}.
Default: `1e7`.
:::

::: {.option}
mode\_compute = INTEGER \| FUNCTION\_NAME

See
`mode_compute <mode_compute = INTEGER | FUNCTION_NAME>`{.interpreted-text
role="opt"}. Default: `13`, i.e. `lsqnonlin` if the Matlab Optimization
Toolbox or the Octave optim-package are present, `4`, i.e. `csminwel`
otherwise.
:::

::: {.option}
additional\_optimizer\_steps = \[INTEGERFUNCTION\_NAME,\...\]

Vector of additional minimization algorithms run after `mode_compute`.
If `verbose`{.interpreted-text role="opt"} option is set, then the
additional estimation results are saved into the `oo_.mom` structure
prefixed with [verbose\_]{.title-ref}. Default: no additional
optimization iterations.
:::

::: {.option}
optim = (NAME, VALUE, \...)

See `optim <optim = (NAME, VALUE, ...)>`{.interpreted-text role="opt"}.
:::

::: {.option}
silent\_optimizer

See `silent_optimizer`{.interpreted-text role="opt"}. Default: not set.
:::

*Numerical algorithms options*

::: {.option}
aim\_solver

See `aim_solver <aim_solver>`{.interpreted-text role="opt"}. Default:
not set.
:::

::: {.option}
k\_order\_solver

See `k_order_solver <k_order_solver>`{.interpreted-text role="opt"}.
Default: disabled for order 1 and 2, enabled for order 3 and above.
:::

::: {.option}
dr = OPTION

See `dr <dr = OPTION>`{.interpreted-text role="opt"}. Default:
`default`, i.e. generalized Schur decomposition.
:::

::: {.option}
dr\_cycle\_reduction\_tol = DOUBLE

See
`dr_cycle_reduction_tol <dr_cycle_reduction_tol = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-7`.
:::

::: {.option}
dr\_logarithmic\_reduction\_tol = DOUBLE

See
`dr_logarithmic_reduction_tol <dr_logarithmic_reduction_tol = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-12`.
:::

::: {.option}
dr\_logarithmic\_reduction\_maxiter = INTEGER

See
`dr_logarithmic_reduction_maxiter <dr_logarithmic_reduction_maxiter = INTEGER>`{.interpreted-text
role="opt"}. Default: `100`.
:::

::: {.option}
lyapunov = OPTION

See `lyapunov <lyapunov = OPTION>`{.interpreted-text role="opt"}.
Default: `default`, i.e. based on Bartlets-Stewart algorithm.
:::

::: {.option}
lyapunov\_complex\_threshold = DOUBLE

See
`lyapunov_complex_threshold <lyapunov_complex_threshold = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-15`.
:::

::: {.option}
lyapunov\_fixed\_point\_tol = DOUBLE

See
`lyapunov_fixed_point_tol <lyapunov_fixed_point_tol = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-10`.
:::

::: {.option}
lyapunov\_doubling\_tol = DOUBLE

See
`lyapunov_doubling_tol <lyapunov_doubling_tol = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-16`.
:::

::: {.option}
sylvester = OPTION

See `sylvester <sylvester = OPTION>`{.interpreted-text role="opt"}.
Default: `default`, i.e. uses `gensylv`.
:::

::: {.option}
sylvester\_fixed\_point\_tol = DOUBLE

See
`sylvester_fixed_point_tol <sylvester_fixed_point_tol = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-12`.
:::

::: {.option}
qz\_criterium = DOUBLE

See `qz_criterium <qz_criterium = DOUBLE>`{.interpreted-text
role="opt"}. Default: `0.999999` as it is assumed that the observables
are weakly stationary.
:::

::: {.option}
qz\_zero\_threshold = DOUBLE

See `qz_zero_threshold <qz_zero_threshold = DOUBLE>`{.interpreted-text
role="opt"}. Default: `1e-6`.
:::

::: {.option}
schur\_vec\_tol = DOUBLE

Tolerance level used to find nonstationary variables in Schur
decomposition of the transition matrix. Default: `1e-11`.
:::

::: {.option}
mode\_check

Plots univariate slices through the moments distance objective function
around the computed minimum for each estimated parameter. This is
helpful to diagnose problems with the optimizer. Default: not set.
:::

::: {.option}
mode\_check\_neighbourhood\_size = DOUBLE

See
`mode_check_neighbourhood_size <mode_check_neighbourhood_size = DOUBLE>`{.interpreted-text
role="opt"}. Default: `0.5`.
:::

::: {.option}
mode\_check\_symmetric\_plots = INTEGER

See
`mode_check_symmetric_plots <mode_check_symmetric_plots = INTEGER>`{.interpreted-text
role="opt"}. Default: `1`.
:::

::: {.option}
mode\_check\_number\_of\_points = INTEGER

See
`mode_check_number_of_points <mode_check_number_of_points = INTEGER>`{.interpreted-text
role="opt"}. Default: `20`.
:::

*Output*

`method_of_moments` stores user options in a structure called
[options\_mom\_]{.title-ref} in the global workspace. After running the
estimation, the parameters `M_.params` and the covariance matrices of
the shocks `M_.Sigma_e` and of the measurement errors `M_.H` are set to
the parameters that minimize the quadratic moments distance objective
function. The estimation results are stored in the `oo_.mom` structure
with the following fields:

::: {.matvar}
[oo]().mom.data\_moments

Variable set by the `method_of_moments` command. Stores the mean of the
selected empirical moments of data. NaN values due to leads/lags or
missing data are omitted when computing the mean. Vector of dimension
equal to the number of orthogonality conditions.
:::

::: {.matvar}
[oo]().mom.m\_data

Variable set by the `method_of_moments` command. Stores the selected
empirical moments at each point in time. NaN values due to leads/lags or
missing data are replaced by the corresponding mean of the moment.
Matrix of dimension time periods times number of orthogonality
conditions.
:::

::: {.matvar}
[oo]().mom.Sw

Variable set by the `method_of_moments` command. Stores the Cholesky
decomposition of the currently used weighting matrix. Square matrix of
dimensions equal to the number of orthogonality conditions.
:::

::: {.matvar}
[oo]().mom.model\_moments

Variable set by the `method_of_moments` command. Stores the implied
selected model moments given the current parameter guess. Model moments
are computed in closed-form from the pruned state-space system for GMM,
whereas for SMM these are based on averages of simulated data. Vector of
dimension equal to the number of orthogonality conditions.
:::

::: {.matvar}
[oo]().mom.Q

Variable set by the `method_of_moments` command. Stores the scalar value
of the quadratic moment\'s distance objective function.
:::

::: {.matvar}
[oo]().mom.model\_moments\_params\_derivs

Variable set by the `method_of_moments` command. Stores the analytically
computed Jacobian matrix of the derivatives of the model moments with
respect to the estimated parameters. Only for GMM with
`analytic_standard_errors`{.interpreted-text role="opt"}. Matrix with
dimension equal to the number of orthogonality conditions times number
of estimated parameters.
:::

::: {.matvar}
[oo]().[mom.gmm\_stage]()\*\_mode
:::

::: {.matvar}
[oo]().[mom.smm\_stage]()\*\_mode
:::

::: {.matvar}
[oo]().[mom.verbose\_gmm\_stage]()\*\_mode
:::

::: {.matvar}
[oo]().[mom.verbose\_smm\_stage]()\*\_mode

Variables set by the `method_of_moments` command when estimating with
GMM or SMM. Stores the estimated values at stages 1, 2,\.... The
structures contain the following fields:

-   `measurement_errors_corr`: estimated correlation between two
    measurement errors
-   `measurement_errors_std`: estimated standard deviation of
    measurement errors
-   `parameters`: estimated model parameters
-   `shocks_corr`: estimated correlation between two structural shocks.
-   `shocks_std`: estimated standard deviation of structural shocks.

If the `verbose`{.interpreted-text role="opt"} option is set, additional
fields prefixed with `verbose_` are saved for all
`additional_optimizer_steps<additional_optimizer_steps = [INTEGER|FUNCTION_NAME,INTEGER|FUNCTION_NAME,...]>`{.interpreted-text
role="opt"}.
:::

::: {.matvar}
[oo]().[mom.gmm\_stage]()\*\_std\_at\_mode
:::

::: {.matvar}
[oo]().[mom.smm\_stage]()\*\_std\_at\_mode
:::

::: {.matvar}
[oo]().[mom.verbose\_gmm\_stage]()\*\_std\_at\_mode
:::

::: {.matvar}
[oo]().[mom.verbose\_smm\_stage]()\*\_std\_at\_mode

Variables set by the `method_of_moments` command when estimating with
GMM or SMM. Stores the estimated standard errors at stages 1, 2,\....
The structures contain the following fields:

-   `measurement_errors_corr`: standard error of estimated correlation
    between two measurement errors
-   `measurement_errors_std`: standard error of estimated standard
    deviation of measurement errors
-   `parameters`: standard error of estimated model parameters
-   `shocks_corr`: standard error of estimated correlation between two
    structural shocks.
-   `shocks_std`: standard error of estimated standard deviation of
    structural shocks.

If the `verbose`{.interpreted-text role="opt"} option is set, additional
fields prefixed with `verbose_` are saved for all
`additional_optimizer_steps<additional_optimizer_steps = [INTEGER|FUNCTION_NAME,INTEGER|FUNCTION_NAME,...]>`{.interpreted-text
role="opt"}.
:::

::: {.matvar}
[oo]().mom.J\_test

Variable set by the `method_of_moments` command. Structure where the
value of the test statistic is saved into a field called `j_stat`, the
degress of freedom into a field called `degrees_freedom` and the p-value
of the test statistic into a field called `p_val`.
:::
:::

Model Comparison
----------------

::: {.command}
model\_comparison FILENAME\[(DOUBLE)\]\...; model\_comparison
(marginal\_density = ESTIMATOR) FILENAME\[(DOUBLE)\]\...;

This command computes odds ratios and estimate a posterior density over
a collection of models (see e.g. *Koop (2003)*, Ch. 1). The priors over
models can be specified as the *DOUBLE* values, otherwise a uniform
prior over all models is assumed. In contrast to frequentist
econometrics, the models to be compared do not need to be nested.
However, as the computation of posterior odds ratios is a Bayesian
technique, the comparison of models estimated with maximum likelihood is
not supported.

It is important to keep in mind that model comparison of this type is
only valid with proper priors. If the prior does not integrate to one
for all compared models, the comparison is not valid. This may be the
case if part of the prior mass is implicitly truncated because Blanchard
and Kahn conditions (instability or indeterminacy of the model) are not
fulfilled, or because for some regions of the parameters space the
deterministic steady state is undefined (or Dynare is unable to find
it). The compared marginal densities should be renormalized by the
effective prior mass, but this not done by Dynare: it is the user's
responsibility to make sure that model comparison is based on proper
priors. Note that, for obvious reasons, this is not an issue if the
compared marginal densities are based on Laplace approximations.

*Options*

::: {.option}
marginal\_density = ESTIMATOR

Specifies the estimator for computing the marginal data density.
*ESTIMATOR* can take one of the following two values: `laplace` for the
Laplace estimator or `modifiedharmonicmean` for the *Geweke (1999)*
Modified Harmonic Mean estimator. Default value: `laplace`
:::

*Output*

The results are stored in `oo_.Model_Comparison`, which is described
below.

*Example*

>     model_comparison my_model(0.7) alt_model(0.3);
>
> This example attributes a 70% prior over `my_model` and 30% prior over
> `alt_model`.
:::

::: {.matvar}
[oo]().Model\_Comparison

Variable set by the `model_comparison` command. Fields are of the form:

    oo_.Model_Comparison.FILENAME.VARIABLE_NAME

where FILENAME is the file name of the model and VARIABLE\_NAME is one
of the following:

> `Prior`
>
> > (Normalized) prior density over the model.
>
> `Log_Marginal_Density`
>
> > Logarithm of the marginal data density.
>
> `Bayes_Ratio`
>
> > Ratio of the marginal data density of the model relative to the one
> > of the first declared model
>
> `Posterior_Model_Probability`
>
> > Posterior probability of the respective model.
:::

Shock Decomposition
-------------------

::: {.command}
shock\_decomposition \[VARIABLE\_NAME\]\...; shock\_decomposition
(OPTIONS\...) \[VARIABLE\_NAME\]\...;

This command computes the historical shock decomposition for a given
sample based on the Kalman smoother, i.e. it decomposes the historical
deviations of the endogenous variables from their respective steady
state values into the contribution coming from the various shocks. The
`variable_names` provided govern for which variables the decomposition
is plotted.

Note that this command must come after either `estimation` (in case of
an estimated model) or `stoch_simul` (in case of a calibrated model).

*Options*

::: {.option}
parameter\_set = OPTION

Specify the parameter set to use for running the smoother. Possible
values for OPTION are:

> -   `calibration`
> -   `prior_mode`
> -   `prior_mean`
> -   `posterior_mode`
> -   `posterior_mean`
> -   `posterior_median`
> -   `mle_mode`

Note that the parameter set used in subsequent commands like
`stoch_simul` will be set to the specified `parameter_set`. Default
value: `posterior_mean` if Metropolis has been run, `mle_mode` if MLE
has been run.
:::

::: {.option}
datafile = FILENAME

See `datafile <dataf>`{.interpreted-text role="ref"}. Useful when
computing the shock decomposition on a calibrated model.
:::

::: {.option}
first\_obs = INTEGER

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
nobs = INTEGER

See `nobs <nobs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
prefilter = INTEGER

See `prefilter <prefilter = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
loglinear

See `loglinear <loglinear>`{.interpreted-text role="opt"}.
:::

::: {.option}
diffuse\_kalman\_tol = DOUBLE

See `diffuse_kalman_tol <diffuse_kalman_tol = DOUBLE>`{.interpreted-text
role="opt"}.
:::

::: {.option}
diffuse\_filter

See `diffuse_filter <diffuse_filter>`{.interpreted-text role="opt"}.
:::

::: {.option}
xls\_sheet = QUOTED\_STRING

See `xls_sheet <xls_sheet = QUOTED_STRING>`{.interpreted-text
role="opt"}.
:::

::: {.option}
xls\_range = RANGE

See `xls_range <xls_range = RANGE>`{.interpreted-text role="opt"}.
:::

::: {.option}
use\_shock\_groups \[= NAME\]

Uses shock grouping defined by the string instead of individual shocks
in the decomposition. The groups of shocks are defined in the
`shock_groups`{.interpreted-text role="bck"} block. If no group name is
given, `default` is assumed.
:::

::: {.option}
colormap = VARIABLE\_NAME

Controls the `colormap` used for the shocks decomposition graphs.
VARIABLE\_NAME must be the name of a MATLAB/Octave variable that has
been declared beforehand and whose value will be passed to the
MATLAB/Octave `colormap` function (see the MATLAB/Octave manual for the
list of acceptable values).
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}. Suppresses the display and
creation only within the `shock_decomposition` command, but does not
affect other commands. See `plot_shock_decomposition`{.interpreted-text
role="comm"} for plotting graphs.
:::

::: {.option}
init\_state = BOOLEAN

If equal to 0, the shock decomposition is computed conditional on the
smoothed state variables in period `0`, i.e. the smoothed shocks
starting in period 1 are used. If equal to `1`, the shock decomposition
is computed conditional on the smoothed state variables in period 1.
Default: `0`.
:::

::: {.option}
with\_epilogue

If set, then also compute the decomposition for variables declared in
the `epilogue` block (see `epilogue`{.interpreted-text role="ref"}).
:::

*Output*

::: {.matvar}
[oo]().shock\_decomposition

The results are stored in the field `oo_.shock_decomposition`, which is
a three dimensional array. The first dimension contains the
`M_.endo_nbr` endogenous variables. The second dimension stores in the
first `M_.exo_nbr` columns the contribution of the respective shocks.
Column `M_.exo_nbr+1` stores the contribution of the initial conditions,
while column `M_.exo_nbr+2` stores the smoothed value of the respective
endogenous variable in deviations from their steady state, i.e. the mean
and trends are subtracted. The third dimension stores the time periods.
Both the variables and shocks are stored in the order of declaration,
i.e. `M_.endo_names` and `M_.exo_names`, respectively.
:::
:::

::: {.block}
shock\_groups ; shock\_groups(OPTIONS\...);

Shocks can be regrouped for the purpose of shock decomposition. The
composition of the shock groups is written in a block delimited by
`shock_groups` and `end`.

Each line defines a group of shocks as a list of exogenous variables:

    SHOCK_GROUP_NAME   = VARIABLE_1 [[,] VARIABLE_2 [,]...];
    'SHOCK GROUP NAME' = VARIABLE_1 [[,] VARIABLE_2 [,]...];

*Options*

::: {.option}
name = NAME

Specifies a name for the following definition of shock groups. It is
possible to use several `shock_groups` blocks in a model file, each
grouping being identified by a different name. This name must in turn be
used in the `shock_decomposition` command. If no name is given,
`default` is used.
:::

*Example*

>     varexo e_a, e_b, e_c, e_d;
>     ...
>
>     shock_groups(name=group1);
>     supply = e_a, e_b;
>     'aggregate demand' = e_c, e_d;
>     end;
>
>     shock_decomposition(use_shock_groups=group1);
>
> This example defines a shock grouping with the name `group1`,
> containing a set of supply and demand shocks and conducts the shock
> decomposition for these two groups.
:::

::: {.command}
realtime\_shock\_decomposition \[VARIABLE\_NAME\]\...;
realtime\_shock\_decomposition (OPTIONS\...) \[VARIABLE\_NAME\]\...;

This command computes the realtime historical shock decomposition for a
given sample based on the Kalman smoother. For each period
$T=[\texttt{presample},\ldots,\texttt{nobs}]$, it recursively computes
three objects:

> -   Real-time historical shock decomposition $Y(t\vert T)$ for
>     $t=[1,\ldots,T]$, i.e. without observing data in
>     $[T+1,\ldots,\texttt{nobs}]$. This results in a standard shock
>     decomposition being computed for each additional datapoint
>     becoming available after `presample`.
> -   Forecast shock decomposition $Y(T+k\vert T)$ for
>     $k=[1,\ldots,forecast]$, i.e. the $k$-step ahead forecast made for
>     every $T$ is decomposed in its shock contributions.
> -   Real-time conditional shock decomposition of the difference
>     between the real-time historical shock decomposition and the
>     forecast shock decomposition. If `vintage <vintage =
>     INTEGER>`{.interpreted-text role="opt"} is equal to `0`, it
>     computes the effect of shocks realizing in period $T$, i.e.
>     decomposes $Y(T\vert T)-Y(T\vert T-1)$. Put differently, it
>     conducts a $1$-period ahead shock decomposition from $T-1$ to $T$,
>     by decomposing the update step of the Kalman filter. If
>     `vintage>0` and smaller than `nobs`, the decomposition is
>     conducted of the forecast revision
>     $Y(T+k\vert T+k)-Y(T+k\vert T)$.

Like `shock_decomposition`{.interpreted-text role="comm"} it decomposes
the historical deviations of the endogenous variables from their
respective steady state values into the contribution coming from the
various shocks. The `variable_names` provided govern for which variables
the decomposition is plotted.

Note that this command must come after either `estimation` (in case of
an estimated model) or `stoch_simul` (in case of a calibrated model).

*Options*

::: {.option}
parameter\_set = OPTION

See `parameter_set <parameter_set = OPTION>`{.interpreted-text
role="opt"} for possible values.
:::

::: {.option}
datafile = FILENAME

See `datafile <dataf>`{.interpreted-text role="ref"}.
:::

::: {.option}
first\_obs = INTEGER

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
nobs = INTEGER

See `nobs <nobs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
use\_shock\_groups \[= NAME\]

See `use_shock_groups <use_shock_groups [= NAME]>`{.interpreted-text
role="opt"}.
:::

::: {.option}
colormap = VARIABLE\_NAME

See `colormap <colormap = VARIABLE_NAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}. Only shock decompositions
are computed and stored in `oo_.realtime_shock_decomposition`,
`oo_.conditional_shock_decomposition` and
`oo_.realtime_forecast_shock_decomposition` but no plot is made (See
`plot_shock_decomposition`{.interpreted-text role="comm"}).
:::

::: {.option}
presample = INTEGER

Data point above which recursive realtime shock decompositions are
computed, *i.e.* for $T=[\texttt{presample+1} \ldots \texttt{nobs}]$.
:::

::: {.option}
forecast = INTEGER

Compute shock decompositions up to $T+k$ periods, i.e. get shock
contributions to k-step ahead forecasts.
:::

::: {.option}
save\_realtime = INTEGER\_VECTOR

Choose for which vintages to save the full realtime shock decomposition.
Default: `0`.
:::

::: {.option}
fast\_realtime = INTEGER fast\_realtime = \[INTEGER1:INTEGER2\]
fast\_realtime = \[INTEGER1 INTEGER2 \...\]

Runs the smoother only for the data vintages provided by the specified
integer (vector).
:::

::: {.option}
with\_epilogue

See `with_epilogue`{.interpreted-text role="opt"}.
:::

*Output*

::: {.matvar}
[oo]().realtime\_shock\_decomposition

Structure storing the results of realtime historical decompositions.
Fields are three-dimensional arrays with the first two dimension equal
to the ones of `oo_.shock_decomposition`{.interpreted-text role="mvar"}.
The third dimension stores the time periods and is therefore of size
`T+forecast`. Fields are of the form:

    oo_.realtime_shock_decomposition.OBJECT

where OBJECT is one of the following:

> `pool`
>
> > Stores the pooled decomposition, i.e. for every real-time shock
> > decomposition terminal period
> > $T=[\texttt{presample},\ldots,\texttt{nobs}]$ it collects the last
> > period's decomposition $Y(T\vert T)$ (see also
> > `plot_shock_decomposition`{.interpreted-text role="comm"}). The
> > third dimension of the array will have size `nobs+forecast`.
>
> `time_*`
>
> > Stores the vintages of realtime historical shock decompositions if
> > `save_realtime` is used. For example, if `save_realtime=[5]` and
> > `forecast=8`, the third dimension will be of size `13`.
:::

::: {.matvar}
[oo]().realtime\_conditional\_shock\_decomposition

Structure storing the results of real-time conditional decompositions.
Fields are of the form:

    oo_.realtime_conditional_shock_decomposition.OBJECT

where OBJECT is one of the following:

> `pool`
>
> > Stores the pooled real-time conditional shock decomposition, i.e.
> > collects the decompositions of $Y(T\vert T)-Y(T\vert T-1)$ for the
> > terminal periods $T=[\texttt{presample},\ldots,\texttt{nobs}]$. The
> > third dimension is of size `nobs`.
>
> `time_*`
>
> > Store the vintages of $k$-step conditional forecast shock
> > decompositions $Y(t\vert T+k)$, for $t=[T \ldots T+k]$. See `vintage
> > <vintage = INTEGER>`{.interpreted-text role="opt"}. The third
> > dimension is of size `1+forecast`.
:::

::: {.matvar}
[oo]().realtime\_forecast\_shock\_decomposition

Structure storing the results of realtime forecast decompositions.
Fields are of the form:

    oo_.realtime_forecast_shock_decomposition.OBJECT

where `OBJECT` is one of the following:

> `pool`
>
> > Stores the pooled real-time forecast decomposition of the $1$-step
> > ahead effect of shocks on the $1$-step ahead prediction, i.e.
> > $Y(T\vert
> > T-1)$.
>
> `time_*`
>
> > Stores the vintages of $k$-step out-of-sample forecast shock
> > decompositions, i.e. $Y(t\vert
> > T)$, for $t=[T \ldots T+k]$. See `vintage
> > <vintage = INTEGER>`{.interpreted-text role="opt"}.
:::
:::

::: {.command}
plot\_shock\_decomposition \[VARIABLE\_NAME\]\...;
plot\_shock\_decomposition (OPTIONS\...) \[VARIABLE\_NAME\]\...;

This command plots the historical shock decomposition already computed
by `shock_decomposition` or `realtime_shock_decomposition`. For that
reason, it must come after one of these commands. The `variable_names`
provided govern which variables the decomposition is plotted for.

Further note that, unlike the majority of Dynare commands, the options
specified below are overwritten with their defaults before every call to
`plot_shock_decomposition`. Hence, if you want to reuse an option in a
subsequent call to `plot_shock_decomposition`, you must pass it to the
command again.

*Options*

::: {.option}
use\_shock\_groups \[= NAME\]

See `use_shock_groups <use_shock_groups [= NAME]>`{.interpreted-text
role="opt"}.
:::

::: {.option}
colormap = VARIABLE\_NAME

See `colormap <colormap = VARIABLE_NAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
nodisplay

See `nodisplay`{.interpreted-text role="opt"}.
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}.
:::

::: {.option}
graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )

See `graph_format <graph_format = FORMAT>`{.interpreted-text
role="opt"}.
:::

::: {.option}
detail\_plot

Plots shock contributions using subplots, one per shock (or group of
shocks). Default: not activated
:::

::: {.option}
interactive

Under MATLAB, add uimenus for detailed group plots. Default: not
activated
:::

::: {.option}
screen\_shocks

For large models (i.e. for models with more than 16 shocks), plots only
the shocks that have the largest historical contribution for chosen
selected `variable_names`. Historical contribution is ranked by the mean
absolute value of all historical contributions.
:::

::: {.option}
steadystate

If passed, the the $y$-axis value of the zero line in the shock
decomposition plot is translated to the steady state level. Default: not
activated
:::

::: {.option}
type = qoq \| yoy \| aoa

For quarterly data, valid arguments are: `qoq` for quarter-on-quarter
plots, `yoy` for year-on-year plots of growth rates, `aoa` for
annualized variables, i.e. the value in the last quarter for each year
is plotted. Default value: empty, i.e. standard period-on-period plots
(`qoq` for quarterly data).
:::

::: {.option}
fig\_name = STRING

Specifies a user-defined keyword to be appended to the default figure
name set by `plot_shock_decomposition`. This can avoid to overwrite
plots in case of sequential calls to `plot_shock_decomposition`.
:::

::: {.option}
write\_xls

Saves shock decompositions to Excel-file in the main directory, named
`FILENAME_shock_decomposition_TYPE_FIG_NAME.xls`. This option requires
your system to be configured to be able to write Excel files.[^8]
:::

::: {.option}
realtime = INTEGER

Which kind of shock decomposition to plot. INTEGER can take the
following values:

> -   `0`: standard historical shock decomposition. See
>     `shock_decomposition`{.interpreted-text role="comm"}.
> -   `1`: realtime historical shock decomposition. See
>     `realtime_shock_decomposition`{.interpreted-text role="comm"}.
> -   `2`: conditional realtime shock decomposition. See
>     `realtime_shock_decomposition`{.interpreted-text role="comm"}.
> -   `3`: realtime forecast shock decomposition. See
>     `realtime_shock_decomposition`{.interpreted-text role="comm"}.

If no vintage is requested, i.e. `vintage=0` then the pooled objects
from `realtime_shock_decomposition`{.interpreted-text role="comm"} will
be plotted and the respective vintage otherwise. Default: `0`.
:::

::: {.option}
vintage = INTEGER

Selects a particular data vintage in $[presample,\ldots,nobs]$ for which
to plot the results from
`realtime_shock_decomposition`{.interpreted-text role="comm"} selected
via the `realtime <realtime = INTEGER>`{.interpreted-text role="opt"}
option. If the standard historical shock decomposition is selected
(`realtime=0`), `vintage` will have no effect. If `vintage=0` the pooled
objects from `realtime_shock_decomposition`{.interpreted-text
role="comm"} will be plotted. If `vintage>0`, it plots the shock
decompositions for vintage $T=\texttt{vintage}$ under the following
scenarios:

> -   `realtime=1`: the full vintage shock decomposition $Y(t\vert T)$
>     for $t=[1,\ldots,T]$
> -   `realtime=2`: the conditional forecast shock decomposition from
>     $T$, i.e. plots $Y(T+j\vert T+j)$ and the shock contributions
>     needed to get to the data $Y(T+j)$ conditional on $T=$ vintage,
>     with $j=[0,\ldots,\texttt{forecast}]$.
> -   `realtime=3`: plots unconditional forecast shock decomposition
>     from $T$, i.e. $Y(T+j\vert
>     T)$, where $T=\texttt{vintage}$ and
>     $j=[0,\ldots,\texttt{forecast}]$.

Default: `0`.
:::

::: {.option}
plot\_init\_date = DATE

If passed, plots decomposition using `plot_init_date` as initial period.
Default: first observation in estimation
:::

::: {.option}
plot\_end\_date = DATE

If passed, plots decomposition using `plot_end_date` as last period.
Default: last observation in estimation
:::

::: {.option}
diff

If passed, plot the decomposition of the first difference of the list of
variables. If used in combination with `flip`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated
:::

::: {.option}
flip

If passed, plot the decomposition of the opposite of the list of
variables. If used in combination with `diff`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated
:::

::: {.option}
max\_nrows

Maximum number of rows in the subplot layout of detailed shock
decomposition graphs. Note that columns are always 3. Default: 6
:::

::: {.option}
with\_epilogue

See `with_epilogue`{.interpreted-text role="opt"}.
:::

::: {.option}
init2shocks init2shocks = NAME

Use the information contained in an `init2shocks`{.interpreted-text
role="bck"} block, in order to attribute initial conditions to shocks.
The name of the block can be explicitly given, otherwise it defaults to
the `default` block.
:::
:::

::: {.block}
init2shocks ; init2shocks (OPTIONS\...);

This blocks gives the possibility of attributing the initial condition
of endogenous variables to the contribution of exogenous variables in
the shock decomposition.

For example, in an AR(1) process, the contribution of the initial
condition on the process variable can naturally be assigned to the
innovation of the process.

Each line of the block should have the syntax:

    VARIABLE_1 [,] VARIABLE_2;

Where VARIABLE\_1 is an endogenous variable whose initial condition will
be attributed to the exogenous VARIABLE\_2.

The information contained in this block is used by the
`plot_shock_decomposition`{.interpreted-text role="comm"} command when
given the `init2shocks` option.

*Options*

::: {.option}
name = NAME

Specifies a name for the block, that can be referenced from
`plot_shock_decomposition`, so that several such blocks can coexist in a
single model file. If the name is unspecified, it defaults to `default`.
:::

*Example*

>     var y y_s R pie dq pie_s de A y_obs pie_obs R_obs;
>     varexo e_R e_q e_ys e_pies e_A;
>     ...
>
>     model;
>       dq = rho_q*dq(-1)+e_q;
>       A = rho_A*A(-1)+e_A;
>       ...
>     end;
>
>     ...
>
>     init2shocks;
>       dq e_q;
>       A e_A;
>     end;
>
>     shock_decomposition(nograph);
>
>     plot_shock_decomposition(init2shocks) y_obs R_obs pie_obs dq de;
>
> In this example, the initial conditions of `dq` and `A` will be
> respectively attributed to `e_q` and `e_A`.
:::

::: {.command}
initial\_condition\_decomposition \[VARIABLE\_NAME\]\...;
initial\_condition\_decomposition (OPTIONS\...) \[VARIABLE\_NAME\]\...;

This command computes and plots the decomposition of the effect of
smoothed initial conditions of state variables. The `variable_names`
provided govern which variables the decomposition is plotted for.

Further note that, unlike the majority of Dynare commands, the options
specified below are overwritten with their defaults before every call to
`initial_condition_decomposition`. Hence, if you want to reuse an option
in a subsequent call to `initial_condition_decomposition`, you must pass
it to the command again.

*Options*

::: {.option}
colormap = VARIABLE\_NAME

See `colormap <colormap = VARIABLE_NAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
nodisplay

See `nodisplay`{.interpreted-text role="opt"}.
:::

::: {.option}
graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )

See `graph_format <graph_format = FORMAT>`{.interpreted-text
role="opt"}.
:::

::: {.option}
detail\_plot

Plots shock contributions using subplots, one per shock (or group of
shocks). Default: not activated
:::

::: {.option}
steadystate

If passed, the the $y$-axis value of the zero line in the shock
decomposition plot is translated to the steady state level. Default: not
activated
:::

::: {.option}
type = qoq \| yoy \| aoa

For quarterly data, valid arguments are: `qoq` for quarter-on-quarter
plots, `yoy` for year-on-year plots of growth rates, `aoa` for
annualized variables, i.e. the value in the last quarter for each year
is plotted. Default value: empty, i.e. standard period-on-period plots
(`qoq` for quarterly data).
:::

::: {.option}
fig\_name = STRING

Specifies a user-defined keyword to be appended to the default figure
name set by `plot_shock_decomposition`. This can avoid to overwrite
plots in case of sequential calls to `plot_shock_decomposition`.
:::

::: {.option}
write\_xls

Saves shock decompositions to Excel-file in the main directory, named
`FILENAME_shock_decomposition_TYPE_FIG_NAME_initval.xls`. This option
requires your system to be configured to be able to write Excel
files.[^9]
:::

::: {.option}
plot\_init\_date = DATE

If passed, plots decomposition using `plot_init_date` as initial period.
Default: first observation in estimation
:::

::: {.option}
plot\_end\_date = DATE

If passed, plots decomposition using `plot_end_date` as last period.
Default: last observation in estimation
:::

::: {.option}
diff

If passed, plot the decomposition of the first difference of the list of
variables. If used in combination with `flip`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated
:::

::: {.option}
flip

If passed, plot the decomposition of the opposite of the list of
variables. If used in combination with `diff`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated
:::
:::

::: {.command}
squeeze\_shock\_decomposition \[VARIABLE\_NAME\]\...;

For large models, the size of the information stored by shock
decompositions (especially various settings of realtime decompositions)
may become huge. This command allows to squeeze this information in two
possible ways:

> -   Automatic (default): only the variables for which plotting has
>     been explicitly required with `plot_shock_decomposition` will have
>     their decomposition left in `oo_` after this command is run;
> -   If a list of variables is passed to the command, then only those
>     variables will have their decomposition left in `oo_` after this
>     command is run.
:::

Calibrated Smoother
-------------------

Dynare can also run the smoother on a calibrated model:

::: {.command}
calib\_smoother \[VARIABLE\_NAME\]\...; calib\_smoother (OPTIONS\...)
\[VARIABLE\_NAME\]\...;

This command computes the smoothed variables (and possible the filtered
variables) on a calibrated model.

A datafile must be provided, and the observable variables declared with
`varobs`. The smoother is based on a first-order approximation of the
model.

By default, the command computes the smoothed variables and shocks and
stores the results in `oo_.SmoothedVariables` and `oo_.SmoothedShocks`.
It also fills `oo_.UpdatedVariables`.

*Options*

::: {.option}
datafile = FILENAME

See `datafile <dataf>`{.interpreted-text role="ref"}.
:::

::: {.option}
filtered\_vars

Triggers the computation of filtered variables. See
`filtered_vars`{.interpreted-text role="opt"}, for more details.
:::

::: {.option}
filter\_step\_ahead = \[INTEGER1:INTEGER2\]

See
`filter_step_ahead <filter_step_ahead = [INTEGER1:INTEGER2]>`{.interpreted-text
role="opt"}.
:::

::: {.option}
prefilter = INTEGER

See `prefilter <prefilter = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
parameter\_set = OPTION

See `parameter_set <parameter_set = OPTION>`{.interpreted-text
role="opt"} for possible values. Default: `calibration`.
:::

::: {.option}
loglinear

See `loglinear <logl>`{.interpreted-text role="ref"}.
:::

::: {.option}
first\_obs = INTEGER

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
filter\_decomposition

See `filter_decomposition`{.interpreted-text role="opt"}.
:::

::: {.option}
filter\_covariance

See `filter_covariance`{.interpreted-text role="opt"}.
:::

::: {.option}
smoother\_redux

See `smoother_redux`{.interpreted-text role="opt"}.
:::

::: {.option}
kalman\_algo = INTEGER

See `kalman_algo <kalman_algo = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
diffuse\_filter = INTEGER

See `diffuse_filter`{.interpreted-text role="opt"}.
:::

::: {.option}
diffuse\_kalman\_tol = DOUBLE

See `diffuse_kalman_tol <diffuse_kalman_tol = DOUBLE>`{.interpreted-text
role="opt"}.
:::

::: {.option}
xls\_sheet = QUOTED\_STRING

See `xls_sheet <xls_sheet = QUOTED_STRING>`{.interpreted-text
role="opt"}.
:::

::: {.option}
xls\_range = RANGE

See `xls_range <xls_range = RANGE>`{.interpreted-text role="opt"}.
:::
:::

Forecasting {#fore}
-----------

On a calibrated model, forecasting is done using the `forecast` command.
On an estimated model, use the `forecast` option of `estimation`
command.

It is also possible to compute forecasts on a calibrated or estimated
model for a given constrained path of the future endogenous variables.
This is done, from the reduced form representation of the DSGE model, by
finding the structural shocks that are needed to match the restricted
paths. Use `conditional_forecast`{.interpreted-text role="comm"},
`conditional_forecast_paths`{.interpreted-text role="bck"} and
`plot_conditional_forecast`{.interpreted-text role="comm"} for that
purpose.

Finally, it is possible to do forecasting with a Bayesian VAR using the
`bvar_forecast`{.interpreted-text role="comm"} command.

::: {.command}
forecast \[VARIABLE\_NAME\...\]; forecast (OPTIONS\...)
\[VARIABLE\_NAME\...\];

This command computes a simulation of a stochastic model from an
arbitrary initial point.

When the model also contains deterministic exogenous shocks, the
simulation is computed conditionally to the agents knowing the future
values of the deterministic exogenous variables.

`forecast` must be called after `stoch_simul`.

`forecast` plots the trajectory of endogenous variables. When a list of
variable names follows the command, only those variables are plotted. A
90% confidence interval is plotted around the mean trajectory. Use
option `conf_sig` to change the level of the confidence interval.

*Options*

::: {.option}
periods = INTEGER

Number of periods of the forecast. Default: `5`.
:::

::: {#confsig}
::: {.option}
conf\_sig = DOUBLE

Level of significance for confidence interval. Default: `0.90`.
:::
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}.
:::

::: {.option}
nodisplay

See `nodisplay`{.interpreted-text role="opt"}.
:::

::: {.option}
graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )

See `graph_format = FORMAT`{.interpreted-text role="opt"}.
:::

*Initial Values*

`forecast` computes the forecast taking as initial values the values
specified in `histval` (see `histval`{.interpreted-text role="bck"}).
When no `histval` block is present, the initial values are the one
stated in `initval`. When `initval` is followed by command `steady`, the
initial values are the steady state (see `steady`{.interpreted-text
role="comm"}).

*Output*

The results are stored in `oo_.forecast`, which is described below.

*Example*

>     varexo_det tau;
>
>     varexo e;
>     ...
>     shocks;
>     var e; stderr 0.01;
>     var tau;
>     periods 1:9;
>     values -0.15;
>     end;
>
>     stoch_simul(irf=0);
>
>     forecast;

::: {.matvar}
[oo]().forecast

Variable set by the `forecast` command, or by the `estimation` command
if used with the `forecast` option and ML or if no Metropolis-Hastings
has been computed (in that case, the forecast is computed for the
posterior mode). Fields are of the form:

    oo_.forecast.FORECAST_MOMENT.VARIABLE_NAME

where `FORECAST_MOMENT` is one of the following:

> `HPDinf`
>
> > Lower bound of a 90% HPD interval[^10] of forecast due to parameter
> > uncertainty, but ignoring the effect of measurement error on
> > observed variables. In case of ML, it stores the lower bound of the
> > confidence interval.
>
> `HPDsup`
>
> > Upper bound of a 90% HPD forecast interval due to parameter
> > uncertainty, but ignoring the effect of measurement error on
> > observed variables. In case of ML, it stores the upper bound of the
> > confidence interval.
>
> `HPDinf_ME`
>
> > Lower bound of a 90% HPD interval[^11] of forecast for observed
> > variables due to parameter uncertainty and measurement error. In
> > case of ML, it stores the lower bound of the confidence interval.
>
> `HPDsup_ME`
>
> > Upper bound of a 90% HPD interval of forecast for observed variables
> > due to parameter uncertainty and measurement error. In case of ML,
> > it stores the upper bound of the confidence interval.
>
> `Mean`
>
> > Mean of the posterior distribution of forecasts.
:::

::: {.matvar}
[oo]().PointForecast

Set by the `estimation` command, if it is used with the `forecast`
option and if either `mh_replic > 0` or the `load_mh_file` option are
used.

Contains the distribution of forecasts taking into account the
uncertainty about both parameters and shocks.

Fields are of the form:

    oo_.PointForecast.MOMENT_NAME.VARIABLE_NAME
:::

::: {.matvar}
[oo]().MeanForecast

Set by the `estimation` command, if it is used with the `forecast`
option and if either `mh_replic > 0` or `load_mh_file` option are used.

Contains the distribution of forecasts where the uncertainty about
shocks is averaged out. The distribution of forecasts therefore only
represents the uncertainty about parameters.

Fields are of the form:

    oo_.MeanForecast.MOMENT_NAME.VARIABLE_NAME
:::
:::

::: {.command}
conditional\_forecast (OPTIONS\...);

This command computes forecasts on an estimated or calibrated model for
a given constrained path of some future endogenous variables. This is
done using the reduced form first order state-space representation of
the DSGE model by finding the structural shocks that are needed to match
the restricted paths. Consider the augmented state space representation
that stacks both predetermined and non-predetermined variables into a
vector $y_{t}$:

> $$y_t=Ty_{t-1}+R\varepsilon_t$$

Both $y_t$ and $\varepsilon_t$ are split up into controlled and
uncontrolled ones, and we assume without loss of generality that the
constrained endogenous variables and the controlled shocks come first :

> $$\begin{aligned}
> \begin{pmatrix}
> y_{c,t}\\
> y_{u,t}
> \end{pmatrix}
> =
> \begin{pmatrix}
> T_{c,c} & T_{c,u}\\
> T_{u,c} & T_{u,u}
> \end{pmatrix}
> \begin{pmatrix}
> y_{c,t-1}\\
> y_{u,t-1}
> \end{pmatrix}
> +
> \begin{pmatrix}
> R_{c,c} & R_{c,u}\\
> R_{u,c} & R_{u,u}
> \end{pmatrix}
> \begin{pmatrix}
> \varepsilon_{c,t}\\
> \varepsilon_{u,t}
> \end{pmatrix}
> \end{aligned}$$

where matrices $T$ and $R$ are partitioned consistently with the vectors
of endogenous variables and innovations. Provided that matrix $R_{c,c}$
is square and full rank (a necessary condition is that the number of
free endogenous variables matches the number of free innovations), given
$y_{c,t}$, $\varepsilon_{u,t}$ and $y_{t-1}$ the first block of
equations can be solved for $\varepsilon_{c,t}$:

> $$\varepsilon_{c,t} = R_{c,c}^{-1}\bigl( y_{c,t} - T_{c,c}y_{c,t} - T_{c,u}y_{u,t}  - R_{c,u}\varepsilon_{u,t}\bigr)$$

and $y_{u,t}$ can be updated by evaluating the second block of
equations:

> $$y_{u,t} = T_{u,c}y_{c,t-1} + T_{u,u}y_{u,t-1} +  R_{u,c}\varepsilon_{c,t} + R_{u,u}\varepsilon_{u,t}$$

By iterating over these two blocks of equations, we can build a forecast
for all the endogenous variables in the system conditional on paths for
a subset of the endogenous variables. If the distribution of the free
innovations $\varepsilon_{u,t}$ is provided (*i.e.* some of them have
positive variances) this exercise is replicated (the number of
replication is controlled by the option `replic`{.interpreted-text
role="opt"} described below) by drawing different sequences of free
innovations. The result is a predictive distribution for the
uncontrolled endogenous variables, $y_{u,t}$, that Dynare will use to
report confidence bands around the point conditional forecast.

A few things need to be noted. First, the controlled exogenous variables
are set to zero for the uncontrolled periods. This implies that there is
no forecast uncertainty arising from these exogenous variables in
uncontrolled periods. Second, by making use of the first order state
space solution, even if a higher-order approximation was performed, the
conditional forecasts will be based on a first order approximation.
Since the controlled exogenous variables are identified on the basis of
the reduced form model (*i.e.* after solving for the expectations), they
are unforeseen shocks from the perspective of the agents in the model.
That is, agents expect the endogenous variables to return to their
respective steady state levels but are surprised in each period by the
realisation of shocks keeping the endogenous variables along a
predefined (unexpected) path. Fourth, if the structural innovations are
correlated, because the calibrated or estimated covariance matrix has
non zero off diagonal elements, the results of the conditional forecasts
will depend on the ordering of the innovations (as declared after
`varexo`). As in VAR models, a Cholesky decomposition is used to
factorise the covariance matrix and identify orthogonal impulses. It is
preferable to declare the correlations in the model block (explicitly
imposing the identification restrictions), unless you are satisfied with
the implicit identification restrictions implied by the Cholesky
decomposition.

This command has to be called after `estimation` or `stoch_simul`.

Use `conditional_forecast_paths`{.interpreted-text role="bck"} block to
give the list of constrained endogenous, and their constrained future
path. Option `controlled_varexo` is used to specify the structural
shocks which will be matched to generate the constrained path.

Use `plot_conditional_forecast`{.interpreted-text role="comm"} to graph
the results.

*Options*

::: {.option}
parameter\_set = OPTION

See `parameter_set <parameter_set = OPTION>`{.interpreted-text
role="opt"} for possible values. No default value, mandatory option.
:::

::: {.option}
controlled\_varexo = (VARIABLE\_NAME\...)

Specify the exogenous variables to use as control variables. No default
value, mandatory option.
:::

::: {.option}
periods = INTEGER

Number of periods of the forecast. Default: `40`. `periods` cannot be
smaller than the number of constrained periods.
:::

::: {.option}
replic = INTEGER

Number of simulations used to compute the conditional forecast
uncertainty. Default: `5000`.
:::

::: {.option}
conf\_sig = DOUBLE

Level of significance for confidence interval. Default: `0.80`.
:::

*Output*

The results are stored in `oo_.conditional_forecast`, which is described
below.

*Example*

>     var y a;
>     varexo e u;
>     ...
>     estimation(...);
>
>     conditional_forecast_paths;
>     var y;
>     periods 1:3, 4:5;
>     values 2, 5;
>     var a;
>     periods 1:5;
>     values 3;
>     end;
>
>     conditional_forecast(parameter_set = calibration, controlled_varexo = (e, u), replic = 3000);
>
>     plot_conditional_forecast(periods = 10) a y;

::: {.matvar}
[oo]().conditional\_forecast.cond

Variable set by the `conditional_forecast` command. It stores the
conditional forecasts. Fields are `periods+1` by `1` vectors storing the
steady state (time 0) and the subsequent `periods` forecasts periods.
Fields are of the form:

    oo_.conditional_forecast.cond.FORECAST_MOMENT.VARIABLE_NAME

where FORECAST\_MOMENT is one of the following:

> `Mean`
>
> > Mean of the conditional forecast distribution.
>
> `ci`
>
> > Confidence interval of the conditional forecast distribution. The
> > size corresponds to `conf_sig`.
:::

::: {.matvar}
[oo]().conditional\_forecast.uncond

Variable set by the `conditional_forecast` command. It stores the
unconditional forecasts. Fields are of the form:

    oo_.conditional_forecast.uncond.FORECAST_MOMENT.VARIABLE_NAME
:::

::: {.matvar}
forecasts.instruments

Variable set by the `conditional_forecast command`. Stores the names of
the exogenous instruments.
:::

::: {.matvar}
[oo]().conditional\_forecast.controlled\_variables

Variable set by the `conditional_forecast` command. Stores the position
of the constrained endogenous variables in declaration order.
:::

::: {.matvar}
[oo]().conditional\_forecast.controlled\_exo\_variables

Variable set by the `conditional_forecast` command. Stores the values of
the controlled exogenous variables underlying the conditional forecasts
to achieve the constrained endogenous variables. Fields are
`[number of constrained periods]` by `1` vectors and are of the form:

    oo_.conditional_forecast.controlled_exo_variables.FORECAST_MOMENT.SHOCK_NAME
:::

::: {.matvar}
[oo]().conditional\_forecast.graphs

Variable set by the `conditional_forecast` command. Stores the
information for generating the conditional forecast plots.
:::
:::

::: {.block}
conditional\_forecast\_paths ;

Describes the path of constrained endogenous, before calling
`conditional_forecast`. The syntax is similar to deterministic shocks in
`shocks`, see `conditional_forecast` for an example.

The syntax of the block is the same as for the deterministic shocks in
the `shocks` blocks (see `shocks-exo`{.interpreted-text role="ref"}).
Note that you need to specify the full path for all constrained
endogenous variables between the first and last specified period. If an
intermediate period is not specified, a value of 0 is assumed. That is,
if you specify only values for periods 1 and 3, the values for period 2
will be 0. Currently, it is not possible to have uncontrolled
intermediate periods.

It is however possible to have different number of controlled periods
for different variables. In that case, the order of declaration of
endogenous controlled variables and of `controlled_varexo` matters: if
the second endogenous variable is controlled for less periods than the
first one, the second `controlled_varexo` isn\'t set for the last
periods.

In case of the presence of `observation_trends`, the specified
controlled path for these variables needs to include the trend
component. When using the `loglinear <logl>`{.interpreted-text
role="ref"} option, it is necessary to specify the logarithm of the
controlled variables.
:::

::: {.block}
filter\_initial\_state ;

This block specifies the initial values of the endogenous states at the
beginning of the Kalman filter recursions. That is, if the Kalman filter
recursion starts with time t=1 being the first observation, this block
provides the state estimate at time 0 given information at time 0,
$E_0(x_0)$. If nothing is specified, the initial condition is assumed to
be at the steady state (which is the unconditional mean for a stationary
model).

This block is terminated by `end;`.

Each line inside of the block should be of the form:

    VARIABLE_NAME(INTEGER)=EXPRESSION;

`EXPRESSION` is any valid expression returning a numerical value and can
contain parameter values. This allows specifying relationships that will
be honored during estimation. `INTEGER` refers to the lag with which a
variable appears. By convention in Dynare, period 1 is the first period.
Going backwards in time, the first period before the start of the
simulation is period 0, then period -1, and so on. Note that the
`filter_initial_state` block does not take non-state variables.

*Example*

>     filter_initial_state;
>     k(0)= ((1/bet-(1-del))/alp)^(1/(alp-1))*l_ss;
>     P(0)=2.5258;
>     m(0)= mst;
>     end;
:::

::: {.command}
plot\_conditional\_forecast \[VARIABLE\_NAME\...\];
plot\_conditional\_forecast (periods = INTEGER) \[VARIABLE\_NAME\...\];

Plots the conditional (plain lines) and unconditional (dashed lines)
forecasts.

To be used after `conditional_forecast`.

*Options*

::: {.option}
periods = INTEGER

Number of periods to be plotted. Default: equal to periods in
`conditional_forecast`. The number of periods declared in
`plot_conditional_forecast` cannot be greater than the one declared in
`conditional_forecast`.
:::
:::

::: {.command}
bvar\_forecast ;

This command computes (out-of-sample) forecasts for an estimated BVAR
model, using Minnesota priors.

See `bvar-a-la-sims.pdf`, which comes with Dynare distribution, for more
information on this command.
:::

If the model contains strong non-linearities or if some perfectly
expected shocks are considered, the forecasts and the conditional
forecasts can be computed using an extended path method. The forecast
scenario describing the shocks and/or the constrained paths on some
endogenous variables should be build. The first step is the forecast
scenario initialization using the function `init_plan`:

::: {.matcomm}
HANDLE = init\_plan (DATES);

Creates a new forecast scenario for a forecast period (indicated as a
dates class, see `dates class members
<dates-members>`{.interpreted-text role="ref"}). This function return a
handle on the new forecast scenario.
:::

The forecast scenario can contain some simple shocks on the exogenous
variables. This shocks are described using the function `basic_plan`:

::: {.matcomm}
HANDLE = basic\_plan (HANDLE, \`VAR\_NAME\', \`SHOCK\_TYPE\', DATES,
MATLAB VECTOR OF DOUBLE \| \[DOUBLE \| EXPR \[DOUBLE \| EXPR\] \] );

Adds to the forecast scenario a shock on the exogenous variable
indicated between quotes in the second argument. The shock type has to
be specified in the third argument between quotes: 'surprise' in case of
an unexpected shock or 'perfect\_foresight' for a perfectly anticipated
shock. The fourth argument indicates the period of the shock using a
dates class (see `dates class
members <dates-members>`{.interpreted-text role="ref"}). The last
argument is the shock path indicated as a MATLAB vector of double. This
function return the handle of the updated forecast scenario.
:::

The forecast scenario can also contain a constrained path on an
endogenous variable. The values of the related exogenous variable
compatible with the constrained path are in this case computed. In other
words, a conditional forecast is performed. This kind of shock is
described with the function `flip_plan`:

::: {.matcomm}
HANDLE = flip\_plan (HANDLE, \`VAR\_NAME\', \`VAR\_NAME\',
\`SHOCK\_TYPE\', DATES, MATLAB VECTOR OF DOUBLE \| \[DOUBLE \| EXPR
\[DOUBLE \| EXPR\] \] );

Adds to the forecast scenario a constrained path on the endogenous
variable specified between quotes in the second argument. The associated
exogenous variable provided in the third argument between quotes, is
considered as an endogenous variable and its values compatible with the
constrained path on the endogenous variable will be computed. The nature
of the expectation on the constrained path has to be specified in the
fourth argument between quotes: 'surprise' in case of an unexpected path
or 'perfect\_foresight' for a perfectly anticipated path. The fifth
argument indicates the period where the path of the endogenous variable
is constrained using a dates class (see `dates class
members <dates-members>`{.interpreted-text role="ref"}). The last
argument contains the constrained path as a MATLAB vector of double.
This function return the handle of the updated forecast scenario.
:::

Once the forecast scenario if fully described, the forecast is computed
with the command `det_cond_forecast`:

::: {.matcomm}
DSERIES = det\_cond\_forecast (HANDLE\[, DSERIES \[, DATES\]\]);

Computes the forecast or the conditional forecast using an extended path
method for the given forecast scenario (first argument). The past values
of the endogenous and exogenous variables provided with a dseries class
(see `dseries class
members <dseries-members>`{.interpreted-text role="ref"}) can be
indicated in the second argument. By default, the past values of the
variables are equal to their steady-state values. The initial date of
the forecast can be provided in the third argument. By default, the
forecast will start at the first date indicated in the
`init_plan command`. This function returns a dset containing the
historical and forecast values for the endogenous and exogenous
variables.
:::

*Example*

>     % conditional forecast using extended path method
>     % with perfect foresight on r path
>
>     var y r;
>     varexo e u;
>     ...
>     smoothed = dseries('smoothed_variables.csv');
>
>     fplan = init_plan(2013Q4:2029Q4);
>     fplan = flip_plan(fplan, 'y', 'u', 'surprise', 2013Q4:2014Q4,  [1 1.1 1.2 1.1 ]);
>     fplan = flip_plan(fplan, 'r', 'e', 'perfect_foresight', 2013Q4:2014Q4,  [2 1.9 1.9 1.9 ]);
>
>     dset_forecast = det_cond_forecast(fplan, smoothed);
>
>     plot(dset_forecast.{'y','u'});
>     plot(dset_forecast.{'r','e'});

::: {.command}
smoother2histval ; smoother2histval(OPTIONS\...);

The purpose of this command is to construct initial conditions (for a
subsequent simulation) that are the smoothed values of a previous
estimation.

More precisely, after an estimation run with the `smoother` option,
`smoother2histval` will extract the smoothed values (from
`oo_.SmoothedVariables`, and possibly from `oo_.SmoothedShocks` if there
are lagged exogenous), and will use these values to construct initial
conditions (as if they had been manually entered through `histval`).

*Options*

::: {.option}
period = INTEGER

Period number to use as the starting point for the subsequent
simulation. It should be between 1 and the number of observations that
were used to produce the smoothed values. Default: the last observation.
:::

::: {.option}
infile = FILENAME

Load the smoothed values from a `_results.mat` file created by a
previous Dynare run. Default: use the smoothed values currently in the
global workspace.
:::

::: {.option}
invars = ( VARIABLE\_NAME \[VARIABLE\_NAME \...\] )

A list of variables to read from the smoothed values. It can contain
state endogenous variables, and also exogenous variables having a lag.
Default: all the state endogenous variables, and all the exogenous
variables with a lag.
:::

::: {.option}
outfile = FILENAME

Write the initial conditions to a file. Default: write the initial
conditions in the current workspace, so that a simulation can be
performed.
:::

::: {.option}
outvars = ( VARIABLE\_NAME \[VARIABLE\_NAME \...\] )

A list of variables which will be given the initial conditions. This
list must have the same length than the list given to `invars`, and
there will be a one-to-one mapping between the two list. Default: same
value as option `invars`.
:::

*Use cases*

There are three possible ways of using this command:

> -   Everything in a single file: run an estimation with a smoother,
>     then run `smoother2histval` (without the `infile` and `outfile`
>     options), then run a stochastic simulation.
> -   In two files: in the first file, run the smoother and then run
>     `smoother2histval` with the `outfile` option; in the second file,
>     run `histval_file` to load the initial conditions, and run a
>     (deterministic or stochastic) simulation.
> -   In two files: in the first file, run the smoother; in the second
>     file, run `smoother2histval` with the `infile` option equal to the
>     `_results.mat` file created by the first file, and then run a
>     (deterministic or stochastic) simulation.
:::

Optimal policy
--------------

Dynare has tools to compute optimal policies for various types of
objectives. You can either solve for optimal policy under commitment
with `ramsey_model`, for optimal policy under discretion with
`discretionary_policy` or for optimal simple rules with `osr` (also
implying commitment).

::: {.command}
planner\_objective MODEL\_EXPRESSION ;

This command declares the policy maker objective, for use with
`ramsey_model` or `discretionary_policy`.

You need to give the one-period objective, not the discounted lifetime
objective. The discount factor is given by the `planner_discount` option
of `ramsey_model` and `discretionary_policy`. The objective function can
only contain current endogenous variables and no exogenous ones. This
limitation is easily circumvented by defining an appropriate auxiliary
variable in the model.

With `ramsey_model`, you are not limited to quadratic objectives: you
can give any arbitrary nonlinear expression.

With `discretionary_policy`, the objective function must be quadratic.
:::

::: {.command}
evaluate\_planner\_objective; evaluate\_planner\_objective
(OPTIONS\...);

This command computes, displays, and stores the value of the planner
objective function under Ramsey policy or discretion in
`oo_.planner_objective_value`. It will provide both unconditional
welfare and welfare conditional on the initial (i.e. period 0) values of
the endogenous and exogenous state variables inherited by the planner.
In a deterministic context, the respective initial values are set using
`initval` or `histval` (depending on the exact context).

In a stochastic context, if no initial state values have been specified
with `histval`, their values are taken to be the steady state values.
Because conditional welfare is computed conditional on optimal policy by
the planner in the first endogenous period (period 1), it is conditional
on the information set in the period 1. This information set includes
both the predetermined states inherited from period 0 (specified via
`histval` for both endogenous and lagged exogenous states) as well as
the period 1 values of the exogenous shocks. The latter are specified
using the perfect foresight syntax of the shocks-block.

At the current stage, the stochastic context does not support the
`pruning` option. At `order>3`, only the computation of conditional
welfare with steady state Lagrange multipliers is supported. Note that
at [order=2]{.title-ref}, the output is based on the second-order
accurate approximation of the variance stored in [oo\_.var]{.title-ref}.

*Options*

::: {.option}
periods = INTEGER

The value of the option specifies the number of periods to use in the
simulations in the computation of unconditional welfare at higher order.

Default: `10000`.
:::

::: {.option}
drop = INTEGER

The number of burn-in draws out of `periods` discarded before computing
the unconditional welfare at higher order. Default: `1000`.
:::

*Example (stochastic context)*

>     var a ...;
>     varexo u;
>
>     model;
>     a = rho*a(-1)+u+u(-1);
>     ...
>     end;
>
>     histval;
>     u(0)=1;
>     a(0)=-1;
>     end;
>
>     shocks;
>     var u; stderr 0.008;
>     var u;
>     periods 1;
>     values 1;
>     end;
>
>     evaluate_planner_objective;

::: {.matvar}
[oo]().planner\_objective\_value.unconditional
:::

Scalar storing the value of unconditional welfare. In a perfect
foresight context, it corresponds to welfare in the long-run,
approximated as welfare in the terminal simulation period.

::: {.matvar}
[oo]().planner\_objective\_value.conditional
:::

In a perfect foresight context, this field will be a scalar storing the
value of welfare conditional on the specified initial condition and zero
initial Lagrange multipliers.

In a stochastic context, it will have two subfields:

::: {.matvar}
[oo]().planner\_objective\_value.conditional.steady\_initial\_multiplier
:::

Stores the value of the planner objective when the initial Lagrange
multipliers associated with the planner's problem are set to their
steady state values (see `ramsey_policy`{.interpreted-text
role="comm"}).

::: {.matvar}
[oo]().planner\_objective\_value.conditional.zero\_initial\_multiplier
:::

Stores the value of the planner objective when the initial Lagrange
multipliers associated with the planner's problem are set to 0, i.e. it
is assumed that the planner exploits its ability to surprise private
agents in the first period of implementing Ramsey policy. This value
corresponds to the planner implementing optimal policy for the first
time and committing not to re-optimize in the future.
:::

### Optimal policy under commitment (Ramsey)

Dynare allows to automatically compute optimal policy choices of a
Ramsey planner who takes the specified private sector equilibrium
conditions into account and commits to future policy choices. Doing so
requires specifying the private sector equilibrium conditions in the
`model`-block and a `planner_objective` as well as potentially some
`instruments` to facilitate computations.

::: {.warning}
::: {.admonition-title}
Warning
:::

Be careful when employing forward-looking auxiliary variables in the
context of timeless perspective Ramsey computations. They may alter the
problem the Ramsey planner will solve for the first period, although
they seemingly leave the private sector equilibrium unaffected. The
reason is the planner optimizes with respect to variables dated `t` and
takes the value of time 0 variables as given, because they are
predetermined. This set of initially predetermined variables will change
with forward-looking definitions. Thus, users are strongly advised to
use model-local variables instead.

*Example*

> Consider a perfect foresight example where the Euler equation for the
> return to capital is given by
>
>     1/C=beta*1/C(+1)*(R(+1)+(1-delta))
>
> The job of the Ramsey planner in period `1` is to choose $C_1$ and
> $R_1$, taking as given $C_0$. The above equation may seemingly
> equivalently be written as
>
>     1/C=beta*1/C(+1)*(R_cap);
>     R_cap=R(+1)+(1-delta);
>
> due to perfect foresight. However, this changes the problem of the
> Ramsey planner in the first period to choosing $C_1$ and $R_1$, taking
> as given both $C_0$ and $R^{cap}_0$. Thus, the relevant return to
> capital in the Euler equation of the first period is not a choice of
> the planner anymore due to the forward-looking nature of the
> definition in the second line!
>
> A correct specification would be to instead define `R_cap` as a
> model-local variable:
>
>     1/C=beta*1/C(+1)*(R_cap);
>     #R_cap=R(+1)+(1-delta);
:::

::: {.command}
ramsey\_model (OPTIONS\...);

This command computes the First Order Conditions for maximizing the
policy maker objective function subject to the constraints provided by
the equilibrium path of the private economy.

The planner objective must be declared with the
`planner_objective`{.interpreted-text role="comm"} command.

This command only creates the expanded model, it doesn't perform any
computations. It needs to be followed by other instructions to actually
perform desired computations. Examples are calls to `steady` to compute
the steady state of the Ramsey economy, to `stoch_simul` with various
approximation orders to conduct stochastic simulations based on
perturbation solutions, to `estimation` in order to estimate models
under optimal policy with commitment, and to perfect foresight
simulation routines.

See `aux-variables`{.interpreted-text role="ref"}, for an explanation of
how Lagrange multipliers are automatically created.

*Options*

This command accepts the following options:

::: {.option}
planner\_discount = EXPRESSION

Declares or reassigns the discount factor of the central planner
`optimal_policy_discount_factor`. Default: `1.0`.
:::

::: {.option}
planner\_discount\_latex\_name = LATEX\_NAME

Sets the LaTeX name of the `optimal_policy_discount_factor` parameter.
:::

::: {.option}
instruments = (VARIABLE\_NAME,\...)

Declares instrument variables for the computation of the steady state
under optimal policy. Requires a `steady_state_model` block or a
`_steadystate.m` file. See below.
:::

*Steady state*

Dynare takes advantage of the fact that the Lagrange multipliers appear
linearly in the equations of the steady state of the model under optimal
policy. Nevertheless, it is in general very difficult to compute the
steady state with simply a numerical guess in `initval` for the
endogenous variables.

It greatly facilitates the computation, if the user provides an
analytical solution for the steady state (in `steady_state_model` block
or in a `_steadystate.m` file). In this case, it is necessary to provide
a steady state solution CONDITIONAL on the value of the instruments in
the optimal policy problem and declared with the option `instruments`.
The initial value of the instrument for steady state finding in this
case is set with `initval`. Note that computing and displaying steady
state values using the `steady`-command or calls to `resid` must come
after the `ramsey_model` statement and the `initval`-block.

Note that choosing the instruments is partly a matter of interpretation
and you can choose instruments that are handy from a mathematical point
of view but different from the instruments you would refer to in the
analysis of the paper. A typical example is choosing inflation or
nominal interest rate as an instrument.
:::

::: {.block}
ramsey\_constraints ;

This block lets you define constraints on the variables in the Ramsey
problem. The constraints take the form of a variable, an inequality
operator (\> or \<) and a constant.

*Example*

>     ramsey_constraints;
>     i > 0;
>     end;
:::

::: {.command}
ramsey\_policy \[VARIABLE\_NAME\...\]; ramsey\_policy (OPTIONS\...)
\[VARIABLE\_NAME\...\];

This command is deprecated and formally equivalent to the calling
sequence

>     ramsey_model;
>     stoch_simul;
>     evaluate_planner_objective;

It computes an approximation of the policy that maximizes the policy
maker's objective function subject to the constraints provided by the
equilibrium path of the private economy and under commitment to this
optimal policy. The Ramsey policy is computed by approximating the
equilibrium system around the perturbation point where the Lagrange
multipliers are at their steady state, i.e. where the Ramsey planner
acts as if the initial multipliers had been set to 0 in the distant
past, giving them time to converge to their steady state value.
Consequently, the optimal decision rules are computed around this steady
state of the endogenous variables and the Lagrange multipliers.

Note that the variables in the list after the `ramsey_policy` or
`stoch_simul`-command can also contain multiplier names, but in a
case-sensititve way (e.g. `MULT_1`). In that case, Dynare will for
example display the IRFs of the respective multipliers when `irf>0`.

The planner objective must be declared with the
`planner_objective`{.interpreted-text role="comm"} command.

*Options*

This command accepts all options of `stoch_simul`, plus:

::: {.option}
planner\_discount = EXPRESSION

See `planner_discount <planner_discount = EXPRESSION>`{.interpreted-text
role="opt"}.
:::

::: {.option}
instruments = (VARIABLE\_NAME,\...)

Declares instrument variables for the computation of the steady state
under optimal policy. Requires a `steady_state_model` block or a
`_steadystate.m` file. See below.
:::

*Output*

This command generates all the output variables of `stoch_simul`. For
specifying the initial values for the endogenous state variables (except
for the Lagrange multipliers), see above.

*Steady state*

See `Ramsey steady state <ramsey_model>`{.interpreted-text role="comm"}.
:::

### Optimal policy under discretion

::: {.command}
discretionary\_policy \[VARIABLE\_NAME\...\]; discretionary\_policy
(OPTIONS\...) \[VARIABLE\_NAME\...\];

This command computes an approximation of the optimal policy under
discretion. The algorithm implemented is essentially an LQ solver, and
is described by *Dennis (2007)*.

You must ensure that your objective is quadratic. Regarding the model,
it must either be linear or solved at first order with an analytical
steady state provided. In the first case, you should set the `linear`
option of the `model` block.

It is possible to use the `estimation`{.interpreted-text role="comm"}
command after the `discretionary_policy` command, in order to estimate
the model with optimal policy under discretion and
`evaluate_planner_objective`{.interpreted-text role="comm"} to compute
welfare.

*Options*

This command accepts the same options as `ramsey_policy`, plus:

::: {.option}
discretionary\_tol = NON-NEGATIVE DOUBLE

Sets the tolerance level used to assess convergence of the solution
algorithm. Default: `1e-7`.
:::

::: {.option}
maxit = INTEGER

Maximum number of iterations. Default: `3000`.
:::
:::

### Optimal Simple Rules (OSR)

::: {.command}
osr \[VARIABLE\_NAME\...\]; osr (OPTIONS\...) \[VARIABLE\_NAME\...\];

This command computes optimal simple policy rules for linear-quadratic
problems of the form:

> $$\min_\gamma E(y'_tWy_t)$$

such that:

> $$A_1 E_ty_{t+1}+A_2 y_t+ A_3 y_{t-1}+C e_t=0$$

where:

> -   $E$ denotes the unconditional expectations operator;
> -   $\gamma$ are parameters to be optimized. They must be elements of
>     the matrices $A_1$, $A_2$, $A_3$, i.e. be specified as parameters
>     in the `params` command and be entered in the `model` block;
> -   $y$ are the endogenous variables, specified in the `var` command,
>     whose (co)-variance enters the loss function;
> -   $e$ are the exogenous stochastic shocks, specified in the
>     `varexo`- ommand;
> -   $W$ is the weighting matrix;

The linear quadratic problem consists of choosing a subset of model
parameters to minimize the weighted (co)-variance of a specified subset
of endogenous variables, subject to a linear law of motion implied by
the first order conditions of the model. A few things are worth
mentioning. First, $y$ denotes the selected endogenous variables'
deviations from their steady state, i.e. in case they are not already
mean 0 the variables entering the loss function are automatically
demeaned so that the centered second moments are minimized. Second,
`osr` only solves linear quadratic problems of the type resulting from
combining the specified quadratic loss function with a first order
approximation to the model's equilibrium conditions. The reason is that
the first order state-space representation is used to compute the
unconditional (co)-variances. Hence, `osr` will automatically select
`order=1`. Third, because the objective involves minimizing a weighted
sum of unconditional second moments, those second moments must be
finite. In particular, unit roots in $y$ are not allowed.

The subset of the model parameters over which the optimal simple rule is
to be optimized, $\gamma$, must be listed with `osr_params`.

The weighting matrix $W$ used for the quadratic objective function is
specified in the `optim_weights` block. By attaching weights to
endogenous variables, the subset of endogenous variables entering the
objective function, $y$, is implicitly specified.

The linear quadratic problem is solved using the numerical optimizer
specified with `opt_algo <opt_algo = INTEGER>`{.interpreted-text
role="opt"}.

*Options*

The `osr` command will subsequently run `stoch_simul` and accepts the
same options, including restricting the endogenous variables by listing
them after the command, as `stoch_simul` (see
`stoch-sol`{.interpreted-text role="ref"}) plus

::: {.option}
opt\_algo = INTEGER

Specifies the optimizer for minimizing the objective function. The same
solvers as for `mode_compute` (see
`mode_compute <mode_compute = INTEGER | FUNCTION_NAME>`{.interpreted-text
role="opt"}) are available, except for `5`, `6`, and `10`.
:::

::: {.option}
optim = (NAME, VALUE, \...)

A list of NAME\`[ and VALUE pairs. Can be used to set options for the
optimization routines. The set of available options depends on the
selected optimization routine (i.e. on the value of option
:opt:\`opt\_algo \<opt\_algo = INTEGER\>]{.title-ref}). See
`optim <optim = (NAME, VALUE, ...)>`{.interpreted-text role="opt"}.
:::

::: {.option}
maxit = INTEGER

Determines the maximum number of iterations used in `opt_algo=4`. This
option is now deprecated and will be removed in a future release of
Dynare. Use `optim` instead to set optimizer-specific values. Default:
`1000`.
:::

::: {.option}
tolf = DOUBLE

Convergence criterion for termination based on the function value used
in `opt_algo=4`. Iteration will cease when it proves impossible to
improve the function value by more than tolf. This option is now
deprecated and will be removed in a future release of Dynare. Use
`optim` instead to set optimizer-specific values. Default: `1e-7`.
:::

::: {.option}
silent\_optimizer

See `silent_optimizer`{.interpreted-text role="opt"}.
:::

::: {.option}
huge\_number = DOUBLE

Value for replacing the infinite bounds on parameters by finite numbers.
Used by some optimizers for numerical reasons (see
`huge_number <huge_number = DOUBLE>`{.interpreted-text role="opt"}).
Users need to make sure that the optimal parameters are not larger than
this value. Default: `1e7`.
:::

The value of the objective is stored in the variable
`oo_.osr.objective_function` and the value of parameters at the optimum
is stored in `oo_.osr.optim_params`. See below for more details.

After running `osr` the parameters entering the simple rule will be set
to their optimal value so that subsequent runs of `stoch_simul` will be
conducted at these values.
:::

::: {.command}
osr\_params PARAMETER\_NAME\...;

This command declares parameters to be optimized by `osr`.
:::

::: {.block}
optim\_weights ;

This block specifies quadratic objectives for optimal policy problems.

More precisely, this block specifies the nonzero elements of the weight
matrix $W$ used in the quadratic form of the objective function in
`osr`.

An element of the diagonal of the weight matrix is given by a line of
the form:

    VARIABLE_NAME EXPRESSION;

An off-the-diagonal element of the weight matrix is given by a line of
the form:

    VARIABLE_NAME,  VARIABLE_NAME EXPRESSION;
:::

*Example*

>     var y inflation r;
>     varexo y_ inf_;
>
>     parameters delta sigma alpha kappa gammarr gammax0 gammac0 gamma_y_ gamma_inf_;
>
>     delta =  0.44;
>     kappa =  0.18;
>     alpha =  0.48;
>     sigma = -0.06;
>
>     gammarr = 0;
>     gammax0 = 0.2;
>     gammac0 = 1.5;
>     gamma_y_ = 8;
>     gamma_inf_ = 3;
>
>     model(linear);
>     y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - inflation(+1)) + y_;
>     inflation  =   alpha * inflation(-1) + (1-alpha) * inflation(+1) + kappa*y + inf_;
>     r = gammax0*y(-1)+gammac0*inflation(-1)+gamma_y_*y_+gamma_inf_*inf_;
>     end;
>
>     shocks;
>     var y_; stderr 0.63;
>     var inf_; stderr 0.4;
>     end;
>
>     optim_weights;
>     inflation 1;
>     y 1;
>     y, inflation 0.5;
>     end;
>
>     osr_params gammax0 gammac0 gamma_y_ gamma_inf_;
>     osr y;

::: {.block}
osr\_params\_bounds ;

This block declares lower and upper bounds for parameters in the optimal
simple rule. If not specified the optimization is unconstrained.

Each line has the following syntax:

    PARAMETER_NAME, LOWER_BOUND, UPPER_BOUND;

Note that the use of this block requires the use of a constrained
optimizer, i.e. setting
`opt_algo <opt_algo = INTEGER>`{.interpreted-text role="opt"} to `1`,
`2`, `5` or `9`.

*Example*

>     osr_params_bounds;
>     gamma_inf_, 0, 2.5;
>     end;
>
>     osr(opt_algo=9) y;
:::

::: {.matvar}
[oo]().osr.objective\_function

After an execution of the `osr` command, this variable contains the
value of the objective under optimal policy.
:::

::: {.matvar}
[oo]().osr.optim\_params

After an execution of the `osr` command, this variable contains the
value of parameters at the optimum, stored in fields of the form
`oo_.osr.optim_params.PARAMETER_NAME`.
:::

::: {.matvar}
[M]().osr.param\_names

After an execution of the `osr` command, this cell contains the names of
the parameters.
:::

::: {.matvar}
[M]().osr.param\_indices

After an execution of the `osr` command, this vector contains the
indices of the OSR parameters in `M_.params`.
:::

::: {.matvar}
[M]().osr.param\_bounds

After an execution of the `osr` command, this two by number of OSR
parameters matrix contains the lower and upper bounds of the parameters
in the first and second column, respectively.
:::

::: {.matvar}
[M]().osr.variable\_weights

After an execution of the `osr` command, this sparse matrix contains the
weighting matrix associated with the variables in the objective
function.
:::

::: {.matvar}
[M]().osr.variable\_indices

After an execution of the `osr` command, this vector contains the
indices of the variables entering the objective function in
`M_.endo_names`.
:::

Sensitivity and identification analysis
---------------------------------------

Dynare provides an interface to the global sensitivity analysis (GSA)
toolbox (developed by the Joint Research Center (JRC) of the European
Commission), which is now part of the official Dynare distribution. The
GSA toolbox can be used to answer the following questions:

> 1.  What is the domain of structural coefficients assuring the
>     stability and determinacy of a DSGE model?
> 2.  Which parameters mostly drive the fit of, e.g., GDP and which the
>     fit of inflation? Is there any conflict between the optimal fit of
>     one observed series versus another?
> 3.  How to represent in a direct, albeit approximated, form the
>     relationship between structural parameters and the reduced form of
>     a rational expectations model?

The discussion of the methodologies and their application is described
in *Ratto (2008)*.

With respect to the previous version of the toolbox, in order to work
properly, the GSA toolbox no longer requires that the Dynare estimation
environment is set up.

### Performing sensitivity analysis

::: {.command}
dynare\_sensitivity ; dynare\_sensitivity(OPTIONS\...);

This command triggers sensitivity analysis on a DSGE model.

::: {#sampl-opt}
*Sampling Options*
:::

::: {.option}
Nsam = INTEGER

Size of the Monte-Carlo sample. Default: `2048`.
:::

::: {.option}
ilptau = INTEGER

If equal to `1`, use $LP_\tau$ quasi-Monte-Carlo. If equal to `0`, use
LHS Monte-Carlo. Default: `1`.
:::

::: {.option}
pprior = INTEGER

If equqal to `1`, sample from the prior distributions. If equal to `0`,
sample from the multivariate normal $N(\bar{\theta},\Sigma)$, where
$\bar{\theta}$ is the posterior mode and $\Sigma=H^{-1}$, $H$ is the
Hessian at the mode. Default: `1`.
:::

::: {.option}
prior\_range = INTEGER

If equal to `1`, sample uniformly from prior ranges. If equal to `0`,
sample from prior distributions. Default: `1`.
:::

::: {.option}
morris = INTEGER

If equal to `0`, ANOVA mapping (Type I error) If equal to `1`, Screening
analysis (Type II error). If equal to `2`, Analytic derivatives (similar
to Type II error, only valid when identification=1). Default: `1` when
`identification=1`, `0` otherwise.
:::

::: {.option}
morris\_nliv = INTEGER

Number of levels in Morris design. Default: `6`.
:::

::: {.option}
morris\_ntra = INTEGER

Number trajectories in Morris design. Default: `20`.
:::

::: {.option}
ppost = INTEGER

If equal to `1`, use Metropolis posterior sample. If equal to `0`, do
not use Metropolis posterior sample. Default: `0`.

NB: This overrides any other sampling option.
:::

::: {.option}
neighborhood\_width = DOUBLE

When `pprior=0` and `ppost=0`, allows for the sampling of parameters
around the value specified in the `mode_file`, in the range
$\texttt{xparam1} \pm \left \vert
\texttt{xparam1} \times \texttt{neighborhood\_width} \right
\vert$. Default: `0`.
:::

*Stability Mapping Options*

::: {.option}
stab = INTEGER

If equal to `1`, perform stability mapping. If equal to `0`, do not
perform stability mapping. Default: `1`.
:::

::: {.option}
load\_stab = INTEGER

If equal to `1`, load a previously created sample. If equal to `0`,
generate a new sample. Default: `0`.
:::

::: {.option}
alpha2\_stab = DOUBLE

Critical value for correlations $\rho$ in filtered samples: plot couples
of parmaters with $\left\vert\rho\right\vert>$ `alpha2_stab`. Default:
`0`.
:::

::: {.option}
pvalue\_ks = DOUBLE

The threshold $pvalue$ for significant Kolmogorov-Smirnov test (i.e.
plot parameters with $pvalue<$ `pvalue_ks`). Default: `0.001`.
:::

::: {.option}
pvalue\_corr = DOUBLE

The threshold $pvalue$ for significant correlation in filtered samples
(i.e. plot bivariate samples when $pvalue<$ `pvalue_corr`). Default:
`1e-5`.
:::

*Reduced Form Mapping Options*

::: {.option}
redform = INTEGER

If equal to `1`, prepare Monte-Carlo sample of reduced form matrices. If
equal to `0`, do not prepare Monte-Carlo sample of reduced form
matrices. Default: `0`.
:::

::: {.option}
load\_redform = INTEGER

If equal to `1`, load previously estimated mapping. If equal to `0`,
estimate the mapping of the reduced form model. Default: `0`.
:::

::: {.option}
logtrans\_redform = INTEGER

If equal to `1`, use log-transformed entries. If equal to `0`, use raw
entries. Default: `0`.
:::

::: {.option}
threshold\_redform = \[DOUBLE DOUBLE\]

The range over which the filtered Monte-Carlo entries of the reduced
form coefficients should be analyzed. The first number is the lower
bound and the second is the upper bound. An empty vector indicates that
these entries will not be filtered. Default: empty.
:::

::: {.option}
ksstat\_redform = DOUBLE

Critical value for Smirnov statistics $d$ when reduced form entries are
filtered. Default: `0.001`.
:::

::: {.option}
alpha2\_redform = DOUBLE

Critical value for correlations $\rho$ when reduced form entries are
filtered. Default: `1e-5`.
:::

::: {.option}
namendo = (VARIABLE\_NAME\...)

List of endogenous variables. ':' indicates all endogenous variables.
Default: empty.
:::

::: {.option}
namlagendo = (VARIABLE\_NAME\...)

List of lagged endogenous variables. ':' indicates all lagged endogenous
variables. Analyze entries \[namendo $\times$ namlagendo\] Default:
empty.
:::

::: {.option}
namexo = (VARIABLE\_NAME\...)

List of exogenous variables. ':' indicates all exogenous variables.
Analyze entries \[namendo $\times$ namexo\]. Default: empty.
:::

*RMSE Options*

::: {.option}
rmse = INTEGER

If equal to `1`, perform RMSE analysis. If equal to `0`, do not perform
RMSE analysis. Default: `0`.
:::

::: {.option}
load\_rmse = INTEGER

If equal to `1`, load previous RMSE analysis. If equal to `0`, make a
new RMSE analysis. Default: `0`.
:::

::: {.option}
lik\_only = INTEGER

If equal to `1`, compute only likelihood and posterior. If equal to `0`,
compute RMSE's for all observed series. Default: `0`.
:::

::: {.option}
var\_rmse = (VARIABLE\_NAME\...)

List of observed series to be considered. ':' indicates all observed
variables. Default: `varobs`.
:::

::: {.option}
pfilt\_rmse = DOUBLE

Filtering threshold for RMSE's. Default: `0.1`.
:::

::: {.option}
istart\_rmse = INTEGER

Value at which to start computing RMSE's (use `2` to avoid big intitial
error). Default: `presample+1`.
:::

::: {.option}
alpha\_rmse = DOUBLE

Critical value for Smirnov statistics $d$: plot parameters with $d>$
`alpha_rmse`. Default: `0.001`.
:::

::: {.option}
alpha2\_rmse = DOUBLE

Critical value for correlation $\rho$: plot couples of parmaters with
$\left\vert\rho\right\vert=$ `alpha2_rmse`. Default: `1e-5`.
:::

::: {.option}
datafile = FILENAME

See `datafile <dataf>`{.interpreted-text role="ref"}.
:::

::: {.option}
nobs = INTEGER nobs = \[INTEGER1:INTEGER2\]

See `nobs <nobs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
first\_obs = INTEGER

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
prefilter = INTEGER

See `prefilter <prefilter = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
presample = INTEGER

See `presample <presample = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
nograph

See `nograph`{.interpreted-text role="opt"}.
:::

::: {.option}
nodisplay

See `nodisplay`{.interpreted-text role="opt"}.
:::

::: {.option}
graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )

See `graph_format <graph_format = FORMAT>`{.interpreted-text
role="opt"}.
:::

::: {.option}
conf\_sig = DOUBLE

See `conf_sig <confsig>`{.interpreted-text role="ref"}.
:::

::: {.option}
loglinear

See `loglinear <logl>`{.interpreted-text role="ref"}.
:::

::: {.option}
mode\_file = FILENAME

See `mode_file <mode_file = FILENAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
kalman\_algo = INTEGER

See `kalman_algo <kalman_algo = INTEGER>`{.interpreted-text role="opt"}.
:::

*Identification Analysis Options*

::: {.option}
identification = INTEGER

If equal to `1`, performs identification analysis (forcing `redform=0`
and `morris=1`) If equal to `0`, no identification analysis. Default:
`0`.
:::

::: {.option}
morris = INTEGER

See `morris <morris = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
morris\_nliv = INTEGER

See `morris_nliv <morris_nliv = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
morris\_ntra = INTEGER

See `morris_ntra <morris_ntra = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
load\_ident\_files = INTEGER

Loads previously performed identification analysis. Default: `0`.
:::

::: {.option}
useautocorr = INTEGER

Use autocorrelation matrices in place of autocovariance matrices in
moments for identification analysis. Default: `0`.
:::

::: {.option}
ar = INTEGER

Maximum number of lags for moments in identification analysis. Default:
`1`.
:::

::: {.option}
diffuse\_filter = INTEGER

See `diffuse_filter`{.interpreted-text role="opt"}.
:::
:::

### IRF/Moment calibration {#irf-momcal}

The `irf_calibration` and `moment_calibration` blocks allow imposing
implicit "endogenous" priors about IRFs and moments on the model. The
way it works internally is that any parameter draw that is inconsistent
with the "calibration" provided in these blocks is discarded, i.e.
assigned a prior density of `0`. In the context of `dynare_sensitivity`,
these restrictions allow tracing out which parameters are driving the
model to satisfy or violate the given restrictions.

IRF and moment calibration can be defined in `irf_calibration` and
`moment_calibration` blocks:

::: {.block}
irf\_calibration ; irf\_calibration (OPTIONS\...);

This block allows defining IRF calibration criteria and is terminated by
`end;`. To set IRF sign restrictions, the following syntax is used:

    VARIABLE_NAME(INTEGER), EXOGENOUS_NAME, -;
    VARIABLE_NAME(INTEGER:INTEGER), EXOGENOUS_NAME, +;

To set IRF restrictions with specific intervals, the following syntax is
used:

    VARIABLE_NAME(INTEGER), EXOGENOUS_NAME, [EXPRESSION, EXPRESSION];
    VARIABLE_NAME(INTEGER:INTEGER), EXOGENOUS_NAME, [EXPRESSION, EXPRESSION];

When `(INTEGER:INTEGER)` is used, the restriction is considered to be
fulfilled by a logical OR. A list of restrictions must always be
fulfilled with logical AND.

*Options*

::: {.option}
relative\_irf

See `relative_irf`{.interpreted-text role="opt"}.
:::

*Example*

>     irf_calibration;
>     y(1:4), e_ys, [-50, 50]; //[first year response with logical OR]
>     @#for ilag in 21:40
>     R_obs(@{ilag}), e_ys, [0, 6]; //[response from 5th to 10th years with logical AND]
>     @#endfor
>     end;
:::

::: {.block}
moment\_calibration ; moment\_calibration (OPTIONS\...);

This block allows defining moment calibration criteria. This block is
terminated by `end;`, and contains lines of the form:

    VARIABLE_NAME1, VARIABLE_NAME2(+/-INTEGER), [EXPRESSION, EXPRESSION];
    VARIABLE_NAME1, VARIABLE_NAME2(+/-INTEGER), +/-;
    VARIABLE_NAME1, VARIABLE_NAME2(+/-(INTEGER:INTEGER)), [EXPRESSION, EXPRESSION];
    VARIABLE_NAME1, VARIABLE_NAME2((-INTEGER:+INTEGER)), [EXPRESSION, EXPRESSION];

When `(INTEGER:INTEGER)` is used, the restriction is considered to be
fulfilled by a logical OR. A list of restrictions must always be
fulfilled with logical AND. The moment restrictions generally apply to
auto- and cross-correlations between variables. The only exception is a
restriction on the unconditional variance of an endogenous variable,
specified as shown in the example below.

*Example*

>     moment_calibration;
>     y_obs,y_obs, [0.5, 1.5]; //[unconditional variance]
>     y_obs,y_obs(-(1:4)), +; //[sign restriction for first year autocorrelation with logical OR]
>     @#for ilag in -2:2
>     y_obs,R_obs(@{ilag}), -; //[-2:2 cross correlation with logical AND]
>     @#endfor
>     @#for ilag in -4:4
>     y_obs,pie_obs(@{ilag}), -; //[-4_4 cross correlation with logical AND]
>     @#endfor
>     end;
:::

### Performing identification analysis {#identification-analysis}

::: {.command}
identification ; identification (OPTIONS\...);

This command triggers:

> 1.  Theoretical identification analysis based on
>
>     -   moments as in *Iskrev (2010)*
>     -   spectral density as in *Qu and Tkachenko (2012)*
>     -   minimal system as in *Komunjer and Ng (2011)*
>     -   reduced-form solution and linear rational expectation model as
>         in *Ratto and Iskrev (2011)*
>
>     Note that for orders 2 and 3, all identification checks are based
>     on the pruned state space system as in *Mutschler (2015)*. That
>     is, theoretical moments and spectrum are computed from the pruned
>     ABCD-system, whereas the minimal system criteria is based on the
>     first-order system, but augmented by the theoretical (pruned) mean
>     at order 2 or 3.
>
> 2.  Identification strength analysis based on (theoretical or
>     simulated) curvature of moment information matrix as in *Ratto and
>     Iskrev (2011)*
> 3.  Parameter checks based on nullspace and multicorrelation
>     coefficients to determine which (combinations of) parameters are
>     involved
:::

*General Options*

> ::: {.option}
> order = 13
>
> Order of approximation. At orders 2 and 3 identification is based on
> the pruned state space system. Note that the order set in other
> functions does not overwrite the default. Default: `1`.
> :::
>
> ::: {.option}
> parameter\_set = OPTION
>
> See `parameter_set <parameter_set = OPTION>`{.interpreted-text
> role="opt"} for possible values. Default: `prior_mean`.
> :::
>
> ::: {.option}
> prior\_mc = INTEGER
>
> Size of Monte-Carlo sample. Default: `1`.
> :::
>
> ::: {.option}
> prior\_range = INTEGER
>
> Triggers uniform sample within the range implied by the prior
> specifications (when `prior_mc>1`). Default: `0`.
> :::
>
> ::: {.option}
> advanced = INTEGER
>
> If set to `1`, shows a more detailed analysis, comprised of an
> analysis for the linearized rational expectation model as well as the
> associated reduced form solution. Further performs a bruteforce search
> of the groups of parameters best reproducing the behavior of each
> single parameter. The maximum dimension of the group searched is
> triggered by `max_dim_cova_group`. Default: `0`.
> :::
>
> ::: {.option}
> max\_dim\_cova\_group = INTEGER
>
> In the brute force search (performed when `advanced=1`) this option
> sets the maximum dimension of groups of parameters that best reproduce
> the behavior of each single model parameter. Default: `2`.
> :::
>
> ::: {.option}
> gsa\_sample\_file = INTEGER\|FILENAME
>
> If equal to `0`, do not use sample file. If equal to `1`, triggers gsa
> prior sample. If equal to `2`, triggers gsa Monte-Carlo sample (i.e.
> loads a sample corresponding to `pprior=0` and `ppost=0` in the
> `dynare_sensitivity` options). If equal to `FILENAME` uses the
> provided path to a specific user defined sample file. Default: `0`.
> :::
>
> ::: {.option}
> diffuse\_filter
>
> Deals with non-stationary cases. See
> `diffuse_filter`{.interpreted-text role="opt"}.
> :::

*Numerical Options*

> ::: {.option}
> analytic\_derivation\_mode = INTEGER
>
> Different ways to compute derivatives either analytically or
> numerically. Possible values are:
>
> > -   `0`: efficient sylvester equation method to compute analytical
> >     derivatives
> > -   `1`: kronecker products method to compute analytical derivatives
> >     (only at order=1)
> > -   `-1`: numerical two-sided finite difference method to compute
> >     all identification Jacobians (numerical tolerance level is equal
> >     to `options_.dynatol.x`)
> > -   `-2`: numerical two-sided finite difference method to compute
> >     derivatives of steady state and dynamic model numerically, the
> >     identification Jacobians are then computed analytically
> >     (numerical tolerance level is equal to `options_.dynatol.x`)
>
> Default: `0`.
> :::
>
> ::: {.option}
> normalize\_jacobians = INTEGER
>
> If set to `1`: Normalize Jacobian matrices by rescaling each row by
> its largest element in absolute value. Normalize Gram (or
> Hessian-type) matrices by transforming into correlation-type matrices.
> Default: `1`
> :::
>
> ::: {.option}
> tol\_rank = DOUBLE
>
> Tolerance level used for rank computations. Default: `1.e-10`.
> :::
>
> ::: {.option}
> tol\_deriv = DOUBLE
>
> Tolerance level for selecting non-zero columns in Jacobians. Default:
> `1.e-8`.
> :::
>
> ::: {.option}
> tol\_sv = DOUBLE
>
> Tolerance level for selecting non-zero singular values. Default:
> `1.e-3`.
> :::
>
> ::: {.option}
> schur\_vec\_tol = DOUBLE
>
> See `schur_vec_tol <schur_vec_tol = DOUBLE>`{.interpreted-text
> role="opt"}.
> :::

*Identification Strength Options*

> ::: {.option}
> no\_identification\_strength
>
> Disables computations of identification strength analysis based on
> sample information matrix.
> :::
>
> ::: {.option}
> periods = INTEGER
>
> When the analytic Hessian is not available (i.e. with missing values
> or diffuse Kalman filter or univariate Kalman filter), this triggers
> the length of stochastic simulation to compute Simulated Moments
> Uncertainty. Default: `300`.
> :::
>
> ::: {.option}
> replic = INTEGER
>
> When the analytic Hessian is not available, this triggers the number
> of replicas to compute Simulated Moments Uncertainty. Default: `100`.
> :::

*Moments Options*

> ::: {.option}
> no\_identification\_moments
>
> Disables computations of identification check based on Iskrev
> (2010)\'s J, i.e. derivative of first two moments.
> :::
>
> ::: {.option}
> ar = INTEGER
>
> Number of lags of computed autocovariances/autocorrelations
> (theoretical moments) in Iskrev (2010)\'s J criteria. Default: `1`.
> :::
>
> ::: {.option}
> useautocorr = INTEGER
>
> If equal to `1`, compute derivatives of autocorrelation. If equal to
> `0`, compute derivatives of autocovariances. Default: `0`.
> :::

*Spectrum Options*

> ::: {.option}
> no\_identification\_spectrum
>
> Disables computations of identification check based on *Qu and
> Tkachenko (2012)*\'s G, i.e. Gram matrix of derivatives of first
> moment plus outer product of derivatives of spectral density.
> :::
>
> ::: {.option}
> grid\_nbr = INTEGER
>
> Number of grid points in \[-pi;pi\] to approximate the integral to
> compute Qu and Tkachenko (2012)\'s G criteria. Default: `5000`.
> :::

*Minimal State Space System Options*

> ::: {.option}
> no\_identification\_minimal
>
> Disables computations of identification check based on *Komunjer and
> Ng (2011)*\'s D, i.e. minimal state space system and observational
> equivalent spectral density transformations.
> :::

*Misc Options*

> ::: {.option}
> nograph
>
> See `nograph`{.interpreted-text role="opt"}.
> :::
>
> ::: {.option}
> nodisplay
>
> See `nodisplay`{.interpreted-text role="opt"}.
> :::
>
> ::: {.option}
> graph\_format = FORMAT graph\_format = ( FORMAT, FORMAT\... )
>
> See `graph_format <graph_format = FORMAT>`{.interpreted-text
> role="opt"}.
> :::
>
> ::: {.option}
> tex
>
> See `tex`{.interpreted-text role="opt"}.
> :::

*Debug Options*

> ::: {.option}
> load\_ident\_files = INTEGER
>
> If equal to `1`, allow Dynare to load previously computed analyzes.
> Default: `0`.
> :::
>
> ::: {.option}
> lik\_init = INTEGER
>
> See `lik_init <lik_init = INTEGER>`{.interpreted-text role="opt"}.
> :::
>
> ::: {.option}
> kalman\_algo = INTEGER
>
> See `kalman_algo <kalman_algo = INTEGER>`{.interpreted-text
> role="opt"}.
> :::
>
> ::: {.option}
> no\_identification\_reducedform
>
> Disables computations of identification check based on steady state
> and reduced-form solution.
> :::
>
> ::: {.option}
> checks\_via\_subsets = INTEGER
>
> If equal to `1`: finds problematic parameters in a bruteforce fashion:
> It computes the rank of the Jacobians for all possible parameter
> combinations. If the rank condition is not fullfilled, these parameter
> sets are flagged as non-identifiable. The maximum dimension of the
> group searched is triggered by `max_dim_subsets_groups`. Default: `0`.
> :::
>
> ::: {.option}
> max\_dim\_subsets\_groups = INTEGER
>
> Sets the maximum dimension of groups of parameters for which the above
> bruteforce search is performed. Default: `4`.
> :::

### Types of analysis and output files

The sensitivity analysis toolbox includes several types of analyses.
Sensitivity analysis results are saved locally in `<mod_file>/gsa`,
where `<mod_file>.mod` is the name of the Dynare model file.

#### Sampling

The following binary files are produced:

> -   `<mod_file>_prior.mat`: this file stores information about the
>     analyses performed sampling from the prior, i.e. `pprior=1` and
>     `ppost=0`;
> -   `<mod_file>_mc.mat`: this file stores information about the
>     analyses performed sampling from multivariate normal, i.e.
>     `pprior=0` and `ppost=0`;
> -   `<mod_file>_post.mat`: this file stores information about analyses
>     performed using the Metropolis posterior sample, i.e. `ppost=1`.

#### Stability Mapping

Figure files produced are of the form `<mod_file>_prior_*.fig` and store
results for stability mapping from prior Monte-Carlo samples:

> -   `<mod_file>_prior_stable.fig`: plots of the Smirnov test and the
>     correlation analyses confronting the cdf of the sample fulfilling
>     Blanchard-Kahn conditions (blue color) with the cdf of the rest of
>     the sample (red color), i.e. either instability or indeterminacy
>     or the solution could not be found (e.g. the steady state solution
>     could not be found by the solver);
> -   `<mod_file>_prior_indeterm.fig`: plots of the Smirnov test and the
>     correlation analyses confronting the cdf of the sample producing
>     indeterminacy (red color) with the cdf of the rest of the sample
>     (blue color);
> -   `<mod_file>_prior_unstable.fig`: plots of the Smirnov test and the
>     correlation analyses confronting the cdf of the sample producing
>     explosive roots (red color) with the cdf of the rest of the sample
>     (blue color);
> -   `<mod_file>_prior_wrong.fig`: plots of the Smirnov test and the
>     correlation analyses confronting the cdf of the sample where the
>     solution could not be found (e.g. the steady state solution could
>     not be found by the solver - red color) with the cdf of the rest
>     of the sample (blue color);
> -   `<mod_file>_prior_calib.fig`: plots of the Smirnov test and the
>     correlation analyses splitting the sample fulfilling
>     Blanchard-Kahn conditions, by confronting the cdf of the sample
>     where IRF/moment restrictions are matched (blue color) with the
>     cdf where IRF/moment restrictions are NOT matched (red color);

Similar conventions apply for `<mod_file>_mc_*.fig` files, obtained when
samples from multivariate normal are used.

#### IRF/Moment restrictions

The following binary files are produced:

> -   `<mod_file>_prior_restrictions.mat`: this file stores information
>     about the IRF/moment restriction analysis performed sampling from
>     the prior ranges, i.e. `pprior=1` and `ppost=0`;
> -   `<mod_file>_mc_restrictions.mat`: this file stores information
>     about the IRF/moment restriction analysis performed sampling from
>     multivariate normal, i.e. `pprior=0` and `ppost=0`;
> -   `<mod_file>_post_restrictions.mat`: this file stores information
>     about IRF/moment restriction analysis performed using the
>     Metropolis posterior sample, i.e. `ppost=1`.

Figure files produced are of the form `<mod_file>_prior_irf_calib_*.fig`
and `<mod_file>_prior_moment_calib_*.fig` and store results for mapping
restrictions from prior Monte-Carlo samples:

> -   `<mod_file>_prior_irf_calib_<ENDO_NAME>_vs_<EXO_NAME>_<PERIOD>.fig`:
>     plots of the Smirnov test and the correlation analyses splitting
>     the sample fulfilling Blanchard-Kahn conditions, by confronting
>     the cdf of the sample where the individual IRF restriction
>     `<ENDO_NAME>` vs. `<EXO_NAME>` at period(s) `<PERIOD>` is matched
>     (blue color) with the cdf where the IRF restriction is NOT matched
>     (red color)
> -   `<mod_file>_prior_irf_calib_<ENDO_NAME>_vs_<EXO_NAME>_ALL.fig`:
>     plots of the Smirnov test and the correlation analyses splitting
>     the sample fulfilling Blanchard-Kahn conditions, by confronting
>     the cdf of the sample where ALL the individual IRF restrictions
>     for the same couple `<ENDO_NAME>` vs. `<EXO_NAME>` are matched
>     (blue color) with the cdf where the IRF restriction is NOT matched
>     (red color)
> -   `<mod_file>_prior_irf_restrictions.fig`: plots visual information
>     on the IRF restrictions compared to the actual Monte Carlo
>     realization from prior sample.
> -   `<mod_file>_prior_moment_calib_<ENDO_NAME1>_vs_<ENDO_NAME2>_<LAG>.fig`:
>     plots of the Smirnov test and the correlation analyses splitting
>     the sample fulfilling Blanchard-Kahn conditions, by confronting
>     the cdf of the sample where the individual acf/ccf moment
>     restriction `<ENDO_NAME1>` vs. `<ENDO_NAME2>` at lag(s) `<LAG>` is
>     matched (blue color) with the cdf where the IRF restriction is NOT
>     matched (red color)
> -   `<mod_file>_prior_moment_calib_<ENDO_NAME>_vs_<EXO_NAME>_ALL.fig`:
>     plots of the Smirnov test and the correlation analyses splitting
>     the sample fulfilling Blanchard-Kahn conditions, by confronting
>     the cdf of the sample where ALL the individual acf/ccf moment
>     restrictions for the same couple `<ENDO_NAME1>` vs. `<ENDO_NAME2>`
>     are matched (blue color) with the cdf where the IRF restriction is
>     NOT matched (red color)
> -   `<mod_file>_prior_moment_restrictions.fig`: plots visual
>     information on the moment restrictions compared to the actual
>     Monte Carlo realization from prior sample.

Similar conventions apply for `<mod_file>_mc_*.fig` and
`<mod_file>_post_*.fig` files, obtained when samples from multivariate
normal or from posterior are used.

#### Reduced Form Mapping

When the option `threshold_redform` is not set, or it is empty (the
default), this analysis estimates a multivariate smoothing spline ANOVA
model (the 'mapping') for the selected entries in the transition matrix
of the shock matrix of the reduce form first order solution of the
model. This mapping is done either with prior samples or with MC samples
with `neighborhood_width`. Unless `neighborhood_width` is set with MC
samples, the mapping of the reduced form solution forces the use of
samples from prior ranges or prior distributions, i.e.: `pprior=1` and
`ppost=0`. It uses 250 samples to optimize smoothing parameters and 1000
samples to compute the fit. The rest of the sample is used for
out-of-sample validation. One can also load a previously estimated
mapping with a new Monte-Carlo sample, to look at the forecast for the
new Monte-Carlo sample.

The following synthetic figures are produced:

> -   `<mod_file>_redform_<endo name>_vs_lags_*.fig`: shows bar charts
>     of the sensitivity indices for the ten most important parameters
>     driving the reduced form coefficients of the selected endogenous
>     variables (`namendo`) versus lagged endogenous variables
>     (`namlagendo`); suffix `log` indicates the results for
>     log-transformed entries;
> -   `<mod_file>_redform_<endo name>_vs_shocks_*.fig`: shows bar charts
>     of the sensitivity indices for the ten most important parameters
>     driving the reduced form coefficients of the selected endogenous
>     variables (`namendo`) versus exogenous variables (`namexo`);
>     suffix `log` indicates the results for log-transformed entries;
> -   `<mod_file>_redform_gsa(_log).fig`: shows bar chart of all
>     sensitivity indices for each parameter: this allows one to notice
>     parameters that have a minor effect for any of the reduced form
>     coefficients.

Detailed results of the analyses are shown in the subfolder
`<mod_file>/gsa/redform_prior` for prior samples and in
`<mod_file>/gsa/redform_mc` for MC samples with option
`neighborhood_width`, where the detailed results of the estimation of
the single functional relationships between parameters $\theta$ and
reduced form coefficient (denoted as $y$ hereafter) are stored in
separate directories named as:

> -   `<namendo>_vs_<namlagendo>`, for the entries of the transition
>     matrix;
> -   `<namendo>_vs_<namexo>`, for entries of the matrix of the shocks.

The following files are stored in each directory (we stick with prior
sample but similar conventions are used for MC samples):

> -   `<mod_file>_prior_<namendo>_vs_<namexo>.fig`: histogram and CDF
>     plot of the MC sample of the individual entry of the shock matrix,
>     in sample and out of sample fit of the ANOVA model;
> -   `<mod_file>_prior_<namendo>_vs_<namexo>_map_SE.fig`: for entries
>     of the shock matrix it shows graphs of the estimated first order
>     ANOVA terms $y = f(\theta_i)$ for each deep parameter $\theta_i$;
> -   `<mod_file>_prior_<namendo>_vs_<namlagendo>.fig`: histogram and
>     CDF plot of the MC sample of the individual entry of the
>     transition matrix, in sample and out of sample fit of the ANOVA
>     model;
> -   `<mod_file>_prior_<namendo>_vs_<namlagendo>_map_SE.fig`: for
>     entries of the transition matrix it shows graphs of the estimated
>     first order ANOVA terms $y = f(\theta_i)$ for each deep parameter
>     $\theta_i$;
> -   `<mod_file>_prior_<namendo>_vs_<namexo>_map.mat`,
>     `<mod_file>_<namendo>_vs_<namlagendo>_map.mat`: these files store
>     info in the estimation;

When option `logtrans_redform` is set, the ANOVA estimation is performed
using a log-transformation of each y. The ANOVA mapping is then
transformed back onto the original scale, to allow comparability with
the baseline estimation. Graphs for this log-transformed case, are
stored in the same folder in files denoted with the `_log` suffix.

When the option `threshold_redform` is set, the analysis is performed
via Monte Carlo filtering, by displaying parameters that drive the
individual entry `y` inside the range specified in `threshold_redform`.
If no entry is found (or all entries are in the range), the MCF
algorithm ignores the range specified in `threshold_redform` and
performs the analysis splitting the MC sample of `y` into deciles.
Setting `threshold_redform=[-inf inf]` triggers this approach for all
`y`'s.

Results are stored in subdirectories of `<mod_file>/gsa/redform_prior`
named

> -   `<mod_file>_prior_<namendo>_vs_<namlagendo>_threshold`, for the
>     entries of the transition matrix;
> -   `<mod_file>_prior_<namendo>_vs_<namexo>_threshold`, for entries of
>     the matrix of the shocks.

The files saved are named:

> -   `<mod_file>_prior_<namendo>_vs_<namexo>_threshold.fig`,
>     `<mod_file>_<namendo>_vs_<namlagendo>_threshold.fig`: graphical
>     outputs;
> -   `<mod_file>_prior_<namendo>_vs_<namexo>_threshold.mat`,
>     `<mod_file>_<namendo>_vs_<namlagendo>_threshold.mat`: info on the
>     analysis;

#### RMSE

The RMSE analysis can be performed with different types of sampling
options:

> 1.  When `pprior=1` and `ppost=0`, the toolbox analyzes the RMSEs for
>     the Monte-Carlo sample obtained by sampling parameters from their
>     prior distributions (or prior ranges): this analysis provides some
>     hints about what parameter drives the fit of which observed
>     series, prior to the full estimation;
> 2.  When `pprior=0` and `ppost=0`, the toolbox analyzes the RMSEs for
>     a multivariate normal Monte-Carlo sample, with covariance matrix
>     based on the inverse Hessian at the optimum: this analysis is
>     useful when maximum likelihood estimation is done (i.e. no
>     Bayesian estimation);
> 3.  When `ppost=1` the toolbox analyzes the RMSEs for the posterior
>     sample obtained by Dynare's Metropolis procedure.

The use of cases 2 and 3 requires an estimation step beforehand. To
facilitate the sensitivity analysis after estimation, the
`dynare_sensitivity` command also allows you to indicate some options of
the `estimation command`. These are:

> -   `datafile`
> -   `nobs`
> -   `first_obs`
> -   `prefilter`
> -   `presample`
> -   `nograph`
> -   `nodisplay`
> -   `graph_format`
> -   `conf_sig`
> -   `loglinear`
> -   `mode_file`

Binary files produced my RMSE analysis are:

> -   `<mod_file>_prior_*.mat`: these files store the filtered and
>     smoothed variables for the prior Monte-Carlo sample, generated
>     when doing RMSE analysis (`pprior=1` and `ppost=0`);
> -   `<mode_file>_mc_*.mat`: these files store the filtered and
>     smoothed variables for the multivariate normal Monte-Carlo sample,
>     generated when doing RMSE analysis (`pprior=0` and `ppost=0`).

Figure files \<mod\_file\>\_[rmse](#rmse)\*.fig store results for the
RMSE analysis.

> -   `<mod_file>_rmse_prior*.fig`: save results for the analysis using
>     prior Monte-Carlo samples;
> -   `<mod_file>_rmse_mc*.fig`: save results for the analysis using
>     multivariate normal Monte-Carlo samples;
> -   `<mod_file>_rmse_post*.fig`: save results for the analysis using
>     Metropolis posterior samples.

The following types of figures are saved (we show prior sample to fix
ideas, but the same conventions are used for multivariate normal and
posterior):

> -   `<mod_file>_rmse_prior_params_*.fig`: for each parameter, plots
>     the cdfs corresponding to the best 10% RMSEs of each observed
>     series (only those cdfs below the significance threshold
>     `alpha_rmse`);
> -   `<mod_file>_rmse_prior_<var_obs>_*.fig`: if a parameter
>     significantly affects the fit of `var_obs`, all possible
>     trade-off's with other observables for same parameter are plotted;
> -   `<mod_file>_rmse_prior_<var_obs>_map.fig`: plots the MCF analysis
>     of parameters significantly driving the fit the observed series
>     `var_obs`;
> -   `<mod_file>_rmse_prior_lnlik*.fig`: for each observed series,
>     plots in BLUE the cdf of the log-likelihood corresponding to the
>     best 10% RMSEs, in RED the cdf of the rest of the sample and in
>     BLACK the cdf of the full sample; this allows one to see the
>     presence of some idiosyncratic behavior;
> -   `<mod_file>_rmse_prior_lnpost*.fig`: for each observed series,
>     plots in BLUE the cdf of the log-posterior corresponding to the
>     best 10% RMSEs, in RED the cdf of the rest of the sample and in
>     BLACK the cdf of the full sample; this allows one to see
>     idiosyncratic behavior;
> -   `<mod_file>_rmse_prior_lnprior*.fig`: for each observed series,
>     plots in BLUE the cdf of the log-prior corresponding to the best
>     10% RMSEs, in RED the cdf of the rest of the sample and in BLACK
>     the cdf of the full sample; this allows one to see idiosyncratic
>     behavior;
> -   `<mod_file>_rmse_prior_lik.fig`: when `lik_only=1`, this shows the
>     MCF tests for the filtering of the best 10% log-likelihood values;
> -   `<mod_file>_rmse_prior_post.fig`: when `lik_only=1`, this shows
>     the MCF tests for the filtering of the best 10% log-posterior
>     values.

#### Screening Analysis

Screening analysis does not require any additional options with respect
to those listed in `Sampling Options <sampl-opt>`{.interpreted-text
role="ref"}. The toolbox performs all the analyses required and displays
results.

The results of the screening analysis with Morris sampling design are
stored in the subfolder `<mod_file>/gsa/screen`. The data file
`<mod_file>_prior` stores all the information of the analysis (Morris
sample, reduced form coefficients, etc.).

Screening analysis merely concerns reduced form coefficients. Similar
synthetic bar charts as for the reduced form analysis with Monte-Carlo
samples are saved:

> -   `<mod_file>_redform_<endo name>_vs_lags_*.fig`: shows bar charts
>     of the elementary effect tests for the ten most important
>     parameters driving the reduced form coefficients of the selected
>     endogenous variables (`namendo`) versus lagged endogenous
>     variables (`namlagendo`);
> -   `<mod_file>_redform_<endo name>_vs_shocks_*.fig`: shows bar charts
>     of the elementary effect tests for the ten most important
>     parameters driving the reduced form coefficients of the selected
>     endogenous variables (`namendo`) versus exogenous variables
>     (`namexo`);
> -   `<mod_file>_redform_screen.fig`: shows bar chart of all elementary
>     effect tests for each parameter: this allows one to identify
>     parameters that have a minor effect for any of the reduced form
>     coefficients.

#### Identification Analysis

Setting the option `identification=1`, an identification analysis based
on theoretical moments is performed. Sensitivity plots are provided that
allow to infer which parameters are most likely to be less identifiable.

Prerequisite for properly running all the identification routines, is
the keyword `identification`; in the Dynare model file. This keyword
triggers the computation of analytic derivatives of the model with
respect to estimated parameters and shocks. This is required for option
`morris=2`, which implements *Iskrev (2010)* identification analysis.

For example, the placing:

    identification;
    dynare_sensitivity(identification=1, morris=2);

in the Dynare model file triggers identification analysis using analytic
derivatives as in *Iskrev (2010)*, jointly with the mapping of the
acceptable region.

The identification analysis with derivatives can also be triggered by
the single command:

    identification;

This does not do the mapping of acceptable regions for the model and
uses the standard random sampler of Dynare. Additionally, using only
`identification;` adds two additional identification checks: namely, of
*Qu and Tkachenko (2012)* based on the spectral density and of *Komunjer
and Ng (2011)* based on the minimal state space system. It completely
offsets any use of the sensitivity analysis toolbox.

Markov-switching SBVAR
----------------------

Given a list of variables, observed variables and a data file, Dynare
can be used to solve a Markov-switching SBVAR model according to *Sims,
Waggoner and Zha (2008)*.[^12] Having done this, you can create
forecasts and compute the marginal data density, regime probabilities,
IRFs, and variance decomposition of the model.

The commands have been modularized, allowing for multiple calls to the
same command within a `<mod_file>.mod` file. The default is to use
`<mod_file>` to tag the input (output) files used (produced) by the
program. Thus, to call any command more than once within a
`<mod_file>.mod` file, you must use the `*_tag` options described below.

::: {.command}
markov\_switching (OPTIONS\...);

Declares the Markov state variable information of a Markov-switching
SBVAR model.

*Options*

::: {.option}
chain = INTEGER

The Markov chain considered. Default: `none`.
:::

::: {.option}
number\_of\_regimes = INTEGER

Specifies the total number of regimes in the Markov Chain. This is a
required option.
:::

::: {.option}
duration = DOUBLE \| \[ROW VECTOR OF DOUBLES\]

The duration of the regimes or regimes. This is a required option. When
passed a scalar real number, it specifies the average duration for all
regimes in this chain. When passed a vector of size equal
`number_of_regimes`, it specifies the average duration of the associated
regimes (`1:number_of_regimes`) in this chain. An absorbing state can be
specified through the `restrictions <restrictions
= [[ROW VECTOR OF 3 DOUBLES],[ROW VECTOR OF 3 DOUBLES],...]>`{.interpreted-text
role="opt"} option.
:::

::: {.option}
restrictions = \[\[ROW VECTOR OF 3 DOUBLES\],\[ROW VECTOR OF 3
DOUBLES\],\...\]

Provides restrictions on this chain's regime transition matrix. Its
vector argument takes three inputs of the form:
`[current_period_regime, next_period_regime, transition_probability]`.

The first two entries are positive integers, and the third is a
non-negative real in the set \[0,1\]. If restrictions are specified for
every transition for a regime, the sum of the probabilities must be 1.
Otherwise, if restrictions are not provided for every transition for a
given regime the sum of the provided transition probabilities msut be
\<1. Regardless of the number of lags, the restrictions are specified
for parameters at time `t` since the transition probability for a
parameter at t is equal to that of the parameter at `t-1`.
:::

In case of estimating a MS-DSGE model,[^13] in addition the following
options are allowed:

::: {.option}
parameters = \[LIST OF PARAMETERS\]

This option specifies which parameters are controlled by this Markov
Chain.
:::

::: {.option}
number\_of\_lags = DOUBLE

Provides the number of lags that each parameter can take within each
regime in this chain.
:::

*Example*

>     markov_switching(chain=1, duration=2.5, restrictions=[[1,3,0],[3,1,0]]);
>
> Specifies a Markov-switching BVAR with a first chain with 3 regimes
> that all have a duration of 2.5 periods. The probability of directly
> going from regime 1 to regime 3 and vice versa is 0.

*Example*

>     markov_switching(chain=2, number_of_regimes=3, duration=[0.5, 2.5, 2.5],
>     parameter=[alpha, rho], number_of_lags=2, restrictions=[[1,3,0],[3,3,1]]);
>
> Specifies a Markov-switching DSGE model with a second chain with 3
> regimes that have durations of 0.5, 2.5, and 2.5 periods,
> respectively. The switching parameters are `alpha` and `rho`. The
> probability of directly going from regime 1 to regime 3 is 0, while
> regime 3 is an absorbing state.
:::

::: {.command}
svar (OPTIONS\...);

Each Markov chain can control the switching of a set of parameters. We
allow the parameters to be divided equation by equation and by variance
or slope and intercept.

*Options*

::: {.option}
coefficients

Specifies that only the slope and intercept in the given equations are
controlled by the given chain. One, but not both, of `coefficients` or
`variances` must appear. Default: `none`.
:::

::: {.option}
variances

Specifies that only variances in the given equations are controlled by
the given chain. One, but not both, of `coefficients` or `variances`
must appear. Default: `none`.
:::

::: {.option}
equations

Defines the equation controlled by the given chain. If not specified,
then all equations are controlled by `chain`. Default: `none`.
:::

::: {.option}
chain = INTEGER

Specifies a Markov chain defined by `markov_switching`{.interpreted-text
role="comm"}. Default: `none`.
:::
:::

::: {.command}
sbvar (OPTIONS\...);

To be documented. For now, see the wiki:
<https://archives.dynare.org/DynareWiki/SbvarOptions>

*Options*

`datafile`, `freq`, `initial_year`, `initial_subperiod`, `final_year`,
`final_subperiod`, `data`, `vlist`, `vlistlog`, `vlistper`,
`restriction_fname`, `nlags`, `cross_restrictions`,
`contemp_reduced_form`, `real_pseudo_forecast`, `no_bayesian_prior`,
`dummy_obs`, `nstates`, `indxscalesstates`, `alpha`, `beta`,
`gsig2_lmdm`, `q_diag`, `flat_prior`, `ncsk`, `nstd`, `ninv`,
`indxparr`, `indxovr`, `aband`, `indxap`, `apband`, `indximf`,
`indxfore`, `foreband`, `indxgforhat`, `indxgimfhat`, `indxestima`,
`indxgdls`, `eq_ms`, `cms`, `ncms`, `eq_cms`, `tlindx`, `tlnumber`,
`cnum`, `forecast`, `coefficients_prior_hyperparameters`
:::

::: {.block}
svar\_identification ;

This block is terminated by `end;` and contains lines of the form:

    UPPER_CHOLESKY;
    LOWER_CHOLESKY;
    EXCLUSION CONSTANTS;
    EXCLUSION LAG INTEGER; EQUATION INTEGER, VARIABLE_NAME [[,] VARIABLE_NAME...];
    RESTRICTION EQUATION INTEGER, EXPRESSION = EXPRESSION;

To be documented. For now, see the wiki:
<https://archives.dynare.org/DynareWiki/MarkovSwitchingInterface>
:::

::: {.command}
ms\_estimation (OPTIONS\...);

Triggers the creation of an initialization file for, and the estimation
of, a Markov-switching SBVAR model. At the end of the run, the $A^0$,
$A^+$, $Q$ and $\zeta$ matrices are contained in the `oo_.ms` structure.

*General Options*

::: {.option}
file\_tag = FILENAME

The portion of the filename associated with this run. This will create
the model initialization file, `init_<file_tag>.dat`. Default:
`<mod_file>`.
:::

::: {.option}
output\_file\_tag = FILENAME

The portion of the output filename that will be assigned to this run.
This will create, among other files, `est_final_<output_file_tag>.out`,
`est_intermediate_<output_file_tag>.out`. Default: `<file_tag>`.
:::

::: {.option}
no\_create\_init

Do not create an initialization file for the model. Passing this option
will cause the *Initialization Options* to be ignored. Further, the
model will be generated from the output files associated with the
previous estimation run (i.e. `est_final_<file_tag>.out`,
`est_intermediate_<file_tag>.out` or `init_<file_tag>.dat`, searched for
in sequential order). This functionality can be useful for continuing a
previous estimation run to ensure convergence was reached or for reusing
an initialization file. NB: If this option is not passed, the files from
the previous estimation run will be overwritten. Default: off (i.e.
create initialization file)
:::

*Initialization Options*

::: {.option}
coefficients\_prior\_hyperparameters = \[DOUBLE1 DOUBLE2 \... DOUBLE6\]

Sets the hyper parameters for the model. The six elements of the
argument vector have the following interpretations:

`1`

> Overall tightness for $A^0$ and $A^+$.

`2`

> Relative tightness for $A^+$.

`3`

> Relative tightness for the constant term.

`4`

> Tightness on lag decay (range: 1.2 - 1.5); a faster decay produces
> better inflation process.

`5`

> Weight on nvar sums of coeffs dummy observations (unit roots).

`6`

> Weight on single dummy initial observation including constant.

Default: `[1.0 1.0 0.1 1.2 1.0 1.0]`
:::

::: {.option}
freq = INTEGER \| monthly \| quarterly \| yearly

Frequency of the data (e.g. `monthly, 12`). Default: `4`.
:::

::: {.option}
initial\_year = INTEGER

The first year of data. Default: `none`.
:::

::: {.option}
initial\_subperiod = INTEGER

The first period of data (i.e. for quarterly data, an integer in
`[1,4]`). Default: `1`.
:::

::: {.option}
final\_year = INTEGER

The last year of data. Default: Set to encompass entire dataset.
:::

::: {.option}
final\_subperiod = INTEGER

The final period of data (i.e. for monthly data, an integer in `[1,12]`.
Default: When final\_year is also missing, set to encompass entire
dataset; when `final_year` is indicated, set to the maximum number of
subperiods given the frequency (i.e. 4 for quarterly data, 12 for
monthly,\...).
:::

::: {.option}
datafile = FILENAME

See `datafile <dataf>`{.interpreted-text role="ref"}.
:::

::: {.option}
xls\_sheet = QUOTED\_STRING

See `xls_sheet <xls_sheet = QUOTED_STRING>`{.interpreted-text
role="opt"}.
:::

::: {.option}
xls\_range = RANGE

See `xls_range <xls_range = RANGE>`{.interpreted-text role="opt"}.
:::

::: {.option}
nlags = INTEGER

The number of lags in the model. Default: `1`.
:::

::: {.option}
cross\_restrictions

Use cross $A^0$ and $A^+$ restrictions. Default: `off`.
:::

::: {.option}
contemp\_reduced\_form

Use contemporaneous recursive reduced form. Default: `off`.
:::

::: {.option}
no\_bayesian\_prior

Do not use Bayesian prior. Default: `off` (i.e. use Bayesian prior).
:::

::: {.option}
alpha = INTEGER

Alpha value for squared time-varying structural shock lambda. Default:
`1`.
:::

::: {.option}
beta = INTEGER

Beta value for squared time-varying structural shock lambda. Default:
`1`.
:::

::: {.option}
gsig2\_lmdm = INTEGER

The variance for each independent $\lambda$ parameter under `SimsZha`
restrictions. Default: `50^2`.
:::

::: {.option}
specification = sims\_zha \| none

This controls how restrictions are imposed to reduce the number of
parameters. Default: `Random Walk`.
:::

*Estimation Options*

::: {.option}
convergence\_starting\_value = DOUBLE

This is the tolerance criterion for convergence and refers to changes in
the objective function value. It should be rather loose since it will
gradually be tightened during estimation. Default: `1e-3`.
:::

::: {.option}
convergence\_ending\_value = DOUBLE

The convergence criterion ending value. Values much smaller than square
root machine epsilon are probably overkill. Default: `1e-6`.
:::

::: {.option}
convergence\_increment\_value = DOUBLE

Determines how quickly the convergence criterion moves from the starting
value to the ending value. Default: `0.1`.
:::

::: {.option}
max\_iterations\_starting\_value = INTEGER

This is the maximum number of iterations allowed in the hill-climbing
optimization routine and should be rather small since it will gradually
be increased during estimation. Default: `50`.
:::

::: {.option}
max\_iterations\_increment\_value = DOUBLE

Determines how quickly the maximum number of iterations is increased.
Default: `2`.
:::

::: {.option}
max\_block\_iterations = INTEGER

The parameters are divided into blocks and optimization proceeds over
each block. After a set of blockwise optimizations are performed, the
convergence criterion is checked and the blockwise optimizations are
repeated if the criterion is violated. This controls the maximum number
of times the blockwise optimization can be performed. Note that after
the blockwise optimizations have converged, a single optimization over
all the parameters is performed before updating the convergence value
and maximum number of iterations. Default: `100`.
:::

::: {.option}
max\_repeated\_optimization\_runs = INTEGER

The entire process described by `max_block_iterations
<max_block_iterations = INTEGER>`{.interpreted-text role="opt"} is
repeated until improvement has stopped. This is the maximum number of
times the process is allowed to repeat. Set this to `0` to not allow
repetitions. Default: `10`.
:::

::: {.option}
function\_convergence\_criterion = DOUBLE

The convergence criterion for the objective function when
`max_repeated_optimizations_runs` is positive. Default: `0.1`.
:::

::: {.option}
parameter\_convergence\_criterion = DOUBLE

The convergence criterion for parameter values when
`max_repeated_optimizations_runs` is positive. Default: `0.1`.
:::

::: {.option}
number\_of\_large\_perturbations = INTEGER

The entire process described by `max_block_iterations
<max_block_iterations = INTEGER>`{.interpreted-text role="opt"} is
repeated with random starting values drawn from the posterior. This
specifies the number of random starting values used. Set this to `0` to
not use random starting values. A larger number should be specified to
ensure that the entire parameter space has been covered. Default: `5`.
:::

::: {.option}
number\_of\_small\_perturbations = INTEGER

The number of small perturbations to make after the large perturbations
have stopped improving. Setting this number much above `10` is probably
overkill. Default: `5`.
:::

::: {.option}
number\_of\_posterior\_draws\_after\_perturbation = INTEGER

The number of consecutive posterior draws to make when producing a small
perturbation. Because the posterior draws are serially correlated, a
small number will result in a small perturbation. Default: `1`.
:::

::: {.option}
max\_number\_of\_stages = INTEGER

The small and large perturbation are repeated until improvement has
stopped. This specifies the maximum number of stages allowed. Default:
`20`.
:::

::: {.option}
random\_function\_convergence\_criterion = DOUBLE

The convergence criterion for the objective function when
`number_of_large_perturbations` is positive. Default: `0.1`.
:::

::: {.option}
random\_parameter\_convergence\_criterion = DOUBLE

The convergence criterion for parameter values when
`number_of_large_perturbations` is positive. Default: `0.1`.
:::

*Example*

>     ms_estimation(datafile=data, initial_year=1959, final_year=2005,
>     nlags=4, max_repeated_optimization_runs=1, max_number_of_stages=0);
>
>     ms_estimation(file_tag=second_run, datafile=data, initial_year=1959,
>     final_year=2005, nlags=4, max_repeated_optimization_runs=1,
>     max_number_of_stages=0);
>
>     ms_estimation(file_tag=second_run, output_file_tag=third_run,
>     no_create_init, max_repeated_optimization_runs=5,
>     number_of_large_perturbations=10);
:::

::: {.command}
ms\_simulation ; ms\_simulation (OPTIONS\...);

Simulates a Markov-switching SBVAR model.

*Options*

::: {.option}
file\_tag = FILENAME

The portion of the filename associated with the `ms_estimation` run.
Default: `<mod_file>`.
:::

::: {.option}
output\_file\_tag = FILENAME

The portion of the output filename that will be assigned to this run.
Default: `<file_tag>`.
:::

::: {.option}
mh\_replic = INTEGER

The number of draws to save. Default: `10,000`.
:::

::: {.option}
drop = INTEGER

The number of burn-in draws. Default: `0.1*mh_replic*thinning_factor`.
:::

::: {.option}
thinning\_factor = INTEGER

The total number of draws is equal to `thinning_factor*mh_replic+drop`.
Default: `1`.
:::

::: {.option}
adaptive\_mh\_draws = INTEGER

Tuning period for Metropolis-Hastings draws. Default: `30,000`.
:::

::: {.option}
save\_draws

Save all elements of $A^0$, $A^+$, $Q$, and $\zeta$, to a file named
`draws_<<file_tag>>.out` with each draw on a separate line. A file that
describes how these matrices are laid out is contained in
`draws_header_<<file_tag>>.out`. A file called `load_flat_file.m` is
provided to simplify loading the saved files into the corresponding
variables `A0`, `Aplus`, `Q`, and `Zeta` in your MATLAB/Octave
workspace. Default: `off`.
:::

*Example*

>     ms_simulation(file_tag=second_run);
>     ms_simulation(file_tag=third_run, mh_replic=5000, thinning_factor=3);
:::

::: {.command}
ms\_compute\_mdd ; ms\_compute\_mdd (OPTIONS\...);

Computes the marginal data density of a Markov-switching SBVAR model
from the posterior draws. At the end of the run, the Muller and Bridged
log marginal densities are contained in the `oo_.ms` structure.

*Options*

::: {.option}
file\_tag = FILENAME

See `file_tag <file_tag = FILENAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
output\_file\_tag = FILENAME

See `output_file_tag <output_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
simulation\_file\_tag = FILENAME

The portion of the filename associated with the simulation run. Default:
`<file_tag>`.
:::

::: {.option}
proposal\_type = INTEGER

The proposal type:

`1`

> Gaussian.

`2`

> Power.

`3`

> Truncated Power.

`4`

> Step.

`5`

> Truncated Gaussian.

Default: `3`
:::

::: {.option}
proposal\_lower\_bound = DOUBLE

The lower cutoff in terms of probability. Not used for `proposal_type`
in `[1,2]`. Required for all other proposal types. Default: `0.1`.
:::

::: {.option}
proposal\_upper\_bound = DOUBLE

The upper cutoff in terms of probability. Not used for `proposal_type`
equal to `1`. Required for all other proposal types. Default: `0.9`.
:::

::: {.option}
mdd\_proposal\_draws = INTEGER

The number of proposal draws. Default: `100,000`.
:::

::: {.option}
mdd\_use\_mean\_center

Use the posterior mean as center. Default: `off`.
:::
:::

::: {.command}
ms\_compute\_probabilities ; ms\_compute\_probabilities (OPTIONS\...);

Computes smoothed regime probabilities of a Markov-switching SBVAR
model. Output `.eps` files are contained in
`<output_file_tag/Output/Probabilities>`.

*Options*

::: {.option}
file\_tag = FILENAME

See `file_tag <file_tag = FILENAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
output\_file\_tag = FILENAME

See `output_file_tag <output_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
filtered\_probabilities

Filtered probabilities are computed instead of smoothed. Default: `off`.
:::

::: {.option}
real\_time\_smoothed

Smoothed probabilities are computed based on time `t` information for
$0\le t\le nobs$. Default: `off`
:::
:::

::: {.command}
ms\_irf ; ms\_irf (OPTIONS\...);

Computes impulse response functions for a Markov-switching SBVAR model.
Output `.eps` files are contained in `<output_file_tag/Output/IRF>`,
while data files are contained in `<output_file_tag/IRF>`.

*Options*

::: {.option}
file\_tag = FILENAME

See `file_tag <file_tag = FILENAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
output\_file\_tag = FILENAME

See `output_file_tag <output_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
simulation\_file\_tag = FILENAME

See
`simulation_file_tag <simulation_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
horizon = INTEGER

The forecast horizon. Default: `12`.
:::

::: {.option}
filtered\_probabilities

Uses filtered probabilities at the end of the sample as initial
conditions for regime probabilities. Only one of
`filtered_probabilities`, `regime` and `regimes` may be passed. Default:
`off`.
:::

::: {.option}
error\_band\_percentiles = \[DOUBLE1 \...\]

The percentiles to compute. Default: `[0.16 0.50 0.84]`. If `median` is
passed, the default is `[0.5]`.
:::

::: {.option}
shock\_draws = INTEGER

The number of regime paths to draw. Default: `10,000`.
:::

::: {.option}
shocks\_per\_parameter = INTEGER

The number of regime paths to draw under parameter uncertainty. Default:
`10`.
:::

::: {.option}
thinning\_factor = INTEGER

Only $1/ \texttt{thinning\_factor}$ of the draws in posterior draws file
are used. Default: `1`.
:::

::: {.option}
free\_parameters = NUMERICAL\_VECTOR

A vector of free parameters to initialize theta of the model. Default:
use estimated parameters
:::

::: {.option}
parameter\_uncertainty

Calculate IRFs under parameter uncertainty. Requires that
`ms_simulation` has been run. Default: `off`.
:::

::: {.option}
regime = INTEGER

Given the data and model parameters, what is the ergodic probability of
being in the specified regime. Only one of `filtered_probabilities`,
`regime` and `regimes` may be passed. Default: `off`.
:::

::: {.option}
regimes

Describes the evolution of regimes. Only one of
`filtered_probabilities`, `regime` and `regimes` may be passed. Default:
`off`.
:::

::: {.option}
median

A shortcut to setting `error_band_percentiles=[0.5]`. Default: `off`.
:::
:::

::: {.command}
ms\_forecast ; ms\_forecast (OPTIONS\...);

Generates forecasts for a Markov-switching SBVAR model. Output `.eps`
files are contained in `<output_file_tag/Output/Forecast>`, while data
files are contained in `<output_file_tag/Forecast>`.

*Options*

::: {.option}
file\_tag = FILENAME

See `file_tag <file_tag = FILENAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
output\_file\_tag = FILENAME

See `output_file_tag <output_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
simulation\_file\_tag = FILENAME

See
`simulation_file_tag <simulation_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
data\_obs\_nbr = INTEGER

The number of data points included in the output. Default: `0`.
:::

::: {.option}
error\_band\_percentiles = \[DOUBLE1 \...\]

See `error_band_percentiles <error_band_percentiles =
[DOUBLE1 ...]>`{.interpreted-text role="opt"}.
:::

::: {.option}
shock\_draws = INTEGER

See `shock_draws <shock_draws = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
shocks\_per\_parameter = INTEGER

See
`shocks_per_parameter <shocks_per_parameter = INTEGER>`{.interpreted-text
role="opt"}.
:::

::: {.option}
thinning\_factor = INTEGER

See `thinning_factor <thinning_factor = INTEGER>`{.interpreted-text
role="opt"}.
:::

::: {.option}
free\_parameters = NUMERICAL\_VECTOR

See
`free_parameters <free_parameters = NUMERICAL_VECTOR>`{.interpreted-text
role="opt"}.
:::

::: {.option}
parameter\_uncertainty

See `parameter_uncertainty`{.interpreted-text role="opt"}.
:::

::: {.option}
regime = INTEGER

See `regime <regime = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
regimes

See `regimes`{.interpreted-text role="opt"}.
:::

::: {.option}
median

See `median`{.interpreted-text role="opt"}.
:::

::: {.option}
horizon = INTEGER

See `horizon <horizon = INTEGER>`{.interpreted-text role="opt"}.
:::
:::

::: {.command}
ms\_variance\_decomposition ; ms\_variance\_decomposition (OPTIONS\...);

Computes the variance decomposition for a Markov-switching SBVAR model.
Output `.eps` files are contained in
`<output_file_tag/Output/Variance_Decomposition>`, while data files are
contained in `<output_file_tag/Variance_Decomposition>`.

*Options*

::: {.option}
file\_tag = FILENAME

See `file_tag <file_tag = FILENAME>`{.interpreted-text role="opt"}.
:::

::: {.option}
output\_file\_tag = FILENAME

See `output_file_tag <output_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
simulation\_file\_tag = FILENAME

See
`simulation_file_tag <simulation_file_tag = FILENAME>`{.interpreted-text
role="opt"}.
:::

::: {.option}
horizon = INTEGER

See `horizon <horizon = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
filtered\_probabilities

See `filtered_probabilities`{.interpreted-text role="opt"}.
:::

::: {.option}
no\_error\_bands

Do not output percentile error bands (i.e. compute mean). Default: `off`
(i.e. output error bands)
:::

::: {.option}
error\_band\_percentiles = \[DOUBLE1 \...\]

See `error_band_percentiles <error_band_percentiles =
[DOUBLE1 ...]>`{.interpreted-text role="opt"}.
:::

::: {.option}
shock\_draws = INTEGER

See `shock_draws <shock_draws = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
shocks\_per\_parameter = INTEGER

See
`shocks_per_parameter <shocks_per_parameter = INTEGER>`{.interpreted-text
role="opt"}.
:::

::: {.option}
thinning\_factor = INTEGER

See `thinning_factor <thinning_factor = INTEGER>`{.interpreted-text
role="opt"}.
:::

::: {.option}
free\_parameters = NUMERICAL\_VECTOR

See
`free_parameters <free_parameters = NUMERICAL_VECTOR>`{.interpreted-text
role="opt"}.
:::

::: {.option}
parameter\_uncertainty

See `parameter_uncertainty`{.interpreted-text role="opt"}.
:::

::: {.option}
regime = INTEGER

See `regime <regime = INTEGER>`{.interpreted-text role="opt"}.
:::

::: {.option}
regimes

See `regimes`{.interpreted-text role="opt"}.
:::
:::

Epilogue Variables {#epilogue}
------------------

::: {.block}
epilogue ;
:::

The epilogue block is useful for computing output variables of interest
that may not be necessarily defined in the model (e.g. various kinds of
real/nominal shares or relative prices, or annualized variables out of a
quarterly model).

It can also provide several advantages in terms of computational
efficiency and flexibility:

-   You can calculate variables in the epilogue block after
    smoothers/simulations have already been run without adding the new
    definitions and equations and rerunning smoothers/simulations. Even
    posterior smoother subdraws can be recycled for computing epilogue
    variables without rerunning subdraws with the new definitions and
    equations.
-   You can also reduce the state space dimension in data
    filtering/smoothing. Assume, for example, you want annualized
    variables as outputs. If you define an annual growth rate in a
    quarterly model, you need lags up to order 7 of the associated
    quarterly variable; in a medium/large scale model this would just
    blow up the state dimension and increase by a huge amount the
    computing time of a smoother.

The `epilogue` block is terminated by `end;` and contains lines of the
form:

> NAME = EXPRESSION;

*Example*

:   epilogue;
        // annualized level of y
        ya = exp(y)+exp(y(-1))+exp(y(-2))+exp(y(-3));
        // annualized growth rate of y
        gya = ya/ya(-4)-1;
        end;

Semi-structural models {#semi-strutural}
----------------------

Dynare provides tools for semi-structural models, in the vain of the
FRB/US model (see *Brayton and Tinsley (1996)*), where expectations are
not necessarily model consistent but based on a VAR auxiliary model. In
the following, it is assumed that each equation is written as
`VARIABLE = EXPRESSION` or `T(VARIABLE) = EXPRESSION` where
`T(VARIABLE)` stands for a transformation of an endogenous variable
(`log` or `diff`). This representation, where each equation determines
the endogenous variable on the LHS, can be exploited when simulating the
model (see algorithms 12 and 14 in
`solve_algo <solvalg>`{.interpreted-text role="ref"}) and is mandatory
to define auxiliary models used for computing expectations (see below).

### Auxiliary models

The two auxiliary models defined in this section are linear
backward-looking models used to form expectations. Both models can be
recast as VAR(1)-processes and therefore offer isomorphic ways of
specifying the expectations process, but differ in their convenience of
specifying features like cointegration and error correction. `var_model`
directly specifies a VAR, while `trend_component_model` allows to define
a trend target to which the endogenous variables may be attracted in the
long-run (i.e. an error correction model).

::: {.command}
var\_model (OPTIONS\...);

Picks equations in the `model` block to form a VAR model. This model can
be used as an auxiliary model in `var_expectation_model` or `pac_model`.
It must be of the following form:

or

if the VAR is structural (see below), where $Y_t$ and $\varepsilon_t$
are $n\times 1$ vectors, $\mathbf{c}$ is a $n\times 1$ vector of
parameters, $A_i$ ($i=0,\ldots,p$) are $n\times n$ matrices of
parameters, and $A_0$ is non singular square matrix. Vector $\mathbf{c}$
and matrices $A_i$ ($i=0,\ldots,p$) are set by Dynare by parsing the
equations in the `model` block. Then, Dynare builds a VAR(1)-companion
form model for $\mathcal{Y}_t = (1, Y_t, \ldots, Y_{t-p+1})'$ as:

assuming that we are dealing with a reduced form VAR (otherwise, the
right-hand side would additionally be premultiplied by $A_0^{-1}.$ to
obtain the reduced for representation). If the VAR does not have a
constant, we remove the first line of the system and the first column of
the companion matrix $\mathcal{C}.$ Dynare only saves the companion
matrix, since that is the only information required to compute the
expectations.

::: {.matvar}
[oo]().var.MODEL\_NAME.CompanionMatrix

Reduced form companion matrix of the `var_model`.
:::

*Options*

> ::: {.option}
> model\_name = STRING
> :::
>
> Name of the VAR model, which will be referenced in
> `var_expectation_model` or `pac_model` as an
> [auxiliary\_model]{.title-ref}. Needs to be a valid Matlab field name.
>
> ::: {.option}
> eqtags = \[QUOTED\_STRING\[, QUOTED\_STRING\[, \...\]\]\]
> :::
>
> List of equations in the `model` block (referenced using the equation
> tag `name`) used to build the VAR model.
>
> ::: {.option}
> structural
> :::
>
> By default the VAR model is not structural, *i.e.* each equation must
> contain exactly one contemporaneous variable (on the LHS). If the
> `structural` option is provided then any variable defined in the
> system can appear at time $t$ in each equation. Internally Dynare will
> rewrite this model as a reduced form VAR (by inverting the implied
> matrix $A_0$).

*Example*

>     var_model(model_name = toto, eqtags = [ 'X', 'Y', 'Z' ]);
>
>     model;
>
>     [ name = 'X' ]
>     x = a*x(-1) + b*x(-2) + c*z(-2) + e_x;
>
>     [ name = 'Z' ]
>     z = f*z(-1) + e_z;
>
>     [ name = 'Y' ]
>     y = d*y(-2) + e*z(-1) + e_y;
>
>     end;
:::

::: {.command}
trend\_component\_model (OPTIONS\...);

Picks equations in the model block to form a trend component model. This
model can be used as an auxiliary model in `var_expectation_model` or
`pac_model`. It must be of the following form:

where $X_t$ and $Z_t$ are $n\times 1$ and $m\times 1$ vectors of
endogenous variables. $Z_t$ defines the trend target to whose linear
combination $C_0 Z_t$ the endogenous variables $X_t$ will be attracted,
provided the implied error correction matrix $A_0$ is negative definite.
$\varepsilon_t$ and $\eta_t$ are $n\times 1$ and $m\times 1$ vectors of
exogenous variables, $A_i$ ($i=0,\ldots,p$) are $n\times n$ matrices of
parameters, and $C_0$ is a $n\times m$ matrix. This model can also be
cast into a VAR(1) model by first rewriting it in levels. Let
$Y_t = (X_t',Z_t')'$ and $\zeta_t = (\varepsilon_t',\eta_t')'$. Then we
have:

with

where $\Lambda = A_0C_0$,

for $i=2,\ldots,p$, and

This VAR(p+1) in levels can again be rewritten as a VAR(1)-companion
model form.

::: {.matvar}
[oo]().trend\_component.MODEL\_NAME.CompanionMatrix

Reduced form companion matrix of the `trend_component_model`.
:::

*Options*

> ::: {.option}
> model\_name = STRING
> :::
>
> Name of the trend component model, will be referenced in
> `var_expectation_model` or `pac_model` as an
> [auxiliary\_model]{.title-ref}. Needs to be a valid Matlab field name.
>
> ::: {.option}
> eqtags = \[QUOTED\_STRING\[, QUOTED\_STRING\[, \...\]\]\]
> :::
>
> List of equations in the `model` block (referenced using the equation
> tag `name`) used to build the trend component model.
>
> ::: {.option}
> targets = \[QUOTED\_STRING\[, QUOTED\_STRING\[, \...\]\]\]
> :::
>
> List of targets, corresponding to the variables in vector $Z_t$,
> referenced using the equation tag `name`) of the associated equation
> in the `model` block. `target` must be a subset of `eqtags`.

*Example*

>     trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);
>
>     model;
>
>     [name='eq:x1']
>     diff(x1) = a_x1_0*(x1(-1)-x1bar(-1))+a_x1_0_*(x2(-1)-x2bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;
>
>     [name='eq:x2']
>     diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;
>
>     [name='eq:x1bar']
>     x1bar = x1bar(-1) + ex1bar;
>
>     [name='eq:x2bar']
>     x2bar = x2bar(-1) + ex2bar;
>
>     end;
:::

### VAR expectations

Suppose we wish to forecast a variable $y_t$ and that $y_t$ is an
element of vector of variables $\mathcal{Y}_t$ whose law of motion is
described by a VAR(1) model $\mathcal{Y}_t =
\mathcal{C}\mathcal{Y}_{t-1}+\epsilon_t$. More generally, $y_t$ may be a
linear combination of the scalar variables in $\mathcal{Y}_t$. Let the
vector $\alpha$ be such that $y_t = \alpha'\mathcal{Y}_t$ ($\alpha$ is a
selection vector if $y_t$ is a variable in $\mathcal{Y}_t$, *i.e.* a
column of an identity matrix, or an arbitrary vector defining the
weights of a linear combination). Then the best prediction, in the sense
of the minimisation of the RMSE, for $y_{t+h}$ given the information set
at $t-\tau$ (which we assume to include all observables up to time
$t-\tau$, $\mathcal{Y}_{\underline{t-\tau}}$) is:

In a semi-structural model, variables appearing in $t+h$ (*e.g.* the
expected output gap in a dynamic IS curve or expected inflation in a
(New Keynesian) Phillips curve) will be replaced by the expectation
implied by an auxiliary VAR model. Another use case is for the
computation of permanent incomes. Typically, consumption will depend on
something like:

Assuming that \$0\<beta\<1\$ and knowing the limit of geometric series,
the conditional expectation of this variable can be evaluated based on
the same auxiliary model:

More generally, it is possible to consider finite discounted sums.

::: {.command}
var\_expectation\_model (OPTIONS\...);

Declares a model used to forecast an endogenous variable or linear
combination of variables in $t+h$. More generally, the same model can be
used to forecast the discounted flow of a variable or a linear
expression of variables:

where $(a,b)\in\mathbb N^2$ with $a<b$, $\beta\in(0,1]$ is a discount
factor, and $\tau$ is a finite positive integer.

*Options*

> ::: {.option}
> model\_name = STRING
> :::
>
> Name of the VAR based expectation model, which will be referenced in
> the `model` block.
>
> ::: {.option}
> auxiliary\_model = STRING
> :::
>
> Name of the associated auxiliary model, defined with `var_model` or
> `trend_component_model`.
>
> ::: {.option}
> expression = STRING \| QUOTED\_STRING
> :::
>
> Name of the variable or expression (linear combination of variables)
> to be expected. The quotes are mandatory if an expression is expected.
>
> ::: {.option}
> discount = PARAMETER\_NAME \| DOUBLE
> :::
>
> Discount factor ($\beta$).
>
> ::: {.option}
> horizon = INTEGER \| \[INTEGER:INTEGER\]
> :::
>
> The upper limit $b$ of the horizon $h$ (in which case $a=0$), or range
> of periods $a:b$ over which the discounted sum is computed (the upper
> bound can be `Inf`).
>
> ::: {.option}
> time\_shift = INTEGER
> :::
>
> Shift of the information set ($\tau$), default value is 0.
:::

::: {.operator}
var\_expectation (NAME\_OF\_VAR\_EXPECTATION\_MODEL);

This operator is used instead of a leaded variable, e.g. `X(1)`, in the
`model` block to substitute a model-consistent forecast with a forecast
based on a VAR model.
:::

> *Example*
>
> >     var_model(model_name=toto, eqtags=['X', 'Y', 'Z']);
> >
> >     var_expectation_model(model_name=varexp, expression=x, auxiliary_model_name=toto, horizon=1, discount=beta);
> >
> >
> >     model;
> >
> >     [name='X']
> >     x = a*x(-1) + b*x(-2) + c*z(-2) + e_x;
> >
> >     [name='Z']
> >     z = f*z(-1) + e_z;
> >
> >     [name='Y']
> >     y = d*y(-2) + e*z(-1) + e_y;
> >
> >     foo = .5*foo(-1) + var_expectation(varexp);
> >
> >     end;
> >
> > In this example `var_expectation(varexp)` stands for the one step
> > ahead expectation of `x`, as a replacement for `x(1)`.

::: {.matcomm}
var\_expectation.initialize(NAME\_OF\_VAR\_EXPECTATION\_MODEL);

Initialise the [var\_expectation\_model]{.title-ref} by building the
companion matrix of the associated auxiliary [var\_model]{.title-ref}.
Needs to be executed before attempts to simulate or estimate the model.
:::

::: {.matcomm}
var\_expectation.update(NAME\_OF\_VAR\_EXPECTATION\_MODEL);

Update/compute the reduced form parameters of
[var\_expectation\_model]{.title-ref}. Needs to be executed before
attempts to simulate or estimate the model and requires the auxiliary
[var\_model]{.title-ref} to have previously been initialized.

*Example (continued)*

>     var_expectation.initialize('varexp');
>
>     var_expectation.update('varexp');
:::

::: {.warning}
::: {.admonition-title}
Warning
:::

Changes to the parameters of the underlying auxiliary
[var\_model]{.title-ref} require calls to
[var\_expectation.initialize]{.title-ref} and
[var\_expectation.update]{.title-ref} to become effective. Changes to
the [var\_expectation\_model]{.title-ref} or its associated parameters
require a call to [var\_expectation.update]{.title-ref}.
:::

### PAC equation

In its simplest form, a PAC equation breaks down changes in a variable
of interest $y$ into three contributions: (*i*) the lagged deviation
from a target $y^{\star}$, (*ii*) the lagged changes in the variable
$y$, and (*iii*) the expected changes in the target $y^{\star}$:

*Brayton et alii (2000)* shows how such an equation can be derived from
the minimisation of a quadratic cost function penalising expected
deviations from the target and non-smoothness of $y$, where future costs
are discounted (with discount factor $\beta$). They also show that the
parameters $(d_i)_{i\in\mathbb N}$ are non-linear functions of the $m$
parameters $a_i$ and the discount factor $\beta$. To simulate or
estimate this equation we need to figure out how to determine the
expected changes of the target. This can be done as in the previous
section using VAR based expectations, or considering model consistent
expectations (MCE).

To ensure that the endogenous variable $y$ is equal to its target
$y^{\star}$ in the (deterministic) long run, *i.e.* that the error
correction term is zero in the long run, we can optionally add a growth
neutrality correction to this equation. Suppose that $g$ is the long run
growth rate, for $y$ and $y^{\star}$, then in the long run (assuming
that the data are in logs) we must have:

Unless additional restrictions are placed on the coefficients
$(a_i)_{i=0}^{m-1}$, i.e. on the form of the minimised cost function,
there is no reason for the right-hand side to be zero. Instead, we can
optionally add the right hand side to the PAC equation, to ensure that
the error correction term is asymptotically zero.

The PAC equations can be generalised by adding exogenous variables. This
can be done in two, non exclusive, manners. We can replace the PAC
equation by a convex combination of the original PAC equation (derived
from an optimisation program) and a linear expression involving
exogenous variables (referred as the rule of thumb part as opposed to
the part derived from the minimisation of a cost function; not to be
confused with exogenous shocks):

where $\lambda\in[0,1]$ is the weight of the pure PAC equation, $\gamma$
is a $k\times 1$ vector of parameters, and $X_t$ a $k\times 1$ vector of
variables in the rule of thumb part. Or we can simply add the exogenous
variables to the PAC equation (without the weight $\lambda$):

::: {.command}
pac\_model (OPTIONS\...);

Declares a PAC model. A [.mod]{.title-ref} file can have more than one
PAC model or PAC equation, but each PAC equation must be associated to a
different PAC model.

*Options*

> ::: {.option}
> model\_name = STRING
> :::
>
> Name of the PAC model, will be referenced in the `model` block.
>
> ::: {.option}
> auxiliary\_model = STRING
> :::
>
> Name of the associated auxiliary model, defined with `var_model` or
> `trend_component_model`, to compute the VAR based expectations for the
> expected changes in the target, *i.e.* to evaluate
> $\sum_{i=0}^{\infty} d_i \Delta y^{\star}_{t+i}$. The infinite sum
> will then be replaced by a linear combination of the variables
> involved in the companion representation of the auxiliary model. The
> weights defining the linear combination are nonlinear functions of the
> $(a_i)_{i=0}^{m-1}$ coefficients in the PAC equation. This option is
> not mandatory, if absent Dynare understands that the expected changes
> of the target have to be computed under the MCE assumption. This is
> done by rewriting recursively the infinite sum as shown in equation 10
> of *Brayton et alii (2000)*.
>
> ::: {.option}
> discount = PARAMETER\_NAME \| DOUBLE
> :::
>
> Discount factor ($\beta$) for future expected costs appearing in the
> definition of the cost function.
>
> ::: {.option}
> growth = PARAMETER\_NAME \| VARIABLE\_NAME \| EXPRESSION \| DOUBLE
> :::
>
> If present a growth neutrality correction is added to the PAC
> equation. The user must ensure that the provided value (or long term
> level if a variable or expression is given) is consistent with the
> asymptotic growth rate of the endogenous variable.
:::

::: {.operator}
pac\_expectation (NAME\_OF\_PAC\_MODEL);

This operator is used instead of the infinite sum,
$\sum_{i=0}^{\infty} d_i \Delta y^{\star}_{t+i}$, in a PAC equation
defined in the `model` block. Depending on the assumption regarding the
formation of expectations, it will be replaced by a linear combination
of the variables involved in the companion representation of the
auxiliary model or by a recursive forward equation.
:::

::: {.matcomm}
pac.initialize(NAME\_OF\_PAC\_MODEL);
:::

::: {.matcomm}
pac.update(NAME\_OF\_PAC\_MODEL);

Same as in the previous section for the VAR expectations, initialise the
PAC model, by building the companion matrix of the auxiliary model, and
computes the reduced form parameters of the PAC equation (the weights in
the linear combination of the variables involved in the companion
representation of the auxiliary model, or the parameters of the
recursive representation of the infinite sum in the MCE case).
:::

*Example*

>     trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);
>
>     pac_model(auxiliary_model_name=toto, discount=beta, growth=diff(x1(-1)), model_name=pacman);
>
>     model;
>
>     [name='eq:y']
>     y = rho_1*y(-1) + rho_2*y(-2) + ey;
>
>     [name='eq:x1']
>     diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;
>
>     [name='eq:x2']
>     diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;
>
>     [name='eq:x1bar']
>     x1bar = x1bar(-1) + ex1bar;
>
>     [name='eq:x2bar']
>     x2bar = x2bar(-1) + ex2bar;
>
>     [name='zpac']
>     diff(z) = e_c_m*(x1(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman) + ez;
>
>     end;
>
>     pac.initialize('pacman');
>
>     pac.update.expectation('pacman');

### Estimation of a PAC equation

The PAC equation, introduced in the previous section, can be estimated.
This equation is nonlinear with respect to the estimated parameters
$(a_i)_{i=0}^{m-1}$, since the reduced form parameters (in the
computation of the infinite sum) are nonlinear functions of the
autoregressive parameters and the error correction parameter. *Brayton
et alii (2000)* shows how to estimate the PAC equation by iterative OLS.
Although this approach is implemented in Dynare, mainly for comparison
purposes, we also propose NLS estimation, which is much preferable
(asymptotic properties of NLS being more solidly grounded).

Note that it is currently not feasible to estimate the PAC equation
jointly with the remaining parameters of the model using e.g. Bayesian
techniques. Thus, estimation of the PAC equation can only be conducted
conditional on the values of the parameters of the auxiliary model.

::: {.warning}
::: {.admonition-title}
Warning
:::

The estimation routines described below require the option
[json=compute]{.title-ref} be passed to the preprocessor (via the
command line or at the top of the [.mod]{.title-ref} file, see
`dyn-invoc`{.interpreted-text role="ref"}).
:::

::: {.matcomm}
pac.estimate.nls(EQNAME, GUESS, DATA, RANGE\[, ALGO\]);
:::

::: {.matcomm}
pac.estimate.iterative\_ols(EQNAME, GUESS, DATA, RANGE);

Trigger the NLS or iterative OLS estimation of a PAC equation. `EQNAME`
is a row char array designating the PAC equation to be estimated (the
PAC equation must have a name specified with an equation tag). `DATA` is
a `dseries` object containing the data required for the estimation
(*i.e.* data for all the endogenous and exogenous variables in the
equation). The residual values of the PAC equation (which correspond to
a defined [varexo]{.title-ref}) must also be a member of `DATA`, but
filled with `NaN` values. `RANGE` is a `dates` object defining the time
span of the sample. `ALGO` is a row char array used to select the method
(or minimisation algorithm) for NLS. Possible values are : `'fmincon'`,
`'fminunc'`, `'fminsearch'`, `'lsqnonlin'`, `'particleswarm'`,
`'csminwel'`, `'simplex'`, `'annealing'`, and `'GaussNewton'`. The first
four algorithms require the Mathworks Optimisation toolbox. The fifth
algorithm requires the Mathworks Global Optimisation toolbox. When the
optimisation algorithm allows it, we impose constraints on the error
correction parameter, which must be positive and smaller than 1 (it the
case for `'fmincon'`, `'lsqnonlin'`, `'particleswarm'`, and
`'annealing'`). The default optimisation algorithm is `'csminwel'`.
`GUESS` is a structure containing the initial guess values for the
estimated parameters. Each field is the name of a parameter in the PAC
equation and holds the initial guess for this parameter. If some
parameters are calibrated, then they should not be members of the
`GUESS` structure (and values have to be provided in the `.mod` file
before the call to the estimation routine).

For the NLS routine the estimation results are displayed in a table
after the estimation. For both the NLS and iterative OLS routines, the
results are saved in `oo_` (under the fields `nls` or `iterative_ols`).
Also, the values of the parameters are updated in `M_.params`.
:::

*Example (continued)*

>     // Set the initial guess for the estimated parameters
>     eparams.e_c_m =  .9;
>     eparams.c_z_1 =  .5;
>     eparams.c_z_2 =  .2;
>
>     // Define the dataset used for estimation
>     edata = TrueData;
>     edata.ez = dseries(NaN); // Set to NaN the residual of the equation.
>
>     pac.estimate.nls('zpac', eparams, edata, 2005Q1:2005Q1+200, 'annealing');

::: {.warning}
::: {.admonition-title}
Warning
:::

The specification of [GUESS]{.title-ref} and [DATA]{.title-ref} involves
the use of structures. As such, their subfields will not be cleared
across Dynare runs as the structures stay in the workspace. Be careful
to clear these structures from the memory (e.g. within the mod-file)
when e.g. changing which parameters are calibrated.
:::

Displaying and saving results
-----------------------------

Dynare has comments to plot the results of a simulation and to save the
results.

::: {.command}
rplot VARIABLE\_NAME\...;

Plots the simulated path of one or several variables, as stored in
`oo_.endo_simul` by either `perfect_foresight_solver`, `simul` (see
`det-simul`{.interpreted-text role="ref"}) or `stoch_simul` with option
`periods` (see `stoch-sol`{.interpreted-text role="ref"}). The variables
are plotted in levels.
:::

::: {.command}
dynatype (FILENAME) \[VARIABLE\_NAME\...\];

This command prints the listed endogenous or exogenous variables in a
text file named FILENAME. If no VARIABLE\_NAME is listed, all endogenous
variables are printed.
:::

::: {.command}
dynasave (FILENAME) \[VARIABLE\_NAME\...\];

This command saves the listed endogenous or exogenous variables in a
binary file named FILENAME. If no VARIABLE\_NAME is listed, all
endogenous variables are saved.

In MATLAB or Octave, variables saved with the `dynasave` command can be
retrieved by the command:

    load(FILENAME,'-mat')
:::

Macro processing language {#macro-proc-lang}
-------------------------

It is possible to use "macro" commands in the `.mod` file for performing
tasks such as: including modular source files, replicating blocks of
equations through loops, conditionally executing some code, writing
indexed sums or products inside equations\...

The Dynare macro-language provides a new set of *macro-commands* which
can be used in `.mod` files. It features:

> -   File inclusion
> -   Loops (`for` structure)
> -   Conditional inclusion (`if/then/else` structures)
> -   Expression substitution

This macro-language is totally independent of the basic Dynare language,
and is processed by a separate component of the Dynare pre-processor.
The macro processor transforms a `.mod` file with macros into a `.mod`
file without macros (doing expansions/inclusions), and then feeds it to
the Dynare parser. The key point to understand is that the macro
processor only does text substitution (like the C preprocessor or the
PHP language). Note that it is possible to see the output of the macro
processor by using the `savemacro` option of the `dynare` command (see
`dyn-invoc`{.interpreted-text role="ref"}).

The macro processor is invoked by placing *macro directives* in the
`.mod` file. Directives begin with an at-sign followed by a pound sign
(`@#`). They produce no output, but give instructions to the macro
processor. In most cases, directives occupy exactly one line of text. If
needed, two backslashes (`\\`) at the end of the line indicate that the
directive is continued on the next line. Macro directives following `//`
are not interpreted by the macro processor. For historical reasons,
directives in commented blocks, *ie* surrounded by `/*` and `*/`, are
interpreted by the macro processor. The user should not rely on this
behavior. The main directives are:

> -   `@#includepath`, paths to search for files that are to be
>     included,
> -   `@#include`, for file inclusion,
> -   `@#define`, for defining a macro processor variable,
> -   `@#if, @#ifdef, @#ifndef, @#elseif, @#else, @#endif` for
>     conditional statements,
> -   `@#for, @#endfor` for constructing loops.

The macro processor maintains its own list of variables (distinct from
model variables and MATLAB/Octave variables). These macro-variables are
assigned using the `@#define` directive and can be of the following
basic types: boolean, real, string, tuple, function, and array (of any
of the previous types).

### Macro expressions {#macro-exp}

Macro-expressions can be used in two places:

> -   Inside macro directives, directly;
> -   In the body of the `.mod` file, between an at-sign and curly
>     braces (like `@{expr}`): the macro processor will substitute the
>     expression with its value

It is possible to construct macro-expressions that can be assigned to
macro-variables or used within a macro-directive. The expressions are
constructed using literals of the basic types (boolean, real, string,
tuple, array), comprehensions, macro-variables, macro-functions, and
standard operators.

::: {.note}
::: {.admonition-title}
Note
:::

Elsewhere in the manual, MACRO\_EXPRESSION designates an expression
constructed as explained in this section.
:::

**Boolean**

The following operators can be used on booleans:

> -   Comparison operators: `==, !=`
> -   Logical operators: `&&, ||, !`

**Real**

The following operators can be used on reals:

> -   Arithmetic operators: `+, -, *, /, ^`
> -   Comparison operators: `<, >, <=, >=, ==, !=`
> -   Logical operators: `&&, ||, !`
> -   Ranges with an increment of `1`: `REAL1:REAL2` (for example, `1:4`
>     is equivalent to real array `[1, 2, 3, 4]`).
>
>     ::: {.versionchanged}
>     4.6 Previously, putting brackets around the arguments to the colon
>     operator (e.g. `[1:4]`) had no effect. Now, `[1:4]` will create an
>     array containing an array (i.e. `[ [1, 2, 3, 4] ]`).
>     :::
>
> -   Ranges with user-defined increment: `REAL1:REAL2:REAL3` (for
>     example, `6:-2.1:-1` is equivalent to real array
>     `[6, 3.9, 1.8, -0.3]`).
> -   Functions:
>     `max, min, mod, exp, log, log10, sin, cos, tan, asin, acos, atan, sqrt, cbrt, sign, floor, ceil, trunc, erf, erfc, gamma, lgamma, round, normpdf, normcdf`.
>     NB `ln` can be used instead of `log`

**String**

String literals have to be enclosed by **double** quotes (like
`"name"`).

The following operators can be used on strings:

> -   Comparison operators: `<, >, <=, >=, ==, !=`
> -   Concatenation of two strings: `+`
> -   Extraction of substrings: if `s` is a string, then `s[3]` is a
>     string containing only the third character of `s`, and `s[4:6]`
>     contains the characters from 4th to 6th
> -   Function: `length`

**Tuple**

Tuples are enclosed by parenthesis and elements separated by commas
(like `(a,b,c)` or `(1,2,3)`).

The following operators can be used on tuples:

> -   Comparison operators: `==, !=`
> -   Functions: `empty, length`

**Array**

Arrays are enclosed by brackets, and their elements are separated by
commas (like `[1,[2,3],4]` or `["US", "FR"]`).

The following operators can be used on arrays:

> -   Comparison operators: `==, !=`
> -   Dereferencing: if `v` is an array, then `v[2]` is its 2nd element
> -   Concatenation of two arrays: `+`
> -   Set union of two arrays: `|`
> -   Set intersection of two arrays: `&`
> -   Difference `-`: returns the first operand from which the elements
>     of the second operand have been removed.
> -   Cartesian product of two arrays: `*`
> -   Cartesian product of one array N times: `^N`
> -   Extraction of sub-arrays: e.g. `v[4:6]`
> -   Testing membership of an array: `in` operator (for example: `"b"`
>     in `["a", "b", "c"]` returns `1`)
> -   Functions: `empty, sum, length`

**Comprehension**

Comprehension syntax is a shorthand way to make arrays from other
arrays. There are three different ways the comprehension syntax can be
employed: [filtering]{.title-ref}, [mapping]{.title-ref}, and [filtering
and mapping]{.title-ref}.

**Filtering**

Filtering allows one to choose those elements from an array for which a
certain condition hold.

> *Example*
>
> Create a new array, choosing the even numbers from the array `1:5`:
>
>     [ i in 1:5 when mod(i,2) == 0 ]
>
> would result in:
>
>     [2, 4]

**Mapping**

Mapping allows you to apply a transformation to every element of an
array.

> *Example*
>
> Create a new array, squaring all elements of the array `1:5`:
>
>     [ i^2 for i in 1:5 ]
>
> would result in:
>
>     [1, 4, 9, 16, 25]

**Filtering and Mapping**

Combining the two preceding ideas would allow one to apply a
transformation to every selected element of an array.

> *Example*
>
> Create a new array, squaring all even elements of the array `1:5`:
>
>     [ i^2 for i in 1:5 when mod(i,2) == 0]
>
> would result in:
>
>     [4, 16]
>
> *Further Examples* :
>
>     [ (j, i+1) for (i,j) in (1:2)^2 ]
>     [ (j, i+1) for (i,j) in (1:2)*(1:2) when i < j ]
>
> would result in:
>
>     [(1, 2), (2, 2), (1, 3), (2, 3)]
>     [(2, 2)]

**Function**

Functions can be defined in the macro processor using the `@#define`
directive (see below). A function is evaluated at the time it is
invoked, not at define time. Functions can be included in expressions
and the operators that can be combined with them depend on their return
type.

**Checking variable type**

Given a variable name or literal, you can check the type it evaluates to
using the following functions: `isboolean`, `isreal`, `isstring`,
`istuple`, and `isarray`.

> *Examples*

  ----------------------------------
  **Code**              **Output**
  --------------------- ------------
  `isboolean(0)`        `false`

  `isboolean(true)`     `true`

  `isreal("str")`       `false`
  ----------------------------------

**Casting between types**

Variables and literals of one type can be cast into another type. Some
type changes are straightforward (e.g. changing a [real]{.title-ref} to
a [string]{.title-ref}) whereas others have certain requirements (e.g.
to cast an [array]{.title-ref} to a [real]{.title-ref} it must be a one
element array containing a type that can be cast to [real]{.title-ref}).

> *Examples*

  -----------------------------------------------
  **Code**                           **Output**
  ---------------------------------- ------------
  `(bool) -1.1`                      `true`

  `(bool) 0`                         `false`

  `(real) "2.2"`                     `2.2`

  `(tuple) [3.3]`                    `(3.3)`

  `(array) 4.4`                      `[4.4]`

  `(real) [5.5]`                     `5.5`

  `(real) [6.6, 7.7]`                `error`

  `(real) "8.8 in a string"`         `error`
  -----------------------------------------------

Casts can be used in expressions:

> *Examples*

  ----------------------------------------
  **Code**                    **Output**
  --------------------------- ------------
  `(bool) 0 && true`          `false`

  `(real) "1" + 2`            `3`

  `(string) (3 + 4)`          `"7"`

  `(array) 5 + (array) 6`     `[5, 6]`
  ----------------------------------------

### Macro directives

::: {.macrodir}
@\#includepath \"PATH\" @\#includepath MACRO\_EXPRESSION

This directive adds the path contained in PATH to the list of those to
search when looking for a `.mod` file specified by `@#include`. If
provided with a MACRO\_EXPRESSION argument, the argument must evaluate
to a string. Note that these paths are added *after* any paths passed
using `-I <-I\<\<path\>\>>`{.interpreted-text role="opt"}.

*Example*

>     @#includepath "/path/to/folder/containing/modfiles"
>     @#includepath folders_containing_mod_files
:::

::: {.macrodir}
@\#include \"FILENAME\" @\#include MACRO\_EXPRESSION

This directive simply includes the content of another file in its place;
it is exactly equivalent to a copy/paste of the content of the included
file. If provided with a MACRO\_EXPRESSION argument, the argument must
evaluate to a string. Note that it is possible to nest includes (i.e. to
include a file from an included file). The file will be searched for in
the current directory. If it is not found, the file will be searched for
in the folders provided by `-I <-I\<\<path\>\>>`{.interpreted-text
role="opt"} and `@#includepath`.

*Example*

>     @#include "modelcomponent.mod"
>     @#include location_of_modfile
:::

::: {.macrodir}
@\#define MACRO\_VARIABLE @\#define MACRO\_VARIABLE = MACRO\_EXPRESSION
@\#define MACRO\_FUNCTION = MACRO\_EXPRESSION

Defines a macro-variable or macro function.

*Example*

>     @#define var                      // Equals true
>     @#define x = 5                    // Real
>     @#define y = "US"                 // String
>     @#define v = [ 1, 2, 4 ]          // Real array
>     @#define w = [ "US", "EA" ]       // String array
>     @#define u = [ 1, ["EA"] ]        // Mixed array
>     @#define z = 3 + v[2]             // Equals 5
>     @#define t = ("US" in w)          // Equals true
>     @#define f(x) = " " + x + y       // Function `f` with argument `x`
>                                       // returns the string ' ' + x + 'US'

*Example*

>     @#define x = 1
>     @#define y = [ "B", "C" ]
>     @#define i = 2
>     @#define f(x) = x + " + " + y[i]
>     @#define i = 1
>
>     model;
>       A = @{y[i] + f("D")};
>     end;
>
> The latter is strictly equivalent to:
>
>     model;
>       A = BD + B;
>     end;
:::

::: {.macrodir}
@\#if MACRO\_EXPRESSION @\#ifdef MACRO\_VARIABLE @\#ifndef
MACRO\_VARIABLE @\#elseif MACRO\_EXPRESSION @\#else @\#endif

Conditional inclusion of some part of the `.mod` file. The lines between
`@#if`, `@#ifdef`, or `@#ifndef` and the next `@#elseif`, `@#else` or
`@#endif` is executed only if the condition evaluates to `true`.
Following the `@#if` body, you can zero or more `@#elseif` branches. An
`@#elseif` condition is only evaluated if the preceding `@#if` or
`@#elseif` condition evaluated to `false`. The `@#else` branch is
optional and is only evaluated if all `@#if` and `@#elseif` statements
evaluate to false.

Note that when using `@#ifdef`, the condition will evaluate to `true` if
the MACRO\_VARIABLE has been previously defined, regardless of its
value. Conversely, `@#ifndef` will evaluate to true if the
MACRO\_VARIABLE has not yet been defined.

Note that when using `@#elseif` you can check whether or not a variable
has been defined by using the `defined` operator. Hence, to enter the
body of an `@#elseif` branch if the variable `X` has not been defined,
you would write: `@#elseif !defined(X)`.

Note that if a real appears as the result of the MACRO\_EXPRESSION, it
will be interpreted as a boolean; a value of `0` is interpreted as
`false`, otherwise it is interpreted as `true`. Further note that
because of the imprecision of reals, extra care must be taken when
testing them in the MACRO\_EXPRESSION. For example, `exp(log(5)) == 5`
will evaluate to `false`. Hence, when comparing real values, you should
generally use a zero tolerance around the value desired, e.g.
`exp(log(5)) > 5-1e-14 && exp(log(5)) < 5+1e-14`

*Example*

> Choose between two alternative monetary policy rules using a
> macro-variable:
>
>     @#define linear_mon_pol = false // 0 would be treated the same
>     ...
>     model;
>     @#if linear_mon_pol
>       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
>     @#else
>       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
>     @#endif
>     ...
>     end;
>
> This would result in:
>
>     ...
>     model;
>       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
>     ...
>     end;

*Example*

> Choose between two alternative monetary policy rules using a
> macro-variable. The only difference between this example and the
> previous one is the use of `@#ifdef` instead of `@#if`. Even though
> `linear_mon_pol` contains the value `false` because `@#ifdef` only
> checks that the variable has been defined, the linear monetary policy
> is output:
>
>     @#define linear_mon_pol = false // 0 would be treated the same
>     ...
>     model;
>     @#ifdef linear_mon_pol
>       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
>     @#else
>       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
>     @#endif
>     ...
>     end;
>
> This would result in:
>
>     ...
>     model;
>       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
>     ...
>     end;
:::

::: {.macrodir}
@\#for MACRO\_VARIABLE in MACRO\_EXPRESSION @\#for MACRO\_VARIABLE in
MACRO\_EXPRESSION when MACRO\_EXPRESSION @\#for MACRO\_TUPLE in
MACRO\_EXPRESSION @\#for MACRO\_TUPLE in MACRO\_EXPRESSION when
MACRO\_EXPRESSION @\#endfor

Loop construction for replicating portions of the `.mod` file. Note that
this construct can enclose variable/parameters declaration,
computational tasks, but not a model declaration.

*Example*

>     model;
>     @#for country in [ "home", "foreign" ]
>       GDP_@{country} = A * K_@{country}^a * L_@{country}^(1-a);
>     @#endfor
>     end;
>
> The latter is equivalent to:
>
>     model;
>       GDP_home = A * K_home^a * L_home^(1-a);
>       GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
>     end;

*Example*

>     model;
>     @#for (i, j) in ["GDP"] * ["home", "foreign"]
>       @{i}_@{j} = A * K_@{j}^a * L_@{j}^(1-a);
>     @#endfor
>     end;
>
> The latter is equivalent to:
>
>     model;
>       GDP_home = A * K_home^a * L_home^(1-a);
>       GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
>     end;

*Example*

>     @#define countries = ["US", "FR", "JA"]
>     @#define nth_co = "US"
>     model;
>     @#for co in countries when co != nth_co
>       (1+i_@{co}) = (1+i_@{nth_co}) * E_@{co}(+1) / E_@{co};
>     @#endfor
>       E_@{nth_co} = 1;
>     end;
>
> The latter is equivalent to:
>
>     model;
>       (1+i_FR) = (1+i_US) * E_FR(+1) / E_FR;
>       (1+i_JA) = (1+i_US) * E_JA(+1) / E_JA;
>       E_US = 1;
>     end;
:::

::: {.macrodir}
@\#echo MACRO\_EXPRESSION

Asks the preprocessor to display some message on standard output. The
argument must evaluate to a string.
:::

::: {.macrodir}
@\#error MACRO\_EXPRESSION

Asks the preprocessor to display some error message on standard output
and to abort. The argument must evaluate to a string.
:::

::: {.macrodir}
@\#echomacrovars @\#echomacrovars MACRO\_VARIABLE\_LIST
@\#echomacrovars(save) MACRO\_VARIABLE\_LIST

Asks the preprocessor to display the value of all macro variables up
until this point. If the `save` option is passed, then values of the
macro variables are saved to `options_.macrovars_line_<<line_numbers>>`.
If `NAME_LIST` is passed, only display/save variables and functions with
that name.

*Example*

>     @#define A = 1
>     @#define B = 2
>     @#define C(x) = x*2
>     @#echomacrovars A C D
>
> The output of the command above is:
>
>     Macro Variables:
>       A = 1
>     Macro Functions:
>       C(x) = (x * 2)
:::

### Typical usages

#### Modularization

The `@#include` directive can be used to split `.mod` files into several
modular components.

Example setup:

`modeldesc.mod`

> Contains variable declarations, model equations, and shocks
> declarations.

`simul.mod`

> Includes `modeldesc.mod`, calibrates parameter,s and runs stochastic
> simulations.

`estim.mod`

> Includes `modeldesc.mod`, declares priors on parameters, and runs
> Bayesian estimation.

Dynare can be called on `simul.mod` and `estim.mod` but it makes no
sense to run it on `modeldesc.mod`.

The main advantage is that you don\'t have to copy/paste the whole model
(at the beginning) or changes to the model (during development).

#### Indexed sums of products

The following example shows how to construct a moving average:

    @#define window = 2

    var x MA_x;
    ...
    model;
    ...
    MA_x = @{1/(2*window+1)}*(
    @#for i in -window:window
            +x(@{i})
    @#endfor
           );
    ...
    end;

After macro processing, this is equivalent to:

    var x MA_x;
    ...
    model;
    ...
    MA_x = 0.2*(
            +x(-2)
            +x(-1)
            +x(0)
            +x(1)
            +x(2)
           );
    ...
    end;

#### Multi-country models

Here is a skeleton example for a multi-country model:

    @#define countries = [ "US", "EA", "AS", "JP", "RC" ]
    @#define nth_co = "US"

    @#for co in countries
    var Y_@{co} K_@{co} L_@{co} i_@{co} E_@{co} ...;
    parameters a_@{co} ...;
    varexo ...;
    @#endfor

    model;
    @#for co in countries
     Y_@{co} = K_@{co}^a_@{co} * L_@{co}^(1-a_@{co});
    ...
    @#if co != nth_co
     (1+i_@{co}) = (1+i_@{nth_co}) * E_@{co}(+1) / E_@{co}; // UIP relation
    @#else
     E_@{co} = 1;
    @#endif
    @#endfor
    end;

#### Endogeneizing parameters

When calibrating the model, it may be useful to consider a parameter as
an endogenous variable (and vice-versa).

For example, suppose production is defined by a CES function:

> $$y = \left(\alpha^{1/\xi} \ell^{1-1/\xi}+(1-\alpha)^{1/\xi}k^{1-1/\xi}\right)^{\xi/(\xi-1)}$$

and the labor share in GDP is defined as:

> $$\textrm{lab\_rat} = (w \ell)/(p y)$$

In the model, $\alpha$ is a (share) parameter and `lab_rat` is an
endogenous variable.

It is clear that calibrating $\alpha$ is not straightforward; on the
contrary, we have real world data for `lab_rat` and it is clear that
these two variables are economically linked.

The solution is to use a method called *variable flipping*, which
consists in changing the way of computing the steady state. During this
computation, $\alpha$ will be made an endogenous variable and `lab_rat`
will be made a parameter. An economically relevant value will be
calibrated for `lab_rat`, and the solution algorithm will deduce the
implied value for $\alpha$.

An implementation could consist of the following files:

`modeqs.mod`

> This file contains variable declarations and model equations. The code
> for the declaration of $\alpha$ and `lab_rat` would look like:
>
>     @#if steady
>       var alpha;
>       parameter lab_rat;
>     @#else
>       parameter alpha;
>       var lab_rat;
>     @#endif

`steady.mod`

> This file computes the steady state. It begins with:
>
>     @#define steady = 1
>     @#include "modeqs.mod"
>
> Then it initializes parameters (including `lab_rat`, excluding
> $\alpha$), computes the steady state (using guess values for
> endogenous, including $\alpha$), then saves values of parameters and
> endogenous at steady state in a file, using the
> `save_params_and_steady_state` command.

`simul.mod`

> This file computes the simulation. It begins with:
>
>     @#define steady = 0
>     @#include "modeqs.mod"
>
> Then it loads values of parameters and endogenous at steady state from
> file, using the `load_params_and_steady_state` command, and computes
> the simulations.

### MATLAB/Octave loops versus macro processor loops

Suppose you have a model with a parameter $\rho$ and you want to run
simulations for three values: $\rho = 0.8, 0.9,
1$. There are several ways of doing this:

*With a MATLAB/Octave loop*

>     rhos = [ 0.8, 0.9, 1];
>     for i = 1:length(rhos)
>       rho = rhos(i);
>       stoch_simul(order=1);
>     end
>
> Here the loop is not unrolled, MATLAB/Octave manages the iterations.
> This is interesting when there are a lot of iterations.

*With a macro processor loop (case 1)*

>     rhos = [ 0.8, 0.9, 1];
>     @#for i in 1:3
>       rho = rhos(@{i});
>       stoch_simul(order=1);
>     @#endfor
>
> This is very similar to the previous example, except that the loop is
> unrolled. The macro processor manages the loop index but not the data
> array (`rhos`).

*With a macro processor loop (case 2)*

>     @#for rho_val in [ 0.8, 0.9, 1]
>       rho = @{rho_val};
>       stoch_simul(order=1);
>     @#endfor
>
> The advantage of this method is that it uses a shorter syntax, since
> the list of values is directly given in the loop construct. The
> inconvenience is that you can not reuse the macro array in
> MATLAB/Octave.

Verbatim inclusion {#verbatim}
------------------

Pass everything contained within the verbatim block to the
`<mod_file>.m` file.

::: {.block}
verbatim ;

By default, whenever Dynare encounters code that is not understood by
the parser, it is directly passed to the preprocessor output.

In order to force this behavior you can use the `verbatim` block. This
is useful when the code you want passed to the `<mod_file>.m` file
contains tokens recognized by the Dynare preprocessor.

*Example*

>     verbatim;
>     % Anything contained in this block will be passed
>     % directly to the <modfile>.m file, including comments
>     var = 1;
>     end;
:::

Misc commands
-------------

::: {.command}
set\_dynare\_seed (INTEGER) set\_dynare\_seed (\`default\')
set\_dynare\_seed (\`clock\') set\_dynare\_seed (\`reset\')
set\_dynare\_seed (\`ALGORITHM\', INTEGER)

Sets the seed used for random number generation. It is possible to set a
given integer value, to use a default value, or to use the clock (by
using the latter, one will therefore get different results across
different Dynare runs). The `reset` option serves to reset the seed to
the value set by the last `set_dynare_seed` command. On MATLAB 7.8 or
above, it is also possible to choose a specific algorithm for random
number generation; accepted values are `mcg16807`, `mlfg6331_64`,
`mrg32k3a`, `mt19937ar` (the default), `shr3cong` and `swb2712`.
:::

::: {.command}
save\_params\_and\_steady\_state (FILENAME);

For all parameters, endogenous and exogenous variables, stores their
value in a text file, using a simple name/value associative table.

> -   for parameters, the value is taken from the last parameter
>     initialization.
> -   for exogenous, the value is taken from the last `initval` block.
> -   for endogenous, the value is taken from the last steady state
>     computation (or, if no steady state has been computed, from the
>     last `initval` block).

Note that no variable type is stored in the file, so that the values can
be reloaded with `load_params_and_steady_state` in a setup where the
variable types are different.

The typical usage of this function is to compute the steady-state of a
model by calibrating the steady-state value of some endogenous variables
(which implies that some parameters must be endogeneized during the
steady-state computation).

You would then write a first `.mod` file which computes the steady state
and saves the result of the computation at the end of the file, using
`save_params_and_steady_state`.

In a second file designed to perform the actual simulations, you would
use `load_params_and_steady_state` just after your variable
declarations, in order to load the steady state previously computed
(including the parameters which had been endogeneized during the steady
state computation).

The need for two separate `.mod` files arises from the fact that the
variable declarations differ between the files for steady state
calibration and for simulation (the set of endogenous and parameters
differ between the two); this leads to different `var` and `parameters`
statements.

Also note that you can take advantage of the `@#include` directive to
share the model equations between the two files (see
`macro-proc-lang`{.interpreted-text role="ref"}).
:::

::: {.command}
load\_params\_and\_steady\_state (FILENAME);

For all parameters, endogenous and exogenous variables, loads their
value from a file created with `save_params_and_steady_state`.

> -   for parameters, their value will be initialized as if they had
>     been calibrated in the `.mod` file.
> -   for endogenous and exogenous variables, their value will be
>     initialized as they would have been from an `initval` block .

This function is used in conjunction with
`save_params_and_steady_state`; see the documentation of that function
for more information.
:::

::: {.command}
compilation\_setup (OPTIONS);

When the `use_dll`{.interpreted-text role="opt"} option is present,
Dynare uses the GCC compiler that was distributed with it to compile the
static and dynamic C files produced by the preprocessor. You can use
this option to change the compiler, flags, and libraries used.

*Options*

> ::: {.option}
> compiler = FILENAME
>
> The path to the compiler.
> :::
>
> ::: {.option}
> substitute\_flags = QUOTED\_STRING
>
> The flags to use instead of the default flags.
> :::
>
> ::: {.option}
> add\_flags = QUOTED\_STRING
>
> The flags to use in addition to the default flags. If
> `substitute_flags` is passed, these flags are added to the flags
> specified there.
> :::
>
> ::: {.option}
> substitute\_libs = QUOTED\_STRING
>
> The libraries to link against instead of the default libraries.
> :::
>
> ::: {.option}
> add\_libs = QUOTED\_STRING
>
> The libraries to link against in addition to the default libraries. If
> `substitute_libs` is passed, these libraries are added to the
> libraries specified there.
> :::
:::

::: {.matcomm}
dynare\_version ;

Output the version of Dynare that is currently being used (i.e. the one
that is highest on the MATLAB/Octave path).
:::

::: {.matcomm}
write\_latex\_definitions ;

Writes the names, LaTeX names and long names of model variables to
tables in a file named `<<M_.fname>>_latex_definitions.tex`. Requires
the following LaTeX packages: `longtable`.
:::

::: {.matcomm}
write\_latex\_parameter\_table ;

Writes the LaTeX names, parameter names, and long names of model
parameters to a table in a file named
`<<M_.fname>>_latex_parameters.tex.` The command writes the values of
the parameters currently stored. Thus, if parameters are set or changed
in the steady state computation, the command should be called after a
steady-command to make sure the parameters were correctly updated. The
long names can be used to add parameter descriptions. Requires the
following LaTeX packages: `longtable, booktabs`.
:::

::: {.matcomm}
write\_latex\_prior\_table ;

Writes descriptive statistics about the prior distribution to a LaTeX
table in a file named `<<M_.fname>>_latex_priors_table.tex`. The command
writes the prior definitions currently stored. Thus, this command must
be invoked after the `estimated_params` block. If priors are defined
over the measurement errors, the command must also be preceeded by the
declaration of the observed variables (with `varobs`). The command
displays a warning if no prior densities are defined (ML estimation) or
if the declaration of the observed variables is missing. Requires the
following LaTeX packages: `longtable, booktabs`.
:::

::: {.matcomm}
collect\_latex\_files ;

Writes a LaTeX file named `<<M_.fname>>_TeX_binder.tex` that collects
all TeX output generated by Dynare into a file. This file can be
compiled using `pdflatex` and automatically tries to load all required
packages. Requires the following LaTeX packages: `breqn`, `psfrag`,
`graphicx`, `epstopdf`, `longtable`, `booktabs`, `caption`, `float,`
`amsmath`, `amsfonts`, and `morefloats`.
:::

**Footnotes**

[^1]: A `.mod` file must have lines that end with a line feed character,
    which is not commonly visible in text editors. Files created on
    Windows and Unix-based systems have always conformed to this
    requirement, as have files created on OS X and macOS. Files created
    on old, pre-OS X Macs used carriage returns as end of line
    characters. If you get a Dynare parsing error of the form
    `ERROR: <<mod file>>: line 1, cols 341-347: syntax error,...` and
    there\'s more than one line in your `.mod` file, know that it uses
    the carriage return as an end of line character. To get more helpful
    error messages, the carriage returns should be changed to line
    feeds.

[^2]: Note that arbitrary MATLAB or Octave expressions can be put in a
    `.mod` file, but those expressions have to be on separate lines,
    generally at the end of the file for post-processing purposes. They
    are not interpreted by Dynare, and are simply passed on unmodified
    to MATLAB or Octave. Those constructions are not addresses in this
    section.

[^3]: In particular, for big models, the compilation step can be very
    time-consuming, and use of this option may be counter-productive in
    those cases.

[^4]: See options `conf_sig <confsig>`{.interpreted-text role="ref"} and
    `mh_conf_sig <mh_conf_sig = DOUBLE>`{.interpreted-text role="opt"}
    to change the size of the HPD interval.

[^5]: See options `conf_sig <confsig>`{.interpreted-text role="ref"} ()
    and `mh_conf_sig <mh_conf_sig = DOUBLE>`{.interpreted-text
    role="opt"} to change the size of the HPD interval.

[^6]: When the shocks are correlated, it is the decomposition of
    orthogonalized shocks via Cholesky decomposition according to the
    order of declaration of shocks (see `var-decl`{.interpreted-text
    role="ref"})

[^7]: See `forecast <forecast = INTEGER>`{.interpreted-text role="opt"}
    for more information.

[^8]: In case of Excel not being installed,
    <https://mathworks.com/matlabcentral/fileexchange/38591-xlwrite--generate-xls-x--files-without-excel-on-mac-linux-win>
    may be helpful.

[^9]: In case of Excel not being installed,
    <https://mathworks.com/matlabcentral/fileexchange/38591-xlwrite--generate-xls-x--files-without-excel-on-mac-linux-win>
    may be helpful.

[^10]: See option `conf_sig <confsig>`{.interpreted-text role="ref"} to
    change the size of the HPD interval.

[^11]: See option `conf_sig <confsig>`{.interpreted-text role="ref"} to
    change the size of the HPD interval.

[^12]: If you want to align the paper with the description herein,
    please note that $A$ is $A^0$ and $F$ is $A^+$.

[^13]: An example can be found at
    <https://git.dynare.org/Dynare/dynare/blob/master/tests/ms-dsge/test_ms_dsge.mod>.
