# The model file

## Conventions

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

    // This is a single line comment

    var x; // This is a comment about x
    
    /* This is another inline comment about alpha */  alpha = 0.3;
    
     /*
      This comment is spanning
      two lines.
     */

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
    '[]';
-   Repeated arguments are indicated by ellipses: "...";
-   Mutually exclusive arguments are separated by vertical bars: '|';
-   INTEGER indicates an integer number;
-   INTEGER_VECTOR indicates a vector of integer numbers separated by
    spaces, enclosed by square brackets;
-   DOUBLE indicates a double precision number. The following syntaxes
    are valid: `1.1e3`, `1.1E3`, `1.1d3`, `1.1D3`. In some places,
    infinite Values `Inf` and `-Inf` are also allowed;
-   NUMERICAL_VECTOR indicates a vector of numbers separated by spaces,
    enclosed by square brackets;
-   EXPRESSION indicates a mathematical expression valid outside the
    model description (see `expr`{.interpreted-text role="ref"});
-   MODEL_EXPRESSION (sometimes MODEL_EXP) indicates a mathematical
    expression valid in the model description (see
    `expr`{.interpreted-text role="ref"} and
    `model-decl`{.interpreted-text role="ref"});
-   MACRO_EXPRESSION designates an expression of the macro processor
    (see `macro-exp`{.interpreted-text role="ref"});
-   VARIABLE_NAME (sometimes VAR_NAME) indicates a variable name
    starting with an alphabetical character and can't contain:
    '()+-*/\^=!;:@#.' or accentuated characters;
-   PARAMETER_NAME (sometimes PARAM_NAME) indicates a parameter name
    starting with an alphabetical character and can't contain:
    '()+-*/\^=!;:@#.' or accentuated characters;
-   LATEX_NAME (sometimes TEX_NAME) indicates a valid LaTeX expression
    in math mode (not including the dollar signs);
-   FUNCTION_NAME indicates a valid MATLAB function name;
-   FILENAME indicates a filename valid in the underlying operating
    system; it is necessary to put it between quotes when specifying the
    extension or if the filename contains a non-alphanumeric character;
-   QUOTED_STRING indicates an arbitrary string enclosed between
    (single) quotes.

## Variable declarations

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

*command*

```
var VAR_NAME [$TEX_NAME$][(long_name=QUOTED_STRINGNAME=QUOTED_STRING)]...;
```

```
var (deflator=MODEL_EXPR) VAR_NAME (... same options apply) var(log,deflator=MODEL_EXPR) VAR_NAME (... same options apply)
```

```
var (log_deflator=MODEL_EXPR) VAR_NAME (... same options apply)
```

This required command declares the endogenous variables in the model.
See `conv`{.interpreted-text role="ref"} for the syntax of *VAR\_NAME*
and *MODEL\_EXPR*. Optionally it is possible to give a LaTeX name to the
variable or, if it is nonstationary, provide information regarding its
deflator. The variables in the list can be separated by spaces or by
commas. `var` commands can appear several times in the file and Dynare
will concatenate them. Dynare stores the list of declared parameters, in
the order of declaration, in a column cell array `M_.endo_names`.

*Options*

If the model is nonstationary and is to be written as such in the
`model` block, Dynare will need the trend deflator for the appropriate
endogenous variables in order to stationarize the model. The trend
deflator must be provided alongside the variables that follow this
trend.

*Options*

- `log`

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


- `deflator = MODEL_EXPR`

    The expression used to detrend an endogenous variable. All trend
    variables, endogenous variables and parameters referenced in MODEL\_EXPR
    must already have been declared by the `trend_var, log_trend_var, var`
    and `parameters` commands. The deflator is assumed to be multiplicative;
    for an additive deflator, use `log_deflator`. This option can be used
    together with the `log` option (the latter must come first).


- `log_deflator = MODEL_EXPR`

    Same as `deflator`, except that the deflator is assumed to be additive
    instead of multiplicative (or, to put it otherwise, the declared
    variable is equal to the log of a variable with a multiplicative trend).
    This option cannot be used together with the `log` option, because it
    would not make much sense from an economic point of view (the
    corresponding auxiliary variable would correspond to the log taken two
    times on a variable with a multiplicative trend).


- `long_name = QUOTED_STRING`

    This is the long version of the variable name. Its value is stored in
    `M_.endo_names_long` (a column cell array, in the same order as
    `M_.endo_names`). In case multiple `long_name` options are provided, the
    last one will be used. Default: `VAR_NAME`.


- `NAME = QUOTED_STRING`

    This is used to create a partitioning of variables. It results in the
    direct output in the `.m` file analogous to:
    `M_.endo_partitions.NAME = QUOTED_STRING`;.
    

*Example (variable partitioning)*

```
var c gnp cva (country=`US', state=`VA')
          cca (country=`US', state=`CA', long_name=`Consumption CA');
var(deflator=A) i b;
var c $C$ (long_name=`Consumption');
```

*Command*:

```
varexo_det VAR_NAME [$TEX_NAME$][(long_name=QUOTED_STRING|NAME=QUOTED_STRING)...];
```

This optional command declares exogenous deterministic variables in a
stochastic model. See `conv`{.interpreted-text role="ref"} for the
syntax of `VAR_NAME`. Optionally it is possible to give a LaTeX name
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

- `long_name = QUOTED_STRING`

Like `long_name <long-name>`{.interpreted-text role="ref"} but value
stored in `M_.exo_det_names_long`.

- NAME = QUOTED_STRING

Like `partitioning <partitioning>`{.interpreted-text role="ref"} but
QUOTED_STRING stored in `M_.exo_det_partitions.NAME`.

*Example*


```
varexo m gov;
varexo_det tau;
```

*Command*:

```
var_remove VAR_NAME | PARAM_NAME...;
```

Removes the listed variables (or parameters) from the model. Removing a
variable that has already been used in a model equation or elsewhere
will lead to an error.

*Command*:

```
predetermined_variables VAR_NAME...;
```

In Dynare, the default convention is that the timing of a variable
reflects when this variable is decided. The typical example is for
capital stock: since the capital stock used at current period is
actually decided at the previous period, then the capital stock entering
the production function is `k(-1)`, and the law of motion of capital
must be written:
```
k = i + (1-delta)*k(-1)
```

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

The following two program snippets are strictly equivalent.

Using default Dynare timing convention:

```
var y, k, i;
...
model;
y = k(-1)^alpha;
k = i + (1-delta)*k(-1);
...
end;
```

Using the alternative timing convention:

```
var y, k, i;
predetermined_variables k;
...
model;
y = k^alpha;
k(+1) = i + (1-delta)*k;
...
end;
```

*command*:

```
trend_var (growth_factor = MODEL_EXPR) VAR_NAME[$LATEX_NAME$]...;
```

This optional command declares the trend variables in the model. See
`conv`{.interpreted-text role="ref"} for the syntax of `MODEL_EXPR` and
`VAR_NAME`. Optionally it is possible to give a LaTeX name to the
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
referenced in `MODEL_EXPR` must already have been declared by the var and
parameters commands.

*Example*

```
trend_var (growth_factor=gA) A;
```

*command*:

```
model_local_variable VARIABLE_NAME [LATEX_NAME]... ;
```

This optional command declares a model local variable. See
`conv`{.interpreted-text role="ref"} for the syntax of `VARIABLE_NAME`.
As you can create model local variables on the fly in the model block
(see `model-decl`{.interpreted-text role="ref"}), the interest of this
command is primarily to assign a `LATEX_NAME` to the model local
variable.

*Example*

```
model_local_variable GDP_US $GDPUS$;
```

## On-the-fly Model Variable Declaration

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

The following two snippets are equivalent:

```
model;
  [endogenous='k',name='law of motion of capital']
  k(+1) = i|e + (1-delta|p)*k;
  y|e = k^alpha|p;
  ...
end;
delta = 0.025;
alpha = 0.36;
```

```
var k, i, y;
parameters delta, alpha;
delta = 0.025;
alpha = 0.36;
...
model;
  [name='law of motion of capital']
  k(1) = i|e + (1-delta|p)*k;
  y|e = k|e^alpha|p;
  ...
end;
```

## Expressions

Dynare distinguishes between two types of mathematical expressions:
those that are used to describe the model, and those that are used
outside the model block (e.g. for initializing parameters or variables,
or as command options). In this manual, those two types of expressions
are respectively denoted by `MODEL_EXPRESSION` and `EXPRESSION`.

Unlike MATLAB or Octave expressions, Dynare expressions are necessarily
scalar ones: they cannot contain matrices or evaluate to matrices.[^2]

Expressions can be constructed using integers (INTEGER), floating point
numbers (`DOUBLE`), parameter names (`PARAMETER_NAME`), variable names
(`VARIABLE_NAME`), operators and functions.

The following special constants are also accepted in some contexts:

*Constant*: `inf`

Represents infinity.

*Constant*: `nan`

"Not a number": represents an undefined or unrepresentable value.

### Parameters and variables

Parameters and variables can be introduced in expressions by simply
typing their names. The semantics of parameters and variables is quite
different whether they are used inside or outside the model block.

#### Inside the model

Parameters used inside the model refer to the value given through
parameter initialization (see `param-init`{.interpreted-text
role="ref"}) or `homotopy_setup` when doing a simulation, or are the
estimated variables when doing an estimation.

Variables used in a `MODEL_EXPRESSION` denote current period values when
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

*Operator*: **STEADY_STATE (MODEL_EXPRESSION)**

This operator is used to take the value of the enclosed expression at
the steady state. A typical usage is in the Taylor rule, where you may
want to use the value of GDP at steady state to compute the output gap.

Exogenous and exogenous deterministic variables may not appear in
`MODEL_EXPRESSION`.

!!! note

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

*Operator*: `EXPECTATION (INTEGER) (MODEL_EXPRESSION`

This operator is used to take the expectation of some expression using a
different information set than the information available at current
period. For example, `EXPECTATION(-1)(x(+1))` is equal to the expected
value of variable x at next period, using the information set available
at the previous period. See `aux-variables`{.interpreted-text
role="ref"} for an explanation of how this operator is handled
internally and how this affects the output.

### Functions

#### Built-in functions

The following standard functions are supported internally for both
`MODEL_EXPRESSION` and `EXPRESSION`:

*Function*: `exp(x)`

Natural exponential.

*Function*: `log(x)`

*Function*: `ln(x)`

Natural logarithm.

*Function*: `log10(x)`

Base 10 logarithm.

*Function*: `sqrt(x)`

Square root.

*Function*: `cbrt(x)`

Cube root.

*Function*: `sign(x)`

Signum function, defined as:

```
$$\begin{aligned}
  \textrm{sign}(x) =
         \begin{cases}
         -1 &\quad\text{if }x<0\\
         0 &\quad\text{if }x=0\\
         1 &\quad\text{if }x>0
         \end{cases}
  \end{aligned}$$
```

Note that this function is not continuous, hence not differentiable, at
$x=0$. However, for facilitating convergence of Newton-type methods,
Dynare assumes that the derivative at $x=0$ is equal to $0$. This
assumption comes from the observation that both the right- and
left-derivatives at this point exist and are equal to $0$, so we can
remove the singularity by postulating that the derivative at $x=0$ is
$0$.

*Function*: `abs(x)`

Absolute value.

Note that this continuous function is not differentiable at $x=0$.
However, for facilitating convergence of Newton-type methods, Dynare
assumes that the derivative at $x=0$ is equal to $0$ (even if the
derivative does not exist). The rational for this mathematically
unfounded definition, rely on the observation that the derivative of
$\mathrm{abs}(x)$ is equal to $\mathrm{sign}(x)$ for any $x\neq 0$ in
$\mathbb R$ and from the convention for the value of $\mathrm{sign}(x)$
at $x=0$).

*Function*: `sin(x)`

*Function*: `cos(x)`

*Function*: `tan(x)`

*Function*: `asin(x)`

*Function*: `acos(x)`

*Function*: `atan(x)`

Trigonometric functions.

*Function*: `sinh(x)`

*Function*: `cosh(x)`

*Function*: `tanh(x)`

*Function*: `asinh(x)`

*Function*: `acosh(x)`

*Function*: `atanh(x)`

Hyperbolic functions.

*Function*: `max(a, b)`

*Function*: `min(a, b)`

Maximum and minimum of two reals.

Note that these functions are differentiable everywhere except on a line
of the 2-dimensional real plane defined by $a=b$. However for
facilitating convergence of Newton-type methods, Dynare assumes that, at
the points of non-differentiability, the partial derivative of these
functions with respect to the first (resp. the second) argument is equal
to $1$ (resp. to $0$) (i.e. the derivatives at the kink are equal to the
derivatives observed on the half-plane where the function is equal to
its first argument).

*Function*: `normcdf(x)` 

*Function*: `normcdf(x, mu, sigma)`

Gaussian cumulative density function, with mean *mu* and standard
deviation *sigma*. Note that `normcdf(x)` is equivalent to
`normcdf(x,0,1)`.

*Function*: `normpdf(x)` 
*Function*: `normpdf(x, mu, sigma)`

Gaussian probability density function, with mean *mu* and standard
deviation *sigma*. Note that `normpdf(x)` is equivalent to
`normpdf(x,0,1)`.

*Function*: `erf(x)`

Gauss error function.

*Function*: `erfc(x)`

Complementary error function, *i.e.*
$\mathrm{erfc}(x) = 1-\mathrm{erf}(x)$.

#### External functions

Any other user-defined (or built-in) MATLAB or Octave function may be
used in both a `MODEL_EXPRESSION` and an `EXPRESSION`, provided that this
function has a scalar argument as a return value.

To use an external function in a `MODEL_EXPRESSION`, one must declare the
function using the `external_function` statement. This is not required
for external functions used in an EXPRESSION outside of a `model` block
or `steady_state_model` block.

*Command*: `external_function (OPTIONS...);`

This command declares the external functions used in the model block. It
is required for every unique function used in the model block.

`external_function` commands can appear several times in the file and
must come before the model block.

*Options*

- `name = NAME`

The name of the function, which must also be the name of the Julia file
implementing it. This option is mandatory.

- `nargs = INTEGER`

The number of arguments of the function. If this option is not provided,
Dynare assumes `nargs = 1`.

- `first_deriv_provided [= NAME]`

If NAME is provided, this tells Dynare that the Jacobian is provided as
the only output of the Julia file given as the option argument. If NAME
is not provided, this tells Dynare that the Julia file specified by the
argument passed to NAME returns the Jacobian as its second output
argument. When this option is not provided, Dynare will use finite
difference approximations for computing the derivatives of the function,
whenever needed.

- `second_deriv_provided [= NAME]`

If NAME is provided, this tells Dynare that the Hessian is provided as
the only output of the Julia file given as the option argument. If NAME
is not provided, this tells Dynare that the Julia file specified by the
argument passed to NAME returns the Hessian as its third output
argument. NB: This option can only be used if the `first_deriv_provided`
option is used in the same `external_function` command. When this option
is not provided, Dynare will use finite difference approximations for
computing the Hessian derivatives of the function, whenever needed.

*Example*

```
external_function(name = funcname);
external_function(name = otherfuncname, nargs = 2, first_deriv_provided, second_deriv_provided);
external_function(name = yetotherfuncname, nargs = 3, first_deriv_provided = funcname_deriv);
```

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

## Parameter initialization

When using Dynare for computing simulations, it is necessary to
calibrate the parameters of the model. This is done through parameter
initialization.

The syntax is the following:

```
PARAMETER_NAME = EXPRESSION;
```

Here is an example of calibration:

```
parameters alpha, beta;

beta = 0.99;
alpha = 0.36;
A = 1-alpha*beta;
```

Internally, the parameter values are stored in `M_.params`:

- *MATLAB/Octave variable*: `M_.params`

Contains the values of model parameters. The parameters are in the order
that was used in the `parameters` command, hence ordered as in
`M_.param_names`.

The parameter names are stored in `M_.param_names`:

- *MATLAB/Octave variable*: `M_.param_names`

Cell array containing the names of the model parameters.

- *MATLAB/Octave variable*: `get_param_by_name ('PARAMETER_NAME')`;

Given the name of a parameter, returns its calibrated value as it is
stored in `M_.params`.

- *MATLAB/Octave variable*: `set_param_value ('PARAMETER_NAME', MATLAB_EXPRESSION);`

Sets the calibrated value of a parameter to the provided expression.
This does essentially the same as the parameter initialization syntax
described above, except that it accepts arbitrary MATLAB/Octave
expressions, and that it works from MATLAB/Octave scripts.

## Model declaration

The model is declared inside a `model` block:

*Block*: `model ;` 

*Block*: `model (OPTIONS...);`

The equations of the model are written in a block delimited by `model`
and `end` keywords.

There must be as many equations as there are endogenous variables in the
model, except when computing the unconstrained optimal policy with
`ramsey_model`, `ramsey_policy` or `discretionary_policy`.

The syntax of equations must follow the conventions for
`MODEL_EXPRESSION` as described in `expr`{.interpreted-text role="ref"}.
Each equation must be terminated by a semicolon (';'). A normal equation
looks like:

```
MODEL_EXPRESSION = MODEL_EXPRESSION;
```

When the equations are written in homogenous form, it is possible to
omit the '=0' part and write only the left hand side of the equation. A
homogenous equation looks like:

``` 
MODEL_EXPRESSION;
```

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

```
#VARIABLE_NAME = MODEL_EXPRESSION;
```

It is possible to tag equations written in the model block. A tag can
serve different purposes by allowing the user to attach arbitrary
informations to each equation and to recover them at runtime. For
instance, it is possible to name the equations with a `name`-tag, using
a syntax like:

```
model;

[name = 'Budget constraint'];
c + k = k^theta*A;

end;
```

Here, `name` is the keyword indicating that the tag names the equation.
If an equation of the model is tagged with a name, the `resid` command
will display the name of the equations (which may be more informative
than the equation numbers) in addition to the equation number. Several
tags for one equation can be separated using a comma:

```   
model;

[name='Taylor rule',mcp = 'r > -1.94478']
r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;

end;
```

More information on tags is available at
[https://git.dynare.org/Dynare/dynare/-/wikis/Equations-Tags](https://git.dynare.org/Dynare/dynare/-/wikis/Equations-Tags).

There can be several `model` blocks, in which case they are simply
concatenated. The set of effective options is also the concatenation of
the options declared in all the blocks, but in that case you may rather
want to use the `model_options`{.interpreted-text role="comm"} command.

*Options*

- `linear`

Declares the model as being linear. It spares oneself from having to
declare initial values for computing the steady state of a stationary
linear model. This option can't be used with non-linear models, it will
NOT trigger linearization of the model.

- `no_static`

Don't create the static model file. This can be useful for models which
don't have a steady state.

- `differentiate_forward_vars` `differentiate_forward_vars = (VARIABLE_NAME [VARIABLE_NAME ...] )`

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

- `parallel_local_files = ( FILENAME [, FILENAME]... )`

Declares a list of extra files that should be transferred to follower
nodes when doing a parallel computation (see
`paral-conf`{.interpreted-text role="ref"}).

- `balanced_growth_test_tol = DOUBLE`

Tolerance used for determining whether cross-derivatives are zero in the
test for balanced growth path (the latter is documented on
[https://archives.dynare.org/DynareWiki/RemovingTrends](https://archives.dynare.org/DynareWiki/RemovingTrends)). Default:
`1e-6`

*Example* (Elementary RBC model)

```
var c k;
varexo x;
parameters aa alph bet delt gam;

model;
c =  - k + aa*x*k(-1)^alph + (1-delt)*k(-1);
c^(-gam) = (aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam)/(1+bet);
end;
```

*Example* (Use of model local variables)

The following program:

```
model;
# gamma = 1 - 1/sigma;
u1 = c1^gamma/gamma;
u2 = c2^gamma/gamma;
end;
```

...is formally equivalent to:

```
model;
u1 = c1^(1-1/sigma)/(1-1/sigma);
u2 = c2^(1-1/sigma)/(1-1/sigma);
end;
```

*Example* (A linear model)

```
model(linear);
x = a*x(-1)+b*y(+1)+e_x;
y = d*y(-1)+e_y;
end;
```

- `model_options (OPTIONS...);`

This command accepts the same options as the `model`{.interpreted-text
role="bck"} block.

The purpose of this statement is to specify the options that apply to
the whole model, when there are several `model` blocks, so as to restore
the symmetry between those blocks (since otherwise one `model` block
would typically bear the options, while the other ones would typically
have no option).

- `model_remove (TAGS...);`

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

```
var c k dummy1 dummy2;

model;
  c + k - aa*x*k(-1)^alph - (1-delt)*k(-1) + dummy1;
  c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
  [ name = 'eq:dummy1', endogenous = 'dummy1' ]
  c*k = dummy1;
  [ foo = 'eq:dummy2' ]
  log(dummy2) = k + 2;
end;

  model_remove('eq:dummy1', foo = 'eq:dummy2');
```
 
In the above example, the last two equations will be removed, `dummy1`
will be turned into an exogenous, and `dummy2` will be removed.

- *block*: `model_replace (TAGS...);`

This block replaces several equations in the model. It removes the
equations given by the tags list (with the same syntax as in
`model_remove`{.interpreted-text role="comm"}), and it adds equations
given within the block (with the same syntax as
`model`{.interpreted-text role="bck"}).

No variable is removed or has its type changed in the process.

*Example*

```
var c k;

model;
  c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
  [ name = 'dummy' ]
  c*k = 1;
end;

model_replace('dummy');
  c^(-gam) = (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;
```
 
In the above example, the dummy equation is replaced by a proper Euler
equation.

Dynare has the ability to output the original list of model equations to
a LaTeX file, using the `write_latex_original_model` command, the list
of transformed model equations using the
`write_latex_dynamic_model command`, and the list of static model
equations using the `write_latex_static_model` command.

- *Command*: `write_latex_original_model (OPTIONS);`

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

- `write\_equation\_tags`

Write the equation tags in the LaTeX output. The equation tags will be
interpreted with LaTeX markups.

- *Command*: `write_latex_dynamic_model ; 

- *Command*: `write_latex_dynamic_model (OPTIONS);`

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

-   The timing convention of predetermined variables (see
    `predetermined_variables`{.interpreted-text role="comm"}) will
    have been changed to the default Dynare timing convention; in
    other words, variables declared as predetermined will be lagged on
    period back,
-   The `EXPECTATION` operators will have been removed, replaced by
    auxiliary variables and new equations (as explained in the
    documentation of
    `EXPECTATION <EXPECTATION (INTEGER) (MODEL_EXPRESSION)>`{.interpreted-text
    role="op"}),
-   Endogenous variables with leads or lags greater or equal than two
    will have been removed, replaced by new auxiliary variables and
    equations,
-   Exogenous variables with leads or lags will also have been
    replaced by new auxiliary variables and equations.

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.

*Options*

- *write\_equation\_tags*

See `write_equation_tags`{.interpreted-text role="opt"}

- *Command*: `write_latex_static_model (OPTIONS);`

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

- *write_equation_tags*: 

See `write_equation_tags`{.interpreted-text role="opt"}.

- *Command*: `write_latex_steady_state_model`

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

## Auxiliary variables

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

See [https://git.dynare.org/Dynare/dynare/-/wikis/Auxiliary-variables](https://git.dynare.org/Dynare/dynare/-/wikis/Auxiliary-variables)
for more technical details on auxiliary variables.

## Initial and terminal conditions

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

- *Block*: `initval ;` 

- *Block*: `initval(OPTIONS...);`

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

```
VARIABLE_NAME = EXPRESSION;
```

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
(`stoch_simul`, `estimation`...).

As such, `initval` allows specifying the initial instrument value for
steady state finding when providing an analytical conditional steady
state file for `ramsey_model`-computations.

It is not necessary to declare 0 as initial value for exogenous
stochastic variables, since it is the only possible value.

The subsequently computed steady state (not the initial values, use
histval for this) will be used as the initial condition at all the
periods preceeding the first simulation period for the three possible
types of simulations in stochastic mode:

-   `stoch_simul`{.interpreted-text role="comm"}, if the `periods`
    option is specified.
-   `forecast`{.interpreted-text role="comm"} as the initial point at
    which the forecasts are computed.
-   `conditional_forecast`{.interpreted-text role="comm"} as the
    initial point at which the conditional forecasts are computed.

To start simulations at a particular set of starting values that are not
a computed steady state, use `histval`{.interpreted-text role="bck"}.

*Options*

- `all_values_required`

Issues an error and stops processing the .mod file if there is at least
one endogenous or exogenous variable that has not been set in the
initval block.

*Example*

```
initval;
   c = 1.2;
   k = 12;
   x = 1;
end;

steady;
```

- *Block*: `endval ;` 

- *Block*: `endval (OPTIONS...);`

This block is terminated by `end;` and contains lines of the form:

```
VARIABLE_NAME = EXPRESSION;
```

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

- `all_values_required`

See `all_values_required`{.interpreted-text role="opt"}.

*Example*

```
var c k;
varexo x;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
c = 1.2;
k = 12;
x = 1;
end;

steady;

endval;
c = 2;
k = 20;
x = 2;
end;

steady;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;
```


In this example, the problem is finding the optimal path for
consumption and capital for the periods $t=1$ to $T=200$, given the
path of the exogenous technology level `x`. `c` is a forward-looking
variable and the exogenous variable `x` appears with a lead in the
expected return of physical capital, while `k` is a purely
backward-looking (state) variable.

The initial equilibrium is computed by `steady` conditional on `x=1`,
and the terminal one conditional on `x=2`. The `initval` block sets
the initial condition for `k` (since it is the only backward-looking
variable), while the `endval` block sets the terminal condition for
`c` (since it is the only forward-looking endogenous variable). The
starting values for the perfect foresight solver are given by the
`endval` block. See below for more details.

*Example*

```
var c k;
varexo x;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
k = 12;
end;

endval;
c = 2;
x = 1.1;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;
```

In this example, there is no [steady]{.title-ref} command, hence the
conditions are exactly those specified in the [initval]{.title-ref}
and [endval]{.title-ref} blocks. We need terminal conditions for `c`
and `x`, since both appear with a lead, and an initial condition for
`k`, since it appears with a lag.

Setting `x=1.1` in the `endval` block without a `shocks` block implies
that technology is at $1.1$ in $t=1$ and stays there forever, because
`endval` is filling all entries of `oo_.endo_simul` and
`oo_.exo_simul` except for the very first one, which stores the
initial conditions and was set to $0$ by the `initval` block when not
explicitly specifying a value for it.

Because the law of motion for capital is backward-looking, we need an
initial condition for `k` at time $0$. Due to the presence of
`endval`, this cannot be done via a `histval` block, but rather must
be specified in the `initval` block. Similarly, because the Euler
equation is forward-looking, we need a terminal condition for `c` at
$t=201$, which is specified in the `endval` block.

As can be seen, it is not necessary to specify `c` and `x` in the
`initval` block and `k` in the `endval` block, because they have no
impact on the results. Due to the optimization problem in the first
period being to choose `c,k` at $t=1$ given the predetermined capital
stock `k` inherited from $t=0$ as well as the current and future
values for technology `x`, the values for `c` and `x` at time $t=0$
play no role. The same applies to the choice of `c,k` at time $t=200$,
which does not depend on `k` at $t=201$. As the Euler equation shows,
that choice only depends on current capital as well as future
consumption `c` and technology `x`, but not on future capital `k`. The
intuitive reason is that those variables are the consequence of
optimization problems taking place in at periods $t=0$ and $t=201$,
respectively, which are not modeled here.

*Example*

```
initval;
c = 1.2;
k = 12;
x = 1;
end;

endval;
c = 2;
k = 20;
x = 1.1;
end;
```

In this example, initial conditions for the forward-looking variables
`x` and `c` are provided, together with a terminal condition for the
backward-looking variable `k`. As shown in the previous example, these
values will not affect the simulation results. Dynare simply takes
them as given and basically assumes that there were realizations of
exogenous variables and states that make those choices equilibrium
values (basically initial/terminal conditions at the unspecified time
periods $t<0$ and $t>201$).

The above example suggests another way of looking at the use of
`steady` after `initval` and `endval`. Instead of saying that the
implicit unspecified conditions before and after the simulation range
have to fit the initial/terminal conditions of the endogenous
variables in those blocks, steady specifies that those conditions at
$t<0$ and $t>201$ are equal to being at the steady state given the
exogenous variables in the `initval` and `endval` blocks. The
endogenous variables at $t=0$ and $t=201$ are then set to the
corresponding steady state equilibrium values.

The fact that `c` at $t=0$ and `k` at $t=201$ specified in `initval`
and `endval` are taken as given has an important implication for
plotting the simulated vector for the endogenous variables, i.e. the
rows of `oo_.endo_simul`: this vector will also contain the initial
and terminal conditions and thus is 202 periods long in the example.
When you specify arbitrary values for the initial and terminal
conditions for forward- and backward-looking variables, respectively,
these values can be very far away from the endogenously determined
values at $t=1$ and $t=200$. While the values at $t=0$ and $t=201$ are
unrelated to the dynamics for $0<t<201$, they may result in
strange-looking large jumps. In the example above, consumption will
display a large jump from $t=0$ to $t=1$ and capital will jump from
$t=200$ to $t=201$ when using `rplot`{.interpreted-text role="comm"}
or manually plotting `oo_.endo_val`.

- *Block* : `histval ;`

- *Block*: `histval (OPTIONS...);`

*In a deterministic perfect foresight context*

In models with lags on more than one period, the `histval` block permits
to specify different historical initial values for different periods of
the state variables. In this case, the `initval` block takes over the
role of specifying terminal conditions and starting values for the
solver. Note that the `histval` block does not take non-state variables.

This block is terminated by `end;` and contains lines of the form:

```
VARIABLE\_NAME(INTEGER) = EXPRESSION;
```

EXPRESSION is any valid expression returning a numerical value and can
contain already initialized variable names.

By convention in Dynare, period 1 is the first period of the simulation.
Going backward in time, the first period before the start of the
simulation is period 0, then period -1, and so on.

State variables not initialized in the `histval` block are assumed to
have a value of zero at period 0 and before. Note that `histval` cannot
be followed by `steady`.

*Example*

```
model;
x=1.5*x(-1)-0.6*x(-2)+epsilon;
log(c)=0.5*x+0.5*log(c(+1));
end;

histval;
x(0)=-1;
x(-1)=0.2;
end;

initval;
c=1;
x=1;
end;
```

In this example, `histval` is used to set the historical conditions
for the two lags of the endogenous variable `x`, stored in the first
column of `oo_.endo_simul`. The `initval` block is used to set the
terminal condition for the forward looking variable `c`, stored in the
last column of `oo_.endo_simul`. Moreover, the `initval` block defines
the starting values for the perfect foresight solver for both
endogenous variables `c` and `x`.

*In a stochastic simulation context*

In the context of stochastic simulations, `histval` allows setting the
starting point of those simulations in the state space. As for the case
of perfect foresight simulations, all not explicitly specified variables
are set to 0. Moreover, as only states enter the recursive policy
functions, all values specified for control variables will be ignored.
This can be used

-   In `stoch_simul`{.interpreted-text role="comm"}, if the `periods`
    option is specified. Note that this only affects the starting
    point for the simulation, but not for the impulse response
    functions. When using the `loglinear <logl>`{.interpreted-text
    role="ref"} option, the `histval` block nevertheless takes the
    unlogged starting values.
-   In `forecast`{.interpreted-text role="comm"} as the initial point
    at which the forecasts are computed. When using the `loglinear
    <logl>`{.interpreted-text role="ref"} option, the `histval` block
    nevertheless takes the unlogged starting values.
-   In `conditional_forecast`{.interpreted-text role="comm"} for a
    calibrated model as the initial point at which the conditional
    forecasts are computed. When using the
    `loglinear <logl>`{.interpreted-text role="ref"} option, the
    histval-block nevertheless takes the unlogged starting values.
-   In `Ramsey policy <ramsey_model>`{.interpreted-text role="comm"},
    where it also specifies the values of the endogenous states
    (including lagged exogenous) at which the objective function of
    the planner is computed. Note that the initial values of the
    Lagrange multipliers associated with the planner's problem cannot
    be set (see `evaluate_planner_objective`{.interpreted-text
    role="comm"}).

*Options*

- `all_values_required`

See `all_values_required`{.interpreted-text role="opt"}.

*Example*

```
var x y;
varexo e;

model;
x = y(-1)^alpha*y(-2)^(1-alpha)+e;

end;

initval;
x = 1;
y = 1;
e = 0.5;
end;

steady;

histval;
y(0) = 1.1;
y(-1) = 0.9;
end;

stoch_simul(periods=100);
```

*Command*: `resid ;`

This command will display the residuals of the static equations of the
model, using the values given for the endogenous in the last `initval`
or `endval` block (or the steady state file if you provided one, see
`st-st`{.interpreted-text role="ref"}).

*Options*

- `non_zero`

Only display non-zero residuals.

*Command*: `initval_file (OPTIONS...);`

In a deterministic setup, this command is used to specify a path for all
endogenous and exogenous variables. The length of these paths must be
equal to the number of simulation periods, plus the number of leads and
the number of lags of the model (for example, with 50 simulation
periods, in a model with 2 lags and 1 lead, the paths must have a length
of 53). Note that these paths cover two different things:

-   The constraints of the problem, which are given by the path for
    exogenous and the initial and terminal values for endogenous
-   The initial guess for the non-linear solver, which is given by the
    path for endogenous variables for the simulation periods
    (excluding initial and terminal conditions)

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

- `datafile = FILENAME filename = FILENAME (deprecated)`

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
    first column may contain the date of each observation.
-   CSV files (extension `.csv`): for each endogenous and
    exogenous variable, the file must contain a column of the same
    name. The first column may contain the date of each
    observation.

- `first_obs = {INTEGER | DATE}`

The observation number or the date (see `dates-members`{.interpreted-text role="ref"}) of 
the first observation to be used in the file

- `first_simulation_period = {INTEGER | DATE}`

The observation number in the file or the date (see
`dates <dates-members>`{.interpreted-text role="ref"}) at which the
simulation (or the forecast) is starting. This option avoids to have to
compute the maximum number of lags in the model. The observation
corresponding to the first period of simulation doesn't need to exist in
the file as the only dates necessary for initialization are before that
date.

- `last_simulation_period = {INTEGER | DATE}`

The observation number in the file or the date (see
`dates <dates-members>`{.interpreted-text role="ref"}) at which the
simulation (or the forecast) is ending. This option avoids to have to
compute the maximum number of leads in the model.

- `last_obs = {INTEGER | DATE}`

The observaton number or the date (see `dates-members`{.interpreted-text role="ref"}) of the last observation
to be used in the file.

- `nobs = INTEGER`

The number of observations to be used in the file (starting with first
of `first_obs` observation).

- `series = DSERIES NAME`

The name of a DSERIES containing the data (see
`dseries-members`{.interpreted-text role="ref"})

*Example 1*

```
var c x;
varexo e;
parameters a b c d;

a = 1.5; b = -0,6; c = 0.5; d = 0.5;
model; x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv);

perfect_foresight_setup(periods=200); 
perfect_foresight_solver;
```

The initial and terminal values are taken from file `mydata.csv`
(nothing guarantees that these vales are the steady state of the model).
The guess value for the trajectories are also taken from the file. The
file must contain at least 203 observations of variables `c`, `x` and
`e`. If there are more than 203 observations available in the file, the
first 203 are used by `perfect_foresight_setup(periods=200)`. Note that
the values for the auxiliary variable corresponding to `x(-2)` are
automatically computed by `initval_file`.

*Example 2*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; b = -0,6; c = 0.5; d = 0.5;

model; x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv,   
             first_obs=10);

perfect_foresight_setup(periods=200); 
perfect_foresight_solver;
```

The initial and terminal values are taken from file `mydata.csv`
starting with the 10th observation in the file. There must be at least
212 observations in the file.

*Example 3*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5;
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

ds = dseries(mydata.csv);
lds = log(ds);

initval_file(series=lds,
            first\_obs=2010Q1);

perfect_foresight_setup(periods=200);
perfect_foresight_solver;
```

The initial and terminal values are taken from dseries `lds`. All
observations are loaded starting with the 1st quarter of 2010 until the
end of the file. There must be data available at least until 2050Q3.

*Example 4*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6;
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv,
            first_simulation_period=2010Q1);

perfect_foresight_setup(periods=200);
perfect_foresight_solver;
```


The initial and terminal values are taken from file `mydata.csv`. The
observations in the file must have dates. All observations are loaded
from the 3rd quarter of 2009 until the end of the file. There must be
data available in the file at least until 2050Q1.

*Example 5*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv,
            last_obs = 212);

perfect_foresight_setup(periods=200);
perfect_foresight_solver;
```

The initial and terminal values are taken from file `mydata.csv`. The
first 212 observations are loaded and the first 203 observations will be
used by `perfect_foresight_setup(periods=200)`.

*Example 6*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6;
c = 0.5; 
d = 0.5;

model; x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv,
            first_obs = 10,nobs = 203);

perfect_foresight_setup(periods=200); 
perfect_foresight_solver;
```

The initial and terminal values are taken from file `mydata.csv`.
Observations 10 to 212 are loaded.

*Example 7*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
  x = a\*x(-1) + b\*x(-2) + e; 
  log(c) = c\*x + d\*log(c(+1));
end;

initval_file(datafile=mydata.csv,
            first\_obs = 10);

steady;
```

The values of the 10th observation of `mydata.csv` are used

as guess value to compute the steady state. The exogenous variables are
set to values found in the file or zero if these variables aren't
present.

*Command*: `histval_file (OPTIONS...);`

This command is equivalent to `histval`, except that it reads its input
from a file, and is typically used in conjunction with
`smoother2histval`.

*Options*

- `datafile = FILENAME` 

   `filename = FILENAME (deprecated)`

The name of the file containing the data. The command accepts the following file formats:

-   M-file (extension `.m`): for each endogenous and exogenous
    variable, the file must contain a row or column vector of the same
    name.
-   MAT-file (extension `.mat`): same as for M-files.
-   Excel file (extension `.xls` or `.xlsx`): for each endogenous and
    exogenous variable, the file must contain a column of the same
    name. NB: Octave only supports the `.xlsx` file extension and must
    have the [io](https://octave.sourceforge.io/io/) package installed
    (easily done via octave by typing '`pkg install -forge io`'). The
    first column may contain the date of each observation.
-   CSV files (extension `.csv`): for each endogenous and exogenous
    variable, the file must contain a column of the same name. The
    first column may contain the date of each observation.

- `first_obs = {INTEGER | DATE}`

The observation number or the date (see
`dates-members`{.interpreted-text role="ref"}) of
the first observation to be used in the file

- `first_simulation_period = {INTEGER | DATE}`

The observation number in the file or the date (see
`dates-members`{.interpreted-text role="ref"}) at which the simulation
(or the forecast) is starting. This option avoids to have to compute the
maximum number of lags in the model. The observation corresponding to
the first period of simulation doesn't need to exist in the file as the
only dates necessary for initialization are before that date.

- `last_simulation_period = {INTEGER | DATE}`

The observation number in the file or the date (see
`dates <dates-members>`{.interpreted-text role="ref"}) at which the
simulation (or the forecast) is ending. This option avoids to have to
compute the maximum number of leads in the model.

- `last_obs = {INTEGER | DATE}`

The observation number or the date (see
`dates-members`{.interpreted-text role="ref"}) of the
last observation to be used in the file.

- `nobs = INTEGER`

The number of observations to be used in the file (starting with first
of `first_obs` observation).

- `series = DSERIES NAME`

The name of a DSERIES containing the data (see
`dseries-members`{.interpreted-text role="ref"})

*Example 1*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a\*x(-1) + b\*x(-2) + e; 
log(c) = c\*x + d\*log(c(+1));
end;

steady_state_model; 
x = 0; 
c = exp(c*x/(1 - d)); 
end;

histval_file(datafile=mydata.csv);

stoch_simul(order=1,periods=100);
```

The initial values for the stochastic simulation are taken from the two
first rows of file `mydata.csv`.

*Example 2*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

histval_file(datafile=mydata.csv,
            first_obs=10);

stoch_simul(order=1,periods=100);
```

The initial values for the stochastic simulation are taken from rows 10
and 11 of file `mydata.csv`.

*Example 3*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

histval_file(datafile=mydata.csv,
            first_obs=2010Q1);

stoch_simul(order=1,periods=100);
```

The initial values for the stochastic simulation are taken from
observations 2010Q1 and 2010Q2 of file `mydata.csv`.

*Example 4*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

histval_file(datafile=mydata.csv,
            first_simulation_period=2010Q1)

stoch_simul(order=1,periods=100);
```

The initial values for the stochastic simulation are taken from
observations 2009Q3 and 2009Q4 of file `mydata.csv`.

*Example 5*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

histval_file(datafile=mydata.csv,
            last_obs = 4);

stoch_simul(order=1,periods=100);
```

The initial values for the stochastic simulation are taken from the two
first rows of file `mydata.csv`.

*Example 6*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv,
            first_obs = 10,
            nobs = 4);

stoch_simul(order=1,periods=100);
```

The initial values for the stochastic simulation are taken from rows 10
and 11 of file `mydata.csv`.

*Example 7*

```
var c x;
varexo e;

parameters a b c d;

a = 1.5; 
b = -0,6; 
c = 0.5; 
d = 0.5;

model; 
x = a*x(-1) + b*x(-2) + e; 
log(c) = c*x + d*log(c(+1));
end;

initval_file(datafile=mydata.csv,
            first_obs=10);

histval_file(datafile=myotherdata.csv);

perfect_foresight_setup(periods=200);
perfect_foresight_solver;
```

Historical initial values for the simulation are taken from the two
first rows of file `myotherdata.csv`.

Terminal values and guess values for the simulation are taken
from file `mydata.csv` starting with the 12th observation in the file.
There must be at least 212 observations in the file.

## Shocks on exogenous variables

In a deterministic context, when one wants to study the transition of
one equilibrium position to another, it is equivalent to analyze the
consequences of a permanent shock and this in done in Dynare through the
proper use of `initval` and `endval`.

Another typical experiment is to study the effects of a temporary shock
after which the system goes back to the original equilibrium (if the
model is stable...). A temporary shock is a temporary change of value
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

- *block*: `shocks ;` 
          
- *block*: `shocks(overwrite);`

See above for the meaning of the `overwrite` option.

*In deterministic context*

For deterministic simulations, the `shocks` block specifies temporary
changes in the value of exogenous variables. For permanent shocks, use
an `endval` block.

The block should contain one or more occurrences of the following group
of three lines:

```
var VARIABLE_NAME;
periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
values DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;
```

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

!!! note

    Note that the first endogenous simulation period is period 1. Thus, a
    shock value specified for the initial period 0 may conflict with (i.e.
    may overwrite or be overwritten by) values for the initial period
    specified with `initval` or `endval` (depending on the exact context).
    Users should always verify the correct setting of `oo_.exo_simul` after
    `perfect_foresight_setup`.

*Example* (with scalar values)

```
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
```

*Example* (with vector values)

```
xx = [1.2; 1.3; 1];

shocks;
var e;
periods 1:3;
values (xx);
end;
```

*In stochastic context*

For stochastic simulations, the `shocks` block specifies the non zero
elements of the covariance matrix of the shocks of exogenous variables.

You can use the following types of entries in the block:

-   Specification of the standard error of an exogenous variable.

```
var VARIABLE_NAME; 
stderr EXPRESSION;
```

-   Specification of the variance of an exogenous variable.

```
var VARIABLE_NAME = EXPRESSION;
```

-   Specification the covariance of two exogenous variables.

```
var VARIABLE_NAME, VARIABLE_NAME = EXPRESSION;
```

-   Specification of the correlation of two exogenous variables.

```
corr VARIABLE_NAME, VARIABLE_NAME = EXPRESSION;
```

In an estimation context, it is also possible to specify variances and
covariances on endogenous variables: in that case, these values are
interpreted as the calibration of the measurement errors on these
variables. This requires the `varobs` command to be specified before the
`shocks` block.

*Example*

```
shocks;
var e = 0.000081;
var u; stderr 0.009;
corr e, u = 0.8;
var v, w = 2;
end;
```

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

```
shocks;
var u; stderr 0.008;
var u;
periods 1;
values 1;
end;
```

*Mixing deterministic and stochastic shocks*

It is possible to mix deterministic and stochastic shocks to build
models where agents know from the start of the simulation about future
exogenous changes. In that case `stoch_simul` will compute the rational
expectation solution adding future information to the state space
(nothing is shown in the output of `stoch_simul`) and `forecast` will
compute a simulation conditional on initial conditions and future
information.

*Example*

```
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
```



*MATLAB/Octave Command*: `get_shock_stderr_by_name ('EXOGENOUS_NAME\');`

Given the name of an exogenous variable, returns its standard deviation,
as set by a previous `shocks` block.

*MATLAB/Octave Command*: `set_shock_stderr_value ('EXOGENOUS_NAME', MATLAB_EXPRESSION);`

Sets the standard deviation of an exgonous variable. This does
essentially the same as setting the standard error via a `shocks` block,
except that it accepts arbitrary Julia expressions, and that it
works from Julia scripts.

## Steady state

There are two ways of computing the steady state (i.e. the static
equilibrium) of a model. The first way is to let Dynare compute the
steady state using a nonlinear Newton-type solver; this should work for
most models, and is relatively simple to use. The second way is to give
more guidance to Dynare, using your knowledge of the model, by providing
it with a method to compute the steady state, either using a
[steady_state_model]{.title-ref} block or writing matlab routine.

### Finding the steady state with Dynare nonlinear solver

- *Command*: `steady ;` 

- *Command*: `steady (OPTIONS...);`

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

- `maxit = INTEGER`

Determines the maximum number of iterations used in the non-linear
solver. The default value of `maxit` is 50.

- `tolf = DOUBLE`

Convergence criterion for termination based on the function value.
Iteration will cease when the residuals are smaller than `tolf`.
Default: `eps^(1/3)`

- `tolx = DOUBLE`

Convergence criterion for termination based on the step tolerance along.
Iteration will cease when the attempted step size is smaller than
`tolx`. Default: `eps^(2/3)`

- `homotopy_mode = INTEGER`

Use a homotopy (or divide-and-conquer) technique to solve for the steady
state. If you use this option, you must specify a `homotopy_setup`
block. This option can take three possible values:

- `1`: In this mode, all the parameters are changed simultaneously, and the
       distance between the boundaries for each parameter is divided in as
       many intervals as there are steps (as defined by the
       `homotopy_steps` option); the problem is solved as many times as
       there are steps.

- `2`: Same as mode `1`, except that only one parameter is changed at a
       time; the problem is solved as many times as steps times number of
       parameters.

- `3`: Dynare tries first the most extreme values. If it fails to compute
       the steady state, the interval between initial and desired values is
       divided by two for all parameters. Every time that it is impossible
       to find a steady state, the previous interval is divided by two.
       When it succeeds to find a steady state, the previous interval is
       multiplied by two. In that last case `homotopy_steps` contains the
       maximum number of computations attempted before giving up.

- `homotopy_steps = INTEGER`

Defines the number of steps when performing a homotopy. See
`homotopy_mode` option for more details.

- `homotopy_force_continue = INTEGER`

This option controls what happens when homotopy fails.

- `0`: `steady` fails with an error message

- `1`: `steady` keeps the values of the last homotopy step that was
        successful and continues. **BE CAREFUL**: parameters and/or
        exogenous variables are NOT at the value expected by the user
        Default is `0`.

- `nocheck`

Don't check the steady state values when they are provided explicitly
either by a steady state file or a `steady_state_model` block. This is
useful for models with unit roots as, in this case, the steady state is
not unique or doesn't exist.


*Example*

See `init-term-cond`{.interpreted-text role="ref"}.

After computation, the steady state is available in the following
variable:

*MATLAB/Octave Command*: `oo.steady_state`

Contains the computed steady state. Endogenous variables are ordered in
the order of declaration used in the `var` command (which is also the
order used in `M_.endo_names`).

- *MATLAB/Octave Command*: `get_mean ('ENDOGENOUS_NAME' [, 'ENDOGENOUS_NAME']... );`

Returns the steady of state of the given endogenous variable(s), as it
is stored in `oo_.steady_state`. Note that, if the steady state has not
yet been computed with `steady`, it will first try to compute it.

*Block*: `homotopy_setup ;`

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

```
VARIABLE_NAME, EXPRESSION, EXPRESSION;
```

This syntax specifies the initial and final values of a given
parameter/exogenous.

There is an alternative syntax:

```
VARIABLE_NAME, EXPRESSION;
```

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

```
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
```

### Providing the steady state to Dynare

If you know how to compute the steady state for your model, you can
provide a MATLAB/Octave function doing the computation instead of using
`steady`. Again, there are two options for doing that:

-   The easiest way is to write a `steady_state_model` block, which is
    described below in more details. See also `fs2000.mod` in the
    `examples` directory for an example. The steady state file
    generated by Dynare will be called `+FILENAME/steadystate.m.`
-   You can write the corresponding MATLAB function by hand. If your
    MOD-file is called `FILENAME.mod`, the steady state file must be
    called `FILENAME_steadystate.m`. See `NK_baseline_steadystate.m`
    in the examples directory for an example. This option gives a bit
    more flexibility (loops and conditional structures can be used),
    at the expense of a heavier programming burden and a lesser
    efficiency.

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

*Block*: `steady_state_model ;`

When the analytical solution of the model is known, this command can be
used to help Dynare find the steady state in a more efficient and
reliable way, especially during estimation where the steady state has to
be recomputed for every point in the parameter space.

Each line of this block consists of a variable (either an endogenous, a
temporary variable or a parameter) which is assigned an expression
(which can contain parameters, exogenous at the steady state, or any
endogenous or temporary variable already declared above). Each line
therefore looks like:

```
VARIABLE_NAME = EXPRESSION;
```

Note that it is also possible to assign several variables at the same
time, if the main function in the right hand side is a MATLAB/Octave
function returning several arguments:

```
[ VARIABLE_NAME, VARIABLE_NAME... ] = EXPRESSION;
```

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

```
var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

...
// parameter calibration, (dynamic) model declaration, shock calibration...
...

steady_state_model;
dA = exp(gam);
gst = 1/dA; // A temporary variable
m = mst;

// Three other temporary variables
khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );

n  = xist/(nust+xist);
P  = xist + nust;
k  = khst*n;

l  = psi*mst*n/( (1-psi)*(1-n) );
c  = mst/P;
d  = l - mst + 1;
y  = k^alp*n^(1-alp)*gst^alp;
R  = mst/bet;

// You can use MATLAB functions which return several arguments
[W, e] = my_function(l, n);

gp_obs = m/dA;
gy_obs = dA;
end;

steady;
```

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

```
var c k;
varexo x;
...
model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
[dynamic] c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
[static] k = ((delt+bet)/(x*aa*alph))^(1/(alph-1));
end;
```

## Getting information about the model

*Command*: `check ;` 

*Command*: `check (OPTIONS...);`

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

- `solve_algo = INTEGER`

See `solve_algo <solvalg>`{.interpreted-text role="ref"}, for the
possible values and their meaning.

- `qz_zero_threshold = DOUBL

Value used to test if a generalized eigenvalue is $0/0$ in the
generalized Schur decomposition (in which case the model does not admit
a unique solution). Default: `1e-6`.

*Output*

`check` returns the eigenvalues in the global variable `oo_.dr.eigval`.

- *MATLAB/Octave Command*: `oo.dr.eigval`

Contains the eigenvalues of the model, as computed by the `check`

*Command*: `model_diagnostics ;`

This command performs various sanity checks on the model, and prints a
message if a problem is detected (missing variables at current period,
invalid steady state, singular Jacobian of static model).

*Command*: `model_info ; 

*Command*: `model_info (OPTIONS...);`

This command provides information about the model.

When used outside the context of the `block` option of the `model`
block, it will provide a list of predetermined state variables,
forward-looking variables, and purely static variables.


## Deterministic simulation

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

!!! note
    Be careful when employing auxiliary variables in the context of perfect
    foresight computations. The same model may work for stochastic
    simulations, but fail for perfect foresight simulations. The issue
    arises when an equation suddenly only contains variables dated `t+1` (or
    `t-1` for that matter). In this case, the derivative in the last (first)
    period with respect to all variables will be 0, rendering the stacked
    Jacobian singular.

*Example*

Consider the following specification of an Euler equation with log
utility:

```
Lambda = beta*C(-1)/C;
Lambda(+1)*R(+1)= 1;
```

Clearly, the derivative of the second equation with respect to all
endogenous variables at time `t` is zero, causing
`perfect_foresight_solver` to generally fail. This is due to the use
of the Lagrange multiplier `Lambda` as an auxiliary variable. Instead,
employing the identical

```
beta*C/C(+1)*R(+1)= 1;
```

will work.

*Command*: `perfect_foresight_setup ;

*Command*: `perfect_foresight_setup (OPTIONS...);`

Prepares a perfect foresight simulation, by extracting the information
in the `initval`, `endval` and `shocks` blocks and converting them into
simulation paths for exogenous and endogenous variables.

This command must always be called before running the simulation with
`perfect_foresight_solver`.

*Options*

- `periods = INTEGER`

Number of periods of the simulation.

- `datafile = FILENAME`

Used to specify path for all endogenous and exogenous variables.
Strictly equivalent to `initval_file`{.interpreted-text role="comm"}.

*Output*

The paths for the exogenous variables are stored into `oo_.exo_simul`.

The initial and terminal conditions for the endogenous variables and the
initial guess for the path of endogenous variables are stored into
`oo_.endo_simul`.

*Command*: `perfect_foresight_solver ;`

*Command*: `perfect_foresight_solver (OPTIONS...);`

Computes the perfect foresight (or deterministic) simulation of the
model.

Note that `perfect_foresight_setup` must be called before this command,
in order to setup the environment for the simulation.

*Options*

- `maxit = INTEGER`

Determines the maximum number of iterations used in the non-linear
solver. The default value of `maxit` is `50`.

- `tolf = DOUBLE`

Convergence criterion for termination based on the function value.
Iteration will cease when it proves impossible to improve the function
value by more than `tolf`. Default: `1e-5`

- `tolx = DOUBLE`

Convergence criterion for termination based on the change in the
function argument. Iteration will cease when the solver attempts to take
a step that is smaller than `tolx`. Default: `1e-5`

- `noprint`

Don't print anything. Useful for loops.

- `print`

Print results (opposite of `noprint`).

- `lmmcp`

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

```
    model;
       ...
       [mcp = 'r > -1.94478']
       r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;
       ...
    end;
```

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

- `endogenous_terminal_period`

The number of periods is not constant across Newton iterations when
solving the perfect foresight model. The size of the nonlinear system of
equations is reduced by removing the portion of the paths (and
associated equations) for which the solution has already been identified
(up to the tolerance parameter). This strategy can be interpreted as a
mix of the shooting and relaxation approaches. Note that round off
errors are more important with this mixed strategy (user should check
the reported value of the maximum absolute error). Only available with
option `stack_solve_algo==0`.

- `linear_approximation`

Solves the linearized version of the perfect foresight model. The model
must be stationary. Only available with option `stack_solve_algo==0` or
`stack_solve_algo==7`.

*Output*

The simulated endogenous variables are available in global matrix
`oo_.endo_simul`.


This variable stores the path of exogenous variables during a simulation
(computed by `perfect_foresight_solver`, `simul`, `stoch_simul` or
`extended_path`). The variables are arranged in columns, in order of
declaration (as in `M_.exo_names`). Periods are in rows. Note that this
convention regarding columns and rows is the opposite of the convention
for `oo_.endo_simul`!


## Stochastic solution and simulation

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

*Command*: `stoch_simul [VARIABLE_NAME...];`

*Command*: `stoch_simul (OPTIONS...)[VARIABLE_NAME...];`

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
[https://archives.dynare.org/DynareWiki/IrFs](https://archives.dynare.org/DynareWiki/IrFs).

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

- `ar = INTEGER`

Order of autocorrelation coefficients to compute and to print. Default:
`5`.

- `drop = INTEGER`

Number of points (burnin) dropped at the beginning of simulation before
computing the summary statistics. Note that this option does not affect
the simulated series stored in `oo_.endo_simul` and the workspace. Here,
no periods are dropped. Default: `100`.

- `hp_filter = DOUBLE`

Uses HP filter with $\lambda =$ `DOUBLE` before computing moments. If
theoretical moments are requested, the spectrum of the model solution is
filtered following the approach outlined in Uhlig (2001). Default: no
filter.

- `one_sided_hp_filter = DOUBLE`

Uses the one-sided HP filter with $\lambda =$ `DOUBLE` described in
*Stock and Watson (1999)* before computing moments. This option is only
available with simulated moments. Default: no filter.

- `bandpass_filter`

Uses a bandpass filter with the default passband before computing
moments. If theoretical moments are requested, the spectrum of the model
solution is filtered using an ideal bandpass filter. If empirical
moments are requested, the *Baxter and King (1999)* filter is used.
Default: no filter.

- `bandpass_filter = [HIGHEST_PERIODICITY LOWEST_PERIODICITY]`

Uses a bandpass filter before computing moments. The passband is set to
a periodicity of to LOWEST_PERIODICITY, e.g. $6$ to $32$ quarters if
the model frequency is quarterly. Default: `[6,32]`.

- `filtered_theoretical_moments_grid = INTEGER`

When computing filtered theoretical moments (with either option
`hp_filter` or option `bandpass_filter`), this option governs the number
of points in the grid for the discrete Inverse Fast Fourier Transform.
It may be necessary to increase it for highly autocorrelated processes.
Default: `512`.

- `irf = INTEGER`

Number of periods on which to compute the IRFs. Setting `irf=0`
suppresses the plotting of IRFs. Default: `40`.

- `irf_shocks = ( VARIABLE_NAME [[,] VARIABLE_NAME ...] )`

The exogenous variables for which to compute IRFs. Default: all.

- `relative_irf`

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

- `irf_plot_threshold = DOUBLE`

Threshold size for plotting IRFs. All IRFs for a particular variable
with a maximum absolute deviation from the steady state smaller than
this value are not displayed. Default: `1e-10`.

- `nocorr`

Don't print the correlation matrix (printing them is the default).

- `nodecomposition`

Don't compute (and don't print) unconditional variance decomposition.

- `nofunctions`

Don't print the coefficients of the approximated solution (printing them
is the default).

- `nomoments`

Don't print moments of the endogenous variables (printing them is the
default).

- `nograph`

Do not create graphs (which implies that they are not saved to the disk
nor displayed). If this option is not used, graphs will be saved to disk
(to the format specified by `graph_format` option, except if
`graph_format=none`) and displayed to screen (unless `nodisplay` option
is used).

- `graph`

Re-enables the generation of graphs previously shut off with `nograph`.

- `nodisplay`

Do not display the graphs, but still save them to disk (unless `nograph`
is used).

- `graph_format = FORMAT graph_format = ( FORMAT, FORMAT... )`

Specify the file format(s) for graphs saved to disk. Possible values are
`eps` (the default), `pdf`, `fig` and `none` (under Octave, `fig` is
unavailable). If the file format is set equal to `none`, the graphs are
displayed but not saved to the disk.

- `noprint`

See `noprint`{.interpreted-text role="opt"}.

- `print`

See `print`{.interpreted-text role="opt"}.

- `order = INTEGER`

Order of Taylor approximation. Note that for third order and above, the
`k_order_solver` option is implied and only empirical moments are
available (you must provide a value for `periods` option). Default: `2`
(except after an `estimation` command, in which case the default is the
value used for the estimation).

- `k_order_solver`

Use a k-order solver (implemented in C++) instead of the default Dynare
solver. This option is not yet compatible with the `bytecode` option
(see `model-decl`{.interpreted-text role="ref"}). Default: disabled for
order 1 and 2, enabled for order 3 and above.

- `periods = INTEGER`

If different from zero, empirical moments will be computed instead of
theoretical moments. The value of the option specifies the number of
periods to use in the simulations. Values of the initval block, possibly
recomputed by `steady`, will be used as starting point for the
simulation. The simulated endogenous variables are made available to the
user in a vector for each variable and in the global matrix
`oo_.endo_simul` (see `oo_.endo_simul`{.interpreted-text role="mvar"}).
The simulated exogenous variables are made available in `oo_.exo_simul`
(see `oo_.exo_simul`{.interpreted-text role="mvar"}). Default: `0`.

- `qz_criterium = DOUBLE`

Value used to split stable from unstable eigenvalues in reordering the
Generalized Schur decomposition used for solving first order problems.
Default: `1.000001` (except when estimating with `lik_init` option equal
to `1`: the default is `0.999999` in that case; see
`estim`{.interpreted-text role="ref"}).

- `qz_zero_threshold = DOUBLE`

See `qz_zero_threshold <qz_zero_threshold = DOUBLE>`{.interpreted-text

- `replic = INTEGER`

Number of simulated series used to compute the IRFs. Default: `1` if
`order=1`, and `50` otherwise.

- `simul_replic = INTEGER`

Number of series to simulate when empirical moments are requested (i.e.
`periods` $>$ 0). Note that if this option is greater than 1, the
additional series will not be used for computing the empirical moments
but will simply be saved in binary form to the file `FILENAME_simul` in
the `FILENAME/Output`-folder. Default: `1`.

- `solve_algo = INTEGER`

See `solve_algo <solvalg>`{.interpreted-text role="ref"}, for the
possible values and their meaning.

- `conditional_variance_decomposition = INTEGER`

- `conditional_variance_decomposition = [INTEGER1:INTEGER2]`

- `conditional_variance_decomposition = [INTEGER1 INTEGER2 ...]`

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
Only available at `order<3` and without `pruning''. 
In case of `order=2`, Dynare provides a second-order accurate approximation to 
the true second moments based on the linear terms of the second-order solution 
(see *Kim, Kim, Schaumburg and Sims (2008)*). Note that the unconditional variance 
decomposition *i.e.* at horizon infinity) is automatically conducted if theoretical moments are requested and if`nodecomposition`[is not set (see :mvar:`oo_.variance_decomposition]{.title-ref}).

- `pruning`

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

- `sylvester = OPTION`

Determines the algorithm used to solve the Sylvester equation for block
decomposed model. Possible values for OPTION are:

   `default`
    
    Uses the default solver for Sylvester equations (`gensylv`) based on
    Ondra Kamenik's algorithm (see
    [here](https://www.dynare.org/assets/dynare++/sylvester.pdf) for
    more information).

  `fixed_point`

    Uses a fixed point algorithm to solve the Sylvester equation
    `gensylv_fp`). This method is faster than the default one for large
    scale models.

Default value is `default`.

- `sylvester_fixed_point_tol = DOUBL

The convergence criterion used in the fixed point Sylvester solver. Its
default value is `1e-12`.

- `dr = OPTION`

Determines the method used to compute the decision rule. Possible values
for OPTION are:

  `default`

    Uses the default method to compute the decision rule based on the
    generalized Schur decomposition (see *Villemot (2011)* for more
    information).

  `cycle_reduction`

    Uses the cycle reduction algorithm to solve the polynomial equation
    for retrieving the coefficients associated to the endogenous
    variables in the decision rule. This method is faster than the
    default one for large scale models.

Default value is `default`.

- `dr_cycle_reduction_tol = DOUBLE`

The convergence criterion used in the cycle reduction algorithm. Its
default value is `1e-7`.

- `tex`

Requests the printing of results and graphs in TeX tables and graphics
that can be later directly included in LaTeX files.

- `dr_display_tol = DOUBLE`

Tolerance for the suppression of small terms in the display of decision
rules. Rows where all terms are smaller than `dr_display_tol` are not
displayed. Default value: `1e-6`.

- `contemporaneous_correlation`

Saves the contemporaneous correlation between the endogenous variables
in `oo_.contemporaneous_correlation`. Requires the `nocorr` option not
to be set.

- `spectral_density`

Triggers the computation and display of the theoretical spectral density
of the (filtered) model variables. Results are stored in
`oo_.SpectralDensity`, defined below. Default: do not request spectral
density estimates.

- `hp_ngrid = INTEGER`

Deprecated option. It has the same effect as
`filtered_theoretical_moments_grid <filtered_theoretical_moments_grid = INTEGER>`{.interpreted-text
role="opt"}.

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

```
shocks;
var e;
stderr 0.0348;
end;

stoch_simul;
```

Performs the simulation of the 2nd-order approximation of a model with
a single stochastic shock `e`, with a standard error of `0.0348`.

*Example*

```
    stoch_simul(irf=60) y k;
```

Performs the simulation of a model and displays impulse response
functions on 60 periods for variables `y` and `k`.

*MATLAB/Octave*: `oo.mean`

After a run of `stoch_simul`, contains the mean of the endogenous
variables. Contains theoretical mean if the `periods` option is not
present, and simulated mean otherwise. The variables are arranged in
declaration order.

*MATLAB/Octave*: `oo.var`

After a run of `stoch_simul`, contains the variance-covariance of the
endogenous variables. Contains theoretical variance if the `periods`
option is not present and simulated variance otherwise. Only available
for `order<4`. At `order=2` it will be be a second-order accurate
approximation (i.e. ignoring terms of order 3 and 4 that would arise
when using the full second-order policy function). At `order=3`,
theoretical moments are only available with `pruning`. The variables are
arranged in declaration order.

*MATLAB/Octave*: `oo.var_list`

The list of variables for which results are displayed.

*MATLAB/Octave*: `oo.skewness`

After a run of `stoch_simul` contains the skewness (standardized third
moment) of the simulated variables if the `periods` option is present.
The variables are arranged in declaration order.

*MATLAB/Octave*: `oo.kurtosis`

After a run of `stoch_simul` contains the excess kurtosis (standardized
fourth moment) of the simulated variables if the `periods` option is
present. The variables are arranged in declaration order.

*MATLAB/Octave*: `oo.autocorr`

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

*MATLAB/Octave*: `oo.gamma_y`

After a run of `stoch_simul`, if theoretical moments have been requested
(i.e. if the `periods` option is not present), this variable contains a
cell array with the following values (where `ar` is the value of the
option of the same name):

  `oo_.gamma{1}`

    Variance/covariance matrix.

   `oo_.gamma{i+1}` (for i=1:ar)

    Autocorrelation function. See `oo_.autocorr`{.interpreted-text
    role="mvar"} for more details. **Beware**, this is the
    autocorrelation function, not the autocovariance function.

  `oo_.gamma{ar+2}`

    Unconditional variance decomposition, see
    `oo_.variance_decomposition`{.interpreted-text role="mvar"}.

   `oo_.gamma{ar+3}`

    If a second order approximation has been requested, contains the
    vector of the mean correction terms.
 
    Only available at `order<4`. In case `order=2`, the theoretical
    second moments are a second order accurate approximation of the true
    second moments. See conditional\_variance\_decomposition. At
    `order=3`, theoretical moments are only available with `pruning`.

*MATLAB/Octave*: `oo.variance_decomposition`

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

*MATLAB/Octave*: `oo.variance_decomposition_ME`

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

*MATLAB/Octave*: `oo.conditional_variance_decomposition`

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

*MATLAB/Octave*: `oo.conditional_variance_decomposition_ME`

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

*MATLAB/Octave*: `oo.contemporaneous_correlation`

After a run of `stoch_simul` with the
`contemporaneous_correlation option`, contains theoretical
contemporaneous correlations if the `periods` option is not present, and
simulated contemporaneous correlations otherwise. Only available for
`order<4`. At `order=2` it will be be a second-order accurate
approximation. At `order=3`, theoretical moments are only available with
`pruning`. The variables are arranged in declaration order.

*MATLAB/Octave*: `oo.SpectralDensity`

After a run of `stoch_simul` with option `spectral_density`, contains
the spectral density of the model variables. There will be a `nvars` by
`nfrequencies` subfield `freqs` storing the respective frequency grid
points ranging from $0$ to $2\pi$ and a same sized subfield `density`
storing the corresponding density.

*MATLAB/Octave*: `oo.irfs`

After a run of `stoch_simul` with option `irf` different from zero,
contains the impulse responses, with the following naming convention:
[VARIABLE_NAME_SHOCK_NAME]{.title-ref}.

For example, `oo_.irfs.gnp_ea` contains the effect on `gnp` of a
one-standard deviation shock on `ea`.

*MATLAB/Octave*: `get_irf ('EXOGENOUS_NAME' [, 'ENDOGENOUS_NAME']... );`

Given the name of an exogenous variables, returns the IRFs for the
requested endogenous variable(s), as they are stored in `oo_.irfs`.

The approximated solution of a model takes the form of a set of decision
rules or transition equations expressing the current value of the
endogenous variables of the model as function of the previous state of
the model and shocks observed at the beginning of the period. The
decision rules are stored in the structure `oo_.dr` which is described
below.

*MATLAB/Octave*: `oo.dr`

Structure storing the decision rules. The subfields for different orders
of approximation are explained below.

*Command*: `extended_path ;`

*Command*: `extended_path (OPTIONS...);`

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

- `periods = INTEGER`

The number of periods for which the simulation is to be computed. No
default value, mandatory option.

- `solver_periods = INTEGER`

The number of periods used to compute the solution of the perfect
foresight at every iteration of the algorithm. Default: `200`.

- `order = INTEGER`

If order is greater than `0` Dynare uses a gaussian quadrature to take
into account the effects of future uncertainty. If `order` $=S$ then the
time series for the endogenous variables are generated by assuming that
the agents believe that there will no more shocks after period $t+S$.
This is an experimental feature and can be quite slow. A non-zero value
is not compatible with either the `bytecode` or the `block` option of
the `model` block. Default: `0`.

- `hybrid`

Use the constant of the second order perturbation reduced form to
correct the paths generated by the (stochastic) extended path algorithm.

- `lmmcp`

Solves the perfect foresight model with a Levenberg-Marquardt mixed
complementarity problem (LMMCP) solver (*Kanzow and Petra (2004)*),
which allows to consider inequality constraints on the endogenous
variables (such as a ZLB on the nominal interest rate or a model with
irreversible investment). For specifying the necessary `mcp`-tag, see
`lmmcp`{.interpreted-text role="opt"}.


### Typology and ordering of variables

Dynare distinguishes four types of endogenous variables:

*Purely backward (or purely predetermined) variables*

Those that appear only at current and past period in the model, but
not at future period (i.e. at $t$ and $t-1$ but not $t+1$). The number
of such variables is equal to `M_.npred`.

*Purely forward variables*

Those that appear only at current and future period in the model, but
not at past period (i.e. at $t$ and $t+1$ but not $t-1$). The number
of such variables is stored in `M_.nfwrd`.

*Mixed variables*

Those that appear at current, past and future period in the model
(i.e. at $t$, $t+1$ and $t-1$). The number of such variables is stored
in `M_.nboth`.

*Static variables*

Those that appear only at current, not past and future period in the
model (i.e. only at $t$, not at $t+1$ or $t-1$). The number of such
variables is stored in `M_.nstatic`.

Note that all endogenous variables fall into one of these four
categories, since after the creation of auxiliary variables (see
`aux-variables`{.interpreted-text role="ref"}), all endogenous have at
most one lead and one lag. We therefore have the following identity:

``` {.sourceCode .matlab}
M_.npred + M_.both + M_.nfwrd + M_.nstatic = M_.endo_nbr
```

*MATLAB/Octave variable*: `M.state_var`

Vector of numerical indices identifying the state variables in the
vector of declared variables. `M_.endo_names(M_.state_var)` therefore
yields the name of all variables that are states in the model
declaration, i.e. that show up with a lag.

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

*MATLAB/Octave variable*: `oo.dr.order_var`

This variables maps DR-order to declaration order.

*MATLAB/Octave variable*: `oo.dr.inv_order_var`

This variable contains the inverse map.

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

*MATLAB/Octave variable*: `oo.dr.state_var`

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


## Estimation based on likelihood

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

*Command*: `varobs VARIABLE_NAME...;`

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

```
     varobs C y rr;
```

Declares endogenous variables `C`, `y` and `rr` as observed variables.

*Example* (with a macro processor loop)

```
   varobs
   @#for co in countries
   GDP_@{co}
   @#endfor
;
```

*Block*: `observation_trends ;`

This block specifies linear trends for observed variables as functions
of model parameters. In case the `loglinear` option is used, this
corresponds to a linear trend in the logged observables, i.e. an
exponential trend in the level of the observables.

Each line inside of the block should be of the form:

```
    VARIABLE_NAME(EXPRESSION);
```

In most cases, variables shouldn't be centered when `observation_trends`
is used.

*Example*

```
     observation_trends;
     Y (eta);
     P (mu/eta);
     end;
```

*Block*: `estimated_params ;`

*Block*: `estimated_params (overwrite) ;`

This block lists all parameters to be estimated and specifies bounds and
priors as necessary.

Each line corresponds to an estimated parameter.

In a maximum likelihood or a method of moments estimation, each line
follows this syntax:

```
    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME
    , INITIAL_VALUE [, LOWER_BOUND, UPPER_BOUND ];
```

In a Bayesian MCMC or a penalized method of moments estimation, each
line follows this syntax:

```
    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME | DSGE_PRIOR_WEIGHT
    [, INITIAL_VALUE [, LOWER_BOUND, UPPER_BOUND]], PRIOR_SHAPE,
    PRIOR_MEAN, PRIOR_STANDARD_ERROR [, PRIOR_3RD_PARAMETER [,
    PRIOR_4TH_PARAMETER [, SCALE_PARAMETER ] ] ];
```

The first part of the line consists of one of the four following
alternatives:

-   `stderr VARIABLE_NAME`

    Indicates that the standard error of the exogenous variable
    VARIABLE_NAME, or of the observation error/measurement errors
    associated with endogenous observed variable VARIABLE_NAME, is to
    be estimated.

-   `corr VARIABLE_NAME1, VARIABLE_NAME2`

    Indicates that the correlation between the exogenous variables
    VARIABLE_NAME1 and VARIABLE_NAME2, or the correlation of the
    observation errors/measurement errors associated with endogenous
    observed variables VARIABLE_NAME1 and VARIABLE_NAME2, is to be
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

- `INITIAL_VALUE`

Specifies a starting value for the posterior mode optimizer or the
maximum likelihood estimation. If unset, defaults to the prior mean.

- `LOWER_BOUND`

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

- `UPPER_BOUND`

Same as `lower_bound`, but specifying an upper bound instead.

- `PRIOR_SHAPE`

A keyword specifying the shape of the prior density. The possible values
are: `beta_pdf`, `gamma_pdf`, `normal_pdf`, `uniform_pdf`,
`inv_gamma_pdf`, `inv_gamma1_pdf`, `inv_gamma2_pdf` and `weibull_pdf`.
Note that `inv_gamma_pdf` is equivalent to `inv_gamma1_pdf`.

- `PRIOR_MEAN`

The mean of the prior distribution.

- `PRIOR_STANDARD_ERROR`

The standard error of the prior distribution.

- `PRIOR_3RD_PARAMETER`

A third parameter of the prior used for generalized beta distribution,
generalized gamma, generalized Weibull and for the uniform distribution.
Default: `0`.

- `PRIOR_4TH_PARAMETER`

A fourth parameter of the prior used for generalized beta distribution
and for the uniform distribution. Default: `1`.

- `SCALE_PARAMETER`

A parameter specific scale parameter for the jumping distribution's
covariance matrix of the Metropolis-Hasting algorithm.

Note that INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND, PRIOR_MEAN,
PRIOR_STANDARD_ERROR, PRIOR_3RD_PARAMETER, PRIOR_4TH_PARAMETER and
SCALE_PARAMETER can be any valid EXPRESSION. Some of them can be empty,
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
SCALE_PARAMETER, you must specify `PRIOR_3RD_PARAMETER` and
`PRIOR_4TH_PARAMETER`. Use empty values, if these parameters don't
apply.

*Example*

```
     corr eps_1, eps_2, 0.5,  ,  , beta_pdf, 0, 0.3, -1, 1;
```

Sets a generalized beta prior for the correlation between `eps_1` and
`eps_2` with mean `0` and variance `0.3`. By setting
`PRIOR_3RD_PARAMETER` to `-1` and `PRIOR_4TH_PARAMETER` to `1` the
standard beta distribution with support `[0,1]` is changed to a
generalized beta with support `[-1,1]`. Note that LOWER_BOUND and
UPPER_BOUND are left empty and thus default to `-1` and `1`,
respectively. The initial value is set to `0.5`.

*Example*

```
     corr eps_1, eps_2, 0.5,  -0.5,  1, beta_pdf, 0, 0.3, -1, 1;
```

Sets the same generalized beta distribution as before, but now
truncates this distribution to `[-0.5,1]` through the use of
LOWER_BOUND and UPPER_BOUND.

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

```
     parameters bet;

     model;
     # sig = 1/bet;
     c = sig*c(+1)*mpk;
     end;

     estimated_params;
     bet, normal_pdf, 1, 0.05;
     end;
```

It is possible to have several `estimated_params` blocks. By default,
subsequent blocks are concatenated with the previous ones; this can be
useful when building models in a modular fashion (see also
`estimated_params_remove`{.interpreted-text role="bck"} for that use
case). However, if an `estimated_params` block has the `overwrite`
option, its contents becomes the new list of estimated parameters,
cancelling previous blocks; this can be useful when doing several
estimations in a single `.mod` file.

*Block*: `estimated_params_init ;`

*Block*: `estimated_params_init (OPTIONS...);`

This block declares numerical initial values for the optimizer when
these ones are different from the prior mean. It should be specified
after the `estimated_params` block as otherwise the specified starting
values are overwritten by the latter.

Each line has the following syntax:

```
    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME, INITIAL_VALUE;
```

*Options*

- use_calibration

For not specifically initialized parameters, use the deep parameters and
the elements of the covariance matrix specified in the `shocks` block
from calibration as starting values for estimation. For components of
the `shocks` block that were not explicitly specified during calibration
or which violate the prior, the prior mean is used.

See `estimated_params`{.interpreted-text role="bck"}, for the meaning
and syntax of the various components.

*Block*: `estimated_params_bounds ;`

This block declares lower and upper bounds for parameters in maximum
likelihood estimation.

Each line has the following syntax:

```
    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME, LOWER_BOUND, UPPER_BOUND;
```

See `estimated_params`{.interpreted-text role="bck"}, for the meaning
and syntax of the various components.

*Block*: `estimated_params_remove ;`

This block partially undoes the effect of a previous
`estimated_params`{.interpreted-text role="bck"} block, by removing some
parameters from the estimation.

Each line has the following syntax:

```
    stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME;
```

*Command*: `estimation [VARIABLE_NAME...];`

*Command*: `estimation (OPTIONS...)[VARIABLE_NAME...];`

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

```
    set_dynare_seed('clock');
```

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

- `datafile = FILENAME`

The datafile: a `.m` file, a `.mat` file, a `.csv` file, or a
`.xls/.xlsx` file (under Octave, the
[io](https://octave.sourceforge.io/io/) package from Octave-Forge is
required for the `.csv` and `.xlsx` formats and the `.xls` file
extension is not supported). Note that the base name (i.e. without
extension) of the datafile has to be different from the base name of the
model file. If there are several files named FILENAME, but with
different file endings, the file name must be included in quoted strings
and provide the file ending like:

```
    estimation(datafile='../fsdat_simul.mat',...);
```

- `dirname = FILENAME`

Directory in which to store `estimation` output. To pass a subdirectory
of a directory, you must quote the argument. Default: `<mod_file>`.

- `xls_sheet = QUOTED_STRING`

The name of the sheet with the data in an Excel file.

- `xls_range = RANGE`

The range with the data in an Excel file. For example,
`xls_range=B2:D200`.

- `nobs = INTEGER`

The number of observations following `first_obs <first_obs
= [INTEGER1:INTEGER2]>`{.interpreted-text role="opt"} to be used.
Default: all observations in the file after `first_obs`.

- `nobs = [INTEGER1:INTEGER2]`

Runs a recursive estimation and forecast for samples of size ranging of
`INTEGER1` to `INTEGER2`. Option `forecast` must also be specified. The
forecasts are stored in the `RecursiveForecast` field of the results
structure (see
`RecursiveForecast <oo_.RecursiveForecast>`{.interpreted-text
role="mvar"}). The respective results structures `oo_` are saved in
`oo_recursive_` (see `oo_recursive_`{.interpreted-text role="mvar"}) and
are indexed with the respective sample length.

- `first_obs = INTEGER`

The number of the first observation to be used. In case of estimating a
DSGE-VAR, `first_obs` needs to be larger than the number of lags.
Default: `1`.

- `first_obs = [INTEGER1:INTEGER2]`

Runs a rolling window estimation and forecast for samples of fixed size
`nobs` starting with the first observation ranging from `INTEGER1` to
`INTEGER2`. Option `forecast` must also be specified. This option is
incompatible with requesting recursive forecasts using an expanding
window (see `nobs
<nobs = [INTEGER1:INTEGER2]>`{.interpreted-text role="opt"}). The
respective results structures `oo_` are saved in `oo_recursive_` (see
`oo_recursive_`{.interpreted-text role="mvar"}) and are indexed with the
respective first observation of the rolling window.

- `prefilter = INTEGER`

A value of 1 means that the estimation procedure will demean each data
series by its empirical mean. If the `loglinear
<logl>`{.interpreted-text role="ref"} option without the
`logdata`{.interpreted-text role="opt"} option is requested, the data
will first be logged and then demeaned. Default: `0`, i.e. no
prefiltering.

- `presample = INTEGER`

The number of observations after `first_obs <first_obs =
[INTEGER1:INTEGER2]>`{.interpreted-text role="opt"} to be skipped before
evaluating the likelihood. These presample observations do not enter the
likelihood, but are used as a training sample for starting the Kalman
filter iterations. This option is incompatible with estimating a
DSGE-VAR. Default: `0`.

- `loglinear`

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

- `logdata`

Dynare applies the $log$ transformation to the provided data if a
log-linearization of the model is requested
(`loglinear`{.interpreted-text role="opt"}) unless `logdata` option is
used. This option is necessary if the user provides data already in
logs, otherwise the $log$ transformation will be applied twice (this may
result in complex data).

- `plot_priors = INTEGER`

Control the plotting of priors.

 `0`

    No prior plot.

 `1`

    Prior density for each estimated parameter is plotted. It is
    important to check that the actual shape of prior densities matches
    what you have in mind. Ill-chosen values for the prior standard
    density can result in absurd prior densities.

Default value is `1`.

- `nograph`

See `nograph`{.interpreted-text role="opt"}.

- `posterior_nograph`

Suppresses the generation of graphs associated with Bayesian IRFs
(`bayesian_irf`{.interpreted-text role="opt"}), posterior smoothed
objects (`smoother`{.interpreted-text role="opt"}), and posterior
forecasts (`forecast`{.interpreted-text role="opt"}).

- `posterior_graph`

Re-enables the generation of graphs previously shut off with
`posterior_nograph`{.interpreted-text role="opt"}.

- `nodisplay`

See `nodisplay`{.interpreted-text role="opt"}.

- `graph_format = FORMAT` 

- `graph_format = ( FORMAT, FORMAT... )`

See
`graph_format <graph_format = ( FORMAT, FORMAT... )>`{.interpreted-text
role="opt"}.

- `no_init_estimation_check_first_obs`

Do not check for stochastic singularity in first period. If used,
[ESTIMATION CHECKS]{.title-ref} does not return an error if the check
fails only in first observation. This should only be used when observing
stock variables (e.g. capital) in first period, on top of their
associated flow (e.g. investment). Using this option may lead to a crash
or provide undesired/wrong results for badly specified problems (e.g.
the additional variable observed in first period is not predetermined).

For advanced use only.

- `lik_init = INTEGER`

Type of initialization of Kalman filter:

 `1`

    For stationary models, the initial matrix of variance of the error
    of forecast is set equal to the unconditional variance of the state
    variables.

 `2`

    For nonstationary models: a wide prior is used with an initial
    matrix of variance of the error of forecast diagonal with 10 on the
    diagonal (follows the suggestion of *Harvey and Phillips(1979)*).

 `3`

    For nonstationary models: use a diffuse filter (use rather the
    `diffuse_filter` option).

 `4`

    The filter is initialized with the fixed point of the Riccati
     equation.

 `5`

    Use i) option 2 for the non-stationary elements by setting their
    initial variance in the forecast error matrix to 10 on the diagonal
    and all covariances to 0 and ii) option 1 for the stationary
    elements.

Default value is 1. For advanced use only.

- `lik_algo = INTEGER`

For internal use and testing only.

- `conf_sig = DOUBLE`

Level of significance of the confidence interval used for classical
forecasting after estimation. Default: 0.9.

- `mh_conf_sig = DOUBLE`

Confidence/HPD interval used for the computation of prior and posterior
statistics like: parameter distributions, prior/posterior moments,
conditional variance decomposition, impulse response functions, Bayesian
forecasting. Default: `0.9`.

- `mh_replic = INTEGER`

Number of replications for each chain of the Metropolis-Hastings
algorithm. The number of draws should be sufficient to achieve
convergence of the MCMC and to meaningfully compute posterior objects.
Default: `20000`.

- `sub_draws = INTEGER`

Number of draws from the MCMC that are used to compute posterior
distribution of various objects (smoothed variable, smoothed shocks,
forecast, moments, IRF). The draws used to compute these posterior
moments are sampled uniformly in the estimated empirical posterior
distribution (i.e. draws of the MCMC). `sub_draws` should be smaller
than the total number of MCMC draws available. Default:
`min(posterior_max_subsample_draws, (Total number of draws)*(number of chains) )`.

- `posterior_max_subsample_draws = INTEGER`

Maximum number of draws from the MCMC used to compute posterior
distribution of various objects (smoothed variable, smoothed shocks,
forecast, moments, IRF), if not overriden by option `sub_draws`.
Default: `1200`.

- `mh_nblocks = INTEGER`

Number of parallel chains for Metropolis-Hastings algorithm. Default:
`2`.

- `mh_drop = DOUBLE`

The fraction of initially generated parameter vectors to be dropped as a
burn-in before using posterior simulations. Default: `0.5`.

- `mh_jscale = DOUBLE`

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

- `mh_init_scale = DOUBLE`

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

- `mh_tune_jscale [= DOUBLE]`

Automatically tunes the scale parameter of the jumping distribution\'s
covariance matrix (Metropolis-Hastings), so that the overall acceptance
ratio is close to the desired level. Default value is `0.33`. It is not
possible to match exactly the desired acceptance ratio because of the
stochastic nature of the algorithm (the proposals and the initial
conditions of the markov chains if `mh_nblocks>1`). This option is only
available for the Random Walk Metropolis Hastings algorithm. Must not be
used in conjunction with `mh_jscale = DOUBLE`{.interpreted-text
role="opt"}.

- `mh_tune_guess = DOUBLE`

Specifies the initial value for the `mh_tune_jscale
<mh_tune_jscale [= DOUBLE]>`{.interpreted-text role="opt"} option.
Default: `0.2`. Must not be set if
`mh_tune_jscale <mh_tune_jscale [= DOUBLE]>`{.interpreted-text
role="opt"} is not used.

- `mh_recover`

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

- `mh_posterior_mode_estimation`

Skip optimizer-based mode-finding and instead compute the mode based on
a run of a MCMC. The MCMC will start at the prior mode and use the prior
variances to compute the inverse Hessian.

- `mode_file = FILENAME`

Name of the file containing previous value for the mode. When computing
the mode, Dynare stores the mode (`xparam1`) and the hessian (`hh`, only
if `cova_compute=1`) in a file called `MODEL_FILENAME_mode.mat` in the
`FILENAME/Output`-folder. After a successful run of the estimation
command, the `mode_file` will be disabled to prevent other function
calls from implicitly using an updated mode-file. Thus, if the `.mod`
file contains subsequent `estimation` commands, the `mode_file` option,
if desired, needs to be specified again.

- `silent_optimizer`

Instructs Dynare to run mode computing/optimization silently without
displaying results or saving files in between. Useful when running
loops.

- `mcmc_jumping_covariance = OPTION`

Tells Dynare which covariance to use for the proposal density of the
MCMC sampler. OPTION can be one of the following:

 `hessian`

    Uses the Hessian matrix computed at the mode.

 `prior_variance`

    Uses the prior variances. No infinite prior variances are allowed in
    this case.

 `identity_matrix`

    Uses an identity matrix.
 
 `FILENAME`

    Loads an arbitrary user-specified covariance matrix from
    `FILENAME.mat`. The covariance matrix must be saved in a variable
    named `jumping_covariance`, must be square, positive definite, and
    have the same dimension as the number of estimated parameters.

Note that the covariance matrices are still scaled with
`mh_jscale <mh_jscale = DOUBLE>`{.interpreted-text role="opt"}. Default
value is `hessian`.

- `mode_check`

Tells Dynare to plot the posterior density for values around the
computed mode for each estimated parameter in turn. This is helpful to
diagnose problems with the optimizer. Note that for `order>1` the
likelihood function resulting from the particle filter is not
differentiable anymore due to the resampling step. For this reason, the
`mode_check` plot may look wiggly.

- `mode_check_neighbourhood_size = DOUBLE`

Used in conjunction with option `mode_check`, gives the width of the
window around the posterior mode to be displayed on the diagnostic
plots. This width is expressed in percentage deviation. The `Inf` value
is allowed, and will trigger a plot over the entire domain (see also
`mode_check_symmetric_plots`). Default:`0.5`.

- `mode_check_symmetric_plots = INTEGER`

Used in conjunction with option `mode_check`, if set to `1`, tells
Dynare to ensure that the check plots are symmetric around the posterior
mode. A value of `0` allows to have asymmetric plots, which can be
useful if the posterior mode is close to a domain boundary, or in
conjunction with `mode_check_neighbourhood_size = Inf` when the domain
in not the entire real line. Default: `1`.

- `mode_check_number_of_points = INTEGER`

Number of points around the posterior mode where the posterior kernel is
evaluated (for each parameter). Default is `20`.

- `prior_trunc = DOUBLE`

Probability of extreme values of the prior density that is ignored when
computing bounds for the parameters. Default: `1e-32`.

- `huge_number = DOUBLE`

Value for replacing infinite values in the definition of (prior) bounds
when finite values are required for computational reasons. Default:
`1e7`.

- `load_mh_file`

Tells Dynare to add to previous Metropolis-Hastings simulations instead
of starting from scratch. Since Dynare 4.5 the proposal density from the
previous run will automatically be loaded. In older versions, to assure
a neat continuation of the chain with the same proposal density, you
should provide the `mode_file` used in the previous run or the same
user-defined `mcmc_jumping_covariance` when using this option. Shouldn't
be used together with `mh_recover`. Note that under Octave, a neat
continuation of the chain with the last random number generator state of
the already present draws is currently not supported.

- `load_results_after_load_mh`

This option is available when loading a previous MCMC run without adding
additional draws, i.e. when `load_mh_file` is specified with
`mh_replic=0`. It tells Dynare to load the previously computed
convergence diagnostics, marginal data density, and posterior statistics
from an existing `_results` file instead of recomputing them.

- `mh_initialize_from_previous_mcmc`

This option allows to pick initial values for new MCMC from a previous
one, where the model specification, the number of estimated parameters,
(some) prior might have changed (so a situation where `load_mh_file`
would not work). If an additional parameter is estimated, it is
automatically initialized from prior\_draw. Note that, if this option is
used to skip the optimization step, you should use a sampling method
which does not require a proposal density, like slice. Otherwise,
optimization should always be done beforehand or a mode file with an
appropriate posterior covariance matrix should be used.

- `mh_initialize_from_previous_mcmc_directory = FILENAME`

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

- `mh_initialize_from_previous_mcmc_record = FILENAME`

If `mh_initialize_from_previous_mcmc` is set, and whenever the standard
file or directory tree is not applicable to load initial values, users
may directly provide here the path to the record file from which to load
values to be used to initialize the new MCMC.

- `mh_initialize_from_previous_mcmc_prior = FILENAME`

If `mh_initialize_from_previous_mcmc` is set, and whenever the standard
file or directory tree is not applicable to load initial values, users
may directly provide here the path to the prior definition file, to get
info in the priors used in previous MCMC.

- `optim = (NAME, VALUE, ...)`

A list of NAME and VALUE pairs. Can be used to set options for the
optimization routines. The set of available options depends on the
selected optimization routine (i.e. on the value of option
`mode_compute <mode_compute = INTEGER |
FUNCTION_NAME>`{.interpreted-text role="opt"}):

 `1, 3, 7, 12, 13`

    Available options are given in the documentation of the MATLAB
    Optimization Toolbox or in Octave's documentation.

 `2`

    Available options are:
  
    `'initial_step_length'`

    Initial step length. Default: `1`.

    `'initial_temperature'`

    Initial temperature. Default: `15`.

    `'MaxIter'`

    Maximum number of function evaluations. Default: `100000`.

    `'neps'`

    Number of final function values used to decide upon termination.
    Default: `10`.

    `'ns'`

    Number of cycles. Default: `10`.
   
    `'nt'`

    Number of iterations before temperature reduction. Default:
    `10`.

    `'step_length_c'`

    Step length adjustment. Default: `0.1`.

    `'TolFun'`

    Stopping criteria. Default: `1e-8`.

    `'rt'`

    Temperature reduction factor. Default: `0.1`.

    `'verbosity'`

    Controls verbosity of display during optimization, ranging from
    `0` (silent) to `3` (each function evaluation). Default: `1` 
    
    `4`

    Available options are:

    `'InitialInverseHessian'`

    Initial approximation for the inverse of the Hessian matrix of
    the posterior kernel (or likelihood). Obviously this
    approximation has to be a square, positive definite and
    symmetric matrix. Default: `'1e-4*eye(nx)'`, where nx is the
    number of parameters to be estimated.

    `'MaxIter'`

    Maximum number of iterations. Default: `1000`.

    `'NumgradAlgorithm'`

    Possible values are `2`, `3` and `5`, respectively,
    corresponding to the two, three and five points formula used to
    compute the gradient of the objective function (see *Abramowitz
    and Stegun (1964)*). Values `13` and `15` are more experimental.
    If perturbations on the right and the left increase the value of
    the objective function (we minimize this function) then we force
    the corresponding element of the gradient to be zero. The idea
    is to temporarily reduce the size of the optimization problem.
    Default: `2`.

    `'NumgradEpsilon'`

    Size of the perturbation used to compute numerically the
    gradient of the objective function. Default: `1e-6`.

    `'TolFun'`

    Stopping criteria. Default: `1e-7`.

    `'verbosity'`

    Controls verbosity of display during optimization. Set to `0` to
    set to silent. Default: `1`.

    `'SaveFiles'`

    Controls saving of intermediate results during optimization. Set
    to `0` to shut off saving. Default: `1`.

  `5`
 
    Available options are:

    `'Hessian'`

    Triggers three types of Hessian computations. `0`: outer product
    gradient; `1`: default Dynare Hessian routine; `2`: 'mixed' outer
    product gradient, where diagonal elements are obtained using
    second order derivation formula and outer product is used for
    correlation structure. Both {0} and {2} options require univariate
    filters, to ensure using maximum number of individual densities
    and a positive definite Hessian. Both {0} and {2} are quicker than
    default Dynare numeric Hessian, but provide decent starting values
    for Metropolis for large models (option {2} being more accurate
    than {0}). Default: `1`.

    `'MaxIter'`

    Maximum number of iterations. Default: `1000`.

    `'TolFun'`

    Stopping criteria. Default: `1e-5` for numerical derivatives,
    `1e-7` for analytic derivatives.

    `'verbosity'`

    Controls verbosity of display during optimization. Set to `0` to
    set to silent. Default: `1`.

    `'SaveFiles'`

    Controls saving of intermediate results during optimization. Set
    to `0` to shut off saving. Default: `1`.

 `6`

    Available options are:

    `'AcceptanceRateTarget'`

    A real number between zero and one. The scale parameter of the
    jumping distribution is adjusted so that the effective
    acceptance rate matches the value of option
    `'AcceptanceRateTarget'`. Default: `1.0/3.0`.

    `'InitialCovarianceMatrix'`

    Initial covariance matrix of the jumping distribution. Default
    is `'previous'` if option `mode_file` is used, `'prior'`
    otherwise.

    `'nclimb-mh'`

    Number of iterations in the last MCMC (climbing mode). Default:
    `200000`.

    `'ncov-mh'`

    Number of iterations used for updating the covariance matrix of
    the jumping distribution. Default: `20000`.

    `'nscale-mh'`

    Maximum number of iterations used for adjusting the scale
    parameter of the jumping distribution. Default: `200000`.

    `'NumberOfMh'`

    Number of MCMC run sequentially. Default: `3`.

  `8`

    Available options are:

    `'InitialSimplexSize'`

    Initial size of the simplex, expressed as percentage deviation
    from the provided initial guess in each direction. Default:
    `.05`.

    `'MaxIter'`

    Maximum number of iterations. Default: `5000`.

    `'MaxFunEvals'`

    Maximum number of objective function evaluations. No default.

    `'MaxFunvEvalFactor'`
 
    Set `MaxFunvEvals` equal to `MaxFunvEvalFactor` times the number
    of estimated parameters. Default: `500`.

    `'TolFun'`

    Tolerance parameter (w.r.t the objective function). Default:
    `1e-4`.

    `'TolX'`

    Tolerance parameter (w.r.t the instruments). Default: `1e-4`.

    `'verbosity'`

    Controls verbosity of display during optimization. Set to `0` to
    set to silent. Default: `1`.

  `9`

    Available options are:
 
    `'CMAESResume'`

    Resume previous run. Requires the `variablescmaes.mat` from the
    last run. Set to `1` to enable. Default: `0`.

    `'MaxIter'`

    Maximum number of iterations.

    `'MaxFunEvals'`

    Maximum number of objective function evaluations. Default:
    `Inf`.

    `'TolFun'`

    Tolerance parameter (w.r.t the objective function). Default:
    `1e-7`.

    `'TolX'`

    Tolerance parameter (w.r.t the instruments). Default: `1e-7`.

    `'verbosity'`

    Controls verbosity of display during optimization. Set to `0` to
    set to silent. Default: `1`.

    `'SaveFiles'`

    Controls saving of intermediate results during optimization. Set
    to `0` to shut off saving. Default: `1`.

  `10`

    Available options are:

    `'EndTemperature'`

    Terminal condition w.r.t the temperature. When the temperature
    reaches `EndTemperature`, the temperature is set to zero and the
    algorithm falls back into a standard simplex algorithm. Default:
    `0.1`.

    `'MaxIter'`

    Maximum number of iterations. Default: `5000`.

    `'MaxFunvEvals'`

    Maximum number of objective function evaluations. No default.

    `'TolFun'`
 
    Tolerance parameter (w.r.t the objective function). Default:
    `1e-4`.

    `'TolX'`

    Tolerance parameter (w.r.t the instruments). Default: `1e-4`.

    `'verbosity'`

    Controls verbosity of display during optimization. Set to `0` to
    set to silent. Default: `1`.

  `101`

    Available options are:

    `'LBGradientStep'`

    Lower bound for the stepsize used for the difference
    approximation of gradients. Default: `1e-11`.

    `'MaxIter'`

    Maximum number of iterations. Default: `15000`

    `'SpaceDilation'`

    Coefficient of space dilation. Default: `2.5`.

    `'TolFun'`

    Tolerance parameter (w.r.t the objective function). Default:
    `1e-6`.

    `'TolX'`

    Tolerance parameter (w.r.t the instruments). Default: `1e-6`.

    `'verbosity'`

    Controls verbosity of display during optimization. Set to `0` to
    set to silent. Default: `1`.

   `102`

    Available options are given in the documentation of the MATLAB
    Global Optimization Toolbox.

*Example*

To change the defaults of `csminwel` (`mode_compute=4`):

```
   estimation(..., mode_compute=4,optim=('NumgradAlgorithm',3,'TolFun',1e-5),...);
```

- `nodiagnostic`

Does not compute the convergence diagnostics for Metropolis-Hastings.
Default: diagnostics are computed and displayed.

- `bayesian_irf`

Triggers the computation of the posterior distribution of IRFs. The
length of the IRFs are controlled by the `irf` option. Results are
stored in `oo_.PosteriorIRF.dsge` (see below for a description of this
variable).

- `relative_irf`

See `relative_irf`{.interpreted-text role="opt"}.

- `posterior_sampling_method = NAME`

Selects the sampler used to sample from the posterior distribution
during Bayesian estimation. Default:`’random_walk_metropolis_hastings’`.

    `'random_walk_metropolis_hastings'`

    Instructs Dynare to use the Random-Walk Metropolis-Hastings. In this
    algorithm, the proposal density is recentered to the previous draw
    in every step.
 
    `'tailored_random_block_metropolis_hastings'`

    Instructs Dynare to use the Tailored randomized block (TaRB)
    Metropolis-Hastings algorithm proposed by *Chib and Ramamurthy
    (2010)* instead of the standard Random-Walk Metropolis-Hastings. In
    this algorithm, at each iteration the estimated parameters are
    randomly assigned to different blocks. For each of these blocks a
    mode-finding step is conducted. The inverse Hessian at this mode is
    then used as the covariance of the proposal density for a
    Random-Walk Metropolis-Hastings step. If the numerical Hessian is
    not positive definite, the generalized Cholesky decomposition of
    *Schnabel and Eskow (1990)* is used, but without pivoting. The
    TaRB-MH algorithm massively reduces the autocorrelation in the MH
    draws and thus reduces the number of draws required to
    representatively sample from the posterior. However, this comes at a
    computational cost as the algorithm takes more time to run.

    `'independent_metropolis_hastings'`

    Use the Independent Metropolis-Hastings algorithm where the proposal
    distribution - in contrast to the Random Walk Metropolis-Hastings
    algorithm - does not depend on the state of the chain.

    `'slice'`

    Instructs Dynare to use the Slice sampler of *Planas, Ratto, and
    Rossi (2015)*. Note that `'slice'` is incompatible with
    `prior_trunc=0`.

- `posterior_sampler_options = (NAME, VALUE, ...)`

A list of NAME and VALUE pairs. Can be used to set options for the
posterior sampling methods. The set of available options depends on the
selected posterior sampling routine (i.e. on the value of option
`posterior_sampling_method
<posterior_sampling_method = NAME>`{.interpreted-text role="opt"}):

    `'random_walk_metropolis_hastings'`

    Available options are:
 
    `'proposal_distribution'`

    Specifies the statistical distribution used for the proposal
    density.
 
    `'rand_multivariate_normal'`

    Use a multivariate normal distribution. This is the default.

    `'rand_multivariate_student'`

    Use a multivariate student distribution.

    `'student_degrees_of_freedom'`

    Specifies the degrees of freedom to be used with the multivariate
    student distribution. Default: `3`.

    `'use_mh_covariance_matrix'`

    Indicates to use the covariance matrix of the draws from a previous
    MCMC run to define the covariance of the proposal distribution.
    Requires the `load_mh_file`{.interpreted-text role="opt"} option to
    be specified. Default: `0`.

    `'scale_file'`

    Provides the name of a `_mh_scale.mat` file storing the tuned scale
    factor from a previous run of `mode_compute=6`.

    `'save_tmp_file'`

    Save the MCMC draws into a `_mh_tmp_blck` file at the refresh rate
    of the status bar instead of just saving the draws when the current
    `_mh*_blck` file is full. Default: `0`

    `'independent_metropolis_hastings'`

    Takes the same options as in the case of
   `random_walk_metropolis_hastings`.

    `'slice'`
 
    `'rotated'`

    Triggers rotated slice iterations using a covariance matrix from
    initial burn-in iterations. Requires either
    `use_mh_covariance_matrix` or `slice_initialize_with_mode`. Default:
    `0`.

    `'mode_files'`

    For multimodal posteriors, provide the name of a file containing a
    `nparam` by `nmodes` variable called `xparams` storing the different
    modes. This array must have one column vector per mode and the
    estimated parameters along the row dimension. With this info, the
    code will automatically trigger the `rotated` and `mode` options.
    Default: `[]`.

    `'slice_initialize_with_mode'`

    The default for slice is to set `mode_compute=0` and start the
    chain(s) from a random location in the prior space. This option
    first runs the mode-finder and then starts the chain from the mode.
    Together with `rotated`, it will use the inverse Hessian from the
    mode to perform rotated slice iterations. Default: `0`.

    `'initial_step_size'`
 
    Sets the initial size of the interval in the stepping-out procedure
    as fraction of the prior support, i.e. the size will be
    `initial_step_size * (UB-LB)`. `initial_step_size` must be a real
    number in the interval `[0,1]`. Default: `0.8`.

    `'use_mh_covariance_matrix'`

    See `use_mh_covariance_matrix <usemhcov>`{.interpreted-text
    role="ref"}. Must be used with `'rotated'`. Default: `0`.

    `'save_tmp_file'`

    See `save_tmp_file <savetmp>`{.interpreted-text role="ref"}.
    Default: `1`.

    `'tailored_random_block_metropolis_hastings'`

    `'proposal_distribution'`

    Specifies the statistical distribution used for the proposal
    density. See
    `proposal_distribution <prop_distrib>`{.interpreted-text
     role="ref"}.

    `new_block_probability = DOUBLE`

    Specifies the probability of the next parameter belonging to a new
    block when the random blocking in the TaRB Metropolis-Hastings
    algorithm is conducted. The higher this number, the smaller is the
    average block size and the more random blocks are formed during each
    parameter sweep. Default: `0.25`.

    `mode_compute = INTEGER`

    Specifies the mode-finder run in every iteration for every block of
    the TaRB Metropolis-Hastings algorithm. See
    `mode_compute <mode_compute =
     INTEGER | FUNCTION_NAME>`{.interpreted-text role="opt"}. Default:
    `4`.

    `optim = (NAME, VALUE,...)`

    Specifies the options for the mode-finder used in the TaRB
    Metropolis-Hastings algorithm. See `optim
    <optim = (NAME, VALUE, ...)>`{.interpreted-text role="opt"}.

    `'scale_file'`

    See `scale_file <scale-file>`{.interpreted-text role="ref"}..

    `'save_tmp_file'`

    See `save_tmp_file <savetmp>`{.interpreted-text role="ref"}.
    Default: `1`.

- `moments_varendo`

Triggers the computation of the posterior distribution of the
theoretical moments of the endogenous variables. Results are stored in
`oo_.PosteriorTheoreticalMoments` (see
`oo_.PosteriorTheoreticalMoments`{.interpreted-text role="mvar"}). The
number of lags in the autocorrelation function is controlled by the `ar`
option.

- `contemporaneous_correlation`

See `contemporaneous_correlation`{.interpreted-text role="opt"}. Results
are stored in `oo_.PosteriorTheoreticalMoments`. Note that the `nocorr`
option has no effect.

- `no_posterior_kernel_density`

Shuts off the computation of the kernel density estimator for the
posterior objects (see `density <dens>`{.interpreted-text role="ref"}
field).

```
conditional_variance_decomposition = INTEGER
conditional_variance_decomposition = [INTEGER1:INTEGER2]
conditional_variance_decomposition = [INTEGER1 INTEGER2 ...]
```

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

- `filtered_vars`

Triggers the computation of the posterior distribution of filtered
endogenous variables/one-step ahead forecasts, i.e. $E_{t}{y_{t+1}}$.
Results are stored in `oo_.FilteredVariables` (see below for a
description of this variable)

- `smoother`

Triggers the computation of the posterior distribution of smoothed
endogenous variables and shocks, i.e. the expected value of variables
and shocks given the information available in all observations up to the
final date ($E_{T}{y_t}$). Results are stored in
`oo_.SmoothedVariables`, `oo_.SmoothedShocks` and
`oo_.SmoothedMeasurementErrors`. Also triggers the computation of
`oo_.UpdatedVariables`, which contains the estimation of the expected
value of variables given the information available at the current date
($E_{t}{y_t}$). See below for a description of all these variables.

- `smoother_redux`

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

    `FilteredVariablesKStepAhead`: will be fully recovered

    `SmoothedVariables`, `FilteredVariables`, `UpdatedVariables`: recovered for all periods beyond period `d+1`,
    where `d` denotes the number of diffuse filtering steps.

    `FilteredVariablesKStepAheadVariances`, `Variance`, and
    `State_uncertainty` cannot be recovered, and ZERO is provided as
    output.

If you need variances for those variables, either do not set the option,
or declare the variable as observed, using NaNs as data points.

- `forecast = INTEGER`

Computes the posterior distribution of a forecast on INTEGER periods
after the end of the sample used in estimation. If no
Metropolis-Hastings is computed, the result is stored in variable
`oo_.forecast` and corresponds to the forecast at the posterior mode. If
a Metropolis-Hastings is computed, the distribution of forecasts is
stored in variables `oo_.PointForecast` and `oo_.MeanForecast`. See
`fore`{.interpreted-text role="ref"}, for a description of these
variables.

- `tex`

See `tex`{.interpreted-text role="opt"}.

- `kalman_algo = INTEGER`

    `0`

    Automatically use the Multivariate Kalman Filter for stationary models
    and the Multivariate Diffuse Kalman Filter for non-stationary models.

    `1`

    Use the Multivariate Kalman Filter.

    `2`
   
    Use the Univariate Kalman Filter.

    `3`

    Use the Multivariate Diffuse Kalman Filter.

    `4`
  
    Use the Univariate Diffuse Kalman Filter.

Default value is `0`. In case of missing observations of single or all
series, Dynare treats those missing values as unobserved states and
uses the Kalman filter to infer their value (see e.g. *Durbin and
Koopman (2012)*, Ch. 4.10) This procedure has the advantage of being
capable of dealing with observations where the forecast error variance
matrix becomes singular for some variable(s). If this happens, the
respective observation enters with a weight of zero in the
log-likelihood, i.e. this observation for the respective variable(s)
is dropped from the likelihood computations (for details see *Durbin
and Koopman (2012)*, Ch. 6.4 and 7.2.5 and *Koopman and Durbin
(2000)*). If the use of a multivariate Kalman filter is specified and
a singularity is encountered, Dynare by default automatically switches
to the univariate Kalman filter for this parameter draw. This behavior
can be changed via the`use_univariate_filters_if_singularity_is_detected
<use_univariate_filters_if_singularity_is_detected = INTEGER>`{.interpreted-text
role="opt"} option.

- `fast_kalman_filter`

Select the fast Kalman filter using Chandrasekhar recursions as
described by `Herbst (2015)`. This setting is only used with
`kalman_algo=1` or `kalman_algo=3`. In case of using the diffuse Kalman
filter (`kalman_algo=3/lik_init=3`), the observables must be stationary.
This option is not yet compatible with
`analytic_derivation`{.interpreted-text role="opt"}.

- `kalman_tol = DOUBLE`

Numerical tolerance for determining the singularity of the covariance
matrix of the prediction errors during the Kalman filter (minimum
allowed reciprocal of the matrix condition number). Default value is
`1e-10`.

- `diffuse_kalman_tol = DOUBLE`

Numerical tolerance for determining the singularity of the covariance
matrix of the prediction errors ($F_{\infty}$) and the rank of the
covariance matrix of the non-stationary state variables ($P_{\infty}$)
during the Diffuse Kalman filter. Default value is `1e-6`.

- `filter_covariance`

Saves the series of one step ahead error of forecast covariance
matrices. With Metropolis, they are saved in
`oo_.FilterCovariance`{.interpreted-text role="mvar"}, otherwise in
`oo_.Smoother.Variance`{.interpreted-text role="mvar"}. Saves also
k-step ahead error of forecast covariance matrices if
`filter_step_ahead` is set.

- `filter_step_ahead = [INTEGER1:INTEGER2] filter_step_ahead = [INTEGER1 INTEGER2 ...]`

Triggers the computation k-step ahead filtered values, i.e.
$E_{t}{y_{t+k}}$. Stores results in `oo_.FilteredVariablesKStepAhead`.
Also stores 1-step ahead values in `oo_.FilteredVariables`.
`oo_.FilteredVariablesKStepAheadVariances` is stored if
`filter_covariance`.

- `filter_decomposition`

Triggers the computation of the shock decomposition of the above k-step
ahead filtered values. Stores results in
`oo_.FilteredVariablesShockDecomposition`.

- `smoothed_state_uncertainty`

Triggers the computation of the variance of smoothed estimates, i.e.
$var_T(y_t)$. Stores results in `oo_.Smoother.State_uncertainty`.

- `diffuse_filter`

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

$$a_t = a_{t-1} + e_t$$

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

- `selected_variables_only`

Only run the classical smoother on the variables listed just after the
`estimation` command. This option is incompatible with requesting
classical frequentist forecasts and will be overridden in this case.
When using Bayesian estimation, the smoother is by default only run on
the declared endogenous variables. Default: run the smoother on all the
declared endogenous variables.

- `cova_compute = INTEGER`

When `0`, the covariance matrix of estimated parameters is not computed
after the computation of posterior mode (or maximum likelihood). This
increases speed of computation in large models during development, when
this information is not always necessary. Of course, it will break all
successive computations that would require this covariance matrix.
Otherwise, if this option is equal to `1`, the covariance matrix is
computed and stored in variable `hh` of `MODEL_FILENAME_mode.mat`.
Default is `1`.

- `solve_algo = INTEGER`

See `solve_algo <solvalg>`{.interpreted-text role="ref"}.

- `order = INTEGER`

Order of approximation around the deterministic steady state. When
greater than 1, the likelihood is evaluated with a particle or nonlinear
filter (see *Fernández-Villaverde and Rubio-Ramírez (2005)*). Default is
`1`, i.e. the likelihood of the linearized model is evaluated using a
standard Kalman filter.

- `irf = INTEGER`

See `irf <irf = INTEGER>`{.interpreted-text role="opt"}. Only used if
`bayesian_irf`{.interpreted-text role="opt"} is passed.

- `irf_shocks = ( VARIABLE_NAME [[,] VARIABLE_NAME ...] )`

See `irf_shocks <irf_shocks = ( VARIABLE_NAME [[,]
VARIABLE_NAME ...] )>`{.interpreted-text role="opt"}. Only used if
`bayesian_irf`{.interpreted-text role="opt"} is passed.

- `irf_plot_threshold = DOUBLE`

See `irf_plot_threshold <irf_plot_threshold =
DOUBLE>`{.interpreted-text role="opt"}. Only used if
`bayesian_irf`{.interpreted-text role="opt"} is passed.

- `aim_solver`

See `aim_solver`{.interpreted-text role="opt"}.

- `sylvester = OPTION`

See `sylvester <sylvester = OPTION>`{.interpreted-text role="opt"}.

- `sylvester_fixed_point_tol = DOUBLE`

See `sylvester_fixed_point_tol <sylvester_fixed_point_tol
= DOUBLE>`{.interpreted-text role="opt"} .

- `lyapunov = OPTION`

Determines the algorithm used to solve the Lyapunov equation to
initialized the variance-covariance matrix of the Kalman filter using
the steady-state value of state variables. Possible values for OPTION
are:

    `default`

    Uses the default solver for Lyapunov equations based on
    Bartels-Stewart algorithm.

    `fixed_point`

    Uses a fixed point algorithm to solve the Lyapunov equation. This
    method is faster than the `default` one for large scale models, but
    it could require a large amount of iterations.

    `doubling`

    Uses a doubling algorithm to solve the Lyapunov equation
    (`disclyap_fast`). This method is faster than the two previous one
    for large scale models.

    `square_root_solver`

    Uses a square-root solver for Lyapunov equations (`dlyapchol`). This
    method is fast for large scale models (available under MATLAB if the
    Control System Toolbox is installed; available under Octave if the
    [control](https://octave.sourceforge.io/control/) package from
    Octave-Forge is installed)

    Default value is `default`.

- `lyapunov_fixed_point_tol = DOUBLE`

This is the convergence criterion used in the fixed point Lyapunov
solver. Its default value is `1e-10`.

- `lyapunov_doubling_tol = DOUBLE`

This is the convergence criterion used in the doubling algorithm to
solve the Lyapunov equation. Its default value is `1e-16`.

- `use_penalized_objective_for_hessian`

Use the penalized objective instead of the objective function to compute
numerically the hessian matrix at the mode. The penalties decrease the
value of the posterior density (or likelihood) when, for some
perturbations, Dynare is not able to solve the model (issues with steady
state existence, Blanchard and Kahn conditions, ...). In pratice, the
penalized and original objectives will only differ if the posterior mode
is found to be near a region where the model is ill-behaved. By default
the original objective function is used.

- `analytic_derivation`

Triggers estimation with analytic gradient at `order=1`. The final
hessian at the mode is also computed analytically. Only works for
stationary models without missing observations, i.e. for
`kalman_algo<3`. Optimizers that rely on analytic gradients are
`mode_compute=1,3,4,5,101`.

- `ar = INTEGER`

See `ar <ar = INTEGER>`{.interpreted-text role="opt"}. Only useful in
conjunction with option `moments_varendo`.

- `endogenous_prior`

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

- `use_univariate_filters_if_singularity_is_detected = INTEGER`

Decide whether Dynare should automatically switch to univariate filter
if a singularity is encountered in the likelihood computation (this is
the behaviour if the option is equal to `1`). Alternatively, if the
option is equal to `0`, Dynare will not automatically change the filter,
but rather use a penalty value for the likelihood when such a
singularity is encountered. Default: `1`.

- `keep_kalman_algo_if_singularity_is_detected`

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

- `rescale_prediction_error_covariance`

Rescales the prediction error covariance in the Kalman filter to avoid
badly scaled matrix and reduce the probability of a switch to univariate
Kalman filters (which are slower). By default no rescaling is done.

- `qz_zero_threshold = DOUBLE`

See `qz_zero_threshold <qz_zero_threshold = DOUBLE>`{.interpreted-text
role="opt"}.

- `taper_steps = [INTEGER1 INTEGER2 ...]`

Percent tapering used for the spectral window in the *Geweke
(1992,1999)* convergence diagnostics (requires
`mh_nblocks=1 <mh_nblocks = INTEGER>`{.interpreted-text role="opt"}).
The tapering is used to take the serial correlation of the posterior
draws into account. Default: `[4 8 15]`.

- `geweke_interval = [DOUBLE DOUBLE]`

Percentage of MCMC draws at the beginning and end of the MCMC chain
taken to compute the *Geweke (1992,1999)* convergence diagnostics
(requires `mh_nblocks=1 <mh_nblocks =
INTEGER>`{.interpreted-text role="opt"}) after discarding the first
`mh_drop = DOUBLE
<mh_drop>`{.interpreted-text role="opt"} percent of draws as a burnin.
Default: \[0.2 0.5\].

- `raftery_lewis_diagnostics`

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

- `raftery_lewis_qrs = [DOUBLE DOUBLE DOUBLE]`

Sets the quantile of the CDF `q` that is estimated with precision `r`
with a probability `s` in the *Raftery and Lewis (1992)* convergence
diagnostics. Default: `[0.025 0.005 0.95]`.

- `consider_all_endogenous`

Compute the posterior moments, smoothed variables, k-step ahead filtered
variables and forecasts (when requested) on all the endogenous
variables. This is equivalent to manually listing all the endogenous
variables after the `estimation` command.

- `consider_all_endogenous_and_auxiliary`

Compute the posterior moments, smoothed variables, k-step ahead filtered
variables and forecasts (when requested) on all the endogenous variables
and the auxiliary variables introduced by the preprocessor. This option
is useful when e.g. running `smoother2histval` on the results of the
Kalman smoother.

- `consider_only_observed`

Compute the posterior moments, smoothed variables, k-step ahead filtered
variables and forecasts (when requested) on all the observed variables.
This is equivalent to manually listing all the observed variables after
the `estimation` command.


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

    `MOMENT_NAME`

    This field can take the following values:

    `HPDinf`

    Lower bound of a 90% HPD interval.[^4]

    `HPDsup`

    Upper bound of a 90% HPD interval.

    `HPDinf_ME`

    Lower bound of a 90% HPD interval[^5] for observables when taking
    measurement error into account (see e.g. *Christoffel et al.
    (2010*), p.17).

    `HPDsup_ME`

    Upper bound of a 90% HPD interval for observables when taking
    measurement error into account.

    `Mean`

    Mean of the posterior distribution.

    `Median`

    Median of the posterior distribution.

    `Std`

    Standard deviation of the posterior distribution.

    `Variance`

    Variance of the posterior distribution.

    `deciles`

    Deciles of the distribution.

    `density`
   
    Non parametric estimate of the posterior density following the
    approach outlined in *Skoeld and Roberts (2003)*. First and second
    columns are respectively abscissa and ordinate coordinates.

    `ESTIMATED_OBJECT`

    This field can take the following values:

    `measurement_errors_corr`

    Correlation between two measurement errors.

    `measurement_errors_std`

    Standard deviation of measurement errors.

    `parameters`

    Parameters.

    `shocks_corr`

    Correlation between two structural shocks.

    `shocks_std`

    Standard deviation of structural shocks.

*MATLAB/Octave Variables*: `oo.MarginalDensity.LaplaceApproximation`

Variable set by the `estimation` command. Stores the marginal data
density based on the Laplace Approximation.

*MATLAB/Octave Variables*: `oo.MarginalDensity.ModifiedHarmonicMean`

Variable set by the `estimation command`, if it is used with
`mh_replic > 0` or `load_mh_file` option. Stores the marginal data
density based on *Geweke (1999)* Modified Harmonic Mean estimator.

*MATLAB/Octave Variables*: `oo.posterior.optimization`

Variable set by the `estimation` command if mode-finding is used. Stores
the results at the mode. Fields are of the form:

```
    oo_.posterior.optimization.OBJECT
```

where OBJECT is one of the following:

    `mode`

    Parameter vector at the mode.

    `Variance`

    Inverse Hessian matrix at the mode or MCMC jumping covariance matrix
    when used with the `MCMC_jumping_covariance <mcmc_jumping_covariance
    = OPTION>`{.interpreted-text role="opt"} option.

    `log_density`

    Log likelihood (ML)/log posterior density (Bayesian) at the mode
    when used with `mode_compute>0`.

*MATLAB/Octave Variables*: `oo.posterior.metropolis`

Variable set by the `estimation` command if `mh_replic>0` is used.
Fields are of the form:

```
    oo_.posterior.metropolis.OBJECT
```

where OBJECT is one of the following:

    `mean`

    Mean parameter vector from the MCMC.

    `Variance`

    Covariance matrix of the parameter draws in the MCMC.

*MATLAB/Octave Variables*: `oo.FilteredVariables`

Variable set by the `estimation` command, if it is used with the
`filtered_vars` option.

After an estimation without Metropolis, fields are of the form:

```
    oo_.FilteredVariables.VARIABLE_NAME
```

After an estimation with Metropolis, fields are of the form:
```
    oo_.FilteredVariables.MOMENT_NAME.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.FilteredVariablesKStepAhead`

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

*MATLAB/Octave Variables*: `oo.FilteredVariablesKStepAheadVariances`

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

*MATLAB/Octave Variables*: `oo.Filtered_Variables_X_step_ahead`

Variable set by the `estimation` command, if it is used with the
`filter_step_ahead option` in the context of Bayesian estimation. Fields
are of the form:

```
    oo_.Filtered_Variables_X_step_ahead.VARIABLE_NAME
```

The n-th entry stores the k-step ahead filtered variable computed at
time n for time n+k.

*MATLAB/Octave Variables*: `oo.FilteredVariablesShockDecomposition`

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

*MATLAB/Octave Variables*: `oo.PosteriorIRF.dsge`

Variable set by the `estimation` command, if it is used with the
`bayesian_irf` option. Fields are of the form:

```
    oo_.PosteriorIRF.dsge.MOMENT_NAME.VARIABLE_NAME_SHOCK_NAME
```

*MATLAB/Octave Variables*: `oo.SmoothedMeasurementErrors`

Variable set by the `estimation` command, if it is used with the
`smoother` option. Fields are of the form:

```
    oo_.SmoothedMeasurementErrors.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.SmoothedShocks`

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command.

After an estimation without Metropolis, or if computed by
`calib_smoother`, fields are of the form:

```
    oo_.SmoothedShocks.VARIABLE_NAME
```

After an estimation with Metropolis, fields are of the form:

```
    oo_.SmoothedShocks.MOMENT_NAME.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.SmoothedVariables`

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command.

After an estimation without Metropolis, or if computed by
`calib_smoother`, fields are of the form:

```
    oo_.SmoothedVariables.VARIABLE_NAME
```

After an estimation with Metropolis, fields are of the form:

```
    oo_.SmoothedVariables.MOMENT_NAME.VARIABLE_NAME
```

*MATLAB/Octave Command*: `get_smooth ('VARIABLE_NAME' [, 'VARIABLE_NAME']...);

Returns the smoothed values of the given endogenous or exogenous
variable(s), as they are stored in the `oo_.SmoothedVariables` and
`oo_.SmoothedShocks` variables.

*MATLAB/Octave Variables*: `oo.UpdatedVariables`

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command. Contains the estimation of
the expected value of variables given the information available at the
current date.

After an estimation without Metropolis, or if computed by
`calib_smoother`, fields are of the form:

```
    oo_.UpdatedVariables.VARIABLE_NAME
```

After an estimation with Metropolis, fields are of the form:

```
    oo_.UpdatedVariables.MOMENT_NAME.VARIABLE_NAME
```

*MATLAB/Octave Command*: `get_update ('VARIABLE_NAME' [, 'VARIABLE_NAME']...);`

Returns the updated values of the given variable(s), as they are stored
in the `oo_.UpdatedVariables` variable.

*MATLAB/Octave Variables*: `oo.FilterCovariance`

Three-dimensional array set by the `estimation` command if used with the
`smoother` and Metropolis, if the `filter_covariance` option has been
requested. Contains the series of one-step ahead forecast error
covariance matrices from the Kalman smoother. The `M_.endo_nbr` times
`M_.endo_nbr` times `T+1` array contains the variables in declaration
order along the first two dimensions. The third dimension of the array
provides the observation for which the forecast has been made. Fields
are of the form:

```
    oo_.FilterCovariance.MOMENT_NAME
```

Note that density estimation is not supported.

*MATLAB/Octave Variables*: `oo.Smoother.Variance`

Three-dimensional array set by the `estimation` command (if used with
the `smoother`) without Metropolis, or by the `calib_smoother` command,
if the `filter_covariance` option has been requested. Contains the
series of one-step ahead forecast error covariance matrices from the
Kalman smoother. The `M_.endo_nbr` times `M_.endo_nbr` times `T+1` array
contains the variables in declaration order along the first two
dimensions. The third dimension of the array provides the observation
for which the forecast has been made.

*MATLAB/Octave Variables*: `oo.Smoother.State_uncertainty`

Three-dimensional array set by the `estimation` command (if used with
the `smoother` option) without Metropolis, or by the `calib_smoother`
command, if the `smoothed_state_uncertainty` option has been requested.
Contains the series of covariance matrices for the state estimate given
the full data from the Kalman smoother. The `M_.endo_nbr` times
`M_.endo_nbr` times `T` array contains the variables in declaration
order along the first two dimensions. The third dimension of the array
provides the observation for which the smoothed estimate has been made.

*MATLAB/Octave Variables*: `oo.Smoother.SteadyState`

Variable set by the `estimation` command (if used with the `smoother`)
without Metropolis, or by the `calib_smoother` command. Contains the
steady state component of the endogenous variables used in the smoother
in order of variable declaration.

*MATLAB/Octave Variables*: `oo.Smoother.TrendCoeffs`

Variable set by the `estimation` command (if used with the `smoother`)
without Metropolis, or by the `calib_smoother` command. Contains the
trend coefficients of the observed variables used in the smoother in
order of declaration of the observed variables.

*MATLAB/Octave Variables*: `oo.Smoother.Trend`

Variable set by the `estimation command` (if used with the `smoother`
option), or by the `calib_smoother` command. Contains the trend
component of the variables used in the smoother.

Fields are of the form:

```
    oo_.Smoother.Trend.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.Smoother.Constant`

Variable set by the `estimation` command (if used with the `smoother`
option), or by the `calib_smoother` command. Contains the constant part
of the endogenous variables used in the smoother, accounting e.g. for
the data mean when using the prefilter option.

Fields are of the form:

```
    oo_.Smoother.Constant.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.Smoother.loglinear`

Indicator keeping track of whether the smoother was run with the
`loglinear <logl>`{.interpreted-text role="ref"} option and thus whether
stored smoothed objects are in logs.

*MATLAB/Octave Variables*: `oo.PosteriorTheoreticalMoments`

Variable set by the `estimation` command, if it is used with the
`moments_varendo` option. Fields are of the form:

```
    oo_.PosteriorTheoreticalMoments.dsge.THEORETICAL_MOMENT.ESTIMATED_OBJECT.MOMENT_NAME.VARIABLE_NAME
```

where *THEORETICAL_MOMENT* is one of the following:

    `covariance`

    Variance-covariance of endogenous variables.

    `contemporaneous_correlation`

    Contemporaneous correlation of endogenous variables when the 
    `contemporaneous_correlation`{.interpreted-text role="opt"} option
    is specified.

    `correlation`

    Auto- and cross-correlation of endogenous variables. Fields are
    vectors with correlations from 1 up to order `options_.ar`.

    `VarianceDecomposition`
   
    Decomposition of variance (unconditional variance, i.e. at horizon
    infinity).[^6]

    `VarianceDecompositionME`

    Same as [VarianceDecomposition](), but contains the decomposition of
    the measured as opposed to the actual variable. The joint
    contribution of the measurement error will be saved in a field named
    `ME`.

    `ConditionalVarianceDecomposition`
     
    Only if the `conditional_variance_decomposition` option has been
    specified. In the presence of measurement error, the field will
    contain the variance contribution after measurement error has been
    taken out, i.e. the decomposition will be conducted of the actual as
    opposed to the measured variables.

    `ConditionalVarianceDecompositionME`

    Only if the `conditional_variance_decomposition` option has been
    specified. Same as [ConditionalVarianceDecomposition](), but
    contains the decomposition of the measured as opposed to the actual
    variable. The joint contribution of the measurement error will be
    saved in a field names `ME`.


*MATLAB/Octave Variables*: `oo.posterior_density`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
    oo_.posterior_density.PARAMETER_NAME
```

*MATLAB/Octave Variables*: `oo.posterior_hpdinf`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
    oo_.posterior_hpdinf.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.posterior_hpdsup`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
oo_.posterior_hpdsup.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.posterior_mean`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
    oo_.posterior_mean.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.posterior_mode`

Variable set by the `estimation` command during mode-finding. Fields are
of the form:

```
    oo_.posterior_mode.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variable*: `oo.posterior_std_at_mode`

Variable set by the `estimation` command during mode-finding. It is
based on the inverse Hessian at `oo_.posterior_mode`. Fields are of the
form:

```
    oo_.posterior_std_at_mode.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variable*: `oo.posterior_std`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
    oo_.posterior_std.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variable*: `oo.posterior_var`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
    oo_.posterior_var.ESTIMATED_OBJECT.VARIABLE_NAME
```

*MATLAB/Octave Variable*: `oo.posterior_median`

Variable set by the `estimation` command, if it is used with
`mh_replic > 0` or `load_mh_file` option. Fields are of the form:

```
    oo_.posterior_median.ESTIMATED_OBJECT.VARIABLE_NAME
```

*Example*

Here are some examples of generated variables:

```
    oo_.posterior_mode.parameters.alp
    oo_.posterior_mean.shocks_std.ex
    oo_.posterior_hpdsup.measurement_errors_corr.gdp_conso
```

*MATLAB/Octave Variable*: `oo.dsge_var.posterior_mode`

Structure set by the `dsge_var` option of the `estimation` command after
mode_compute.

The following fields are saved:

    `PHI_tilde`

    Stacked posterior DSGE-BVAR autoregressive matrices at the mode
    (equation (28) of *Del Negro and Schorfheide (2004)*).

    `SIGMA_u_tilde`

    Posterior covariance matrix of the DSGE-BVAR at the mode (equation
    (29) of *Del Negro and Schorfheide (2004)*).

    `iXX`

    Posterior population moments in the DSGE-BVAR at the mode (
    $inv(\lambda T \Gamma_{XX}^*+ X'X)$).

    `prior`

    Structure storing the DSGE-BVAR prior.

    `PHI_star`

    Stacked prior DSGE-BVAR autoregressive matrices at the mode
    (equation (22) of *Del Negro and Schorfheide (2004)*).

    `SIGMA_star`
 
    Prior covariance matrix of the DSGE-BVAR at the mode (equation (23)
    of *Del Negro and Schorfheide (2004)*).

    `ArtificialSampleSize`

    Size of the artifical prior sample ( $inv(\lambda T)$).

    `DF`

    Prior degrees of freedom ( $inv(\lambda T-k-n)$).

    `iGXX_star`

    Inverse of the theoretical prior "covariance" between X and X
    ($\Gamma_{xx}^*$ in *Del Negro and Schorfheide (2004)*).

*MATLAB/Octave Variable*: `oo.RecursiveForecast`

Variable set by the `forecast` option of the `estimation` command when
used with the nobs = [INTEGER1:INTEGER2] option (see
`nobs <nobs = [INTEGER1:INTEGER2]>`{.interpreted-text role="opt"}).

Fields are of the form:

```
    oo_.RecursiveForecast.FORECAST_OBJECT.VARIABLE_NAME
```

where `FORECAST_OBJECT` is one of the following[^7] :

    `Mean`

    Mean of the posterior forecast distribution.

    `HPDinf/HPDsup`

    Upper/lower bound of the 90% HPD interval taking into account only
    parameter uncertainty (corresponding to
    `oo_.MeanForecast`{.interpreted-text role="mvar"}).

    `HPDTotalinf/HPDTotalsup`.

    Upper/lower bound of the 90% HPD interval taking into account both
    parameter and future shock uncertainty (corresponding to
    `oo_.PointForecast`{.interpreted-text role="mvar"})

    `VARIABLE_NAME` contains a matrix of the following size: number of time
    periods for which forecasts are requested using the
    `nobs = [INTEGER1:INTEGER2]` option times the number of forecast
    horizons requested by the forecast option. i.e., the row indicates the
    period at which the forecast is performed and the column the respective
    k-step ahead forecast. The starting periods are sorted in ascending
    order, not in declaration order.

*MATLAB/Octave Variable*: `oo.convergence.geweke`

Variable set by the convergence diagnostics of the `estimation` command
when used with `mh_nblocks=1` option (see
`mh_nblocks <mh_nblocks = INTEGER>`{.interpreted-text role="opt"}).

Fields are of the form:

```
    oo_.convergence.geweke.VARIABLE_NAME.DIAGNOSTIC_OBJECT
```

where *DIAGNOSTIC_OBJECT* is one of the following:

    `posteriormean`

    Mean of the posterior parameter distribution.

    `posteriorstd`

    Standard deviation of the posterior parameter distribution.

    `nse_iid`

    Numerical standard error (NSE) under the assumption of iid draws.

    `rne_iid`

    Relative numerical efficiency (RNE) under the assumption of iid draws.

    `nse_x`

    Numerical standard error (NSE) when using an x% taper.

    `rne_x`

    Relative numerical efficiency (RNE) when using an x% taper.

    `pooled_mean`

    Mean of the parameter when pooling the beginning and end parts of the
    chain specified in `geweke_interval
    <geweke_interval = [DOUBLE DOUBLE]>`{.interpreted-text role="opt"} and
    weighting them with their relative precision. It is a vector
    containing the results under the iid assumption followed by the ones
    using the `taper_steps` option (see `taper_steps <taper_steps
    = [INTEGER1 INTEGER2 ...]>`{.interpreted-text role="opt"}).

    `pooled_nse`

    NSE of the parameter when pooling the beginning and end parts of the
    chain and weighting them with their relative precision. See
    `pooled_mean`.

    `prob_chi2_test`

    p-value of a chi-squared test for equality of means in the beginning
    and the end of the MCMC chain. See `pooled_mean`. A value above 0.05
    indicates that the null hypothesis of equal means and thus convergence
    cannot be rejected at the 5 percent level. Differing values along the
    `taper_steps` signal the presence of significant autocorrelation in
    draws. In this case, the estimates using a higher tapering are usually
    more reliable.

*Command*: `unit_root_vars VARIABLE_NAME...;`

This command is deprecated. Use `estimation` option `diffuse_filter`
instead for estimating a model with non-stationary observed variables or
`steady` option `nocheck` to prevent `steady` to check the steady state
returned by your steady state file.

Dynare also has the ability to estimate Bayesian VARs:

*Command*: `bvar_density ;`

Computes the marginal density of an estimated BVAR model, using
Minnesota priors.

See `bvar-a-la-sims.pdf`, which comes with Dynare distribution, for more
information on this command.


## Shock Decomposition

*Command*: `shock_decomposition [VARIABLE_NAME]...; 

*Command*: `shock_decomposition (OPTIONS...) [VARIABLE_NAME]...;`

This command computes the historical shock decomposition for a given
sample based on the Kalman smoother, i.e. it decomposes the historical
deviations of the endogenous variables from their respective steady
state values into the contribution coming from the various shocks. The
`variable_names` provided govern for which variables the decomposition
is plotted.

Note that this command must come after either `estimation` (in case of
an estimated model) or `stoch_simul` (in case of a calibrated model).

*Options*

- `parameter_set = OPTION`

Specify the parameter set to use for running the smoother. Possible
values for OPTION are:

-   `calibration`
-   `prior_mode`
-   `prior_mean`
-   `posterior_mode`
-   `posterior_mean`
-   `posterior_median`
-   `mle_mode`

Note that the parameter set used in subsequent commands like
`stoch_simul` will be set to the specified `parameter_set`. Default
value: `posterior_mean` if Metropolis has been run, `mle_mode` if MLE
has been run.

- `datafile = FILENAME`

See `datafile <dataf>`{.interpreted-text role="ref"}. Useful when
computing the shock decomposition on a calibrated model.

- `first_obs = INTEGER`

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.

- `nobs = INTEGER`

See `nobs <nobs = INTEGER>`{.interpreted-text role="opt"}.

- `prefilter = INTEGER`

See `prefilter <prefilter = INTEGER>`{.interpreted-text role="opt"}.

- `loglinear`

See `loglinear <loglinear>`{.interpreted-text role="opt"}.

- `diffuse_kalma_tol = DOUBLE`

See `diffuse_kalman_tol <diffuse_kalman_tol = DOUBLE>`{.interpreted-text
role="opt"}.

- `diffuse_filter`

See `diffuse_filter <diffuse_filter>`{.interpreted-text role="opt"}.

- `xls_sheet = QUOTED_STRING`

See `xls_sheet <xls_sheet = QUOTED_STRING>`{.interpreted-text
role="opt"}.

- `xls_range = RANGE`

See `xls_range <xls_range = RANGE>`{.interpreted-text role="opt"}.

- `use_shock_groups [= NAME]`

Uses shock grouping defined by the string instead of individual shocks
in the decomposition. The groups of shocks are defined in the
`shock_groups`{.interpreted-text role="bck"} block. If no group name is
given, `default` is assumed.

- `colormap = VARIABLE_NAME`

Controls the `colormap` used for the shocks decomposition graphs.
VARIABLE_NAME must be the name of a MATLAB/Octave variable that has
been declared beforehand and whose value will be passed to the
MATLAB/Octave `colormap` function (see the MATLAB/Octave manual for the
list of acceptable values).

- `nograph`

See `nograph`{.interpreted-text role="opt"}. Suppresses the display and
creation only within the `shock_decomposition` command, but does not
affect other commands. See `plot_shock_decomposition`{.interpreted-text
role="comm"} for plotting graphs.

- `init_state = BOOLEAN`

If equal to 0, the shock decomposition is computed conditional on the
smoothed state variables in period `0`, i.e. the smoothed shocks
starting in period 1 are used. If equal to `1`, the shock decomposition
is computed conditional on the smoothed state variables in period 1.
Default: `0`.

- `wit_epilogue`

If set, then also compute the decomposition for variables declared in
the `epilogue` block (see `epilogue`{.interpreted-text role="ref"}).

*Output*

*MATLAB/Octave Variable*: `oo.shock_decomposition`

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

*Block*: `shock_groups ;` 

*Block*: `shock_groups(OPTIONS...);`

Shocks can be regrouped for the purpose of shock decomposition. The
composition of the shock groups is written in a block delimited by
`shock_groups` and `end`.

Each line defines a group of shocks as a list of exogenous variables:

```
    SHOCK_GROUP_NAME   = VARIABLE_1 [[,] VARIABLE_2 [,]...];
    'SHOCK GROUP NAME' = VARIABLE_1 [[,] VARIABLE_2 [,]...];
```

*Options*

- `name = NAME`

Specifies a name for the following definition of shock groups. It is
possible to use several `shock_groups` blocks in a model file, each
grouping being identified by a different name. This name must in turn be
used in the `shock_decomposition` command. If no name is given,
`default` is used.

*Example*

```
     varexo e_a, e_b, e_c, e_d;
     ...

     shock_groups(name=group1);
     supply = e_a, e_b;
     'aggregate demand' = e_c, e_d;
     end;

     shock_decomposition(use_shock_groups=group1);
```

This example defines a shock grouping with the name `group1`,
containing a set of supply and demand shocks and conducts the shock
decomposition for these two groups.

*Command*: `realtime_shock_decomposition [VARIABLE_NAME]...;`

*Command*: `realtime_shock_decomposition (OPTIONS...) [VARIABLE_NAME]...;`

This command computes the realtime historical shock decomposition for a
given sample based on the Kalman smoother. For each period
$T=[\texttt{presample},\ldots,\texttt{nobs}]$, it recursively computes
three objects:

-   Real-time historical shock decomposition $Y(t\vert T)$ for
    $t=[1,\ldots,T]$, i.e. without observing data in
    $[T+1,\ldots,\texttt{nobs}]$. This results in a standard shock
    decomposition being computed for each additional datapoint
    becoming available after `presample`.
-   Forecast shock decomposition $Y(T+k\vert T)$ for
    $k=[1,\ldots,forecast]$, i.e. the $k$-step ahead forecast made for
    every $T$ is decomposed in its shock contributions.
-   Real-time conditional shock decomposition of the difference
    between the real-time historical shock decomposition and the
    forecast shock decomposition. If `vintage <vintage =
    INTEGER>`{.interpreted-text role="opt"} is equal to `0`, it
    computes the effect of shocks realizing in period $T$, i.e.
    decomposes $Y(T\vert T)-Y(T\vert T-1)$. Put differently, it
    conducts a $1$-period ahead shock decomposition from $T-1$ to $T$,
    by decomposing the update step of the Kalman filter. If
    `vintage>0` and smaller than `nobs`, the decomposition is
    conducted of the forecast revision
    $Y(T+k\vert T+k)-Y(T+k\vert T)$.

Like `shock_decomposition`{.interpreted-text role="comm"} it decomposes
the historical deviations of the endogenous variables from their
respective steady state values into the contribution coming from the
various shocks. The `variable_names` provided govern for which variables
the decomposition is plotted.

Note that this command must come after either `estimation` (in case of
an estimated model) or `stoch_simul` (in case of a calibrated model).

*Options*

- `parameter_set = OPTION`

See `parameter_set <parameter_set = OPTION>`{.interpreted-text
role="opt"} for possible values.

- `datafile = FILENAME`

See `datafile <dataf>`{.interpreted-text role="ref"}.

- `first_obs = INTEGER`

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.

- `nobs = INTEGER`

See `nobs <nobs = INTEGER>`{.interpreted-text role="opt"}.

- `use_shock_groups [= NAME]`

See `use_shock_groups <use_shock_groups [= NAME]>`{.interpreted-text
role="opt"}.

- `colormap = VARIABLE_NAME`

See `colormap <colormap = VARIABLE_NAME>`{.interpreted-text role="opt"}.

- `nograph`

See `nograph`{.interpreted-text role="opt"}. Only shock decompositions
are computed and stored in `oo_.realtime_shock_decomposition`,
`oo_.conditional_shock_decomposition` and
`oo_.realtime_forecast_shock_decomposition` but no plot is made (See
`plot_shock_decomposition`{.interpreted-text role="comm"}).

- `presample = INTEGER`

Data point above which recursive realtime shock decompositions are
computed, *i.e.* for $T=[\texttt{presample+1} \ldots \texttt{nobs}]$.

- `forecast = INTEGER`

Compute shock decompositions up to $T+k$ periods, i.e. get shock
contributions to k-step ahead forecasts.

- `save_realtime = INTEGER_VECTOR`

Choose for which vintages to save the full realtime shock decomposition.
Default: `0`.

- `fast_realtime = INTEGER fast_realtime = [INTEGER1:INTEGER2]`

- `fast_realtime = [INTEGER1 INTEGER2 ...]`

Runs the smoother only for the data vintages provided by the specified
integer (vector).

- `with_epilogue`

See `with_epilogue`{.interpreted-text role="opt"}.

*Output*

*MATLAB/Octave Variable*: `oo.realtime_shock_decomposition`

Structure storing the results of realtime historical decompositions.
Fields are three-dimensional arrays with the first two dimension equal
to the ones of `oo_.shock_decomposition`{.interpreted-text role="mvar"}.
The third dimension stores the time periods and is therefore of size
`T+forecast`. Fields are of the form:

```
    oo_.realtime_shock_decomposition.OBJECT
```

where OBJECT is one of the following:

    `pool`

    Stores the pooled decomposition, i.e. for every real-time shock
    decomposition terminal period
    $T=[\texttt{presample},\ldots,\texttt{nobs}]$ it collects the last
    period's decomposition $Y(T\vert T)$ (see also
    `plot_shock_decomposition`{.interpreted-text role="comm"}). The
    third dimension of the array will have size `nobs+forecast`.

    `time_*`

    Stores the vintages of realtime historical shock decompositions if
    `save_realtime` is used. For example, if `save_realtime=[5]` and
    `forecast=8`, the third dimension will be of size `13`.

*MATLAB/Octave Variable*: `oo.realtime_conditional_shock_decomposition`

Structure storing the results of real-time conditional decompositions.
Fields are of the form:

```
    oo_.realtime_conditional_shock_decomposition.OBJECT
```

where OBJECT is one of the following:

    `pool`

    Stores the pooled real-time conditional shock decomposition, i.e.
    collects the decompositions of $Y(T\vert T)-Y(T\vert T-1)$ for the
    terminal periods $T=[\texttt{presample},\ldots,\texttt{nobs}]$. The
    third dimension is of size `nobs`.

    `time_*`

    Store the vintages of $k$-step conditional forecast shock
    decompositions $Y(t\vert T+k)$, for $t=[T \ldots T+k]$. See `vintage
    <vintage = INTEGER>`{.interpreted-text role="opt"}. The third
    dimension is of size `1+forecast`.

*MATLAB/Octave Variable*: `oo.realtime_forecast_shock_decomposition`

Structure storing the results of realtime forecast decompositions.
Fields are of the form:

```
    oo_.realtime_forecast_shock_decomposition.OBJECT
```

where `OBJECT` is one of the following:
   
    `pool`

    Stores the pooled real-time forecast decomposition of the $1$-step
    ahead effect of shocks on the $1$-step ahead prediction, i.e.
    $Y(T\vert
    T-1)$.

    `time_*`

    Stores the vintages of $k$-step out-of-sample forecast shock
    decompositions, i.e. $Y(t\vert
    T)$, for $t=[T \ldots T+k]$. See `vintage
    <vintage = INTEGER>`{.interpreted-text role="opt"}.

- `plot_shock_decomposition [VARIABLE_NAME]...;`

- `plot_shock_decomposition (OPTIONS...) [VARIABLE_NAME]...;`

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

- `use_shock_groups [= NAME]`

See `use_shock_groups <use_shock_groups [= NAME]>`{.interpreted-text
role="opt"}.

- `colormap = VARIABLE_NAME`

See `colormap <colormap = VARIABLE_NAME>`{.interpreted-text role="opt"}.

- `nodisplay`

See `nodisplay`{.interpreted-text role="opt"}.

- `nograph`

See `nograph`{.interpreted-text role="opt"}.

- `graph_format = FORMAT graph_format = ( FORMAT, FORMAT... )`

See `graph_format <graph_format = FORMAT>`{.interpreted-text
role="opt"}.

- `detail_plot`

Plots shock contributions using subplots, one per shock (or group of
shocks). Default: not activated

- `interactive`

Under MATLAB, add uimenus for detailed group plots. Default: not
activated

- `screen_shocks`

For large models (i.e. for models with more than 16 shocks), plots only
the shocks that have the largest historical contribution for chosen
selected `variable_names`. Historical contribution is ranked by the mean
absolute value of all historical contributions.

- `steadystate`

If passed, the the $y$-axis value of the zero line in the shock
decomposition plot is translated to the steady state level. Default: not
activated

- `type = qoq | yoy | aoa`

For quarterly data, valid arguments are: `qoq` for quarter-on-quarter
plots, `yoy` for year-on-year plots of growth rates, `aoa` for
annualized variables, i.e. the value in the last quarter for each year
is plotted. Default value: empty, i.e. standard period-on-period plots
(`qoq` for quarterly data).

- `fig_name = STRING`

Specifies a user-defined keyword to be appended to the default figure
name set by `plot_shock_decomposition`. This can avoid to overwrite
plots in case of sequential calls to `plot_shock_decomposition`.

- `write_xls`

Saves shock decompositions to Excel-file in the main directory, named
`FILENAME_shock_decomposition_TYPE_FIG_NAME.xls`. This option requires
your system to be configured to be able to write Excel files.[^8]

- `realtime = INTEGER`

Which kind of shock decomposition to plot. INTEGER can take the
following values:

-   `0`: standard historical shock decomposition. See
    `shock_decomposition`{.interpreted-text role="comm"}.
-   `1`: realtime historical shock decomposition. See
    `realtime_shock_decomposition`{.interpreted-text role="comm"}.
-   `2`: conditional realtime shock decomposition. See
    `realtime_shock_decomposition`{.interpreted-text role="comm"}.
-   `3`: realtime forecast shock decomposition. See
    `realtime_shock_decomposition`{.interpreted-text role="comm"}.

If no vintage is requested, i.e. `vintage=0` then the pooled objects
from `realtime_shock_decomposition`{.interpreted-text role="comm"} will
be plotted and the respective vintage otherwise. Default: `0`.

- `vintage = INTEGER`

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

-   `realtime=1`: the full vintage shock decomposition $Y(t\vert T)$
    for $t=[1,\ldots,T]$
-   `realtime=2`: the conditional forecast shock decomposition from
    $T$, i.e. plots $Y(T+j\vert T+j)$ and the shock contributions
    needed to get to the data $Y(T+j)$ conditional on $T=$ vintage,
    with $j=[0,\ldots,\texttt{forecast}]$.
-   `realtime=3`: plots unconditional forecast shock decomposition
    from $T$, i.e. $Y(T+j\vert
    T)$, where $T=\texttt{vintage}$ and
    $j=[0,\ldots,\texttt{forecast}]$.

Default: `0`.

- `plot_init_date = DATE`

If passed, plots decomposition using `plot_init_date` as initial period.
Default: first observation in estimation

- `plot_end_date = DATE`

If passed, plots decomposition using `plot_end_date` as last period.
Default: last observation in estimation

- `diff`

If passed, plot the decomposition of the first difference of the list of
variables. If used in combination with `flip`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated

- `flip`

If passed, plot the decomposition of the opposite of the list of
variables. If used in combination with `diff`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated

- `max_nrows`

Maximum number of rows in the subplot layout of detailed shock
decomposition graphs. Note that columns are always 3. Default: 6

- `with_epilogue`

See `with_epilogue`{.interpreted-text role="opt"}.

- `init2shocks init2shocks = NAME`

Use the information contained in an `init2shocks`{.interpreted-text
role="bck"} block, in order to attribute initial conditions to shocks.
The name of the block can be explicitly given, otherwise it defaults to
the `default` block.

*Block*: `init2shocks ;`

*Block*: `init2shocks (OPTIONS...);`

This blocks gives the possibility of attributing the initial condition
of endogenous variables to the contribution of exogenous variables in
the shock decomposition.

For example, in an AR(1) process, the contribution of the initial
condition on the process variable can naturally be assigned to the
innovation of the process.

Each line of the block should have the syntax:

```
    VARIABLE_1 [,] VARIABLE_2;
```

Where VARIABLE_1 is an endogenous variable whose initial condition will
be attributed to the exogenous VARIABLE_2.

The information contained in this block is used by the
`plot_shock_decomposition`{.interpreted-text role="comm"} command when
given the `init2shocks` option.

*Options*

- `name = NAME`

Specifies a name for the block, that can be referenced from
`plot_shock_decomposition`, so that several such blocks can coexist in a
single model file. If the name is unspecified, it defaults to `default`.

*Example*

```
     var y y_s R pie dq pie_s de A y_obs pie_obs R_obs;
     varexo e_R e_q e_ys e_pies e_A;
     ...

     model;
       dq = rho_q*dq(-1)+e_q;
       A = rho_A*A(-1)+e_A;
       ...
     end;

     ...

     init2shocks;
       dq e_q;
       A e_A;
     end;

     shock_decomposition(nograph);

     plot_shock_decomposition(init2shocks) y_obs R_obs pie_obs dq de;
```

In this example, the initial conditions of `dq` and `A` will be
respectively attributed to `e_q` and `e_A`.

*Command*: `initial_condition_decomposition [VARIABLE_NAME]...;`

*Command*: `initial_condition_decomposition (OPTIONS...) [VARIABLE_NAME]...;`

This command computes and plots the decomposition of the effect of
smoothed initial conditions of state variables. The `variable_names`
provided govern which variables the decomposition is plotted for.

Further note that, unlike the majority of Dynare commands, the options
specified below are overwritten with their defaults before every call to
`initial_condition_decomposition`. Hence, if you want to reuse an option
in a subsequent call to `initial_condition_decomposition`, you must pass
it to the command again.

*Options*

- `colormap = VARIABLE_NAME`

See `colormap <colormap = VARIABLE_NAME>`{.interpreted-text role="opt"}.

- `nodisplay`

See `nodisplay`{.interpreted-text role="opt"}.

- `graph_format = FORMAT graph_format = ( FORMAT, FORMAT... )`

See `graph_format <graph_format = FORMAT>`{.interpreted-text
role="opt"}.

- `detail_plot`

Plots shock contributions using subplots, one per shock (or group of
shocks). Default: not activated

- `steadystate`

If passed, the the $y$-axis value of the zero line in the shock
decomposition plot is translated to the steady state level. Default: not
activated

- `type = qoq | yoy | aoa`

For quarterly data, valid arguments are: `qoq` for quarter-on-quarter
plots, `yoy` for year-on-year plots of growth rates, `aoa` for
annualized variables, i.e. the value in the last quarter for each year
is plotted. Default value: empty, i.e. standard period-on-period plots
(`qoq` for quarterly data).

- `fig_name = STRING`

Specifies a user-defined keyword to be appended to the default figure
name set by `plot_shock_decomposition`. This can avoid to overwrite
plots in case of sequential calls to `plot_shock_decomposition`.

- `write_xls`

Saves shock decompositions to Excel-file in the main directory, named
`FILENAME_shock_decomposition_TYPE_FIG_NAME_initval.xls`. This option
requires your system to be configured to be able to write Excel
files.[^9]

- `plot_init_date = DATE`

If passed, plots decomposition using `plot_init_date` as initial period.
Default: first observation in estimation

- `plot_end_date = DATE`

If passed, plots decomposition using `plot_end_date` as last period.
Default: last observation in estimation

- `diff`

If passed, plot the decomposition of the first difference of the list of
variables. If used in combination with `flip`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated

- `flip`

If passed, plot the decomposition of the opposite of the list of
variables. If used in combination with `diff`{.interpreted-text
role="opt"}, the `diff` operator is first applied. Default: not
activated

*Command*: `squeeze_shock_decomposition [VARIABLE_NAME]...;`

For large models, the size of the information stored by shock
decompositions (especially various settings of realtime decompositions)
may become huge. This command allows to squeeze this information in two
possible ways:

-   Automatic (default): only the variables for which plotting has
    been explicitly required with `plot_shock_decomposition` will have
    their decomposition left in `oo_` after this command is run;
-   If a list of variables is passed to the command, then only those
    variables will have their decomposition left in `oo_` after this
    command is run.

## Calibrated Smoother

Dynare can also run the smoother on a calibrated model:

*Command*: `calib_smoother [VARIABLE_NAME]...;` 

*Command*: `calib_smoother (OPTIONS...)[VARIABLE_NAME]...;`

This command computes the smoothed variables (and possible the filtered
variables) on a calibrated model.

A datafile must be provided, and the observable variables declared with
`varobs`. The smoother is based on a first-order approximation of the
model.

By default, the command computes the smoothed variables and shocks and
stores the results in `oo_.SmoothedVariables` and `oo_.SmoothedShocks`.
It also fills `oo_.UpdatedVariables`.

*Options*

- `datafile = FILENAME`

See `datafile <dataf>`{.interpreted-text role="ref"}.

- `filtered_vars`

Triggers the computation of filtered variables. See
`filtered_vars`{.interpreted-text role="opt"}, for more details.

- `filter_step_ahead = [INTEGER1:INTEGER2]`

See
`filter_step_ahead <filter_step_ahead = [INTEGER1:INTEGER2]>`{.interpreted-text
role="opt"}.

- `prefilter = INTEGER`

See `prefilter <prefilter = INTEGER>`{.interpreted-text role="opt"}.

- `parameter_set = OPTION`

See `parameter_set <parameter_set = OPTION>`{.interpreted-text
role="opt"} for possible values. Default: `calibration`.

- `loglinear`

See `loglinear <logl>`{.interpreted-text role="ref"}.

- `first_obs = INTEGER`

See `first_obs <first_obs = INTEGER>`{.interpreted-text role="opt"}.

- `filter_decomposition`

See `filter_decomposition`{.interpreted-text role="opt"}.

- `filter_covariance`

See `filter_covariance`{.interpreted-text role="opt"}.

- `smoother_redux`

See `smoother_redux`{.interpreted-text role="opt"}.

- `kalman_algo = INTEGER`

See `kalman_algo <kalman_algo = INTEGER>`{.interpreted-text role="opt"}.

- `diffuse_filter = INTEGER`

See `diffuse_filter`{.interpreted-text role="opt"}.

- `diffuse_kalman_tol = DOUBLE`

See `diffuse_kalman_tol <diffuse_kalman_tol = DOUBLE>`{.interpreted-text
role="opt"}.

- xls_sheet = QUOTED_STRING`

See `xls_sheet <xls_sheet = QUOTED_STRING>`{.interpreted-text
role="opt"}.

- `xls_range = RANGE`

See `xls_range <xls_range = RANGE>`{.interpreted-text role="opt"}.

## Forecasting

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

*Command*: `forecast [VARIABLE_NAME...];`

*Command*: `forecast (OPTIONS...)[VARIABLE_NAME...];`

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

- `periods = INTEGER`

Number of periods of the forecast. Default: `5`.

- `conf_sig = DOUBLE`

Level of significance for confidence interval. Default: `0.90`.

- `nograph`

See `nograph`{.interpreted-text role="opt"}.

- `nodisplay`

See `nodisplay`{.interpreted-text role="opt"}.

- `graph_format = FORMAT graph_format = ( FORMAT, FORMAT... )`

See `graph_format = FORMAT`{.interpreted-text role="opt"}.

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

```
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
```

*MATLAB/Octave Variables*: `oo.forecast`

Variable set by the `forecast` command, or by the `estimation` command
if used with the `forecast` option and ML or if no Metropolis-Hastings
has been computed (in that case, the forecast is computed for the
posterior mode). Fields are of the form:

```
    oo_.forecast.FORECAST_MOMENT.VARIABLE_NAME
```

where `FORECAST_MOMENT` is one of the following:

    `HPDinf`

    Lower bound of a 90% HPD interval[^10] of forecast due to parameter
    uncertainty, but ignoring the effect of measurement error on
    observed variables. In case of ML, it stores the lower bound of the
    confidence interval.

    `HPDsup`

    Upper bound of a 90% HPD forecast interval due to parameter
    uncertainty, but ignoring the effect of measurement error on
    observed variables. In case of ML, it stores the upper bound of the
    confidence interval.

   `HPDinf_ME`

    Lower bound of a 90% HPD interval[^11] of forecast for observed
    variables due to parameter uncertainty and measurement error. In
    case of ML, it stores the lower bound of the confidence interval.

    `HPDsup_ME`

    Upper bound of a 90% HPD interval of forecast for observed variables
    due to parameter uncertainty and measurement error. In case of ML,
    it stores the upper bound of the confidence interval.

    `Mean`

    Mean of the posterior distribution of forecasts.

*MATLAB/Octave Variables*: `oo.PointForecast`

Set by the `estimation` command, if it is used with the `forecast`
option and if either `mh_replic > 0` or the `load_mh_file` option are
used.

Contains the distribution of forecasts taking into account the
uncertainty about both parameters and shocks.

Fields are of the form:

```
    oo_.PointForecast.MOMENT_NAME.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `oo.MeanForecast`

Set by the `estimation` command, if it is used with the `forecast`
option and if either `mh_replic > 0` or `load_mh_file` option are used.

Contains the distribution of forecasts where the uncertainty about
shocks is averaged out. The distribution of forecasts therefore only
represents the uncertainty about parameters.

Fields are of the form:

```
    oo_.MeanForecast.MOMENT_NAME.VARIABLE_NAME
```

*Command*: `conditional_forecast (OPTIONS...);`

This command computes forecasts on an estimated or calibrated model for
a given constrained path of some future endogenous variables. This is
done using the reduced form first order state-space representation of
the DSGE model by finding the structural shocks that are needed to match
the restricted paths. Consider the augmented state space representation
that stacks both predetermined and non-predetermined variables into a
vector $y_{t}$:

$$y_t=Ty_{t-1}+R\varepsilon_t$$

Both $y_t$ and $\varepsilon_t$ are split up into controlled and
uncontrolled ones, and we assume without loss of generality that the
constrained endogenous variables and the controlled shocks come first :

$$\begin{aligned}
\begin{pmatrix}
 y_{c,t}\\
 y_{u,t}
 \end{pmatrix}
 =
 \begin{pmatrix}
 T_{c,c} & T_{c,u}\\
 T_{u,c} & T_{u,u}
 \end{pmatrix}
 \begin{pmatrix}
 y_{c,t-1}\\
 y_{u,t-1}
 \end{pmatrix}
 +
 \begin{pmatrix}
 R_{c,c} & R_{c,u}\\
 R_{u,c} & R_{u,u}
 \end{pmatrix}
 \begin{pmatrix}
 \varepsilon_{c,t}\\
 \varepsilon_{u,t}
 \end{pmatrix}
 \end{aligned}$$

where matrices $T$ and $R$ are partitioned consistently with the vectors
of endogenous variables and innovations. Provided that matrix $R_{c,c}$
is square and full rank (a necessary condition is that the number of
free endogenous variables matches the number of free innovations), given
$y_{c,t}$, $\varepsilon_{u,t}$ and $y_{t-1}$ the first block of
equations can be solved for $\varepsilon_{c,t}$:

$$\varepsilon_{c,t} = R_{c,c}^{-1}\bigl( y_{c,t} - T_{c,c}y_{c,t} - T_{c,u}y_{u,t}  - R_{c,u}\varepsilon_{u,t}\bigr)$$

and $y_{u,t}$ can be updated by evaluating the second block of
equations:

$$y_{u,t} = T_{u,c}y_{c,t-1} + T_{u,u}y_{u,t-1} +  R_{u,c}\varepsilon_{c,t} + R_{u,u}\varepsilon_{u,t}$$

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

- `parameter_set = OPTION`

See `parameter_set <parameter_set = OPTION>`{.interpreted-text
role="opt"} for possible values. No default value, mandatory option.

- `controlled_varexo = (VARIABLE_NAME...)`

Specify the exogenous variables to use as control variables. No default
value, mandatory option.

- `periods = INTEGER`

Number of periods of the forecast. Default: `40`. `periods` cannot be
smaller than the number of constrained periods.

- `replic = INTEGER`

Number of simulations used to compute the conditional forecast
uncertainty. Default: `5000`.

- `conf_sig = DOUBLE`

Level of significance for confidence interval. Default: `0.80`.

*Output*

The results are stored in `oo_.conditional_forecast`, which is described
below.

*Example*

```
     var y a;
     varexo e u;
     ...
     estimation(...);

     conditional_forecast_paths;
     var y;
     periods 1:3, 4:5;
     values 2, 5;
     var a;
     periods 1:5;
     values 3;
     end;

     conditional_forecast(parameter_set = calibration, controlled_varexo = (e, u), replic = 3000);
     plot_conditional_forecast(periods = 10) a y;
```

*MATLAB/Octave Variables*: `oo.conditional_forecast.cond`

Variable set by the `conditional_forecast` command. It stores the
conditional forecasts. Fields are `periods+1` by `1` vectors storing the
steady state (time 0) and the subsequent `periods` forecasts periods.
Fields are of the form:

```
    oo_.conditional_forecast.cond.FORECAST_MOMENT.VARIABLE_NAME
```

where FORECAST_MOMENT is one of the following:

    `Mean`

    Mean of the conditional forecast distribution.

    `ci`

    Confidence interval of the conditional forecast distribution. The
    size corresponds to `conf_sig`.

*MATLAB/Octave Variables*: `oo.conditional_forecast.uncond`

Variable set by the `conditional_forecast` command. It stores the
unconditional forecasts. Fields are of the form:

```
    oo_.conditional_forecast.uncond.FORECAST_MOMENT.VARIABLE_NAME
```

*MATLAB/Octave Variables*: `forecasts.instruments`

Variable set by the `conditional_forecast command`. Stores the names of
the exogenous instruments.

*MATLAB/Octave Variables*: `oo.conditional_forecast.controlled_variables`

Variable set by the `conditional_forecast` command. Stores the position
of the constrained endogenous variables in declaration order.

*MATLAB/Octave Variables*: `oo.conditional_forecast.controlled_exo_variables`

Variable set by the `conditional_forecast` command. Stores the values of
the controlled exogenous variables underlying the conditional forecasts
to achieve the constrained endogenous variables. Fields are
`[number of constrained periods]` by `1` vectors and are of the form:

```
    oo_.conditional_forecast.controlled_exo_variables.FORECAST_MOMENT.SHOCK_NAME
```

*MATLAB/Octave Variables*: `oo.conditional_forecast.graphs`

Variable set by the `conditional_forecast` command. Stores the
information for generating the conditional forecast plots.

*Block*: `conditional_forecast_paths ;`

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

*Block*: `filter_initial_state ;`

This block specifies the initial values of the endogenous states at the
beginning of the Kalman filter recursions. That is, if the Kalman filter
recursion starts with time t=1 being the first observation, this block
provides the state estimate at time 0 given information at time 0,
$E_0(x_0)$. If nothing is specified, the initial condition is assumed to
be at the steady state (which is the unconditional mean for a stationary
model).

This block is terminated by `end;`.

Each line inside of the block should be of the form:

```
    VARIABLE_NAME(INTEGER)=EXPRESSION;
```

`EXPRESSION` is any valid expression returning a numerical value and can
contain parameter values. This allows specifying relationships that will
be honored during estimation. `INTEGER` refers to the lag with which a
variable appears. By convention in Dynare, period 1 is the first period.
Going backwards in time, the first period before the start of the
simulation is period 0, then period -1, and so on. Note that the
`filter_initial_state` block does not take non-state variables.

*Example*

```
     filter_initial_state;
     k(0)= ((1/bet-(1-del))/alp)^(1/(alp-1))*l_ss;
     P(0)=2.5258;
     m(0)= mst;
     end;
```

*Command*: `plot_conditional_forecast [VARIABLE_NAME...];`

*Command*: `plot_conditional_forecast (periods = INTEGER) [VARIABLE_NAME...];`

Plots the conditional (plain lines) and unconditional (dashed lines)
forecasts.

To be used after `conditional_forecast`.

*Options*

- `periods = INTEGER`

Number of periods to be plotted. Default: equal to periods in
`conditional_forecast`. The number of periods declared in
`plot_conditional_forecast` cannot be greater than the one declared in
`conditional_forecast`.

- `bvar_forecast ;`

This command computes (out-of-sample) forecasts for an estimated BVAR
model, using Minnesota priors.

See `bvar-a-la-sims.pdf`, which comes with Dynare distribution, for more
information on this command.

If the model contains strong non-linearities or if some perfectly
expected shocks are considered, the forecasts and the conditional
forecasts can be computed using an extended path method. The forecast
scenario describing the shocks and/or the constrained paths on some
endogenous variables should be build. The first step is the forecast
scenario initialization using the function `init_plan`:

- `HANDLE = init_plan (DATES);`

Creates a new forecast scenario for a forecast period (indicated as a
dates class, see `dates class members
<dates-members>`{.interpreted-text role="ref"}). This function return a
handle on the new forecast scenario.

The forecast scenario can contain some simple shocks on the exogenous
variables. This shocks are described using the function `basic_plan`:

*MATLAB/Octave Command*: 

```
HANDLE = basic_plan (HANDLE, `VAR_NAME', `SHOCK_TYPE', DATES, MATLAB VECTOR OF DOUBLE | [DOUBLE | EXPR [DOUBLE | EXPR] ] );
```

Adds to the forecast scenario a shock on the exogenous variable
indicated between quotes in the second argument. The shock type has to
be specified in the third argument between quotes: 'surprise' in case of
an unexpected shock or 'perfect_foresight' for a perfectly anticipated
shock. The fourth argument indicates the period of the shock using a
dates class (see `dates class
members <dates-members>`{.interpreted-text role="ref"}). The last
argument is the shock path indicated as a MATLAB vector of double. This
function return the handle of the updated forecast scenario.

The forecast scenario can also contain a constrained path on an
endogenous variable. The values of the related exogenous variable
compatible with the constrained path are in this case computed. In other
words, a conditional forecast is performed. This kind of shock is
described with the function `flip_plan`:

*MATLAB/Octave Command*: 

```
HANDLE = flip_plan (HANDLE, `VAR_NAME', `VAR_NAME', `SHOCK_TYPE', DATES, MATLAB VECTOR OF DOUBLE | [DOUBLE | EXPR [DOUBLE | EXPR] ] );
```

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

Once the forecast scenario if fully described, the forecast is computed
with the command `det_cond_forecast`:

*MATLAB/Octave Command*:
```
DSERIES = det_cond_forecast (HANDLE[, DSERIES [, DATES]]);
```

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

*Example*

```
     % conditional forecast using extended path method
     % with perfect foresight on r path

     var y r;
     varexo e u;
     ...
     smoothed = dseries('smoothed_variables.csv');

     fplan = init_plan(2013Q4:2029Q4);
     fplan = flip_plan(fplan, 'y', 'u', 'surprise', 2013Q4:2014Q4,  [1 1.1 1.2 1.1 ]);
     fplan = flip_plan(fplan, 'r', 'e', 'perfect_foresight', 2013Q4:2014Q4,  [2 1.9 1.9 1.9 ]);

     dset_forecast = det_cond_forecast(fplan, smoothed);

     plot(dset_forecast.{'y','u'});
     plot(dset_forecast.{'r','e'});
```

*Command*: `smoother2histval ;` 

*Command*: `smoother2histval(OPTIONS...);`

The purpose of this command is to construct initial conditions (for a
subsequent simulation) that are the smoothed values of a previous
estimation.

More precisely, after an estimation run with the `smoother` option,
`smoother2histval` will extract the smoothed values (from
`oo_.SmoothedVariables`, and possibly from `oo_.SmoothedShocks` if there
are lagged exogenous), and will use these values to construct initial
conditions (as if they had been manually entered through `histval`).

*Options*

- `period = INTEGER`

Period number to use as the starting point for the subsequent
simulation. It should be between 1 and the number of observations that
were used to produce the smoothed values. Default: the last observation.

- `infile = FILENAME`

Load the smoothed values from a `_results.mat` file created by a
previous Dynare run. Default: use the smoothed values currently in the
global workspace.

- `invars = ( VARIABLE_NAME [VARIABLE_NAME ...] )`

A list of variables to read from the smoothed values. It can contain
state endogenous variables, and also exogenous variables having a lag.
Default: all the state endogenous variables, and all the exogenous
variables with a lag.

- `outfile = FILENAME`

Write the initial conditions to a file. Default: write the initial
conditions in the current workspace, so that a simulation can be
performed.

- `outvars = ( VARIABLE_NAME [VARIABLE_NAME ...] )`

A list of variables which will be given the initial conditions. This
list must have the same length than the list given to `invars`, and
there will be a one-to-one mapping between the two list. Default: same
value as option `invars`.

*Use cases*

There are three possible ways of using this command:

-   Everything in a single file: run an estimation with a smoother,
    then run `smoother2histval` (without the `infile` and `outfile`
    options), then run a stochastic simulation.
-   In two files: in the first file, run the smoother and then run
    `smoother2histval` with the `outfile` option; in the second file,
    run `histval_file` to load the initial conditions, and run a
    (deterministic or stochastic) simulation.
-   In two files: in the first file, run the smoother; in the second
    file, run `smoother2histval` with the `infile` option equal to the
    `_results.mat` file created by the first file, and then run a
    (deterministic or stochastic) simulation.

## Optimal policy

Dynare has tools to compute optimal policies for various types of
objectives. You can either solve for optimal policy under commitment
with `ramsey_model`, for optimal policy under discretion with
`discretionary_policy` or for optimal simple rules with `osr` (also
implying commitment).

*Command*: `planner_objective MODEL_EXPRESSION ;`

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

*Command*: `evaluate_planner_objective;` 

*Command*: `evaluate_planner_objective (OPTIONS...);`

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

- `periods = INTEGER`

The value of the option specifies the number of periods to use in the
simulations in the computation of unconditional welfare at higher order.

Default: `10000`.

- `drop = INTEGER`

The number of burn-in draws out of `periods` discarded before computing
the unconditional welfare at higher order. Default: `1000`.

*Example (stochastic context)*

```
     var a ...;
     varexo u;

     model;
     a = rho*a(-1)+u+u(-1);
     ...
     end;

     histval;
     u(0)=1;
     a(0)=-1;
     end;

     shocks;
     var u; stderr 0.008;
     var u;
     periods 1;
     values 1;
     end;

     evaluate_planner_objective;
```

*MATLAB/Octave Variables*: `oo.planner_objective_value.unconditional`

Scalar storing the value of unconditional welfare. In a perfect
foresight context, it corresponds to welfare in the long-run,
approximated as welfare in the terminal simulation period.

*MATLAB/Octave Variables*: `oo.planner_objective_value.conditional`

In a perfect foresight context, this field will be a scalar storing the
value of welfare conditional on the specified initial condition and zero
initial Lagrange multipliers.

In a stochastic context, it will have two subfields:

*MATLAB/Octave Variables*: `oo.planner_objective_value.conditional.steady_initial_multiplier`

Stores the value of the planner objective when the initial Lagrange
multipliers associated with the planner's problem are set to their
steady state values (see `ramsey_policy`{.interpreted-text
role="comm"}).

*MATLAB/Octave Variables*: `oo.planner_objective_value.conditional.zero_initial_multiplier`

Stores the value of the planner objective when the initial Lagrange
multipliers associated with the planner's problem are set to 0, i.e. it
is assumed that the planner exploits its ability to surprise private
agents in the first period of implementing Ramsey policy. This value
corresponds to the planner implementing optimal policy for the first
time and committing not to re-optimize in the future.

### Optimal policy under commitment (Ramsey)

Dynare allows to automatically compute optimal policy choices of a
Ramsey planner who takes the specified private sector equilibrium
conditions into account and commits to future policy choices. Doing so
requires specifying the private sector equilibrium conditions in the
`model`-block and a `planner_objective` as well as potentially some
`instruments` to facilitate computations.

!!! note

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

Consider a perfect foresight example where the Euler equation for the
return to capital is given by

```
     1/C=beta*1/C(+1)*(R(+1)+(1-delta))
```

The job of the Ramsey planner in period `1` is to choose $C_1$ and
$R_1$, taking as given $C_0$. The above equation may seemingly
equivalently be written as

```
     1/C=beta*1/C(+1)*(R_cap);
     R_cap=R(+1)+(1-delta);
```

due to perfect foresight. However, this changes the problem of the
Ramsey planner in the first period to choosing $C_1$ and $R_1$, taking
as given both $C_0$ and $R^{cap}_0$. Thus, the relevant return to
capital in the Euler equation of the first period is not a choice of
the planner anymore due to the forward-looking nature of the
definition in the second line!

A correct specification would be to instead define `R_cap` as a
model-local variable:

```
     1/C=beta*1/C(+1)*(R_cap);
     #R_cap=R(+1)+(1-delta);
```

*Command*: `ramsey_model (OPTIONS...);`  

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

- `planner_discount = EXPRESSION`

Declares or reassigns the discount factor of the central planner
`optimal_policy_discount_factor`. Default: `1.0`.

- `planner_discount_latex_name = LATEX_NAME`

Sets the LaTeX name of the `optimal_policy_discount_factor` parameter.

- `instruments = (VARIABLE_NAME,...)`

Declares instrument variables for the computation of the steady state
under optimal policy. Requires a `steady_state_model` block or a
`_steadystate.m` file. See below.

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

*Block*: `ramsey_constraints ;`

This block lets you define constraints on the variables in the Ramsey
problem. The constraints take the form of a variable, an inequality
operator (> or <) and a constant.

*Example*

```
     ramsey_constraints;
     i > 0;
     end;
```

*Command*: `ramsey_policy [VARIABLE_NAME...];`

*Command*: `ramsey_policy (OPTIONS...)[VARIABLE_NAME...];`

This command is deprecated and formally equivalent to the calling
sequence

```
     ramsey_model;
     stoch_simul;
     evaluate_planner_objective;
```

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

- planner_discount = EXPRESSION`

See `planner_discount <planner_discount = EXPRESSION>`{.interpreted-text
role="opt"}.

- `instruments = (VARIABLE_NAME,...)`

Declares instrument variables for the computation of the steady state
under optimal policy. Requires a `steady_state_model` block or a
`_steadystate.m` file. See below.

*Output*

This command generates all the output variables of `stoch_simul`. For
specifying the initial values for the endogenous state variables (except
for the Lagrange multipliers), see above.

*Steady state*

See `Ramsey steady state <ramsey_model>`{.interpreted-text role="comm"}.

### Optimal policy under discretion

*Command*: `discretionary_policy [VARIABLE_NAME...];` 

*Command*: `discretionary_policy (OPTIONS...) [VARIABLE_NAME...];`

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

- `discretionary_tol = NON-NEGATIVE DOUBLE`

Sets the tolerance level used to assess convergence of the solution
algorithm. Default: `1e-7`.

- `maxit = INTEGER`

Maximum number of iterations. Default: `3000`.

### Optimal Simple Rules (OSR)

*Command*: `osr [VARIABLE_NAME...];`

*Command*: `osr (OPTIONS...) [VARIABLE_NAME...];`

This command computes optimal simple policy rules for linear-quadratic
problems of the form:

$$\min_\gamma E(y'_tWy_t)$$

such that:

$$A_1 E_ty_{t+1}+A_2 y_t+ A_3 y_{t-1}+C e_t=0$$

where:

-   $E$ denotes the unconditional expectations operator;
-   $\gamma$ are parameters to be optimized. They must be elements of
    the matrices $A_1$, $A_2$, $A_3$, i.e. be specified as parameters
    in the `params` command and be entered in the `model` block;
-   $y$ are the endogenous variables, specified in the `var` command,
    whose (co)-variance enters the loss function;
-   $e$ are the exogenous stochastic shocks, specified in the
    `varexo`- ommand;
-   $W$ is the weighting matrix;

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

- `opt_algo = INTEGER`

Specifies the optimizer for minimizing the objective function. The same
solvers as for `mode_compute` (see
`mode_compute <mode_compute = INTEGER | FUNCTION_NAME>`{.interpreted-text
role="opt"}) are available, except for `5`, `6`, and `10`.

- `optim = (NAME, VALUE, ...)`

A list of NAME`[ and VALUE pairs. Can be used to set options for the
optimization routines. The set of available options depends on the
selected optimization routine (i.e. on the value of option
:opt:`opt_algo \<opt\_algo = INTEGER\>]{.title-ref}). See
`optim <optim = (NAME, VALUE, ...)>`{.interpreted-text role="opt"}.

- `maxit = INTEGER`

Determines the maximum number of iterations used in `opt_algo=4`. This
option is now deprecated and will be removed in a future release of
Dynare. Use `optim` instead to set optimizer-specific values. Default:
`1000`.

- `tolf = DOUBLE`

Convergence criterion for termination based on the function value used
in `opt_algo=4`. Iteration will cease when it proves impossible to
improve the function value by more than tolf. This option is now
deprecated and will be removed in a future release of Dynare. Use
`optim` instead to set optimizer-specific values. Default: `1e-7`.

- `silent_optimizer`

See `silent_optimizer`{.interpreted-text role="opt"}.

- `huge_number = DOUBLE`

Value for replacing the infinite bounds on parameters by finite numbers.
Used by some optimizers for numerical reasons (see
`huge_number <huge_number = DOUBLE>`{.interpreted-text role="opt"}).
Users need to make sure that the optimal parameters are not larger than
this value. Default: `1e7`.

The value of the objective is stored in the variable
`oo_.osr.objective_function` and the value of parameters at the optimum
is stored in `oo_.osr.optim_params`. See below for more details.

After running `osr` the parameters entering the simple rule will be set
to their optimal value so that subsequent runs of `stoch_simul` will be
conducted at these values.

*Command*: `osr_params PARAMETER_NAME...;`

This command declares parameters to be optimized by `osr`.

*Block*: `optim_weights ;`

This block specifies quadratic objectives for optimal policy problems.

More precisely, this block specifies the nonzero elements of the weight
matrix $W$ used in the quadratic form of the objective function in
`osr`.

An element of the diagonal of the weight matrix is given by a line of
the form:

```
    VARIABLE_NAME EXPRESSION;
```

An off-the-diagonal element of the weight matrix is given by a line of
the form:

```
    VARIABLE_NAME,  VARIABLE_NAME EXPRESSION;
```

*Example*

```
     var y inflation r;
     varexo y_ inf_;

     parameters delta sigma alpha kappa gammarr gammax0 gammac0 gamma_y_ gamma_inf_;

     delta =  0.44;
     kappa =  0.18;
     alpha =  0.48;
     sigma = -0.06;

     gammarr = 0;
     gammax0 = 0.2;
     gammac0 = 1.5;
     gamma_y_ = 8;
     gamma_inf_ = 3;

     model(linear);
     y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - inflation(+1)) + y_;
     inflation  =   alpha * inflation(-1) + (1-alpha) * inflation(+1) + kappa*y + inf_;
     r = gammax0*y(-1)+gammac0*inflation(-1)+gamma_y_*y_+gamma_inf_*inf_;
     end;

     shocks;
     var y_; stderr 0.63;
     var inf_; stderr 0.4;
     end;

     optim_weights;
     inflation 1;
     y 1;
     y, inflation 0.5;
     end;

     osr_params gammax0 gammac0 gamma_y_ gamma_inf_;
     osr y;
```

*Block*: `osr_params_bounds ;`

This block declares lower and upper bounds for parameters in the optimal
simple rule. If not specified the optimization is unconstrained.

Each line has the following syntax:

```
    PARAMETER_NAME, LOWER_BOUND, UPPER_BOUND;
```

Note that the use of this block requires the use of a constrained
optimizer, i.e. setting
`opt_algo <opt_algo = INTEGER>`{.interpreted-text role="opt"} to `1`,
`2`, `5` or `9`.

*Example*

```
     osr_params_bounds;
     gamma_inf_, 0, 2.5;
     end;

     osr(opt_algo=9) y;
```

*MATLAB/Octave Variables*: `oo.osr.objective_function`

After an execution of the `osr` command, this variable contains the
value of the objective under optimal policy.

*MATLAB/Octave Variables*: `oo.osr.optim_params`

After an execution of the `osr` command, this variable contains the
value of parameters at the optimum, stored in fields of the form
`oo_.osr.optim_params.PARAMETER_NAME`.

*MATLAB/Octave Variables*: `M.osr.param_names`

After an execution of the `osr` command, this cell contains the names of
the parameters.

*MATLAB/Octave Variables*: `M.osr.param_indices`

After an execution of the `osr` command, this vector contains the
indices of the OSR parameters in `M_.params`.

*MATLAB/Octave Variables*: `M.osr.param_bounds`

After an execution of the `osr` command, this two by number of OSR
parameters matrix contains the lower and upper bounds of the parameters
in the first and second column, respectively.

*MATLAB/Octave Variables*: `M.osr.variable_weights`

After an execution of the `osr` command, this sparse matrix contains the
weighting matrix associated with the variables in the objective
function.

*MATLAB/Octave Variables*: `M.osr.variable_indices`

After an execution of the `osr` command, this vector contains the
indices of the variables entering the objective function in
`M_.endo_names`.


## Displaying and saving results

Dynare has comments to plot the results of a simulation and to save the
results.

*Command*: `rplot VARIABLE_NAME...;`

Plots the simulated path of one or several variables, as stored in
`oo_.endo_simul` by either `perfect_foresight_solver`, `simul` (see
`det-simul`{.interpreted-text role="ref"}) or `stoch_simul` with option
`periods` (see `stoch-sol`{.interpreted-text role="ref"}). The variables
are plotted in levels.

*Command*: `dynatype (FILENAME) [VARIABLE_NAME...];`

This command prints the listed endogenous or exogenous variables in a
text file named FILENAME. If no VARIABLE\_NAME is listed, all endogenous
variables are printed.

*Command*: `dynasave (FILENAME) [VARIABLE_NAME...];`

This command saves the listed endogenous or exogenous variables in a
binary file named FILENAME. If no VARIABLE\_NAME is listed, all
endogenous variables are saved.

In MATLAB or Octave, variables saved with the `dynasave` command can be
retrieved by the command:

```
    load(FILENAME,'-mat')
```

## Macro processing language

It is possible to use "macro" commands in the `.mod` file for performing
tasks such as: including modular source files, replicating blocks of
equations through loops, conditionally executing some code, writing
indexed sums or products inside equations...

The Dynare macro-language provides a new set of *macro-commands* which
can be used in `.mod` files. It features:

-   File inclusion
-   Loops (`for` structure)
-   Conditional inclusion (`if/then/else` structures)
-   Expression substitution

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

-   `@#includepath`, paths to search for files that are to be
    included,
-   `@#include`, for file inclusion,
-   `@#define`, for defining a macro processor variable,
-   `@#if, @#ifdef, @#ifndef, @#elseif, @#else, @#endif` for
    conditional statements,
-   `@#for, @#endfor` for constructing loops.

The macro processor maintains its own list of variables (distinct from
model variables and MATLAB/Octave variables). These macro-variables are
assigned using the `@#define` directive and can be of the following
basic types: boolean, real, string, tuple, function, and array (of any
of the previous types).

### Macro expressions {#macro-exp}

Macro-expressions can be used in two places:

-   Inside macro directives, directly;
-   In the body of the `.mod` file, between an at-sign and curly
    braces (like `@{expr}`): the macro processor will substitute the
    expression with its value

It is possible to construct macro-expressions that can be assigned to
macro-variables or used within a macro-directive. The expressions are
constructed using literals of the basic types (boolean, real, string,
tuple, array), comprehensions, macro-variables, macro-functions, and
standard operators.

!!! note

    Elsewhere in the manual, MACRO_EXPRESSION designates an expression
    constructed as explained in this section.

**Boolean**

The following operators can be used on booleans:

-   Comparison operators: `==, !=`
-   Logical operators: `&&, ||, !`

**Real**

The following operators can be used on reals:

-   Arithmetic operators: `+, -, *, /, ^`
-   Comparison operators: `<, >, <=, >=, ==, !=`
-   Logical operators: `&&, ||, !`
-   Ranges with an increment of `1`: `REAL1:REAL2` (for example, `1:4`
    is equivalent to real array `[1, 2, 3, 4]`).

    4.6 Previously, putting brackets around the arguments to the colon
    operator (e.g. `[1:4]`) had no effect. Now, `[1:4]` will create an
    array containing an array (i.e. `[ [1, 2, 3, 4] ]`).

-   Ranges with user-defined increment: `REAL1:REAL2:REAL3` (for
    example, `6:-2.1:-1` is equivalent to real array
    `[6, 3.9, 1.8, -0.3]`).
-   Functions:
    `max, min, mod, exp, log, log10, sin, cos, tan, asin, acos, atan, sqrt, cbrt, sign, floor, ceil, trunc, erf, erfc, gamma, lgamma, round, normpdf, normcdf`.
    NB `ln` can be used instead of `log`

**String**

String literals have to be enclosed by **double** quotes (like
`"name"`).

The following operators can be used on strings:

-   Comparison operators: `<, >, <=, >=, ==, !=`
-   Concatenation of two strings: `+`
-   Extraction of substrings: if `s` is a string, then `s[3]` is a
    string containing only the third character of `s`, and `s[4:6]`
    contains the characters from 4th to 6th
-   Function: `length`

**Tuple**

Tuples are enclosed by parenthesis and elements separated by commas
(like `(a,b,c)` or `(1,2,3)`).

The following operators can be used on tuples:

-   Comparison operators: `==, !=`
-   Functions: `empty, length`

**Array**

Arrays are enclosed by brackets, and their elements are separated by
commas (like `[1,[2,3],4]` or `["US", "FR"]`).

The following operators can be used on arrays:

-   Comparison operators: `==, !=`
-   Dereferencing: if `v` is an array, then `v[2]` is its 2nd element
-   Concatenation of two arrays: `+`
-   Set union of two arrays: `|`
-   Set intersection of two arrays: `&`
-   Difference `-`: returns the first operand from which the elements
    of the second operand have been removed.
-   Cartesian product of two arrays: `*`
-   Cartesian product of one array N times: `^N`
-   Extraction of sub-arrays: e.g. `v[4:6]`
-   Testing membership of an array: `in` operator (for example: `"b"`
    in `["a", "b", "c"]` returns `1`)
-   Functions: `empty, sum, length`

**Comprehension**

Comprehension syntax is a shorthand way to make arrays from other
arrays. There are three different ways the comprehension syntax can be
employed: [filtering]{.title-ref}, [mapping]{.title-ref}, and [filtering
and mapping]{.title-ref}.

**Filtering**

Filtering allows one to choose those elements from an array for which a
certain condition hold.

*Example*

Create a new array, choosing the even numbers from the array `1:5`:

```
     [ i in 1:5 when mod(i,2) == 0 ]
```

would result in:

```
   [2, 4]
```

**Mapping**

Mapping allows you to apply a transformation to every element of an
array.

*Example*

Create a new array, squaring all elements of the array `1:5`:

```
   [ i^2 for i in 1:5 ]
```

would result in:

```
    [1, 4, 9, 16, 25]
```

**Filtering and Mapping**

Combining the two preceding ideas would allow one to apply a
transformation to every selected element of an array.

*Example*

Create a new array, squaring all even elements of the array `1:5`:

```
   [ i^2 for i in 1:5 when mod(i,2) == 0]
```

would result in:

```
    [4, 16]
```

*Further Examples* :

```
   [ (j, i+1) for (i,j) in (1:2)^2 ]
   [ (j, i+1) for (i,j) in (1:2)*(1:2) when i < j ]
```

would result in:

```
     [(1, 2), (2, 2), (1, 3), (2, 3)]
     [(2, 2)]
```

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

*Examples*

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

*Examples*

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

*Examples*

  ----------------------------------------
  **Code**                    **Output**
  --------------------------- ------------
  `(bool) 0 && true`          `false`

  `(real) "1" + 2`            `3`

  `(string) (3 + 4)`          `"7"`

  `(array) 5 + (array) 6`     `[5, 6]`
  ----------------------------------------

### Macro directives

*Macro Directive*: `@#includepath "PATH"`

*Macro Directive* `@#includepath MACRO_EXPRESSION`

This directive adds the path contained in PATH to the list of those to
search when looking for a `.mod` file specified by `@#include`. If
provided with a MACRO_EXPRESSION argument, the argument must evaluate
to a string. Note that these paths are added *after* any paths passed
using `-I <-I\<\<path\>\>>`{.interpreted-text role="opt"}.

*Example*

```
     @#includepath "/path/to/folder/containing/modfiles"
     @#includepath folders_containing_mod_files
```

*Macro Directive*: `@#include "FILENAME"` 

*Macro Directive*: `@#include MACRO_EXPRESSION`

This directive simply includes the content of another file in its place;
it is exactly equivalent to a copy/paste of the content of the included
file. If provided with a MACRO_EXPRESSION argument, the argument must
evaluate to a string. Note that it is possible to nest includes (i.e. to
include a file from an included file). The file will be searched for in
the current directory. If it is not found, the file will be searched for
in the folders provided by `-I <-I\<\<path\>\>>`{.interpreted-text
role="opt"} and `@#includepath`.

*Example*

```
     @#include "modelcomponent.mod"
     @#include location_of_modfile
```

*Macro Directive*: `@#define MACRO_VARIABLE` 

*Macro Directive*: `@#define MACRO_VARIABLE = MACRO_EXPRESSION`

*Macro Directive*: `@#define MACRO_FUNCTION = MACRO_EXPRESSION`

Defines a macro-variable or macro function.

*Example*

```
     @#define var                      // Equals true
     @#define x = 5                    // Real
     @#define y = "US"                 // String
     @#define v = [ 1, 2, 4 ]          // Real array
     @#define w = [ "US", "EA" ]       // String array
     @#define u = [ 1, ["EA"] ]        // Mixed array
     @#define z = 3 + v[2]             // Equals 5
     @#define t = ("US" in w)          // Equals true
     @#define f(x) = " " + x + y       // Function `f` with argument `x`
                                       // returns the string ' ' + x + 'US'
```

*Example*

```
     @#define x = 1
     @#define y = [ "B", "C" ]
     @#define i = 2
     @#define f(x) = x + " + " + y[i]
     @#define i = 1

     model;
       A = @{y[i] + f("D")};
     end;
```

The latter is strictly equivalent to:

```
     model;
       A = BD + B;
     end;
```

*Macro Directive*: `@#if MACRO_EXPRESSION`

*Macro Directive*: `@#ifdef MACRO_VARIABLE`

*Macro Directive*: `@#ifndef MACRO_VARIABLE` 

*Macro Directive*: `@#elseif MACRO_EXPRESSION` 

*Macro Directive*: `@#else @#endif`

Conditional inclusion of some part of the `.mod` file. The lines between
`@#if`, `@#ifdef`, or `@#ifndef` and the next `@#elseif`, `@#else` or
`@#endif` is executed only if the condition evaluates to `true`.
Following the `@#if` body, you can zero or more `@#elseif` branches. An
`@#elseif` condition is only evaluated if the preceding `@#if` or
`@#elseif` condition evaluated to `false`. The `@#else` branch is
optional and is only evaluated if all `@#if` and `@#elseif` statements
evaluate to false.

Note that when using `@#ifdef`, the condition will evaluate to `true` if
the MACRO_VARIABLE has been previously defined, regardless of its
value. Conversely, `@#ifndef` will evaluate to true if the
MACRO_VARIABLE has not yet been defined.

Note that when using `@#elseif` you can check whether or not a variable
has been defined by using the `defined` operator. Hence, to enter the
body of an `@#elseif` branch if the variable `X` has not been defined,
you would write: `@#elseif !defined(X)`.

Note that if a real appears as the result of the MACRO_EXPRESSION, it
will be interpreted as a boolean; a value of `0` is interpreted as
`false`, otherwise it is interpreted as `true`. Further note that
because of the imprecision of reals, extra care must be taken when
testing them in the MACRO_EXPRESSION. For example, `exp(log(5)) == 5`
will evaluate to `false`. Hence, when comparing real values, you should
generally use a zero tolerance around the value desired, e.g.
`exp(log(5)) > 5-1e-14 && exp(log(5)) < 5+1e-14`

*Example*

Choose between two alternative monetary policy rules using a
macro-variable:

```
     @#define linear_mon_pol = false // 0 would be treated the same
     ...
     model;
     @#if linear_mon_pol
       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
     @#else
       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
     @#endif
     ...
     end;
```

This would result in:

```
     ...
     model;
       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
     ...
     end;
```

*Example*

Choose between two alternative monetary policy rules using a
macro-variable. The only difference between this example and the
previous one is the use of `@#ifdef` instead of `@#if`. Even though
`linear_mon_pol` contains the value `false` because `@#ifdef` only
checks that the variable has been defined, the linear monetary policy
is output:

```
     @#define linear_mon_pol = false // 0 would be treated the same
     ...
     model;
     @#ifdef linear_mon_pol
       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
     @#else
       i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
     @#endif
     ...
     end;
```

This would result in:

```
     ...
     model;
       i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
     ...
     end;
```

*Macro Directive*: `@#for MACRO_VARIABLE in MACRO_EXPRESSION` 

*Macro Directive*: `@#for MACRO_VARIABLE in MACRO_EXPRESSION when MACRO_EXPRESSION`

*Macro Directive*: `@#for MACRO_TUPLE in MACRO_EXPRESSION` 

*Macro Directive*: `@#for MACRO_TUPLE in MACRO_EXPRESSION when MACRO\_EXPRESSION` 

*Macro Directive*: `@#endfor`

Loop construction for replicating portions of the `.mod` file. Note that
this construct can enclose variable/parameters declaration,
computational tasks, but not a model declaration.

*Example*

```
     model;
     @#for country in [ "home", "foreign" ]
       GDP_@{country} = A * K_@{country}^a * L_@{country}^(1-a);
     @#endfor
     end;
```

The latter is equivalent to:

```
     model;
       GDP_home = A * K_home^a * L_home^(1-a);
       GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
     end;
```

*Example*

```
     model;
     @#for (i, j) in ["GDP"] * ["home", "foreign"]
       @{i}_@{j} = A * K_@{j}^a * L_@{j}^(1-a);
     @#endfor
     end;
```

The latter is equivalent to:

```
     model;
       GDP_home = A * K_home^a * L_home^(1-a);
       GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
     end;
```

*Example*

```
     @#define countries = ["US", "FR", "JA"]
     @#define nth_co = "US"
     model;
     @#for co in countries when co != nth_co
       (1+i_@{co}) = (1+i_@{nth_co}) * E_@{co}(+1) / E_@{co};
     @#endfor
       E_@{nth_co} = 1;
     end;
```

The latter is equivalent to:

```
     model;
       (1+i_FR) = (1+i_US) * E_FR(+1) / E_FR;
       (1+i_JA) = (1+i_US) * E_JA(+1) / E_JA;
       E_US = 1;
     end;
```

*Macro Directive*: `@#echo MACRO_EXPRESSION`

Asks the preprocessor to display some message on standard output. The
argument must evaluate to a string.

*Macro Directive*: `@#error MACRO_EXPRESSION`

Asks the preprocessor to display some error message on standard output
and to abort. The argument must evaluate to a string.

*Macro Directive*: `@#echomacrovars`

*Macro Directive*: `@#echomacrovars MACRO_VARIABLE_LIST`

*Macro Directive*: `@#echomacrovars(save) MACRO_VARIABLE_LIST`

Asks the preprocessor to display the value of all macro variables up
until this point. If the `save` option is passed, then values of the
macro variables are saved to `options_.macrovars_line_<<line_numbers>>`.
If `NAME_LIST` is passed, only display/save variables and functions with
that name.

*Example*

```
     @#define A = 1
     @#define B = 2
     @#define C(x) = x*2
     @#echomacrovars A C D
```

The output of the command above is:

```
     Macro Variables:
       A = 1
     Macro Functions:
       C(x) = (x * 2)
```

### Typical usages

#### Modularization

The `@#include` directive can be used to split `.mod` files into several
modular components.

Example setup:

`modeldesc.mod`

Contains variable declarations, model equations, and shocks
declarations.

`simul.mod`

Includes `modeldesc.mod`, calibrates parameter,s and runs stochastic
simulations.

`estim.mod`

Includes `modeldesc.mod`, declares priors on parameters, and runs
Bayesian estimation.

Dynare can be called on `simul.mod` and `estim.mod` but it makes no
sense to run it on `modeldesc.mod`.

The main advantage is that you don't have to copy/paste the whole model
(at the beginning) or changes to the model (during development).

#### Indexed sums of products

The following example shows how to construct a moving average:

```
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
```

After macro processing, this is equivalent to:

```
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
```

#### Multi-country models

Here is a skeleton example for a multi-country model:

```
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
```

#### Endogeneizing parameters

When calibrating the model, it may be useful to consider a parameter as
an endogenous variable (and vice-versa).

For example, suppose production is defined by a CES function:

$$y = \left(\alpha^{1/\xi} \ell^{1-1/\xi}+(1-\alpha)^{1/\xi}k^{1-1/\xi}\right)^{\xi/(\xi-1)}$$

and the labor share in GDP is defined as:

$$\textrm{lab\_rat} = (w \ell)/(p y)$$

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

This file contains variable declarations and model equations. The code
for the declaration of $\alpha$ and `lab_rat` would look like:

```
     @#if steady
       var alpha;
       parameter lab_rat;
     @#else
       parameter alpha;
       var lab_rat;
     @#endif
```

`steady.mod`

This file computes the steady state. It begins with:

```
     @#define steady = 1
     @#include "modeqs.mod"
```

Then it initializes parameters (including `lab_rat`, excluding
$\alpha$), computes the steady state (using guess values for
endogenous, including $\alpha$), then saves values of parameters and
endogenous at steady state in a file, using the
`save_params_and_steady_state` command.

`simul.mod`

This file computes the simulation. It begins with:

```
     @#define steady = 0
     @#include "modeqs.mod"
```

Then it loads values of parameters and endogenous at steady state from
file, using the `load_params_and_steady_state` command, and computes
the simulations.

### MATLAB/Octave loops versus macro processor loops

Suppose you have a model with a parameter $\rho$ and you want to run
simulations for three values: $\rho = 0.8, 0.9,
1$. There are several ways of doing this:

*With a MATLAB/Octave loop*

```
     rhos = [ 0.8, 0.9, 1];
     for i = 1:length(rhos)
       rho = rhos(i);
       stoch_simul(order=1);
     end
```

Here the loop is not unrolled, MATLAB/Octave manages the iterations.
This is interesting when there are a lot of iterations.

*With a macro processor loop (case 1)*

```
     rhos = [ 0.8, 0.9, 1];
     @#for i in 1:3
       rho = rhos(@{i});
       stoch_simul(order=1);
     @#endfor
```

This is very similar to the previous example, except that the loop is
unrolled. The macro processor manages the loop index but not the data
array (`rhos`).

*With a macro processor loop (case 2)*

```
     @#for rho_val in [ 0.8, 0.9, 1]
       rho = @{rho_val};
       stoch_simul(order=1);
     @#endfor
```

The advantage of this method is that it uses a shorter syntax, since
the list of values is directly given in the loop construct. The
inconvenience is that you can not reuse the macro array in
MATLAB/Octave.

## Verbatim inclusion

Pass everything contained within the verbatim block to the
`<mod_file>.m` file.

*Block*: `verbatim ;`

By default, whenever Dynare encounters code that is not understood by
the parser, it is directly passed to the preprocessor output.

In order to force this behavior you can use the `verbatim` block. This
is useful when the code you want passed to the `<mod_file>.m` file
contains tokens recognized by the Dynare preprocessor.

*Example*

```
     verbatim;
     % Anything contained in this block will be passed
     % directly to the <modfile>.m file, including comments
     var = 1;
     end;
```

## Misc commands

*Command*: `set_dynare_seed (INTEGER) set_dynare_seed ('default')`

*Command*: `set_dynare_seed ('clock') set_dynare_seed ('reset')`

*Command*: `set_dynare_seed ('ALGORITHM', INTEGER)`

Sets the seed used for random number generation. It is possible to set a
given integer value, to use a default value, or to use the clock (by
using the latter, one will therefore get different results across
different Dynare runs). The `reset` option serves to reset the seed to
the value set by the last `set_dynare_seed` command. On MATLAB 7.8 or
above, it is also possible to choose a specific algorithm for random
number generation; accepted values are `mcg16807`, `mlfg6331_64`,
`mrg32k3a`, `mt19937ar` (the default), `shr3cong` and `swb2712`.

*Command*: `save_params_and_steady_state (FILENAME);    

For all parameters, endogenous and exogenous variables, stores their
value in a text file, using a simple name/value associative table.

-   for parameters, the value is taken from the last parameter
    initialization.
-   for exogenous, the value is taken from the last `initval` block.
-   for endogenous, the value is taken from the last steady state
    computation (or, if no steady state has been computed, from the
    last `initval` block).

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

- `load_params_and_steady_state (FILENAME);`

For all parameters, endogenous and exogenous variables, loads their
value from a file created with `save_params_and_steady_state`.

-   for parameters, their value will be initialized as if they had
    been calibrated in the `.mod` file.
-   for endogenous and exogenous variables, their value will be
    initialized as they would have been from an `initval` block .

This function is used in conjunction with
`save_params_and_steady_state`; see the documentation of that function
for more information.

*Command*: `compilation_setup (OPTIONS);`

When the `use_dll`{.interpreted-text role="opt"} option is present,
Dynare uses the GCC compiler that was distributed with it to compile the
static and dynamic C files produced by the preprocessor. You can use
this option to change the compiler, flags, and libraries used.

*Options*


- `compiler = FILENAME`

The path to the compiler.

- `substitute_flags = QUOTED_STRING`

The flags to use instead of the default flags.

- `add_flags = QUOTED_STRING`

The flags to use in addition to the default flags. If
`substitute_flags` is passed, these flags are added to the flags
specified there.

- `substitute_libs = QUOTED_STRING`
 
The libraries to link against instead of the default libraries.

- `add_libs = QUOTED_STRING`

The libraries to link against in addition to the default libraries. If
`substitute_libs` is passed, these libraries are added to the
libraries specified there.

*MATLAB/Ocatve Command*: `dynare_version ;`

Output the version of Dynare that is currently being used (i.e. the one
that is highest on the MATLAB/Octave path).

*MATLAB/Ocatve Command*: `write_latex_definitions ;`

Writes the names, LaTeX names and long names of model variables to
tables in a file named `<<M_.fname>>_latex_definitions.tex`. Requires
the following LaTeX packages: `longtable`.

*MATLAB/Ocatve Command*: `write_latex_parameter_table ;`

Writes the LaTeX names, parameter names, and long names of model
parameters to a table in a file named
`<<M_.fname>>_latex_parameters.tex.` The command writes the values of
the parameters currently stored. Thus, if parameters are set or changed
in the steady state computation, the command should be called after a
steady-command to make sure the parameters were correctly updated. The
long names can be used to add parameter descriptions. Requires the
following LaTeX packages: `longtable, booktabs`.

*MATLAB/Ocatve Command*: `write_latex_prior_table ;`

Writes descriptive statistics about the prior distribution to a LaTeX
table in a file named `<<M_.fname>>_latex_priors_table.tex`. The command
writes the prior definitions currently stored. Thus, this command must
be invoked after the `estimated_params` block. If priors are defined
over the measurement errors, the command must also be preceeded by the
declaration of the observed variables (with `varobs`). The command
displays a warning if no prior densities are defined (ML estimation) or
if the declaration of the observed variables is missing. Requires the
following LaTeX packages: `longtable, booktabs`.

*MATLAB/Ocatve Command*: `collect_latex_files ;`

Writes a LaTeX file named `<<M_.fname>>_TeX_binder.tex` that collects
all TeX output generated by Dynare into a file. This file can be
compiled using `pdflatex` and automatically tries to load all required
packages. Requires the following LaTeX packages: `breqn`, `psfrag`,
`graphicx`, `epstopdf`, `longtable`, `booktabs`, `caption`, `float,`
`amsmath`, `amsfonts`, and `morefloats`.

**Footnotes**

[^1]: A `.mod` file must have lines that end with a line feed character,
    which is not commonly visible in text editors. Files created on
    Windows and Unix-based systems have always conformed to this
    requirement, as have files created on OS X and macOS. Files created
    on old, pre-OS X Macs used carriage returns as end of line
    characters. If you get a Dynare parsing error of the form
    `ERROR: <<mod file>>: line 1, cols 341-347: syntax error,...` and
    there's more than one line in your `.mod` file, know that it uses
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