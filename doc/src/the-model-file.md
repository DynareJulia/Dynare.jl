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