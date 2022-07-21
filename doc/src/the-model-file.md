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