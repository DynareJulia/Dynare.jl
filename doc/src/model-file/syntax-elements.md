## Conventions

A model file contains a list of commands and of blocks. Each command and
each element of a block is terminated by a semicolon (;). Blocks are
terminated by `end;`.

If Dynare encounters an unknown expression at the beginning of a line or
after a semicolon, it will parse the rest of that line as native Julia
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

Note that these comment marks should not be used in native 
Julia code regions where the [#] should be preferred instead to
introduce a comment. In a `verbatim` block, see
`verbatim`, this would result in a crash
since `//` is not a valid Julia statement).

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
    model description (see `expr`);
-   MODEL_EXPRESSION (sometimes MODEL_EXP) indicates a mathematical
    expression valid in the model description (see
    `expr` and
    `model-decl`);
-   MACRO_EXPRESSION designates an expression of the macro processor
    (see `macro-exp`);
-   VARIABLE_NAME (sometimes VAR_NAME) indicates a variable name
    starting with an alphabetical character and can't contain:
    '()+-*/\^=!;:@#.' or accentuated characters;
-   PARAMETER_NAME (sometimes PARAM_NAME) indicates a parameter name
    starting with an alphabetical character and can't contain:
    '()+-*/\^=!;:@#.' or accentuated characters;
-   LATEX_NAME (sometimes TEX_NAME) indicates a valid LaTeX expression
    in math mode (not including the dollar signs);
-   FUNCTION_NAME indicates a valid Julia function name;
-   FILENAME indicates a filename valid in the underlying operating
    system; it is necessary to put it between quotes when specifying the
    extension or if the filename contains a non-alphanumeric character;
-   QUOTED_STRING indicates an arbitrary string enclosed between
    (single) quotes. Note that Dynare commands call for single quotes
    around a string while in Julia strings are enclosed between double quotes.

## Variable declarations

While Dynare allows the user to choose their own variable names, there
are some restrictions to be kept in mind. First, variables and
parameters must not have the same name as Dynare commands or built-in
functions. In this respect, Dynare is not case-sensitive. For example,
do not use `Ln` or `Sigma_e` to name your variable. Not conforming to
this rule might yield hard-to-debug error messages or crashes. 

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
Optionally it is possible to give a LaTeX name to the
variable or, if it is nonstationary, provide information regarding its
deflator. The variables in the list can be separated by spaces or by
commas. `var` commands can appear several times in the file and Dynare
will concatenate them. 

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

    

*Example*

```
var c gnp cva 
          cca (long_name=`Consumption CA');
var(deflator=A) i b;
var c $C$ (long_name=`Consumption');
```

 *Command*:

``` 
 varexo VAR_NAME [$TEX_NAME$] [(long_name=QUOTED_STRING|NAME=QUOTED_STRING)...];
```

    This optional command declares the exogenous variables in the model. Optionally it is possible to give a LaTeX name to the variable. Exogenous variables are required if the user wants to be able to apply shocks to her model. The variables in the list can be separated by spaces or by commas. varexo commands can appear several times in the file and Dynare will concatenate them.

*Options*

- `long_name = QUOTED_STRING`

*Example*

        varexo m gov;

*Remarks*

    An exogenous variable is an innovation, in the sense that this
    variable cannot be predicted from the knowledge of the current
    state of the economy. For instance, if logged TFP is a first order
    autoregressive process:
	```
    a_t=ρa_{t−1}+ε_t
    ```
then logged TFP is an endogenous variable to be declared with var, while the innovation `ε_t` is to be declared with varexo.

*Command*:

```
varexo_det VAR_NAME [$TEX_NAME$][(long_name=QUOTED_STRING|NAME=QUOTED_STRING)...];
```

This optional command declares exogenous deterministic variables in a
stochastic model. Optionally it is possible to give a LaTeX name
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

- `NAME = QUOTED_STRING`


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
`conv` for the syntax of `MODEL_EXPR` and
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
`conv` for the syntax of `VARIABLE_NAME`.
As you can create model local variables on the fly in the model block
(see `model-decl`), the interest of this
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

Unlike Julia expressions, Dynare expressions are necessarily
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
parameter initialization (see `param-init`) or `homotopy_setup` when doing a simulation, or are the
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
`param-init`); for an endogenous or
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
at the previous period. See `aux-variables` for an explanation of how this operator is handled
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

```math
\begin{aligned}
  \textrm{sign}(x) =
         \begin{cases}
         -1 &\quad\text{if }x<0\\
         0 &\quad\text{if }x=0\\
         1 &\quad\text{if }x>0
         \end{cases}
  \end{aligned}
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

Internally, the parameter values are stored in `context.work.params`

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

[^2]: Note that arbitrary Julia expressions can be put in a
    `.mod` file, but those expressions have to be on separate lines,
    generally at the end of the file for post-processing purposes. They
    are not interpreted by Dynare, and are simply passed on unmodified
    to Julia. Those constructions are not addresses in this
    section.
