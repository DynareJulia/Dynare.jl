# Model File
## Syntax elements
### Conventions

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
-   MODEL\_EXPRESSION (sometimes MODEL\_EXP) indicates a mathematical
    expression valid in the model description (see
    `expr` and
    `model-decl`);
-   MACRO_EXPRESSION designates an expression of the macro processor
    (see `macro-exp`);
-   VARIABLE\_NAME (sometimes VAR\_NAME) indicates a variable name
    starting with an alphabetical character and can't contain:
    '()+-*/\^=!;:@#.' or accentuated characters;
-   PARAMETER\_NAME (sometimes PARAM\_NAME) indicates a parameter name
    starting with an alphabetical character and can't contain:
    '()+-*/\^=!;:@#.' or accentuated characters;
-   LATEX\_NAME (sometimes TEX\_NAME) indicates a valid LaTeX expression
    in math mode (not including the dollar signs);
-   FUNCTION_NAME indicates a valid Julia function name;
-   FILENAME indicates a filename valid in the underlying operating
    system; it is necessary to put it between quotes when specifying the
    extension or if the filename contains a non-alphanumeric character;
-   QUOTED_STRING indicates an arbitrary string enclosed between
    (single) quotes. Note that Dynare commands call for single quotes
    around a string while in Julia strings are enclosed between double quotes.


### Expressions

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
