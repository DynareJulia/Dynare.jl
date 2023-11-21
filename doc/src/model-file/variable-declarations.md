
While Dynare allows the user to choose their own variable names, there
are some restrictions to be kept in mind. First, variables and
parameters must not have the same name as Dynare commands or built-in
functions. In this respect, Dynare is not case-sensitive. For example,
do not use `Ln` or `Sigma_e` to name your variable. Not conforming to
this rule might yield hard-to-debug error messages or crashes. 

## var

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

## varexo
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
	```math
    a_t=ρa_{t−1}+ε_t
    ```
then logged TFP is an endogenous variable to be declared with var, while the innovation `ε_t` is to be declared with varexo.

## varexo_det

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

## parameters
*Command*:

```
parameters PARAM_NAME [$TEX_NAME$] [(long_name=QUOTED_STRING|NAME=QUOTED_STRING)...];
```
    This command declares parameters used in the model, in variable
    initialization or in shocks declaration. Optionally it is possible to give a
    LaTeX name to the parameter.

    The parameters must subsequently be assigned values.

    The parameters in the list can be separated by spaces or by commas.
    ``parameters`` commands can appear several times in the file and Dynare
    will concatenate them.

*Options*

- `long_name = QUOTED_STRING`

- `NAME = QUOTED_STRING`

*Example*

```
parameters alpha, bet;
```

*Parameter initialization*

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

## Advanced commands

### change_type

*Command*

```
change_type (var|varexo|varexo_det|parameters) VAR_NAME | PARAM_NAME...;
```

Changes the types of the specified variables/parameters to another
type: endogenous, exogenous, exogenous deterministic or
parameter. It is important to understand that this command has a
global effect on the ``.mod`` file: the type change is effective
after, but also before, the ``change_type`` command. This command
is typically used when flipping some variables for steady state
calibration: typically a separate model file is used for
calibration, which includes the list of variable declarations with
the macro processor, and flips some variable.

*Example*

```
var y, w;
parameters alpha, beta;
...
change_type(var) alpha, beta;
change_type(parameters) y, w;
```

Here, in the whole model file, ``alpha`` and ``beta`` will be
endogenous and ``y`` and ``w`` will be parameters.

### var_remove

*Command*:

```
var_remove VAR_NAME | PARAM_NAME...;
```

Removes the listed variables (or parameters) from the model. Removing a
variable that has already been used in a model equation or elsewhere
will lead to an error.

### predetermined_variables

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

### trend_var

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

### make\_local\_variables

*command*:

```
model_local_variable VARIABLE_NAME [LATEX_NAME]... ;
```

This optional command declares a model local variable. See
`conv` for the syntax of `VARIABLE_NAME`.
As you can create model local variables on the fly in the model block, the interest of this
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
