To study the effects of a temporary shock
after which the system goes back to the original equilibrium (if the
model is stable...) one uses a temporary shock. A temporary shock is a temporary change of value
of one or several exogenous variables in the model. Temporary shocks are
specified with the command `shocks`.

In a deterministic context, when one wants to study the transition of
one equilibrium position to another, it is equivalent to analyze the
consequences of a permanent shock. In Dynare this is done with `initval`,
`endval` and `steady`.

In a stochastic framework, the exogenous variables take random values in
each period. In Dynare, these random values follow a normal distribution
with zero mean, but it belongs to the user to specify the variability of
these shocks. The non-zero elements of the matrix of variance-covariance
of the shocks can be entered with the `shocks` command.

### Dynare commands

#### shocks
- *block*: `shocks ;` 
- *block*: `shocks(overwrite);`

##### Options

- `overwrite`: By default, if there are several `shocks` blocks
in the same `.mod` file, then they are cumulative: all the shocks
declared in all the blocks are considered; however, if a `shocks` block is declared with the `overwrite` option, then it
replaces all the previous `shocks` blocks.

#### In a deterministic context

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

##### Example 1

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

##### Example 2

```
xx = [1.2; 1.3; 1];

shocks;
var e;
periods 1:3;
values (xx);
end;
```

#### In a stochastic context

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

##### Example

```
shocks;
var e = 0.000081;
var u; stderr 0.009;
corr e, u = 0.8;
var v, w = 2;
end;
```

##### Remark

If the variance of an exogenous variable is set to zero, this variable
will appear in the report on policy and transition functions, but isn't
used in the computation of moments and of Impulse Response Functions.
Setting a variance to zero is an easy way of removing an exogenous
shock.


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
### Julia function

#### scenario!()
The Julia function `scenario!()` lets you
- declare shocks on exogenous variables as the `shocks` block
- set the future value of endogenous variables (for conditional
  forecasts)
- add the date at which the above information is made available to the
  agents in the model
  
```@docs
scenario!
```

##### Examples
```
scenario!(name = :e, value = 0.1, period = 2)
```
Exogenous variable `e`, takes value 0.1 in period 2.

```
scenario!(name = :y, value = 0.2, period=2, exogenous = :u)
```
Endogenous variable `y` is set to 0.2 in period 2 and exogenous
variable `u` is treated as endogenous in the same period. Agents in
the model know at the beginning of period 1 that this will happen.

```
scenario!(infoperiod = 2, name = :y, value = 0.2, period = 2,
                exogenous = :u)
```
Endogenous variable `y` is set to 0.2 in period 2 and exogenous
variable `u` is treated as endogenous in the same period. Agents in
the model only learn at the beginning of period 2 that this will happen.

