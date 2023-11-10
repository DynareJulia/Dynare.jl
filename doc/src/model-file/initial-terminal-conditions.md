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
