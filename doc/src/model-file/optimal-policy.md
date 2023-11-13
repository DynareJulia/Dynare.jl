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
the `shocks` blocks (see `shocks-exo`).
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
<dates-members>`). This function return a
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
members <dates-members>`). The last
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
members <dates-members>`). The last
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
members <dseries-members>`) can be
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
