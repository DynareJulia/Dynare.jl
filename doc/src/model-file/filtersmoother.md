From a statistical point of view, DSGE models are unobserved
components models: only a few variables are observed. Filtering or
smoothing provide estimate of the unobserved variables given the
observations.

Filtering provides estimates conditional only on past observations:
```math
\mathbb{E}(y^{no}_t|Y^o_{t-1})
```
where $$y^{no}_t$$ are unobserved variables at period `t` and
$$Y^o_{t-1}$$ represent the set observations until period `t-1` included.

Smoothing provides estimates of unobserved variables conditional on
the entire sample of observations:
```math
\mathbb{E}(y^{no}_t|Y^o_T)
```
where $$Y^o_T$$ represents the all observations in the sample.

### Dynare command

*Command*: `varobs VARIABLE_NAME...;`

This command lists the name of observed endogenous variables for the
estimation procedure. These variables must be available in the data file
(see `estimation_cmd <estim-comm>`).

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

See `datafile <dataf>`.

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

See `loglinear <logl>`.

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
