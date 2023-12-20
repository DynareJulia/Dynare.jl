From a statistical point of view, DSGE models are unobserved
components models: only a few variables are observed. Filtering or
smoothing provide estimate of the unobserved variables given the
observations.

The model is put in state space form
```math
\begin{align*}
y^o_t &= M s_t + N\epsilon_t\\
s_t &= Ts_{t-1} + R\eta_t
\end{align*}
```
where $$y^o_t$$ represents observed variable at period `t`. The coefficient matrices of the transition equation, `T` and `R` are provided by the solution of the linear(-isze) rational expectation model. $$\epsilon_t$$ are possible measurement errors and $$\eta_t$$ the structural shocks. Most often matrix `M` is a selection matrix.

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

#### varobs

Observed variables are declared with the `varobs` command

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

##### Example

```
     varobs C y rr;
```

#### observation\_trends

It is possible to declare a deterministic linear trend that is removed for the computations and added back in the results

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

##### Example

```
     observation_trends;
     Y (eta);
     P (mu/eta);
     end;
```

#### calib\_smoother

This command triggers the computation of the filter and smoother for calibrated models

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

- `datafile = FILENAME`: file containing the observation in CSV format.

- `filtered_vars`: triggers the computation of filtered variables.

- `first_obs = INTEGER`: first observation

- `diffuse_filter = INTEGER`: use a diffuse filter for nonstationary models.

### Julia functions

```@docs
calibsmoother!
```

