
Provided that you have observations on some endogenous variables, it is
possible to use Dynare to estimate some or all parameters. Bayesian techniques (as in
*Fernández-Villaverde and Rubio-Ramírez (2004)*, *Rabanal and
Rubio-Ramirez (2003)*, *Schorfheide (2000)* or *Smets and Wouters
(2003)*) are available. Using Bayesian methods, it is possible to
estimate DSGE models.

Note that in order to avoid stochastic singularity, you must have at
least as many shocks or measurement errors in your model as you have
observed variables.

Before using the estimation commands described below, you need to define some elements
of the state space representation of the model. At the minimum,
you need to declare the observed variables with `var_obs` and, possibly, deterministic trends
with `observation_trends` (see the previous section: `State space, filtering and smoothing`)

### Dynare commands

#### estimated_params

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
`prior_3rd_parameter <PRIOR_3RD_PARAMETER>`. In this case the prior density is effectively modified
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
`PRIOR_3RD_PARAMETER` and
`PRIOR_4TH_PARAMETER` or via mean and
standard deviation using `PRIOR_MEAN`,
`PRIOR_STANDARD_ERROR`. The other two will
automatically be filled out. Note that providing both sets of
hyperparameters will yield an error message.

As one uses options more towards the end of the list, all previous
options must be filled: for example, if you want to specify
SCALE_PARAMETER, you must specify `PRIOR_3RD_PARAMETER` and
`PRIOR_4TH_PARAMETER`. Use empty values, if these parameters don't
apply.

##### Example

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

##### Example

```
     corr eps_1, eps_2, 0.5,  -0.5,  1, beta_pdf, 0, 0.3, -1, 1;
```

Sets the same generalized beta distribution as before, but now
truncates this distribution to `[-0.5,1]` through the use of
LOWER_BOUND and UPPER_BOUND.

#### Parameter transformation

Sometimes, it is desirable to estimate a transformation of a parameter
appearing in the model, rather than the parameter itself. It is of
course possible to replace the original parameter by a function of the
estimated parameter everywhere is the model, but it is often
unpractical.

In such a case, it is possible to declare the parameter to be estimated
in the parameters statement and to define the transformation, using a
pound sign (\#) expression.

##### Example

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
`estimated_params_remove` for that use
case). However, if an `estimated_params` block has the `overwrite`
option, its contents becomes the new list of estimated parameters,
cancelling previous blocks; this can be useful when doing several
estimations in a single `.mod` file.

#### estimated\_params\_init

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

#### estimation

*Command*: `estimation [VARIABLE_NAME...];`

*Command*: `estimation (OPTIONS...)[VARIABLE_NAME...];`

This command runs Bayesian  estimation.


##### Options

- `datafile = FILENAME`

The datafile must be in CSV format

```
    estimation(datafile='../fsdat_simul.csv',...);
```

- `nobs = INTEGER`

The number of observations following `first_obs` to be used.
Default: all observations in the file after `first_obs`.

- `first_obs = INTEGER`

The number of the first observation to be used. In case of estimating a
DSGE-VAR, `first_obs` needs to be larger than the number of lags.
Default: `1`.

- `plot_priors = INTEGER`: Control the plotting of priors, `0`, no prior plot, `1`, pPrior density for each estimated parameter is plotted. It is
    important to check that the actual shape of prior densities matches
    what you have in mind. Ill-chosen values for the prior standard
    density can result in absurd prior densities (default valueL `1`).

- `mh_replic = INTEGER`

Number of replications for each chain of the Metropolis-Hastings
algorithm. The number of draws should be sufficient to achieve
convergence of the MCMC and to meaningfully compute posterior objects.
Default: `20000`.

- `mh_nblocks = INTEGER`

Number of parallel chains for Metropolis-Hastings algorithm. Default:
`2`.

- `mh_jscale = DOUBLE`

The scale parameter of the jumping distribution's covariance matrix. The default value is rarely
satisfactory. This option must be tuned to obtain, ideally, an
acceptance ratio of 25%-33%. Basically, the idea is to increase the
variance of the jumping distribution if the acceptance ratio is too
high, and decrease the same variance if the acceptance ratio is too low.
In some situations it may help to consider parameter-specific values for
this scale parameter. This can be done in the
`estimated_params` block.
 Default: `0.2`.


### Julia functions

```@docs
covariance
```

```@docs
mode_compute!
```

```@docs
plot_priors 
```

```@docs
prior!
```

```@docs
rwmh_compute!
```
