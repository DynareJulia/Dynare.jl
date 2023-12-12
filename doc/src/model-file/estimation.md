
Provided that you have observations on some endogenous variables, it is
possible to use Dynare to estimate some or all parameters. Bayesian techniques (as in
*Fernández-Villaverde and Rubio-Ramírez (2004)*, *Rabanal and
Rubio-Ramirez (2003)*, *Schorfheide (2000)* or *Smets and Wouters
(2003)*) are available. Using Bayesian methods, it is possible to
estimate DSGE models.

Note that in order to avoid stochastic singularity, you must have at
least as many shocks or measurement errors in your model as you have
observed variables.

### Dynare commands

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

Declares endogenous variables `C`, `y` and `rr` as observed variables.

*Example* (with a macro processor loop)

```
   varobs
   @#for co in countries
   GDP_@{co}
   @#endfor
;
```

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

*Example*

```
     observation_trends;
     Y (eta);
     P (mu/eta);
     end;
```

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
`prior_3rd_parameter <PRIOR_3RD_PARAMETER>`{.interpreted-text
role="opt"}. In this case the prior density is effectively modified
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
`PRIOR_3RD_PARAMETER`{.interpreted-text role="opt"} and
`PRIOR_4TH_PARAMETER`{.interpreted-text role="opt"} or via mean and
standard deviation using `PRIOR_MEAN`{.interpreted-text role="opt"},
`PRIOR_STANDARD_ERROR`{.interpreted-text role="opt"}. The other two will
automatically be filled out. Note that providing both sets of
hyperparameters will yield an error message.

As one uses options more towards the end of the list, all previous
options must be filled: for example, if you want to specify
SCALE_PARAMETER, you must specify `PRIOR_3RD_PARAMETER` and
`PRIOR_4TH_PARAMETER`. Use empty values, if these parameters don't
apply.

*Example*

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

*Example*

```
     corr eps_1, eps_2, 0.5,  -0.5,  1, beta_pdf, 0, 0.3, -1, 1;
```

Sets the same generalized beta distribution as before, but now
truncates this distribution to `[-0.5,1]` through the use of
LOWER_BOUND and UPPER_BOUND.

*Parameter transformation*

Sometimes, it is desirable to estimate a transformation of a parameter
appearing in the model, rather than the parameter itself. It is of
course possible to replace the original parameter by a function of the
estimated parameter everywhere is the model, but it is often
unpractical.

In such a case, it is possible to declare the parameter to be estimated
in the parameters statement and to define the transformation, using a
pound sign (\#) expression (see `model-decl`{.interpreted-text
role="ref"}).

*Example*

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
`estimated_params_remove`{.interpreted-text role="bck"} for that use
case). However, if an `estimated_params` block has the `overwrite`
option, its contents becomes the new list of estimated parameters,
cancelling previous blocks; this can be useful when doing several
estimations in a single `.mod` file.

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

*Command*: `estimation [VARIABLE_NAME...];`

*Command*: `estimation (OPTIONS...)[VARIABLE_NAME...];`

This command runs Bayesian or maximum likelihood estimation.

The following information will be displayed by the command:

-   Results from posterior optimization (also for maximum likelihood)
-   Marginal log data density
-   Posterior mean and highest posterior density interval (shortest
    credible set) from posterior simulation
-   Convergence diagnostic table when only one MCM chain is used or
    Metropolis-Hastings convergence graphs documented in
    *Pfeifer (2014)* in case of multiple MCM chains
-   Table with numerical inefficiency factors of the MCMC
-   Graphs with prior, posterior, and mode
-   Graphs of smoothed shocks, smoothed observation errors, smoothed and
    historical variables

Note that the posterior moments, smoothed variables, k-step ahead
filtered variables and forecasts (when requested) will only be computed
on the variables listed after the `estimation` command. Alternatively,
one can choose to compute these quantities on all endogenous or on all
observed variables (see `consider_all_endogenous`,
`consider_all_endogenous_and_auxiliary`, and `consider_only_observed`
options below). If no variable is listed after the estimation command,
then Dynare will interactively ask which variable set to use.

Also, during the MCMC (Bayesian estimation with `mh_replic` $>0$) a
(graphical or text) waiting bar is displayed showing the progress of the
Monte-Carlo and the current value of the acceptance ratio. Note that if
the `load_mh_file` option is used (see below) the reported acceptance
ratio does not take into account the draws from the previous MCMC. In
the literature there is a general agreement for saying that the
acceptance ratio should be close to one third or one quarter. If this
not the case, you can stop the MCMC (`Ctrl-C`) and change the value of
option `mh_jscale` (see below).

Note that by default Dynare generates random numbers using the algorithm
`mt199937ar` (i.e. Mersenne Twister method) with a seed set equal to
`0`. Consequently the MCMCs in Dynare are deterministic: one will get
exactly the same results across different Dynare runs (*ceteris
paribus*). For instance, the posterior moments or posterior densities
will be exactly the same. This behaviour allows to easily identify the
consequences of a change on the model, the priors or the estimation
options. But one may also want to check that across multiple runs, with
different sequences of proposals, the returned results are almost
identical. This should be true if the number of iterations (i.e. the
value of `mh_replic`) is important enough to ensure the convergence of
the MCMC to its ergodic distribution. In this case the default behaviour
of the random number generators in not wanted, and the user should set
the seed according to the system clock before the estimation command
using the following command:

```
    set_dynare_seed('clock');
```

so that the sequence of proposals will be different across different
runs.

Finally, Dynare does not always properly distinguish between maximum
likelihood and Bayesian estimation in its field names. While there is an
important conceptual distinction between frequentist confidence
intervals and Bayesian highest posterior density intervals (HPDI) as
well as between posterior density and likelilhood, Dynare sometimes uses
the Bayesian terms as a stand-in in its display of maximum likelihood
results. An example is the storage of the output of the
`forecast`-option of `estimation` with ML, which will use
`HPDinf/HPDsup` to denote the confidence interval.

*Algorithms*

The Monte Carlo Markov Chain (MCMC) diagnostics are generated by the
estimation command if
`mh_replic <mh_replic = INTEGER>`{.interpreted-text role="opt"} is
larger than 2000 and if option `nodiagnostic`{.interpreted-text
role="opt"} is not used. If
`mh_nblocks <mh_nblocks = INTEGER>`{.interpreted-text role="opt"} is
equal to one, the convergence diagnostics of *Geweke (1992,1999)* is
computed. It uses a chi-square test to compare the means of the first
and last draws specified by `geweke_interval
<geweke_interval = [DOUBLE DOUBLE]>`{.interpreted-text role="opt"} after
discarding the burn-in of `mh_drop <mh_drop = DOUBLE>`{.interpreted-text
role="opt"}. The test is computed using variance estimates under the
assumption of no serial correlation as well as using tapering windows
specified in `taper_steps
<taper_steps = [INTEGER1 INTEGER2 ...]>`{.interpreted-text role="opt"}.
If `mh_nblocks
<mh_nblocks = INTEGER>`{.interpreted-text role="opt"} is larger than 1,
the convergence diagnostics of *Brooks and Gelman (1998)* are used
instead. As described in section 3 of *Brooks and Gelman (1998)* the
univariate convergence diagnostics are based on comparing pooled and
within MCMC moments (Dynare displays the second and third order moments,
and the length of the Highest Probability Density interval covering 80%
of the posterior distribution). Due to computational reasons, the
multivariate convergence diagnostic does not follow *Brooks and Gelman
(1998)* strictly, but rather applies their idea for univariate
convergence diagnostics to the range of the posterior likelihood
function instead of the individual parameters. The posterior kernel is
used to aggregate the parameters into a scalar statistic whose
convergence is then checked using the *Brooks and Gelman (1998)*
univariate convergence diagnostic.

The inefficiency factors are computed as in *Giordano et al.(2011)*
based on Parzen windows as in e.g. *Andrews (1991)*.

*Options*

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

- `plot_priors = INTEGER`

Control the plotting of priors.

 `0`

    No prior plot.

 `1`

    Prior density for each estimated parameter is plotted. It is
    important to check that the actual shape of prior densities matches
    what you have in mind. Ill-chosen values for the prior standard
    density can result in absurd prior densities.

Default value is `1`.


- `mh_replic = INTEGER`

Number of replications for each chain of the Metropolis-Hastings
algorithm. The number of draws should be sufficient to achieve
convergence of the MCMC and to meaningfully compute posterior objects.
Default: `20000`.

- `mh_nblocks = INTEGER`

Number of parallel chains for Metropolis-Hastings algorithm. Default:
`2`.

- `mh_jscale = DOUBLE`

The scale parameter of the jumping distribution's covariance matrix
(Metropolis-Hastings or TaRB-algorithm). The default value is rarely
satisfactory. This option must be tuned to obtain, ideally, an
acceptance ratio of 25%-33%. Basically, the idea is to increase the
variance of the jumping distribution if the acceptance ratio is too
high, and decrease the same variance if the acceptance ratio is too low.
In some situations it may help to consider parameter-specific values for
this scale parameter. This can be done in the
`estimated_params`{.interpreted-text role="bck"} block.
 Default: `0.2`.


- `diffuse_filter`

Uses the diffuse Kalman filter (as described in *Durbin and Koopman
(2012)* and *Koopman and Durbin (2003)* for the multivariate and
*Koopman and Durbin (2000)* for the univariate filter) to estimate
models with non-stationary observed variables. This option will also
reset the `qz_criterium` to count unit root variables towards the stable
variables. Trying to estimate a model with unit roots will otherwise
result in a Blanchard-Kahn error.

When `diffuse_filter` is used the `lik_init` option of `estimation` has
no effect.

When there are nonstationary exogenous variables in a model, there is no
unique deterministic steady state. For instance, if productivity is a
pure random walk:

$$a_t = a_{t-1} + e_t$$

any value of $\bar a$ of $a$ is a deterministic steady state for
productivity. Consequently, the model admits an infinity of steady
states. In this situation, the user must help Dynare in selecting one
steady state, except if zero is a trivial model's steady state, which
happens when the `linear` option is used in the model declaration. The
user can either provide the steady state to Dynare using a
`steady_state_model` block (or writing a steady state file) if a closed
form solution is available, see `steady_state_model`{.interpreted-text
role="bck"}, or specify some constraints on the steady state, see
`equation_tag_for_conditional_steady_state <eq-tag-ss>`{.interpreted-text
role="ref"}, so that Dynare computes the steady state conditionally on
some predefined levels for the non stationary variables. In both cases,
the idea is to use dummy values for the steady state level of the
exogenous non stationary variables.

Note that the nonstationary variables in the model must be integrated
processes (their first difference or k-difference must be stationary).


## Calibrated Smoother

Dynare can also run the smoother on a calibrated model:


### Julia functions

```@docs
calibsmoother!
```

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
