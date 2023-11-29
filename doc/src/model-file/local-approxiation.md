

In a stochastic context, Dynare computes one or several simulations
corresponding to a random draw of the shocks.

The main algorithm for solving stochastic models relies on a Taylor
approximation, up to second order, of the solution function (see
*Judd (1996)*, *Collard and Juillard (2001a, 2001b)*, and *Schmitt-Grohé
and Uríbe (2004)*). The details of the Dynare implementation of the
first order solution are given in *Villemot (2011)*. Such a solution is
computed using the `stoch_simul` command.

### Dynare commands

#### stoch\_simul

*Command*: `stoch_simul;
*Command*: `stoch_simul (OPTIONS...);

Solves a stochastic (i.e. rational expectations) model, using
perturbation techniques.

More precisely, `stoch_simul` computes a Taylor approximation of the
model around the deterministic steady state and solves of the the
decision and transition functions for the approximated model. Using
this, it computes impulse response functions and various descriptive
statistics (moments, variance decomposition, correlation and
autocorrelation coefficients). For correlated shocks, the variance
decomposition is computed as in the VAR literature through a Cholesky
decomposition of the covariance matrix of the exogenous variables. When
the shocks are correlated, the variance decomposition depends upon the
order of the variables in the `varexo` command.

The IRFs are computed as the difference between the trajectory of a
variable following a shock at the beginning of period `1` and its steady
state value. More details on the computation of IRFs can be found at
[https://archives.dynare.org/DynareWiki/IrFs](https://archives.dynare.org/DynareWiki/IrFs).

Variance decomposition, correlation, autocorrelation are only displayed
for variables with strictly positive variance. Impulse response
functions are only plotted for variables with response larger than
$10^{-10}$.

Variance decomposition is computed relative to the sum of the
contribution of each shock. Normally, this is of course equal to
aggregate variance, but if a model generates very large variances, it
may happen that, due to numerical error, the two differ by a significant
amount. Dynare issues a warning if the maximum relative difference
between the sum of the contribution of each shock and aggregate variance
is larger than `0.01%`.

The covariance matrix of the shocks is specified with the `shocks`
command (see `shocks-exo`).

##### Options

- ar = INTEGER

Order of autocorrelation coefficients to compute. Default: 5

- `irf = INTEGER`

Number of periods on which to compute the IRFs. Setting `irf=0`
suppresses the plotting of IRFs. Default: `40`.

- `nonstationary`: declares the model as nonstationary.

- `noprint`: don't print the results

- `order = INTEGER`

Order of Taylor approximation. Note that for third order and above, the
`k_order_solver` option is implied and only empirical moments are
available (you must provide a value for `periods` option). Default: `2`

- `periods = INTEGER`

If different from zero, empirical moments will be computed instead of
theoretical moments. The value of the option specifies the number of
periods to use in the simulations. Values of the initval block, possibly
recomputed by `steady`, will be used as starting point for the
simulation. The simulated endogenous variables are made available to the
user in Julia variable `context.results.model_results[1].simulation`. Default: `0`.


- `dr = OPTION`

Determines the method used to compute the decision rule. Possible values
for OPTION are:

  `default`

  Uses the default method to compute the decision rule based on the
    generalized Schur decomposition (see *Villemot (2011)* for more
    information).

  `cycle_reduction`

  Uses the cycle reduction algorithm to solve the polynomial equation
    for retrieving the coefficients associated to the endogenous
    variables in the decision rule. This method is faster than the
    default one for large scale models.

Default value is `default`.

##### Output

The derivatives of the approximated solution function are availabe in the
vector of matrices `context.results.model_results[1].solution_derivatives`.
The first element contains the matrix of first order derivatives. The second
element, the matrix of second order derivatives.

The matrix of first order derivatives is a ``n x (n_s + n_x + 1)`` matrix where `n` is the number
of endogenous variables, ``n_s``, the number of state variables (variables appearing in the model with a lag), and ``n_x``, the number of exogenous variables. An element of this matrix is
```math
\begin{align*}
X_{i,j} &= \frac{\partial g_i}{\partial y_j},\;\;j=1,\ldots,n_s\\
X_{i,n_s+j} &= \frac{\partial g_i}{\partial x_j},\;\;j=1,\ldots,n_x\\
X_{i,n_s+n_k+1} &= \frac{\partial g_i}{\partial \sigma} = 0
\end{align*}
```
where ``g_i`` is the solution function for variable `i`, `y`, the vector of endogenous variables, `x`, the vector en exogenous variables and ``\sigma`` the stochastic scale of the model. Note that at first order, this derivative is alwasy equal to zero. 

The matrix of second order derivatives is ``n x n^2`` matrix where each column contains derivatives with respect to a pair of endogenous variab
eigenvalues::Vector{Complex{Float64}}
    g1::Matrix{Float64}  # full approximation
    gs1::Matrix{Float64} # state transition matrices: states x states
    hs1::Matrix{Float64} # states x shocks
    gns1::Matrix{Float64} # non states x states
    hns1::Matrix{Float64} # non states x shocsks
    # solution first order derivatives w.r. to state variables
    g1_1::SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int}}, UnitRange{Int}}, true}
    # solution first order derivatives w.r. to current exogenous variables
    g1_2::SubArray{Float64, 2, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int}}, UnitRange{Int}}, true}
    endogenous_variance::Matrix{Float64}
    #    g1_3::SubArray # solution first order derivatives w.r. to lagged exogenous variables
    stationary_variables::Vector{Bool}

This command sets `oo_.dr`, `oo_.mean`, `oo_.var`, `oo_.var_list`, and
`oo_.autocorr`, which are described below.

If the `periods` option is present, sets `oo_.skewness`, `oo_.kurtosis`,
and `oo_.endo_simul` (see `oo_.endo_simul`{.interpreted-text
role="mvar"}), and also saves the simulated variables in MATLAB/Octave
vectors of the global workspace with the same name as the endogenous
variables.

If option `irf` is different from zero, sets `oo_.irfs` (see below) and
also saves the IRFs in MATLAB/Octave vectors of the global workspace
(this latter way of accessing the IRFs is deprecated and will disappear
in a future version).

If the option `contemporaneous_correlation` is different from `0`, sets
`oo_.contemporaneous_correlation`, which is described below.

*Example*

```
shocks;
var e;
stderr 0.0348;
end;

stoch_simul;
```

Performs the simulation of the 2nd-order approximation of a model with
a single stochastic shock `e`, with a standard error of `0.0348`.

*Example*

```
    stoch_simul(irf=60) y k;
```

Performs the simulation of a model and displays impulse response
functions on 60 periods for variables `y` and `k`.


*MATLAB/Octave*: `oo.irfs`

After a run of `stoch_simul` with option `irf` different from zero,
contains the impulse responses, with the following naming convention:
[VARIABLE_NAME_SHOCK_NAME]{.title-ref}.

For example, `oo_.irfs.gnp_ea` contains the effect on `gnp` of a
one-standard deviation shock on `ea`.

*MATLAB/Octave*: `get_irf ('EXOGENOUS_NAME' [, 'ENDOGENOUS_NAME']... );`


### Typology and ordering of variables

Dynare distinguishes four types of endogenous variables:

*Purely backward (or purely predetermined) variables*

Those that appear only at current and past period in the model, but
not at future period (i.e. at $t$ and $t-1$ but not $t+1$). The number
of such variables is equal to `M_.npred`.

*Purely forward variables*

Those that appear only at current and future period in the model, but
not at past period (i.e. at $t$ and $t+1$ but not $t-1$). The number
of such variables is stored in `M_.nfwrd`.

*Mixed variables*

Those that appear at current, past and future period in the model
(i.e. at $t$, $t+1$ and $t-1$). The number of such variables is stored
in `M_.nboth`.

*Static variables*

Those that appear only at current, not past and future period in the
model (i.e. only at $t$, not at $t+1$ or $t-1$). The number of such
variables is stored in `M_.nstatic`.

Note that all endogenous variables fall into one of these four
categories, since after the creation of auxiliary variables (see
`aux-variables`), all endogenous have at
most one lead and one lag. We therefore have the following identity:

```
M_.npred + M_.both + M_.nfwrd + M_.nstatic = M_.endo_nbr
```


### First-order approximation

The approximation has the stylized form:

$$y_t = y^s + A y^h_{t-1} + B u_t$$

where $y^s$ is the steady state value of $y$ and $y^h_t=y_t-y^s$.

*MATLAB/Octave variable*: `oo.dr.state_var`

Vector of numerical indices identifying the state variables in the
vector of declared variables, *given the current parameter values* for
which the decision rules have been computed. It may differ from
`M_.state_var` in case a state variable drops from the model given the
current parameterization, because it only gets 0 coefficients in the
decision rules. See `M_.state_var`.
:::

The coefficients of the decision rules are stored as follows:

-   $y^s$ is stored in `oo_.dr.ys`. The vector rows correspond to all
    endogenous in the declaration order.
-   $A$ is stored in `oo_.dr.ghx`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to state
    variables in DR-order, as given by `oo_.dr.state_var`. (N.B.: if the
    `block` option to the `model` block has been specified, then rows
    are in declaration order, and columns are ordered according to
    `oo_.dr.state_var` which may differ from DR-order.)
-   $B$ is stored `oo_.dr.ghu`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to exogenous
    variables in declaration order. (N.B.: if the `block` option to the
    `model` block has been specified, then rows are in declaration
    order.)

Of course, the shown form of the approximation is only stylized, because
it neglects the required different ordering in $y^s$ and $y^h_t$. The
precise form of the approximation that shows the way Dynare deals with
differences between declaration and DR-order, is

$$y_t(\mathrm{oo\_.dr.order\_var}) =
y^s(\mathrm{oo\_.dr.order\_var}) + A \cdot
y_{t-1}(\mathrm{oo\_.dr.order\_var(k2)}) -
y^s(\mathrm{oo\_.dr.order\_var(k2)}) + B\cdot u_t$$

where $\mathrm{k2}$ selects the state variables, $y_t$ and $y^s$ are in
declaration order and the coefficient matrices are in DR-order.
Effectively, all variables on the right hand side are brought into DR
order for computations and then assigned to $y_t$ in declaration order.

### Second-order approximation

The approximation has the form:

$$y_t = y^s + 0.5 \Delta^2 + A y^h_{t-1} + B u_t + 0.5 C (y^h_{t-1}\otimes y^h_{t-1}) + 0.5 D (u_t \otimes u_t) + E (y^h_{t-1} \otimes u_t)$$

where $y^s$ is the steady state value of $y$, $y^h_t=y_t-y^s$, and
$\Delta^2$ is the shift effect of the variance of future shocks. For the
reordering required due to differences in declaration and DR order, see
the first order approximation.

The coefficients of the decision rules are stored in the variables
described for first order approximation, plus the following variables:

-   $\Delta^2$ is stored in `oo_.dr.ghs2`. The vector rows correspond to
    all endogenous in DR-order.
-   $C$ is stored in `oo_.dr.ghxx`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of the vector of state variables in DR-order.
-   $D$ is stored in `oo_.dr.ghuu`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of exogenous variables in declaration order.
-   $E$ is stored in `oo_.dr.ghxu`. The matrix rows correspond to all
    endogenous in DR-order. The matrix columns correspond to the
    Kronecker product of the vector of state variables (in DR-order) by
    the vector of exogenous variables (in declaration order).

