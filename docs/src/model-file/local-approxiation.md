

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

The matrix of second order derivatives is ``n \times n^2`` matrix where each column contains derivatives with respect to a pair of endogenous variables
  -  eigenvalues::Vector{Complex{Float64}}
  -  g1::Matrix{Float64}  # full approximation
  -  gs1::Matrix{Float64} # state transition matrices: states x states
  -  hs1::Matrix{Float64} # states x shocks
  -  gns1::Matrix{Float64} # non states x states
  -  hns1::Matrix{Float64} # non states x shocsks
  -  g1_1::SubArray{Float64, 2, Matrix{Float64},
     Tuple{Base.Slice{Base.OneTo{Int}}, UnitRange{Int}}, true}	 # solution first order derivatives w.r. to state variables
  -  g1_2::SubArray{Float64, 2, Matrix{Float64},
     Tuple{Base.Slice{Base.OneTo{Int}}, UnitRange{Int}}, true}   # solution first order derivatives w.r. to current exogenous variables
  -  endogenous_variance::Matrix{Float64}
  -  stationary_variables::Vector{Bool}


##### Example 1

```
shocks;
var e;
stderr 0.0348;
end;

stoch_simul;
```

Performs the simulation of the 1st-order approximation of a model with
a single stochastic shock `e`, with a standard error of `0.0348`.

###### Example 2

```
    stoch_simul(irf=60);
```

Performs the simulation of a model and displays impulse response
functions on 60 periods.

### Julia function
```@docs
localapproximation!
```

### First-order approximation

The approximation has the stylized form:
```math
y_t = y^s + A \phi(y_{t-1}) + B u_t
```

where $y^s$ is the steady state value of $y$ and
$\phi(y_{t-1})=y_{t-1}-y^s$. Matrices of coefficients $A$ and $B$ are
computed by Dynare.


### Second-order approximation

The approximation has the form:
```math
y_t = y^s + 0.5 \Delta^2 + A \phi(y_{t-1}) + B u_t + 0.5 C
(\phi(y_{t-1})\otimes \phi(y_{t-1})) + 0.5 D (u_t \otimes u_t) + E
(\phi(y_{t-1}) \otimes u_t)
```

where $y^s$ is the steady state value of $y$, $\phi(y_{t-1})=y_{t-1}-y^s$, and
$\Delta^2$ is the shift effect of the variance of future
shocks. Matrices of coefficients $A$, $B$, $C$, $D$ and $E$ are
computed by Dynare. 


