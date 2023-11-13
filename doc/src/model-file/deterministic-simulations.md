## Deterministic simulation

### Perfect foresight

When the framework is deterministic, Dynare can be used for models with
the assumption of perfect foresight. Typically, the system is supposed
to be in a state of equilibrium before a period `1` when the news of a
contemporaneous or of a future shock is learned by the agents in the
model. The purpose of the simulation is to describe the reaction in
anticipation of, then in reaction to the shock, until the system returns
to the old or to a new state of equilibrium. In most models, this return
to equilibrium is only an asymptotic phenomenon, which one must
approximate by an horizon of simulation far enough in the future.
Another exercise for which Dynare is well suited is to study the
transition path to a new equilibrium following a permanent shock. For
deterministic simulations, the numerical problem consists of solving a
nonlinear system of simultaneous equations in `n` endogenous variables
in `T` periods. Dynare offers several algorithms for solving this
problem, which can be chosen via the `stack_solve_algo` option. By
default (`stack_solve_algo=0`), Dynare uses a Newton-type method to
solve the simultaneous equation system. Because the resulting Jacobian
is in the order of `n` by `T` and hence will be very large for long
simulations with many variables, Dynare makes use of the sparse matrix
capacities of MATLAB/Octave. A slower but potentially less memory
consuming alternative (`stack_solve_algo=1`) is based on a Newton-type
algorithm first proposed by *Laffargue (1990)* and *Boucekkine (1995)*,
which avoids ever storing the full Jacobian. The details of the
algorithm can be found in *Juillard (1996)*. The third type of
algorithms makes use of block decomposition techniques
(divide-and-conquer methods) that exploit the structure of the model.
The principle is to identify recursive and simultaneous blocks in the
model structure and use this information to aid the solution process.
These solution algorithms can provide a significant speed-up on large
models.

!!! note
    Be careful when employing auxiliary variables in the context of perfect
    foresight computations. The same model may work for stochastic
    simulations, but fail for perfect foresight simulations. The issue
    arises when an equation suddenly only contains variables dated `t+1` (or
    `t-1` for that matter). In this case, the derivative in the last (first)
    period with respect to all variables will be 0, rendering the stacked
    Jacobian singular.

*Example*

Consider the following specification of an Euler equation with log
utility:

```
Lambda = beta*C(-1)/C;
Lambda(+1)*R(+1)= 1;
```

Clearly, the derivative of the second equation with respect to all
endogenous variables at time `t` is zero, causing
`perfect_foresight_solver` to generally fail. This is due to the use
of the Lagrange multiplier `Lambda` as an auxiliary variable. Instead,
employing the identical

```
beta*C/C(+1)*R(+1)= 1;
```

will work.

*Command*: `perfect_foresight_setup ;

*Command*: `perfect_foresight_setup (OPTIONS...);`

Prepares a perfect foresight simulation, by extracting the information
in the `initval`, `endval` and `shocks` blocks and converting them into
simulation paths for exogenous and endogenous variables.

This command must always be called before running the simulation with
`perfect_foresight_solver`.

*Options*

- `periods = INTEGER`

Number of periods of the simulation.

- `datafile = FILENAME`

Used to specify path for all endogenous and exogenous variables.
Strictly equivalent to `initval_file`{.interpreted-text role="comm"}.

*Output*

The paths for the exogenous variables are stored into `oo_.exo_simul`.

The initial and terminal conditions for the endogenous variables and the
initial guess for the path of endogenous variables are stored into
`oo_.endo_simul`.

*Command*: `perfect_foresight_solver ;`

*Command*: `perfect_foresight_solver (OPTIONS...);`

Computes the perfect foresight (or deterministic) simulation of the
model.

Note that `perfect_foresight_setup` must be called before this command,
in order to setup the environment for the simulation.

*Options*

- `maxit = INTEGER`

Determines the maximum number of iterations used in the non-linear
solver. The default value of `maxit` is `50`.

- `tolf = DOUBLE`

Convergence criterion for termination based on the function value.
Iteration will cease when it proves impossible to improve the function
value by more than `tolf`. Default: `1e-5`

- `tolx = DOUBLE`

Convergence criterion for termination based on the change in the
function argument. Iteration will cease when the solver attempts to take
a step that is smaller than `tolx`. Default: `1e-5`

- `noprint`

Don't print anything. Useful for loops.

- `print`

Print results (opposite of `noprint`).

- `lmmcp`

Solves the perfect foresight model with a Levenberg-Marquardt mixed
complementarity problem (LMMCP) solver (*Kanzow and Petra (2004)*),
which allows to consider inequality constraints on the endogenous
variables (such as a ZLB on the nominal interest rate or a model with
irreversible investment). This option is equivalent to
`stack_solve_algo=7` **and** `solve_algo=10`. Using the LMMCP solver
requires a particular model setup as the goal is to get rid of any
min/max operators and complementary slackness conditions that might
introduce a singularity into the Jacobian. This is done by attaching an
equation tag (see `model-decl`) with the
`mcp` keyword to affected equations. This tag states that the equation
to which the tag is attached has to hold unless the expression within
the tag is binding. For instance, a ZLB on the nominal interest rate
would be specified as follows in the model block:

```
    model;
       ...
       [mcp = 'r > -1.94478']
       r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;
       ...
    end;
```

where `1.94478` is the steady state level of the nominal interest rate
and `r` is the nominal interest rate in deviation from the steady state.
This construct implies that the Taylor rule is operative, unless the
implied interest rate `r<=-1.94478`, in which case the `r` is fixed at
`-1.94478` (thereby being equivalent to a complementary slackness
condition). By restricting the value of `r` coming out of this equation,
the `mcp` tag also avoids using `max(r,-1.94478)` for other occurrences
of `r` in the rest of the model. It is important to keep in mind that,
because the `mcp` tag effectively replaces a complementary slackness
condition, it cannot be simply attached to any equation. Rather, it must
be attached to the correct affected equation as otherwise the solver
will solve a different problem than originally intended. Also, since the
problem to be solved is nonlinear, the sign of the residuals of the
dynamic equation matters. In the previous example, for the nominal
interest rate rule, if the LHS and RHS are reversed the sign of the
residuals (the difference between the LHS and the RHS) will change and
it may happen that solver fails to identify the solution path. More
generally, convergence of the nonlinear solver is not guaranteed when
using mathematically equivalent representations of the same equation.

Note that in the current implementation, the content of the `mcp`
equation tag is not parsed by the preprocessor. The inequalities must
therefore be as simple as possible: an endogenous variable, followed by
a relational operator, followed by a number (not a variable, parameter
or expression).

- `endogenous_terminal_period`

The number of periods is not constant across Newton iterations when
solving the perfect foresight model. The size of the nonlinear system of
equations is reduced by removing the portion of the paths (and
associated equations) for which the solution has already been identified
(up to the tolerance parameter). This strategy can be interpreted as a
mix of the shooting and relaxation approaches. Note that round off
errors are more important with this mixed strategy (user should check
the reported value of the maximum absolute error). Only available with
option `stack_solve_algo==0`.

- `linear_approximation`

Solves the linearized version of the perfect foresight model. The model
must be stationary. Only available with option `stack_solve_algo==0` or
`stack_solve_algo==7`.

### Julia function
```@docs
perfect_foresight!
```

*Output*

The simulated endogenous variables are available in
`context.results.model_results[1].simulations`. This is a vector of
`AxisArrayTable`, one for each simulations stored in `context`. Each
`AxisArrayTable` contains the trajectories for endogenous and
exogenous variables


