When the framework is deterministic, Dynare can be used for models with
the assumption of perfect foresight. The system is supposed
to be in a given state before a period `1` (often a steady state) 
when the news of a
contemporaneous or of a future shock is learned by the agents in the
model. The purpose of the simulation is to describe the reaction in
anticipation of, then in reaction to the shock, until the system returns
to equilibrium. This return
to equilibrium is only an asymptotic phenomenon, which one must
approximate by an horizon of simulation far enough in the future.
Another exercise for which Dynare is well suited is to study the
transition path to a new equilibrium following a permanent shock. For
deterministic simulations, the numerical problem consists of solving a
nonlinear system of simultaneous equations in `n` endogenous variables
in `T` periods. Dynare uses a Newton-type method to
solve the simultaneous equation system. Because the resulting Jacobian
is in the order of `n` by `T` and hence will be very large for long
simulations with many variables, Dynare makes use of the sparse matrix
code .

### Dynare commands

#### perfect\_foresight\_setup
*Command*: `perfect\_foresight\_setup ;

*Command*: `perfect\_foresight\_setup (OPTIONS...);`

Prepares a perfect foresight simulation, by extracting the information
in the `initval`, `endval` and `shocks` blocks and converting them into
simulation paths for exogenous and endogenous variables.

This command must always be called before running the simulation with
`perfect\_foresight\_solver`.

##### Options

- `periods = INTEGER`

Number of periods of the simulation.

- `datafile = FILENAME`

Used to specify path for all endogenous and exogenous variables.
Strictly equivalent to `initval_file`.

##### Output

The paths for the exogenous variables are stored into `context.results.model_resultst[1].simulations`.

The initial and terminal conditions for the endogenous variables and the
initial guess for the path of endogenous variables are stored into
`context.results.model_results[1].simulations`.

#### perfect\_foresight\_solver

*Command*: `perfect\_foresight\_solver ;`

*Command*: `perfect\_foresight\_solver (OPTIONS...);`

Computes the perfect foresight (or deterministic) simulation of the
model.

Note that `perfect\_foresight\_setup` must be called before this command,
in order to setup the environment for the simulation.

##### Options

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

Solves mixed complementarity problems (the term refers to the LMMCP solver (_Kanzow and Petra, 2004_),
that is used by DynareMatlab.  DynareJulia uses the PATHSovler package)

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

#### Remark
    Be careful when employing auxiliary variables in the context of perfect
    foresight computations. The same model may work for stochastic
    simulations, but fail for perfect foresight simulations. The issue
    arises when an equation suddenly only contains variables dated `t+1` (or
    `t-1` for that matter). In this case, the derivative in the last (first)
    period with respect to all variables will be 0, rendering the stacked
    Jacobian singular.

##### Example

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


### Julia function
```@docs
perfect_foresight!
```

##### Output

The simulated endogenous variables are available in
`context.results.model_results[1].simulations`. This is a vector of
`AxisArrayTable`, one for each simulations stored in `context`. Each
`AxisArrayTable` contains the trajectories for endogenous and
exogenous variables

### Solving mixed complementarity problems

requires a particular model setup as the goal is to get rid of any
min/max operators and complementary slackness conditions that might
introduce a singularity into the Jacobian. This is done by attaching an
equation tag (see `model-decl`) with the
`mcp` keyword to affected equations. The format of the `mcp` tag is
```
[mcp = 'VARIABBLENAME OP CONSTANT']
```
where VARIABLENAME is an endogenous variable and OP is either `>` or `<`.
For complicated occasionally binding constraints, it may be necessary
to declare a new endogenous variable.

This tag states that the equation
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

where  `r` is the nominal interest rate in deviation from the steady state.
This construct implies that the Taylor rule is operative, unless the
implied interest rate `r<=-1.94478`, in which case the `r` is fixed at
`-1.94478`. This is equavalant to

```math
(r_t > -1.94478)\;\; \bot\;\; r_t = \rho r_{t-1} + (1-\rho) (g_\pi Infl_t+g_y YGap_t) + e_t
```

By restricting the value of `r` coming out of this equation,
the `mcp` tag also avoids using `max(r,-1.94478)` for other occurrences
of `r` in the rest of the model. It is important to keep in mind that,
because the `mcp` tag effectively replaces a complementary slackness
condition, it cannot be simply attached to any equation.

Note that in the current implementation, the content of the `mcp`
equation tag is not parsed by the preprocessor. The inequalities must
therefore be as simple as possible: an endogenous variable, followed by
a relational operator, followed by a number (not a variable, parameter
or expression).


