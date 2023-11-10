## Steady state

There are three ways of computing the steady state (i.e. the static
equilibrium) of a model. The first way is to provide the equations of
the steady state in a [steady_state_model]{.title-ref} block. When it
is possible to derive the steady state by hand, this is the
recommended way as it faster and more accurate.

The second way is to provide only a partial solution in the
[steady_state_model]{.title-ref} block and to compute the solution for
the other variables numerically. Guess values for these other
variables must be declared in a [initval]{.title-ref} block. The less
variables the better.

The third way is to compute the steady state value of all variables
numerically. There is no [steady_state_model] block and a guess value
must be declared for all variables. A guess value of 0 can be omitted,
but be careful with variables appearing at the denominator of a fraction. 

### Finding the steady state with Dynare nonlinear solver

#### Computing the steady state
- *Command*: `steady ;` 

- *Command*: `steady (OPTIONS...);`

This command computes the steady state of a model using a nonlinear
Newton-type solver and displays it.

More precisely, it computes the equilibrium value of the endogenous
variables for the value of the exogenous variables specified in the
previous `initval` or `endval` block.

`steady` uses an iterative procedure and takes as initial guess the
value of the endogenous variables set in the previous `initval` or
`endval` block.

For complicated models, finding good numerical initial values for the
endogenous variables is the trickiest part of finding the equilibrium of
that model. Often, it is better to start with a smaller model and add
new variables one by one.

#### Options

- `maxit = INTEGER`

Determines the maximum number of iterations used in the non-linear
solver. The default value of `maxit` is 50.

- `tolf = DOUBLE`

Convergence criterion for termination based on the function value.
Iteration will cease when the residuals are smaller than `tolf`.
Default: `eps^(1/3)`

- `tolx = DOUBLE`

Convergence criterion for termination based on the step tolerance along.
Iteration will cease when the attempted step size is smaller than
`tolx`. Default: `eps^(2/3)`

- `homotopy_mode = INTEGER`

Use a homotopy (or divide-and-conquer) technique to solve for the steady
state. If you use this option, you must specify a `homotopy_setup`
block. This option can take three possible values:

- `1`: In this mode, all the parameters are changed simultaneously, and the
       distance between the boundaries for each parameter is divided in as
       many intervals as there are steps (as defined by the
       `homotopy_steps` option); the problem is solved as many times as
       there are steps.

- `2`: Same as mode `1`, except that only one parameter is changed at a
       time; the problem is solved as many times as steps times number of
       parameters. [Not yet implemented]

- `3`: Dynare tries first the most extreme values. If it fails to compute
       the steady state, the interval between initial and desired values is
       divided by two for all parameters. Every time that it is impossible
       to find a steady state, the previous interval is divided by two.
       When it succeeds to find a steady state, the previous interval is
       multiplied by two. In that last case `homotopy_steps` contains the
       maximum number of computations attempted before giving up. [Not
       yet implemented]

- `homotopy_steps = INTEGER`

Defines the number of steps when performing a homotopy. See
`homotopy_mode` option for more details.

#### Output

After computation, the steady state is available in the following
variable:

*Julia variable*: `context.results.model_results[1].trends.endogenous_steady_state`

Contains the computed steady state. Endogenous variables are ordered in
the order of declaration used in the `var` command,

#### Provinding guess values
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

```

#### Homotopy setup

*Block*: `homotopy_setup ;`

This block is used to declare initial and final values for the
parameters and exogenous variables when using a
homotopy method. It is used in conjunction with the option
`homotopy_mode` of the steady command.

The idea of homotopy is
to subdivide the problem of finding the steady state into smaller
problems. It assumes that you know how to compute the steady state for a
given set of parameters, and it helps you finding the steady state for
another set of parameters, by incrementally moving from one to another
set of parameters.

The purpose of the `homotopy_setup` block is to declare the final (and
possibly also the initial) values for the parameters or exogenous that
will be changed during the homotopy. It should contain lines of the
form:

```
VARIABLE_NAME, EXPRESSION, EXPRESSION;
```

This syntax specifies the initial and final values of a given
parameter/exogenous.

There is an alternative syntax:

```
VARIABLE_NAME, EXPRESSION;
```

Here only the final value is specified for a given parameter/exogenous;
the initial value is taken from the preceeding `initval` block.

A necessary condition for a successful homotopy is that Dynare must be
able to solve the steady state for the initial parameters/exogenous
without additional help (using the guess values given in the `initval`
block).

If the homotopy fails, a possible solution is to increase the number of
steps (given in `homotopy_steps` option of `steady`).

*Example*

In the following example, Dynare will first compute the steady state for
the initial values (`gam=0.5` and `x=1`), and then subdivide the problem
into 50 smaller problems to find the steady state for the final values
(`gam=2` and `x=2`):

```
var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
delt=0.02;
aa=0.5;
bet=0.05;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(aa*x*alph))^(1/(alph-1));
c = aa*x*k^alph-delt*k;
end;

homotopy_setup;
gam, 0.5, 2;
x, 2;
end;

steady(homotopy_mode = 1, homotopy_steps = 50);
```

### Providing the steady state to Dynare

If you know how to compute the steady state for your model, you can
provide a MATLAB/Octave function doing the computation instead of using
`steady`. Again, there are two options for doing that:

-   The easiest way is to write a `steady_state_model` block, which is
    described below in more details. See also `fs2000.mod` in the
    `examples` directory for an example. The steady state file
    generated by Dynare will be called `+FILENAME/steadystate.m.`
-   You can write the corresponding MATLAB function by hand. If your
    MOD-file is called `FILENAME.mod`, the steady state file must be
    called `FILENAME_steadystate.m`. See `NK_baseline_steadystate.m`
    in the examples directory for an example. This option gives a bit
    more flexibility (loops and conditional structures can be used),
    at the expense of a heavier programming burden and a lesser
    efficiency.

Note that both files allow to update parameters in each call of the
function. This allows for example to calibrate a model to a labor supply
of 0.2 in steady state by setting the labor disutility parameter to a
corresponding value (see `NK_baseline_steadystate.m` in the `examples`
directory). They can also be used in estimation where some parameter may
be a function of an estimated parameter and needs to be updated for
every parameter draw. For example, one might want to set the capital
utilization cost parameter as a function of the discount rate to ensure
that capacity utilization is 1 in steady state. Treating both parameters
as independent or not updating one as a function of the other would lead
to wrong results. But this also means that care is required. Do not
accidentally overwrite your parameters with new values as it will lead
to wrong results.

*Block*: `steady_state_model ;`

When the analytical solution of the model is known, this command can be
used to help Dynare find the steady state in a more efficient and
reliable way, especially during estimation where the steady state has to
be recomputed for every point in the parameter space.

Each line of this block consists of a variable (either an endogenous, a
temporary variable or a parameter) which is assigned an expression
(which can contain parameters, exogenous at the steady state, or any
endogenous or temporary variable already declared above). Each line
therefore looks like:

```
VARIABLE_NAME = EXPRESSION;
```

Note that it is also possible to assign several variables at the same
time, if the main function in the right hand side is a MATLAB/Octave
function returning several arguments:

```
[ VARIABLE_NAME, VARIABLE_NAME... ] = EXPRESSION;
```

Dynare will automatically generate a steady state file (of the form
`+FILENAME/steadystate.m`) using the information provided in this block.

*Steady state file for deterministic models*

The `steady_state_model` block also works with deterministic models. An
`initval` block and, when necessary, an `endval` block, is used to set
the value of the exogenous variables. Each `initval` or `endval` block
must be followed by `steady` to execute the function created by
`steady_state_model` and set the initial, respectively terminal, steady
state.

*Example*

```
var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

...
// parameter calibration, (dynamic) model declaration, shock calibration...
...

steady_state_model;
dA = exp(gam);
gst = 1/dA; // A temporary variable
m = mst;

// Three other temporary variables
khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );

n  = xist/(nust+xist);
P  = xist + nust;
k  = khst*n;

l  = psi*mst*n/( (1-psi)*(1-n) );
c  = mst/P;
d  = l - mst + 1;
y  = k^alp*n^(1-alp)*gst^alp;
R  = mst/bet;

// You can use MATLAB functions which return several arguments
[W, e] = my_function(l, n);

gp_obs = m/dA;
gy_obs = dA;
end;

steady;
```

### Replace some equations during steady state computations {#eq-tag-ss}

When there is no steady state file, Dynare computes the steady state by
solving the static model, i.e. the model from the `.mod` file from which
leads and lags have been removed.

In some specific cases, one may want to have more control over the way
this static model is created. Dynare therefore offers the possibility to
explicitly give the form of equations that should be in the static
model.

More precisely, if an equation is prepended by a `[static]` tag, then it
will appear in the static model used for steady state computation, but
that equation will not be used for other computations. For every
equation tagged in this way, you must tag another equation with
`[dynamic]`: that equation will not be used for steady state
computation, but will be used for other computations.

This functionality can be useful on models with a unit root, where there
is an infinity of steady states. An equation (tagged `[dynamic]`) would
give the law of motion of the nonstationary variable (like a random
walk). To pin down one specific steady state, an equation tagged
`[static]` would affect a constant value to the nonstationary variable.
Another situation where the `[static]` tag can be useful is when one has
only a partial closed form solution for the steady state.

*Example*

This is a trivial example with two endogenous variables. The second
equation takes a different form in the static model:

```
var c k;
varexo x;
...
model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
[dynamic] c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
[static] k = ((delt+bet)/(x*aa*alph))^(1/(alph-1));
end;
```
