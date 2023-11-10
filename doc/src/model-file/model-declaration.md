## Model declaration

The model is declared inside a `model` block:

*Block*: `model ;` 

*Block*: `model (OPTIONS...);`

The equations of the model are written in a block delimited by `model`
and `end` keywords.

There must be as many equations as there are endogenous variables in the
model, except when computing the unconstrained optimal policy with
`ramsey_model`, `ramsey_policy` or `discretionary_policy`.

The syntax of equations must follow the conventions for
`MODEL_EXPRESSION` as described in `expr`{.interpreted-text role="ref"}.
Each equation must be terminated by a semicolon (';'). A normal equation
looks like:

```
MODEL_EXPRESSION = MODEL_EXPRESSION;
```

When the equations are written in homogenous form, it is possible to
omit the '=0' part and write only the left hand side of the equation. A
homogenous equation looks like:

``` 
MODEL_EXPRESSION;
```

Inside the model block, Dynare allows the creation of *model-local
variables*, which constitute a simple way to share a common expression
between several equations. The syntax consists of a pound sign (\#)
followed by the name of the new model local variable (which must **not**
be declared as in `var-decl`{.interpreted-text role="ref"}, but may have
been declared by `model_local_variable`{.interpreted-text role="comm"}),
an equal sign, and the expression for which this new variable will
stand. Later on, every time this variable appears in the model, Dynare
will substitute it by the expression assigned to the variable. Note that
the scope of this variable is restricted to the model block; it cannot
be used outside. To assign a LaTeX name to the model local variable, use
the declaration syntax outlined by
`model_local_variable`{.interpreted-text role="comm"}. A model local
variable declaration looks like:

```
#VARIABLE_NAME = MODEL_EXPRESSION;
```

It is possible to tag equations written in the model block. A tag can
serve different purposes by allowing the user to attach arbitrary
informations to each equation and to recover them at runtime. For
instance, it is possible to name the equations with a `name`-tag, using
a syntax like:

```
model;

[name = 'Budget constraint'];
c + k = k^theta*A;

end;
```

Here, `name` is the keyword indicating that the tag names the equation.
If an equation of the model is tagged with a name, the `resid` command
will display the name of the equations (which may be more informative
than the equation numbers) in addition to the equation number. Several
tags for one equation can be separated using a comma:

```   
model;

[name='Taylor rule',mcp = 'r > -1.94478']
r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;

end;
```

More information on tags is available at
[https://git.dynare.org/Dynare/dynare/-/wikis/Equations-Tags](https://git.dynare.org/Dynare/dynare/-/wikis/Equations-Tags).

There can be several `model` blocks, in which case they are simply
concatenated. The set of effective options is also the concatenation of
the options declared in all the blocks, but in that case you may rather
want to use the `model_options`{.interpreted-text role="comm"} command.

*Options*

- `linear`

Declares the model as being linear. It spares oneself from having to
declare initial values for computing the steady state of a stationary
linear model. This option can't be used with non-linear models, it will
NOT trigger linearization of the model.

- `no_static`

Don't create the static model file. This can be useful for models which
don't have a steady state.

- `differentiate_forward_vars` `differentiate_forward_vars = (VARIABLE_NAME [VARIABLE_NAME ...] )`

Tells Dynare to create a new auxiliary variable for each endogenous
variable that appears with a lead, such that the new variable is the
time differentiate of the original one. More precisely, if the model
contains `x(+1)`, then a variable `AUX_DIFF_VAR` will be created such
that `AUX_DIFF_VAR=x-x(-1)`, and `x(+1)` will be replaced with
`x+AUX_DIFF_VAR(+1)`.

The transformation is applied to all endogenous variables with a lead if
the option is given without a list of variables. If there is a list, the
transformation is restricted to endogenous with a lead that also appear
in the list.

This option can useful for some deterministic simulations where
convergence is hard to obtain. Bad values for terminal conditions in the
case of very persistent dynamics or permanent shocks can hinder correct
solutions or any convergence. The new differentiated variables have
obvious zero terminal conditions (if the terminal condition is a steady
state) and this in many cases helps convergence of simulations.

- `parallel_local_files = ( FILENAME [, FILENAME]... )`

Declares a list of extra files that should be transferred to follower
nodes when doing a parallel computation (see
`paral-conf`{.interpreted-text role="ref"}).

- `balanced_growth_test_tol = DOUBLE`

Tolerance used for determining whether cross-derivatives are zero in the
test for balanced growth path (the latter is documented on
[https://archives.dynare.org/DynareWiki/RemovingTrends](https://archives.dynare.org/DynareWiki/RemovingTrends)). Default:
`1e-6`

*Example* (Elementary RBC model)

```
var c k;
varexo x;
parameters aa alph bet delt gam;

model;
c =  - k + aa*x*k(-1)^alph + (1-delt)*k(-1);
c^(-gam) = (aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam)/(1+bet);
end;
```

*Example* (Use of model local variables)

The following program:

```
model;
# gamma = 1 - 1/sigma;
u1 = c1^gamma/gamma;
u2 = c2^gamma/gamma;
end;
```

...is formally equivalent to:

```
model;
u1 = c1^(1-1/sigma)/(1-1/sigma);
u2 = c2^(1-1/sigma)/(1-1/sigma);
end;
```

*Example* (A linear model)

```
model(linear);
x = a*x(-1)+b*y(+1)+e_x;
y = d*y(-1)+e_y;
end;
```

- `model_options (OPTIONS...);`

This command accepts the same options as the `model`{.interpreted-text
role="bck"} block.

The purpose of this statement is to specify the options that apply to
the whole model, when there are several `model` blocks, so as to restore
the symmetry between those blocks (since otherwise one `model` block
would typically bear the options, while the other ones would typically
have no option).

- `model_remove (TAGS...);`

This command removes equations that appeared in a previous
`model`{.interpreted-text role="bck"} block.

The equations must be specified by a list of tag values, separated by
commas. Each element of the list is either a simple quoted string, in
which case it designates an equation by its `name` tag; or a tag name
(without quotes), followed by an equal sign, then by the tag value
(within quotes).

Each removed equation must either have an `endogenous` tag, or have a
left hand side containing a single endogenous variable. The
corresponding endogenous variable will be either turned into an
exogenous (if it is still used in somewhere in the model at that point),
otherwise it will be removed from the model.

*Example*

```
var c k dummy1 dummy2;

model;
  c + k - aa*x*k(-1)^alph - (1-delt)*k(-1) + dummy1;
  c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
  [ name = 'eq:dummy1', endogenous = 'dummy1' ]
  c*k = dummy1;
  [ foo = 'eq:dummy2' ]
  log(dummy2) = k + 2;
end;

  model_remove('eq:dummy1', foo = 'eq:dummy2');
```
 
In the above example, the last two equations will be removed, `dummy1`
will be turned into an exogenous, and `dummy2` will be removed.

- *block*: `model_replace (TAGS...);`

This block replaces several equations in the model. It removes the
equations given by the tags list (with the same syntax as in
`model_remove`{.interpreted-text role="comm"}), and it adds equations
given within the block (with the same syntax as
`model`{.interpreted-text role="bck"}).

No variable is removed or has its type changed in the process.

*Example*

```
var c k;

model;
  c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
  [ name = 'dummy' ]
  c*k = 1;
end;

model_replace('dummy');
  c^(-gam) = (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;
```
 
In the above example, the dummy equation is replaced by a proper Euler
equation.

Dynare has the ability to output the original list of model equations to
a LaTeX file, using the `write_latex_original_model` command, the list
of transformed model equations using the
`write_latex_dynamic_model command`, and the list of static model
equations using the `write_latex_static_model` command.

- *Command*: `write_latex_original_model (OPTIONS);`

This command creates two LaTeX files: one containing the model as
defined in the model block and one containing the LaTeX document header
information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/original.tex`, which includes a file called
`FILENAME/latex/original_content.tex` (also created by Dynare)
containing the list of all the original model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Time subscripts (`t`, `t+1`, `t-1`, \...) will be appended to the
variable names, as LaTeX subscripts.

Compiling the TeX file requires the following LaTeX packages:
`geometry, fullpage, breqn`.

*Options*

- `write\_equation\_tags`

Write the equation tags in the LaTeX output. The equation tags will be
interpreted with LaTeX markups.

- *Command*: `write_latex_dynamic_model ; 

- *Command*: `write_latex_dynamic_model (OPTIONS);`

This command creates two LaTeX files: one containing the dynamic model
and one containing the LaTeX document header information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/dynamic.tex`, which includes a file called
`FILENAME/latex/dynamic_content.tex` (also created by Dynare) containing
the list of all the dynamic model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Time subscripts (`t`, `t+1`, `t-1`, \...) will be appended to the
variable names, as LaTeX subscripts.

Note that the model written in the TeX file will differ from the model
declared by the user in the following dimensions:

-   The timing convention of predetermined variables (see
    `predetermined_variables`{.interpreted-text role="comm"}) will
    have been changed to the default Dynare timing convention; in
    other words, variables declared as predetermined will be lagged on
    period back,
-   The `EXPECTATION` operators will have been removed, replaced by
    auxiliary variables and new equations (as explained in the
    documentation of
    `EXPECTATION <EXPECTATION (INTEGER) (MODEL_EXPRESSION)>`{.interpreted-text
    role="op"}),
-   Endogenous variables with leads or lags greater or equal than two
    will have been removed, replaced by new auxiliary variables and
    equations,
-   Exogenous variables with leads or lags will also have been
    replaced by new auxiliary variables and equations.

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.

*Options*

- *write\_equation\_tags*

See `write_equation_tags`{.interpreted-text role="opt"}

- *Command*: `write_latex_static_model (OPTIONS);`

This command creates two LaTeX files: one containing the static model
and one containing the LaTeX document header information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/static.tex`, which includes a file called
`FILENAME/latex/static_content.tex` (also created by Dynare) containing
the list of all the steady state model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Note that the model written in the TeX file will differ from the model
declared by the user in the some dimensions (see
`write_latex_dynamic_model`{.interpreted-text role="comm"} for details).

Also note that this command will not output the contents of the optional
`steady_state_model` block (see `steady_state_model`{.interpreted-text
role="bck"}); it will rather output a static version (i.e. without leads
and lags) of the dynamic `model` declared in the model block. To write
the LaTeX contents of the `steady_state_model` see
`write_latex_steady_state_model`{.interpreted-text role="comm"}.

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.

*Options*

- *write_equation_tags*: 

See `write_equation_tags`{.interpreted-text role="opt"}.

- *Command*: `write_latex_steady_state_model`

This command creates two LaTeX files: one containing the steady state
model and one containing the LaTeX document header information.

If your `.mod` file is `FILENAME.mod`, then Dynare will create a file
called `FILENAME/latex/steady_state.tex`, which includes a file called
`FILENAME/latex/steady_state_content.tex` (also created by Dynare)
containing the list of all the steady state model equations.

If LaTeX names were given for variables and parameters (see
`var-decl`{.interpreted-text role="ref"}), then those will be used;
otherwise, the plain text names will be used.

Note that the model written in the `.tex` file will differ from the
model declared by the user in some dimensions (see
`write_latex_dynamic_model`{.interpreted-text role="comm"} for details).

For the required LaTeX packages, see
`write_latex_original_model`{.interpreted-text role="comm"}.

## Auxiliary variables

The model which is solved internally by Dynare is not exactly the model
declared by the user. In some cases, Dynare will introduce auxiliary
endogenous variables---along with corresponding auxiliary
equations---which will appear in the final output.

The main transformation concerns leads and lags. Dynare will perform a
transformation of the model so that there is only one lead and one lag
on endogenous variables and no leads/lags on exogenous variables.

This transformation is achieved by the creation of auxiliary variables
and corresponding equations. For example, if `x(+2)` exists in the
model, Dynare will create one auxiliary variable
`AUX_ENDO_LEAD = x(+1)`, and replace `x(+2)` by `AUX_ENDO_LEAD(+1)`.

A similar transformation is done for lags greater than 2 on endogenous
(auxiliary variables will have a name beginning with `AUX_ENDO_LAG`),
and for exogenous with leads and lags (auxiliary variables will have a
name beginning with `AUX_EXO_LEAD` or `AUX_EXO_LAG` respectively).

Another transformation is done for the `EXPECTATION` operator. For each
occurrence of this operator, Dynare creates an auxiliary variable
defined by a new equation, and replaces the expectation operator by a
reference to the new auxiliary variable. For example, the expression
`EXPECTATION(-1)(x(+1))` is replaced by `AUX_EXPECT_LAG_1(-1)`, and the
new auxiliary variable is declared as `AUX_EXPECT_LAG_1 = x(+2)`.

Auxiliary variables are also introduced by the preprocessor for the
`ramsey_model` and `ramsey_policy` commands. In this case, they are used
to represent the Lagrange multipliers when first order conditions of the
Ramsey problem are computed. The new variables take the form `MULT_i`,
where *i* represents the constraint with which the multiplier is
associated (counted from the order of declaration in the model block).

Auxiliary variables are also introduced by the
`differentiate_forward_vars` option of the model block. The new
variables take the form `AUX_DIFF_FWRD_i`, and are equal to `x-x(-1)`
for some endogenous variable `x`.

Finally, auxiliary variables will arise in the context of employing the
`diff`-operator.

Once created, all auxiliary variables are included in the set of
endogenous variables. The output of decision rules (see below) is such
that auxiliary variable names are replaced by the original variables
they refer to.

The number of endogenous variables before the creation of auxiliary
variables is stored in `M_.orig_endo_nbr`, and the number of endogenous
variables after the creation of auxiliary variables is stored in
`M_.endo_nbr`.

See [https://git.dynare.org/Dynare/dynare/-/wikis/Auxiliary-variables](https://git.dynare.org/Dynare/dynare/-/wikis/Auxiliary-variables)
for more technical details on auxiliary variables.
