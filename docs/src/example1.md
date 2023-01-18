# Internals: example1.mod

This section documents the inner working of DynareJulia when running
[example1.mod](../../test/models/example1/example1.mod).

## 1. Running Dynare
The ``example1.mod`` example is run with the following command:

```
context = @dynare "./test/models/example1/example1";
```
Please, run the command now: it will create the files discussed below.

## 2. Setting things up
The ``@dynare`` macro is defined in [Dynare.jl](../../src/Dynare.jl)
It mainly performs two functions:

### a. Calling  the Dynare preprocessor
The preprocessor creates the following files in
``./test/models/example1/example1/model/julia``

- [DynamicParamsDerivs.jl](../../test/models/example1/example1/model/julia/DynamicParamsDerivs.jl)
- [SparseDynamicG1!.jl](../../test/models/example1/example1/model/julia/SparseDynamicG1!.jl)
- [SparseDynamicG1TT!.jl](../../test/models/example1/example1/model/julia/SparseDynamicG1TT!.jl)
- [SparseDynamicResid!.jl](../../test/models/example1/example1/model/julia/SparseDynamicResid!.jl)
- [SparseDynamicResidTT!.jl](../../test/models/example1/example1/model/julia/SparseDynamicResidTT!.jl)
- [SparseStaticG1!.jl](../../test/models/example1/example1/model/julia/SparseStaticG1!.jl)
- [SparseStaticG1TT!.jl](../../test/models/example1/example1/model/julia/SparseStaticG1TT!.jl)
- [SparseStaticResid!.jl](../../test/models/example1/example1/model/julia/SparseStaticResid!.jl)
- [SparseStaticResidTT.jl](../../test/models/example1/example1/model/julia/SparseStaticResidTT!.jl)
- [StaticParamsDerivs.jl](../../test/models/example1/example1/model/julia/StaticParamsDerivs.jl)
- [SteadyState2.jl](../../test/models/example1/example1/model/julia/SteadyState2.jl)

that contains various version of the equations of the model as Julia
functions.

The preprocessor also  creates the directory ``./test/models/example1/example1`` with the
following structure:

```
graphs
model
  json
  julia
output
```
The declarations, initializations and computing tasks contained in
[example1.mod](../../test/models/example1/example1.mod) is translated
in a JSON file [modfile.json](../../test/models/example1/example1/model/json/modfile.json)
that is stored in the ``./test/models/example1/example1/model/json``
directory.

The Julia functions representing dynamic and static versions of the
equations of model and their derivatives are store in directory the
``./test/models/example1/example1/model/julia``:

 - [SparseDynamicG1!.jl](../../test/models/example1/example1/model/julia/SparseDynamicG1!.jl)
 - [SparseDynamicG1TT!.jl](../../test/models/example1/example1/model/julia/SparseDynamicG1TT!.jl)
 - [SparseDynamicResid!.jl](../../test/models/example1/example1/model/julia/SparseDynamicResid!.jl)
 - [SparseDynamicResidTT!.jl](../../test/models/example1/example1/model/julia/SparseDynamicResidTT!.jl)
 - [SparseStaticG1!.jl](../../test/models/example1/example1/model/julia/SparseStaticG1!.jl)
 - [SparseStaticG1TT!.jl](../../test/models/example1/example1/model/julia/SparseStaticG1TT!.jl)
 - [SparseStaticResid!.jl](../../test/models/example1/example1/model/julia/SparseStaticResid!.jl)
 - [SparseStaticResidTT!.jl](../../test/models/example1/example1/model/julia/SparseStaticResidTT!.jl)

The functions containing ```TT`` in their name compute temporary terms
that are factorized. The functions containing ``Resid`` in their name
compute the residuals of the equations when the value of the
endogenous variables are not solution of the model. The functions
containing ``G1`` in their name compute the first order derivatives
with respect to endogenous and exogenous variables.


### b. Calling Dynare parser

The Dynare parser is defined in
[DynareParser.jl](../../src/DynareParser.jl). It performs the
following tasks:

 - Initialize the structure ``context`` that contains  model inputs and
  results stored in the following fields:
   - symboltable::SymbolTable
   - models::Vector{Model}
   - modfileinfo::ModFileInfo
   - results::Results
   - work::Work
 - Parse relevant information from
  [modfile.json](../../test/models/example1/example1/model/json/modfile.json):
   - ```
  var y, c, k, a, h, b;
  varexo e, u;
  parameters beta, rho, alpha, delta, theta, psi, tau;
  ```
  are stored in ``context.symboltable`` from fields ``endogenous``,
  ``exogenous`` and ``parameters`` in ``modfile.json``
  - parameters describing the model (mostly from field ``model_info``
    in ``modfile.json``) are stored in field ``models[1]`` (``models``
    is a vector to support future multiple models in the same
    ``*.mod``)
  - ``modfileinfo`` contains summary features of elements contained in
    the ``*.mod`` file.
  - ``results`` stores the results of the various computations
  - ``work`` stores temporary computations

The main loop on lines 225-302 of function ``parse_statements!`` in
[DynareParser.jl](../../src/DynareParser.jl) triggers actions for the
following statement:

 - ```
 shocks;
   var e; stderr 0.009;
   var u; stderr 0.009;
   var e, u = 0.009*0.009;
 end;
 ```
 is further parsed by ``shocks!`` in
 [initialization.jl](../../src/initialization.jl). The information is
 stored in ``context.models[1].Sigma_e``
 - ```
  check;
 ```
 is not currently implemented
 - ```
  stoch_simul(dr=cycle_reduction, order=1);
  ```
  calls the computation of first order approximation of the model with
  the cyclic reduction algorithm and of
  related results with function ``stoch_simul!`` in [perturbations.jl](../../src/perturbations.jl)

## 3. Solving stochastic models
The model is solved by local approximation at 1st order.

- ``stoch_simul!`` in [perturbations.jl](../../src/perturbations.jl)
  calls un cascade of functions in the same file:
  - ``stoch_simul_core!``
  - ``compute_stoch_simul!``
  - ``compute_first_order_solution!`` that in turns calls
    ``LinearRationalExpectations.first_order_solver!``

## 4. Package ``LinearRationalExpectations.jl`` solves

  $$E_t \{ A y_{t+1} + B y_t + C y_{t-1} + D u_t + e\} = 0$$ 

  by computing $G_y$ and $G_u$ such that

  $$y_t - \bar y= G_y (y_{t-1} - \bar y) + G_u u_t$$

  is solution of the above equation.
- $G_y$ is solution of the matrix polynomial equation 

  $$ AG_yG_y + B_G_y + C = 0$$
  
  that is computed in package ``PolynomialMatrixEquations.jl``
- $G_u$ is then solution of the linear system

  $$ $$G_u = -(A G_y + B)^{-1}Du_t$$
  
- in addition the package provides a function to compute the
unconditional variance of $y_$.

## 5. Package ``PolynomialMatrixEquations.jl``
- The package provides two algorithms to solve

  $$ AG_yG_y + B_G_y + C = 0$$
  
  - cyclic reduction
  - generalized Schur decomposition
  
- the package uses package ``FastLapackInterface`` for non--allocating
  Lapack functions  




