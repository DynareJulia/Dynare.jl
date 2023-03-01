# Computing recursive blocks

This notes present the identification of recursive blocks and the
construction of function to compute residuals and Jacobian by blocks

## Identification of recrusive blocks

Assuming that we have the sparse representation of the Jacobian matrix
of the model where the variables are in columns and the equations in rows

### Getting the static and dynamic Jacobian matrices
```
using Dynare

context = @dynare "models/example1/example1.mod";

model = context.models[1]
results = context.results.model_results[1]
wsd = Dynare.DynamicWs(context)
wss = Dynare.StaticWs(context)
params = context.work.params
steadystate = results.trends.endogenous_steady_state
endogenous = repeat(steadystate, 3)
exogenous = results.trends.exogenous_steady_state 
Jdynamic = Dynare.get_dynamic_jacobian!(wsd, params, endogenous, exogenous, steadystate, model, 200)
Jstatic = Dynare.get_static_jacobian!(wss, params, steadystate, exogenous, model) 
```

### The incidence matrix

The incidence matrix of the model is a boolean matrix indicating which
variable appears in which equation. 

#### Incidence matrix in the static model
```
using SparseArrays

ONES = Vector{Bool}(undef, length(Jstatic.rowval)) # all ONES elements must be true
fill!(ONES, true)
IncidenceStatic = SparseMatrixCSC(size(Jstatic, 1), size(Jstatic,2), Jstatic.colptr, Jstatic.rowval, ONES)
```
The incidence matrix of the static model is a square matrix and its
order corresponds to the number of endogenous variables.

#### Incidence matrix in the dynamic model
```
using SparseArrays

ONES = Vector{Bool}(undef, length(Jdynamic.rowval)) # all ONES elements must be true
fill!(ONES, true)
IncidenceDynamic = SparseMatrixCSC(size(Jdynamic, 1), size(Jdynamic,2), Jdynamic.colptr, Jdynamic.rowval, ONES)
```
The incidence matrix of the dynamic model is a rectangular matrix
$n\times 3n$ as variables in period $t-1$, $t$, and $t+1$ are treated
as different variables. Depending on how we are going to use the
blocks, we can consider different kind of unknowns.

When we treat variables identically at all leads and lags, we compute
the incidence matrix in the static model. This is appropriate when we
want to compute dynamic simulations block by block.

If we want to compute solution functions block by block, we consider
that variables in previous periods are predetermined and we consider
only variables appearing in an equation at the current or future periods.


### Normalization of the model

We assign each variable to one and only one equation that "determines"
this variable. Several assignements are possible. For our purpose any
one is fine. If no normalization is available, that means that the
model is ill specified, the Jacobian is singular, and two or more
equations are linearly dependent.

Normalization is obtained by computing the maximum cardinality
matching of the graph of the model

```
using BipartiteMatching
using Graphs

n = model.endogenous_nbr

U1 = BitMatrix(IncidenceStatic) #bipartite graph
matching1, matched1 = findmaxcardinalitybipartitematching(U1) #maximum cardinality matching of the graph
!all(matched1) && error("Model can't be normalized")

U2 = BitMatrix(IncidenceDynamic[:,n+1:2*n] .| IncidenceDynamic[:,2n+1:3*n]) #bipartite graph
matching2, matched2 = findmaxcardinalitybipartitematching(U2) #maximum cardinality matching of the graph
!all(matched2) && error("Model can't be normalized")
```



### Recursive blocks

Once we have the normalization of the model, we can assign each
equation to a different variable and obtain a simple graph describing
which variable enters directly in the determination of which variable.
Equations must be in the same order as variables. As the purpose of
this procedure is to partition equations, we keep equations in the
original order and reorder the variables (columns) of the incidence matrix.

```
#Reorder columns of incidence matrix in the static model
iorder1 = [p[2] for p in sort(collect(pairs(matching1)), by=x -> x[1])]

#Reorder columns of incidence matrix in the dynamic model
iorder2 = [p[2] for p in sort(collect(pairs(matching2)), by=x -> x[1])]
```

The strongly connected components of this simple graph provide the
subsets of variables that need to be evaluated or solved for
simultaneously. These are the recursive blocks.

```
# Strongly connected components in the static model
g1 = SimpleDiGraph(U1[:, iorder1]) 
sccStatic = strongly_connected_components(g1)

# Strongly connected components in the dynamic model
g2 = SimpleDiGraph(U2[:, iorder2]) 
sccDynamic = strongly_connected_components(g2)
```


## Writing block functions

The general idea is to take advantage of the already parsed version of
the functions for the model a a whole

### Functions as expressions

```
e = :(function f(x); x=1; return x; end) #:(function f(x)
```
Expression `e` is made of two fields: `head` and `args`:
```
e.head #:function
e.args #2-element Vector{Any}:

```
- `e.head` indicates that expression `e` represents a function
- `e.args[1]` is the calling sequence
- `e.args[2]` contains the body of the function
```
e.args[2].head #:block
e.args[2].args #5-element Vector{Any}:
```
- `e.args[2].args[[1, 2,4]]` indicates where the lines are coming from
- `e.args[2].args[3]` is `:(x[1] = 1]`
- `e.args[2].args[5]` is `:(return x)`

and 

```
e.args[2].args[3].head #:(=)
e.args[2].args[3].args #2-element Vector{Any}:
e.args[2].args[3].args[1].head #ERROR: type Symbol has no field head
e.args[2].args[3].args[1].args #ERROR: type Symbol has no field args
```
Observe that the element of `x` that is modified by the function is
indicated in `e.args[2].args[3].args[1].args[2]`

### The SparseDynamicResid! function as an expression

```
julia> Dynare.DFunctions.SparseDynamicResid!
RuntimeGeneratedFunction(#=in Dynare.DFunctions=#, #=using Dynare.DFunctions=#, :((T, residual, y, x, params, steady_state)->begin
          #= none:1 =#
          #= none:2 =#
          #= none:2 =# @assert length(T) >= 5
          #= none:3 =#
          #= none:3 =# @assert length(residual) == 6
          #= none:4 =#
          #= none:4 =# @assert length(y) == 18
          #= none:5 =#
          #= none:5 =# @assert length(x) == 2
          #= none:6 =#
          #= none:6 =# @assert length(params) == 7
          #= none:7 =#
          #= none:7 =# @inbounds begin
                  #= none:8 =#
                  residual[1] = y[8] * params[5] * T[1] - (1 - params[3]) * y[7]
                  #= none:9 =#
                  residual[2] = y[9] - params[1] * T[2] * (params[3] * exp(y[18]) * y[13] + y[9] * (1 - params[4]))
                  #= none:10 =#
                  residual[3] = y[7] - T[5]
                  #= none:11 =#
                  residual[4] = y[9] - (exp(y[12]) * (y[7] - y[8]) + (1 - params[4]) * y[3])
                  #= none:12 =#
                  residual[5] = y[10] - (params[2] * y[4] + params[7] * y[6] + x[1])
                  #= none:13 =#
                  residual[6] = y[12] - (y[4] * params[7] + params[2] * y[6] + x[2])
              end
          #= none:15 =#
          return nothing
      end))
```
Using the same logic as above, we get the residual of the first
equation of the model as
```
julia> Dynare.DFunctions.SparseDynamicResid!.body.args[13].args[3].args[2]
:(residual[1] = y[8] * params[5] * T[1] - (1 - params[3]) * y[7])
```
and the number of the first equation is revealed by
```
julia> Dynare.DFunctions.SparseDynamicResid!.body.args[13].args[3].args[2].args[1].args[2]
1
```
It is possible to exploit this technique to dispatch the computation
of the residuals belonging to each block in different functions
