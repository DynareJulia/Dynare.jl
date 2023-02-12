using Dynare
using SparseArrays
using BipartiteMatching

context = @dynare "models/example1/example1.mod";

#Get Jacobian
model = context.models[1]
results = context.results.model_results[1]
ws = Dynare.DynamicWs(context)
params = context.work.params
steadystate = results.trends.endogenous_steady_state
endogenous = repeat(steadystate, 3)
exogenous = results.trends.exogenous_steady_state # exogenous = repeat(exo_steadystate', 3)
J = Dynare.get_dynamic_jacobian!(ws, params, endogenous, exogenous, steadystate, model, 200)

#Incidence matrix
using SparseArrays
T=Vector{Bool}(undef, length(J.rowval))
# all T elements must be true
fill!(T, true)
Incidence = SparseMatrixCSC(size(J, 1), size(J,2), J.colptr, J.rowval, T)

#Model normalization
n = model.endogenous_nbr
U = BitMatrix(Incidence[:,n+1:2*n] .| Incidence[:,2n+1:3*n]) #bipartite graph
matching, matched = findmaxcardinalitybipartitematching(U) #maximum cardinality matching of the graph
# (Dict(5 => 4, 4 => 3, 6 => 6, 2 => 2, 3 => 5, 1 => 1)
# Equation 1 determines variable 1
# Equation 2 determines variable 3
# Equation 3 determines variable 5
# Equation 4 determines variable 3
# Equation 5 determines variable 4
# Equation 6 determines variable 6
!all(matched) && error("Model can't be normalized")
