using Dynare
using SparseArrays
using BipartiteMatching
using Graphs

context = @dynare "models/example1/example1.mod";

#Get Jacobian
model = context.models[1]
results = context.results.model_results[1]
ws = Dynare.DynamicWs(context)
params = context.work.params
steadystate = results.trends.endogenous_steady_state
endogenous = repeat(steadystate, 3)
exogenous = results.trends.exogenous_steady_state 
J = Dynare.get_dynamic_jacobian!(ws, params, endogenous, exogenous, steadystate, model, 200)

#Incidence matrix
T=Vector{Bool}(undef, length(J.rowval)) # all T elements must be true
fill!(T, true)
Incidence = SparseMatrixCSC(size(J, 1), size(J,2), J.colptr, J.rowval, T)

#Model normalization
n = model.endogenous_nbr
U = BitMatrix(Incidence[:,n+1:2*n] .| Incidence[:,2n+1:3*n]) #bipartite graph
matching, matched = findmaxcardinalitybipartitematching(U) #maximum cardinality matching of the graph
!all(matched) && error("Model can't be normalized")

# Strongly connected components
g = SimpleGraph(n) 
edge_list = collect(matching)
g = SimpleDiGraph(Edge.(edge_list))
scc = strongly_connected_components(g)

