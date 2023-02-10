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
Incidence = SparseMatrixCSC(size(J, 1), size(J,2), J.colptr, J.rowval, T)

#Model normalization
for n in 1:length(steadystate)
    U = BitMatrix(hcat(Incidence[:,n+1:2*n], Incidence[:,2n+1:3*n])) #bipartite graph
    matching, matched = findmaxcardinalitybipartitematching(U) #maximum cardinality matching of the graph
end