module TestBlocks

using Dynare
using SparseArrays
using BipartiteMatching
using Graphs
using RuntimeGeneratedFunctions

RuntimeGeneratedFunctions.init(@__MODULE__)

context = @dynare "models/example1/example1.mod" "notmpterms";

#Get Jacobian
model = context.models[1]
results = context.results.model_results[1]
params = context.work.params
steadystate = results.trends.endogenous_steady_state
endogenous = repeat(steadystate, 3)
exogenous = results.trends.exogenous_steady_state 


function get_strongly_connected_components(J::SparseMatrixCSC)
    #Incidence matrix
    T=Vector{Bool}(undef, length(J.rowval)) # all T elements must be true
    fill!(T, true)
    Incidence = SparseMatrixCSC(size(J, 1), size(J,2), J.colptr, J.rowval, T)

    #Model normalization
    n = model.endogenous_nbr
    U = BitMatrix(Incidence[:,n+1:2*n] .| Incidence[:,2n+1:3*n]) #bipartite graph
    matching, matched = findmaxcardinalitybipartitematching(U) #maximum cardinality matching of the graph
    !all(matched) && error("Model can't be normalized")

    #Reorder columns of incidence matrix
    iorder = [p[2] for p in sort(collect(pairs(matching)), by=x -> x[1])]

    # Strongly connected components
    g = SimpleDiGraph(U[:, iorder]) 
    scc = strongly_connected_components(g)
    return (matching, scc)
end

function make_block_function(equs, block_component)
    # make function template
    func =  :((T, residual, y, x, params, steady_state) -> @inbounds begin; end)
    # add equations belonging to this block
    for i = 2:2:length(equs)
        # in sdr each equations is the second of two elements
        eq = copy(equs[i])
        @show eq
        # change the index of the residual
        eq.args[1].args[2] = i
        push!(func.args[2].args[2].args[3].args, eq)
    end
    return @RuntimeGeneratedFunction(func)
end

function make_blocks(func0, block_components)
    block_functions = Function[]
    equs = func0.body.args[13].args[3].args
    @show equs
    for block in block_components
        push!(block_functions, make_block_function(equs, block))
    end
    return block_functions
end

wsd = Dynare.DynamicWs(context)
Jdynamic = Dynare.get_dynamic_jacobian!(wsd, params, endogenous, exogenous, steadystate, model, 200)
wss = Dynare.StaticWs(context)
Jstatic = Dynare.get_static_jacobian!(wss, params, steadystate, exogenous, model)

scc_dynamic = get_strongly_connected_components(Jdynamic) 
dynamic_residual_block_functions = make_blocks(Dynare.DFunctions.SparseDynamicResid!, scc_dynamic) 
dynamic_jacobian_block_functions = make_blocks(Dynare.DFunctions.SparseDynamicG1!, scc_dynamic) 

scc_static = get_strongly_connected_components(Jstatic) 
dynamic_residual_block_functions = make_blocks(Dynare.DFunctions.SparseStaticResid!, scc_dynamic) 
dynamic_jacobian_block_functions = make_blocks(Dynare.DFunctions.SparseDynamcG1!, scc_dynamic) 


end # end module          
