using Dynare

using BipartiteMatching
using Graphs
using SparseArrays

function get_dynamic_incidence_matrix(context)
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
    return Incidence
end 

function get_incidence_bitmatrix_current_forward(context)
    model = context.models[1]
    incidence = get_dynamic_incidence_matrix(context)
    #Model normalization
    n = model.endogenous_nbr
    U = BitMatrix(incidence[:,n+1:2*n] .| incidence[:,2n+1:3*n]) #bipartite graph
    return U
end 

function get_maximum_cardinality_matching(context)
    U = get_incidence_bitmatrix_current_forward(context)
    matching, matched = findmaxcardinalitybipartitematching(U) #maximum cardinality matching of the graph
    !all(matched) && error("Model can't be normalized")
    return matching
end

function get_recursive_blocks(context, matching, U)
    #Reorder columns of incidence matrix
    iorder = [p[2] for p in sort(collect(pairs(matching)), by=x -> x[1])]

    # Strongly connected components
    g = SimpleDiGraph(U[:, iorder]) 
    scc = strongly_connected_components(g)
    return scc
end

function contains(ex::Union{Expr, Symbol}, subex::Union{Expr, Symbol})
    ex == subex && return true
    for a in ex.args
        typeof(a) == Expr && contains(a, subex) && return true
        typeof(a) == Symbol && a == subex && return true
    end
    # expression without args
    return false
end  

function contains_forwardvariable(ex::Union{Expr,Symbol}, context)
    forwardoffset = 2*context.models[1].endogenous_nbr
    typeof(ex) == Symbol && return false
    ex.head == :ref && ex.args[1] == :y && ex.args[2] > forwardoffset && return true
    for a in ex.args
        if typeof(a) == Expr
            contains_forwardvariable(a, context)
        end 
    end
    return  false
end 

function find_inbounds(f::F) where F <: Function
    for (i, a) in enumerate(f.body.args)
        typeof(a) == Expr && a.head == :macrocall && a.args[1] == Symbol("@inbounds") && return i
    end
    return nothing
end

function analyze_SparseDynamicResid(eq_nbr, context, matching)
    endo_nbr = context.models[1].endogenous_nbr
    f = Dynare.DFunctions.SparseDynamicResid!

    eq_offset = find_inbounds(f)

    eq = f.body.args[eq_offset].args[3].args[2*eq_nbr]
    @show eq.args[2]
    @show contains_forwardvariable(eq.args[2], context)
    e = @eval :(y[$(endo_nbr + matching[eq_nbr])])
    @show eq.args[2].args[2]
    @show e
    @show contains(eq.args[2].args[2], e)
end
    

context = @dynare "models/example3/example3" "notmpterms"

matching = get_maximum_cardinality_matching(context)
U = get_incidence_bitmatrix_current_forward(context)

rb = get_recursive_blocks(context, matching, U)

@show matching
@show rb

f = Dynare.DFunctions.SparseDynamicResid!

analyze_SparseDynamicResid(5, context, matching)

#=
m = context.models[1]
results = context.results.model_results[1]
trends = results.trends

dynamic_resid! = Dynare.DFunctions.dynamic_resid!

ws = Dynare.DynamicWs(context)

steadystate =  trends.endogenous_steady_state 
endogenous = steadystate
exogenous = trends.exogenous_steady_state
params = context.work.params

f = Dynare.DFunctions.SparseDynamicResid!

using Symbolics

@variables y[1:14], T[1:10], params[1:7], y[1:14], x[1:2], residuals[1:7]

ex = f.body.args[13].args[3].args[10].args[2]

Symbolics.solve_for(:(ex.args[2] ~ ex.args[3], :y[11], check=true))
=#
