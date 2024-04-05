#module Blocks
using BipartiteMatching: findmaxcardinalitybipartitematching
using Graphs: SimpleDiGraph, strongly_connected_components
using RuntimeGeneratedFunctions
using SparseArrays

RuntimeGeneratedFunctions.init(@__MODULE__)

struct BlockIndices_1
    equation_pointers::Vector{Int}
    variable_pointers::Vector{Int}
    expressions::Vector{Expr}
end

BlockIndices = BlockIndices_1

struct Block_1{F1 <: Function, F2 <: Function, F3 <: Function}
    assignment::Bool
    forward::Bool
    jacobian::SparseMatrixCSC{Float64, Int}
    assigment_fcn::F1
    jacobian_fcn::F2
    residual_fcn::F3
    indices::BlockIndices
end

abstract type AbstractBlock end
abstract type AbstractPreambleBlock <: AbstractBlock end

struct AssignmentBlock_1 <: AbstractPreambleBlock
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
end

struct PreambleBlock_1 <: AbstractPreambleBlock
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
end

struct ForwardBlock_1 <: AbstractBlock
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
end

struct BackwardBlock_1 <: AbstractBlock    
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
end

AssignmentBlock = AssignmentBlock_1
PreambleBlock = PreambleBlock_1
ForwardBlock = ForwardBlock_1
BackwardBlock = BackwardBlock_1


function get_dynamic_incidence_matrix(context)
    #Get Jacobian
    model = context.models[1]
    results = context.results.model_results[1]
    ws = DynamicWs(context)
    params = context.work.params
    steadystate = results.trends.endogenous_steady_state
    endogenous = repeat(steadystate, 3)
    exogenous = results.trends.exogenous_steady_state 
    J = get_dynamic_jacobian!(ws, params, endogenous, exogenous, steadystate, model, 200)

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

function get_maximum_cardinality_matching(context, U)
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

function contains_forwardvariable(ex::Union{Expr,Symbol}, forwardoffset)
    ex.head == :ref && ex.args[1] == :y && ex.args[2] > forwardoffset && return true
    for a in ex.args
        if typeof(a) == Expr
            contains_forwardvariable(a, forwardoffset) && return true
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

function is_block_normalized(b, matching, endogenous_nbr, eq_offset)
    f = DFunctions.SparseDynamicResid!
    for eq_no in b
        eq = f.body.args[eq_offset].args[3].args[2*eq_no]
        k = endogenous_nbr + matching[eq_no]
        e = Expr(:ref, :y, k) 
        is_normalized =  eq.args[2].args[2].args[2] == k && !contains(eq.args[2].args[3], e)
        !is_normalized && return false
    end
    return true
end

function is_block_forward_looking(b, equation_xref_table, endogenous_nbr)
    for eq_no in b
        for v in equation_xref_table[eq_no]
            v > 2*endogenous_nbr && return false
        end
    end
    return true
end

function is_block_linear(b, context)
    index = context.models[1].dynamic_g2_sparse_indices
    imax = maximum(b)
    for i in index
        i[1] in b && return false
        i[1] > imax && break
    end
    return true
end

function add_block_to_preamble!(preamble_expressions, b, eq_offset)
    f = DFunctions.SparseDynamicResid!
    for eq_no in b
        eq = f.body.args[eq_offset].args[3].args[2*eq_no]
        push!(preamble_expressions, Expr(:(=), eq.args[2].args[2], eq.args[2].args[3]))
    end
end

function  nothing_fcn()
end

function analyze_SparseDynamicResid!(blocks, forward_expressions, preamble_expressions, system_expressions, blocks_eqs, matching, context)
    f = DFunctions.SparseDynamicResid!
    endogenous_nbr = context.models[1].endogenous_nbr

    eq_offset = find_inbounds(f)

    predetermined_variables = Vector{Int}(undef, 0)
    preamble_eqs = similar(predetermined_variables)
    system_equations = similar(predetermined_variables)
    forward_expressions_eqs = similar(predetermined_variables)
    system_expressions_eqs = similar(predetermined_variables)
    k1 = k2 = 1
    preamble = true
    for b in blocks_eqs
        if is_block_normalized(b, matching, endogenous_nbr, eq_offset) &&
            is_block_linear(b, context) && preamble
            add_block_to_preamble!(preamble_expressions, b, eq_offset)
            for eq in b
                # predetermined variables at the current period
                push!(predetermined_variables, matching[eq] + endogenous_nbr)
                push!(preamble_eqs, eq)
            end
        else
            preamble = false
            for eq_no in b
                push!(system_equations, eq_no)
                eq = f.body.args[eq_offset].args[3].args[2*eq_no]
                if contains_forwardvariable(eq.args[2], 2*endogenous_nbr)
                    # residuals[k1] = eq.args[2]
                    push!(forward_expressions, Expr(:(=), :(residuals[$k1]), eq.args[2]))
                    k1 += 1
                    push!(forward_expressions_eqs, eq_no)
                else
                    # residuals[k2] = eq.args[2]
                    push!(system_expressions, Expr(:(=), :(residuals[$k2]), eq.args[2]))
                    k2 += 1
                    push!(system_expressions_eqs, eq_no)
                end
            end
        end
    end
    return (sort!(predetermined_variables), sort!(system_equations), preamble_eqs, forward_expressions_eqs, system_expressions_eqs)
end

function make_jacobian_submatrix(colptr::AbstractVector{I},
                               rowval::AbstractVector{I}, colval::AbstractVector{I},
                               rows::AbstractVector{I}, cols::AbstractVector{I}) where I <: Integer
    scolptr = Vector{Int}(undef, length(cols) + 1)
    srowval = Vector{Int}(undef, 0)
    f = Dynare.DFunctions.SparseDynamicG1!
    eq_offset = Dynare.find_inbounds(f)
    feqs = f.body.args[eq_offset].args[3]
    deriv_expressions = Vector{Expr}(undef, 0)
    scolptr[1] = 1
    oldptr = 1
    kr = 1
    kc = 1
    for (i, r) in enumerate(rowval)
        row_small = findfirst(r .== rows)
        col_small = findfirst(colval[i] .== cols)
        if !isnothing(row_small) && !isnothing(col_small)
            push!(srowval, row_small)
            # renumber output vector elements
            e = feqs.args[2*i]
            e.args[1].args[2] = kr
            push!(deriv_expressions, feqs.args[2*i])
            if col_small != kc
                while col_small > kc + 1
                    kc += 1
                    scolptr[kc] = oldptr
                end
                if col_small == kc + 1
                    kc += 1
                    scolptr[kc] = kr
                    oldptr = kr
                end 
            end
            kr += 1
        end
    end
    while kc <= length(cols)
        kc += 1
        scolptr[kc] = kr
    end
    snzval = similar(srowval, Float64)

    return (SparseMatrixCSC(length(rows), length(cols), scolptr, srowval, snzval),
            deriv_expressions)
end

function make_assignment_jacobian(colptr::AbstractVector{I},
                                    rowval::AbstractVector{I}, colval::AbstractVector{I},
                                    rows::AbstractVector{I}) where I <: Integer
    scolptr = Vector{Int}(undef, length(colptr))
    srowval = Vector{Int}(undef, 0)
    f = Dynare.DFunctions.SparseDynamicG1!
    eq_offset = Dynare.find_inbounds(f)
    feqs = f.body.args[eq_offset].args[3]
    deriv_expressions = Vector{Expr}(undef, 0)
    scolptr[1] = 1
    oldptr = 1
    kr = 1
    kc = 1
    for (i, r) in enumerate(rowval)
        row_small = findfirst(r .== rows)
        # in assignment block ignore derivative of assigned variable
        if !isnothing(row_small) && r + endogenous_nbr != kc
            push!(srowval, row_small)
            # renumber output vector elements
            e = feqs.args[2*i]
            e.args[1].args[2] = kr
            push!(deriv_expressions, feqs.args[2*i])
            while r == colptr[kc + 1]
                kc += 1
                scolptr[kc] = oldptr
            end
            if r < colptr[kc + 1]
                scolptr[kc] = kr
                oldptr = kr
            end 
            kr += 1
        end
    end
    while kc < length(colptr)
        kc += 1
        scolptr[kc] = kr
    end
    snzval = similar(srowval, Float64)

    return (SparseMatrixCSC(length(rows), length(colptr) - 1, scolptr, srowval, snzval),
            deriv_expressions)
end

function make_residual_jacobian(colptr::AbstractVector{I},
                                rowval::AbstractVector{I}, colval::AbstractVector{I},
                                rows::AbstractVector{I}) where I <: Integer
    scolptr = Vector{Int}(undef, length(colptr))
    srowval = Vector{Int}(undef, 0)
    f = Dynare.DFunctions.SparseDynamicG1!
    eq_offset = Dynare.find_inbounds(f)
    feqs = f.body.args[eq_offset].args[3]
    deriv_expressions = Vector{Expr}(undef, 0)
    scolptr[1] = 1
    oldptr = 1
    kr = 1
    kc = 1
    for (i, r) in enumerate(rowval)
        row_small = findfirst(r .== rows)
        # in assignment block ignore derivative of assigned variable
        if !isnothing(row_small)
            push!(srowval, row_small)
            # renumber output vector elements
            e = feqs.args[2*i]
            e.args[1].args[2] = kr
            push!(deriv_expressions, feqs.args[2*i])
            while r == colptr[kc + 1]
                kc += 1
                scolptr[kc] = oldptr
            end
            if r < colptr[kc + 1]
                scolptr[kc] = kr
                oldptr = kr
            end 
            kr += 1
        end
    end
    while kc < length(colptr)
        kc += 1
        scolptr[kc] = kr
    end
    snzval = similar(srowval, Float64)

    @show length(rows), length(colptr)
    @show scolptr
    @show length(scolptr)
    @show srowval
    @show snzval
    return (SparseMatrixCSC(length(rows), length(colptr) - 1, scolptr, srowval, snzval),
            deriv_expressions)
end

function make_assignment_function(fname, expressions)
    f_call = Expr(:call,
                  Symbol(fname),
                  Expr(:(::), Symbol(:T), Vector{Float64}),
                  Expr(:(::), Symbol(:y), Vector{Float64}),
                  Expr(:(::), Symbol(:x), Vector{Float64}),
                  Expr(:(::), Symbol(:params), Vector{Float64}),
                  Expr(:(::), Symbol(:steadystate), Vector{Float64})
                  )
    f_body = Expr(:block, expressions...)
    return @RuntimeGeneratedFunction(Expr(:function, f_call, f_body))
end

function make_system_function(fname, expressions)
    f_call = Expr(:call,
                  Symbol(fname),
                  Expr(:(::), Symbol(:T), Vector{Float64}),
                  Expr(:(::), Symbol(:residuals), Vector{Float64}),
                  Expr(:(::), Symbol(:y), Vector{Float64}),
                  Expr(:(::), Symbol(:x), Vector{Float64}),
                  Expr(:(::), Symbol(:params), Vector{Float64}),
                  Expr(:(::), Symbol(:steadystate), Vector{Float64})
                  )
    f_body = Expr(:block, expressions...)
    return @RuntimeGeneratedFunction(Expr(:function, f_call, f_body))
end

function make_system_jacobian(fname, expressions)
    f_call = Expr(:call,
                  Symbol(fname),
                  Expr(:(::), Symbol(:T), Vector{Float64}),
                  Expr(:(::), Symbol(:g1_v), Vector{Float64}),
                  Expr(:(::), Symbol(:y), Vector{Float64}),
                  Expr(:(::), Symbol(:x), Vector{Float64}),
                  Expr(:(::), Symbol(:params), Vector{Float64}),
                  Expr(:(::), Symbol(:steadystate), Vector{Float64})
                  )
    f_body = Expr(:block, expressions...)
    return @RuntimeGeneratedFunction(Expr(:function, f_call, f_body))
end

function get_state_variables(predetermined_variables,
                             system_equations,
                             equation_xref_list,
                             variable_xref_list,
                             endogenous_nbr)
    states_ = Set{Int}()
    for v in 1:endogenous_nbr
        for e in variable_xref_list[v]
            e in system_equations && push!(states_, v)
        end
    end
    for v in predetermined_variables
        for e in variable_xref_list[v]
            e in system_equations && push!(states_, v)
        end
    end
    return sort(collect(states_))
end

function make_block_functions(context)
    model = context.models[1]
    endogenous_nbr = model.endogenous_nbr
    U = get_incidence_bitmatrix_current_forward(context)
    matching = get_maximum_cardinality_matching(context, U)
    
    rb = get_recursive_blocks(context, matching, U)
    
    f = DFunctions.SparseDynamicResid!
    
    forward_expressions = Vector{Expr}(undef, 0)
    preamble_expressions = Vector{Expr}(undef, 0)
    other_expressions = Vector{Expr}(undef, 0)
    blocks = []
    
    (predetermined_variables, system_equations, preamble_eqs, forward_expressions_eqs, system_expressions_eqs) =
        analyze_SparseDynamicResid!(blocks,
                                    forward_expressions,
                                    preamble_expressions,
                                    other_expressions, 
                                    rb,
                                    matching,
                                    context)
                                                            
    forward_equations_nbr = length(forward_expressions)
    other_equations_nbr = length(other_expressions)
    
    global forward_block_ = make_system_function(:forward_block, forward_expressions)
    global preamble_block_ = make_assignment_function(:preamble_block, preamble_expressions)
    global other_block_ = make_system_function(:system_block, other_expressions)
#    global system_jacobian = make_system_jacobian(:system_jacobian, jacobian_expressions)
    
    equation_xref_list, variable_xref_list = xref_lists(context)
    states = get_state_variables(predetermined_variables,
                                 system_equations,
                                 equation_xref_list,
                                 variable_xref_list,
                                 endogenous_nbr)
    system_variables = [matching[e] + endogenous_nbr for e in system_equations]
    sort!(system_variables)

    #rows = union(forward_expressions_eqs, system_expressions_eqs)
    #make_sparse_submatrix(SM, colval, rows , system_variables)


    assignment = true
    forward = false
    jacobian = make_assignment_jacobian(context.models[1].dynamic_g1_sparse_colptr,
                                        context.models[1].dynamic_g1_sparse_colval,
                                        context.models[1].dynamic_g1_sparse_rowval,
                                        preamble_eqs)
    
    assignment_fcn = make_assignment_function(:preamble_block, preamble_expressions)
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    indices = BlockIndices(preamble_eqs, predetermined_variables .- endogenous_nbr, preamble_expressions)
    push!(blocks, Block(assignment,
                        forward,
                        jacobian[1],
                        assignment_fcn,
                        nothing_fcn,
                        nothing_fcn,
                        indices))
    @show blocks
    @show indices

    assignment = false
    forward = true
    jacobian = make_residual_jacobian(context.models[1].dynamic_g1_sparse_colptr,
                                        context.models[1].dynamic_g1_sparse_colval,
                                        context.models[1].dynamic_g1_sparse_rowval,
                                        forward_expressions_eqs)
    
    assignment_fcn = nothing_fcn
    jacobian_fnc = nothing_fcn
    residual_fcn = make_system_function(:forward_block, forward_expressions)
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    indices = BlockIndices(forward_expressions_eqs, [], forward_expressions)
    push!(blocks, Block(assignment,
                        forward,
                        jacobian[1],
                        assignment_fcn,
                        nothing_fcn,
                        nothing_fcn,
                        indices))
    @show blocks
    @show indices

    assignment = false
    forward = false
    jacobian = make_residual_jacobian(context.models[1].dynamic_g1_sparse_colptr,
                                        context.models[1].dynamic_g1_sparse_colval,
                                        context.models[1].dynamic_g1_sparse_rowval,
                                        system_expressions_eqs)
    
    assignment_fcn = nothing_fcn
    jacobian_fnc = nothing_fcn
    residual_fcn = make_system_function(:other_block, other_expressions)
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    indices = BlockIndices(system_expressions_eqs, [], other_expressions)
    push!(blocks, Block(assignment,
                        forward,
                        jacobian[1],
                        assignment_fcn,
                        nothing_fcn,
                        nothing_fcn,
                        indices))
    @show blocks
    @show indices
error()
    return (states, predetermined_variables, system_variables,
            forward_equations_nbr, other_equations_nbr, preamble_eqs,
            forward_expressions_eqs, system_expressions_eqs)
end

function forward_block!(T, residuals, y, x, params, steadystate)
    DFunctions.SparseDynamicResidTT!(T, y, x, params, steadystate)
    forward_block_(T, residuals, y, x, params, steadystate)
    return nothing
end

function preamble_block!(T, y, x, params, steadystate)
    DFunctions.SparseDynamicResidTT!(T, y, x, params, steadystate)
    preamble_block_(T, y, x, params, steadystate)
    return nothing
end

function other_block!(T, residuals, y, x, params, steadystate)
    DFunctions.SparseDynamicResidTT!(T, y, x, params, steadystate)
    other_block_(T, residuals, y, x, params, steadystate)
    return nothing
end

"""
    xref_lists(context::Context) -> (equation_xref_list, variable_xref_list)

computes equation and variable cross-reference lists

## Output
- equation_xref_list::Vector{Vector{Int}}: list of variables apearing in each equation
- variable_xref_list::Vector{Vector{Int}}:: list of equations where a give variable is appearing
"""
function xref_lists(context)
    trends = context.results.model_results[1].trends
    steadystate = trends.endogenous_steady_state
    endogenous = repeat(steadystate, 3)
    exogenous = trends.exogenous_steady_state
    model = context.models[1]
    params = context.work.params
    ws = DynamicWs(context) 
    jacobian = get_dynamic_jacobian!(
        ws,
        params,
        endogenous,
        exogenous,
        steadystate,
        model,
        2,
    )

    variable_xref_list = [Int[] for i in axes(jacobian, 2)]
    equation_xref_list = [Int[] for i in axes(jacobian, 1)]
    nc = length(endogenous) + length(exogenous)
    colptr = jacobian.colptr
    rowval = jacobian.rowval
    for ic = 1:nc
        for j in colptr[ic]:colptr[ic+1] - 1
            push!(variable_xref_list[ic], rowval[j]) 
            push!(equation_xref_list[rowval[j]], ic)
        end
    end
    return(equation_xref_list, variable_xref_list)
end


    
    
#end #module Blocks
