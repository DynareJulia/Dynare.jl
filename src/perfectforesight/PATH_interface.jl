using PATHSolver

function mcp_perfectforesight_core!(
    perfect_foresight_ws::PerfectForesightWs,
    context::Context,
    periods::Int64,
    guess_values::Vector{Float64},
    initialvalues::Vector{Float64},
    terminalvalues::Vector{Float64},
    dynamic_ws::DynamicWs,
)
    m = context.models[1]
    results = context.results.model_results[1]
    work = context.work
    residuals = zeros(periods * m.endogenous_nbr)
    dynamic_variables = dynamic_ws.dynamic_variables
    temp_vec = dynamic_ws.temporary_values
    steadystate = results.trends.endogenous_steady_state
    params = work.params
    JJ = perfect_foresight_ws.J
    lb = perfect_foresight_ws.lb
    ub = perfect_foresight_ws.ub
    permutationsR = perfect_foresight_ws.permutationsR
    permutationsJ = perfect_foresight_ws.permutationsJ
    mcp!(lb, ub, permutationsJ, m.mcps, context, periods)

    exogenous = perfect_foresight_ws.x

    ntmp =  sum(m.dynamic_tmp_nbr[1:2])
    ws_threaded = [
        Dynare.DynamicWs(
            m.endogenous_nbr,
            m.exogenous_nbr,
            ntmp,
            m.dynamic_g1_sparse_colptr,
            m.dynamic_g1_sparse_rowval,
        ) for i = 1:Threads.nthreads()
    ]
    nzval = similar(m.dynamic_g1_sparse_rowval, Float64)
    nzval1 = similar(nzval)

    function f!(residuals, y)
        @debug "$(now()): start f!"
        get_residuals!(
            residuals,
            y,
            initialvalues,
            terminalvalues,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
            periods,
            temp_vec,
            permutations = permutationsR,
        )
        @debug "$(now()): end f!"
    end

    function JA!(
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        y::AbstractVector{Float64},
    )
        @debug "$(now()): start J!"
        A = updateJacobian!(
            JJ,
            DFunctions.dynamic_derivatives!,
            y,
            initialvalues,
            terminalvalues,
            dynamic_variables,
            exogenous,
            periods,
            temp_vec,
            params,
            steadystate,
            m.dynamic_g1_sparse_colptr,
            nzval,
            m.endogenous_nbr,
            m.exogenous_nbr,
            perfect_foresight_ws.permutationsJ,
            nzval1
        )
        
        @debug count(!iszero, A) / prod(size(A))
        @debug "$(now()): end J!"
        return A
    end

    function fj!(residuals, JJ, y)
        f!(residuals, y)
        JA!(JJ, y)
    end

    @debug "$(now()): start makeJacobian"
    
    @debug "$(now()): end makeJacobian"
    @debug "$(now()): start f!"
    f!(residuals, guess_values)
    @debug "$(now()): end f!"
    @debug "$(now()): start J!"
    
    function J!(y::AbstractVecOrMat{Float64})
        JA!(JJ, y)
        return JA!(JJ, y)
    end

    @debug "$(now()): end J!"
    df = OnceDifferentiable(f!, J!, vec(guess_values), residuals, JJ)
    @debug "$(now()): start nlsolve"

    A = J!(guess_values)
    (status, results, info) = solve_path!(f!, J!, lb, ub, vec(guess_values))
    @debug "$(now()): end nlsolve"
    endogenous_names = get_endogenous_longname(context.symboltable)
    push!(
        context.results.model_results[1].simulations,
        Simulation(
            "Sim1",
            "",
            TimeDataFrame(
                DataFrame(
                    transpose(reshape(results, m.endogenous_nbr, periods)),
                    endogenous_names,
                ),
                UndatedDate(1),
            ),
        ),
    )
end

function make_pf1_residuals(
            initialvalues::AbstractVector{T},
            terminalvalues::AbstractVector{T},
            exogenous::AbstractVector{T},
            dynamic_variables::AbstractVector{T},
            steadystate::AbstractVector{T},
            params::AbstractVector{T},
            m::Model,
            periods::Int,
            temp_vec::AbstractVector{T},
            permutations::Vector{Tuple{Int64,Int64}},
            residuals::AbstractVector{Float64},
        ) where T <: Real
    function f!(y::AbstractVector{T})
        get_residuals!(
            residuals,
            vec(y),
            initialvalues,
            terminalvalues,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
            periods,
            temp_vec,
            permutations = permutations
        )
        return residuals
    end
    return f!
end

function make_pf1_jacobian(
    dynamic_derivatives!::Function,
    initialvalues::AbstractVector{T},
    terminalvalues::AbstractVector{T},
    dynamic_variables::AbstractVector{T},
    exogenous::AbstractVector{T},
    periods::N,
    temp_vec::AbstractVector{T},
    params::AbstractVector{T},
    steadystate::AbstractVector{T},
    dynamic_g1_sparse_colptr::AbstractVector{N},
    nzval::AbstractVector{T},
    endogenous_nbr::N,
    exogenous_nbr::N,
    permutations::Vector{Tuple{Int64,Int64}},
    nzval1::AbstractVector{T},
    A::SparseArrays.SparseMatrixCSC{Float64,Int64},
) where {T <: Real, N <: Integer}
    function J!(
        y::AbstractVector{Float64},
    )
        updateJacobian!(
            A,
            dynamic_derivatives!,
            y,
            initialvalues,
            terminalvalues,
            dynamic_variables,
            exogenous,
            periods,
            temp_vec,
            params,
            steadystate,
            dynamic_g1_sparse_colptr,
            nzval,
            endogenous_nbr,
            exogenous_nbr,
            permutations,
            nzval1
        )
        return A
    end
end

# Using PATHSolver
function solve_path!(f!, J!, JJ, lb, ub, initial_values, nnz; kwargs...)
    n = length(initial_values)
    @assert n == length(lb) == length(ub)

    function function_callback(n::Cint, z::Vector{Cdouble}, F_val::Vector{Cdouble})
        @assert n == length(z) == length(F_val)
        f!(F_val, z)
        return Cint(0)
    end

    function jacobian_callback(
        n::Cint,
        nnz::Cint,
        z::Vector{Cdouble},
        col_start::Vector{Cint},
        col_len::Vector{Cint},
        row::Vector{Cint},
        data::Vector{Cdouble},
    )
        @assert n == length(z) == length(col_start) == length(col_len)
        @assert nnz == length(row) == length(data)

        J!(JJ, z)  # SparseMatrixCSC

        # Pouring Jac::SparseMatrixCSC to the format used in PATH
        copyto!(col_start, 1, JJ.colptr, 1, n)
        copyto!(row, JJ.rowval)
        copyto!(data, JJ.nzval)
           
        for i = 1:n
            col_len[i] = JJ.colptr[i+1] - JJ.colptr[i]
        end

        return Cint(0)
    end

    # var_name, F_name in RawIndex
    var_name = Vector{String}(undef, 0)
    F_name = Vector{String}(undef, 0)

    nnz = SparseArrays.nnz(J(initial_values))


    F_val = zeros(n)
    # Solve the MCP using PATHSolver
    status, z, info = PATHSolver.solve_mcp(
        function_callback,
        jacobian_callback,
        lb,
        ub,
        initial_values;
        nnz = nnz,
        variable_names = var_name,
        constraint_names = F_name,
        kwargs...,
    )

    return Int(status), z, info
end


function ramsey_constraints!(context, field)
    for c in field["ramsey_model_constraints"]
        p1, p2, p3 = split(c["constraint"], limit = 3)
        eq = context.symboltable[p1].orderintype
        push!(context.models[1].mcps, (eq, eq, p2, p3))
    end
end

function get_mcps!(mcps, model)
    for (i, eq) in enumerate(model)
        tags = get(eq, "tags", "")
        if !isempty(tags)
            tag = get(eq["tags"], "mcp", "")
            if !isempty(tag)
                p1, p2, p3 = split(tag, limit = 3)
                iv = context.symboltable[p1].orderintype
                push!(mcps, (i, iv, p2, p3))
            end
        end
    end
end

function mcp!(lb, ub, permutations, mcps, context, periods)
    n = context.models[1].endogenous_nbr
    resize!(lb, periods * n)
    resize!(ub, periods * n)
    fill!(lb, -Inf)
    fill!(ub, Inf)
    for m in mcps
        (eqn, var, op, expr) = m
        boundary = Dynare.dynare_parse_eval(String(expr), context)
        if op[1] == '<'
            for i = var:n:n*periods
                ub[i] = boundary
            end
        elseif op[1] == '>'
            for i = var:n:n*periods
                lb[i] = boundary
            end
        else
            error("MCP operator must be <, <=, > or >=")
        end
#        if var != eqn
#            push!(permutations, (var, eqn))
#        end
    end
end
