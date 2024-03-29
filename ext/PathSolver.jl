module PathSolver

using Dates
using Dynare, PATHSolver
using LinearAlgebra
using SparseArrays

function Dynare.mcp_perfectforesight_core!(
    ::Dynare.PathNLS,
    perfect_foresight_ws::Dynare.PerfectForesightWs,
    context::Dynare.Context,
    periods::Int64,
    guess_values::Vector{Float64},
    initialvalues::Vector{Float64},
    terminalvalues::Vector{Float64},
    dynamic_ws::Dynare.DynamicWs;
    maxit = 50,
    tolf = 1e-5,
    tolx = 1e-5
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
        Dynare.get_residuals!(
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
        A = Dynare.updateJacobian!(
            JJ,
            Dynare.DFunctions.dynamic_derivatives!,
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

    @debug "$(now()): start nlsolve"
    show_iterations = (("JULIA_DEBUG" => "PathSolver") in ENV) ? "yes" : "no"
    (status, results, info) = solve_path!(f!,
        JA!,
        JJ,
        lb,
        ub,
        vec(guess_values),
        output_major_iterarions = show_iterations,
        output_minor_iterations = show_iterations,
        convergence_tolerance = tolf,
    )

    @debug "$(now()): end nlsolve"
    Dynare.make_simulation_results!(context::Context, results, exogenous, terminalvalues, periods)
end

function Dynare.mcp_perfectforesight_core!(
    ::Dynare.PathNLS,
    perfect_foresight_ws::Dynare.PerfectForesightWs,
    context::Context,
    periods::Int64,
    guess_values::Vector{Float64},
    initialvalues::Vector{Float64},
    terminalvalues::Vector{Float64},
    dynamic_ws::Dynare.DynamicWs,
    flipinfo::Dynare.FlipInformation,
    infoperiod;
    maxit = 50,
    tolf = 1e-5,
    tolx = 1e-5
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

    ix_stack = flipinfo.ix_stack
    function f!(residuals, y)
        @debug "$(now()): start f!"
        Dynare.flip!(y, exogenous, ix_stack)
        Dynare.get_residuals!(
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
        Dynare.flip!(y, exogenous, ix_stack)
        @debug "$(now()): end f!"
        return residuals
    end

    function JA!(
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        y::AbstractVector{Float64},
    )
        @debug "$(now()): start J!"
        A = Dynare.updateJacobian!(
            JJ,
            Dynare.DFunctions.dynamic_derivatives!,
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
            flipinfo,
            nzval1
        )
        
        @debug count(!iszero, A) / prod(size(A))
        @debug "$(now()): end J!"
        return A
    end

    Dynare.set_future_information!(guess_values, exogenous, context, periods, infoperiod)
    Dynare.flip!(guess_values, exogenous, flipinfo.ix_stack)

    @debug "$(now()): start nlsolve"
    show_iterations = (("JULIA_DEBUG" => "PathSolver") in ENV) ? "yes" : "no"
    (status, results, info) = solve_path!(f!,
        JA!,
        JJ,
        lb,
        ub,
        vec(guess_values),
        output_major_iterarions = show_iterations,
        output_minor_iterations = show_iterations,
        convergence_tolerance = tolf,
    )
    @debug "$(now()): end nlsolve"
    Dynare.make_simulation_results!(context::Context, results, exogenous, terminalvalues, periods)
end

function Dynare.mcp_sg_core!(::Dynare.PathNLS, grid, pol_guess, pol, gridZero, nPols, nodes, weights, exogenous, params, steadystate, forward_equations_nbr, 
                      endogenous_nbr, exogenous_nbr, state_variables, system_variables, y, bmcps, dynamicws, model,
                      preamble_eqs, predetermined_variables, preamble_lagged_variables, forward_system_variables,
                      forward_expressions_eqs, other_expressions_eqs, ids, aNumAdd, aPoints1, JJ, lb, ub)
    let state
        function f1!(fx, x)
            fx .= Dynare.sysOfEqs(x, y, exogenous, state, gridZero, nPols, nodes, weights,
                                  params, steadystate, forward_equations_nbr, endogenous_nbr,
                                  exogenous_nbr, state_variables, system_variables, bmcps)
        end
        
        function JA1!(J, x)
            Dynare.sysOfEqs_derivatives!(J, x, y, state, gridZero,
                                         nPols, exogenous, nodes, weights, params, steadystate,
                                         forward_equations_nbr, endogenous_nbr, exogenous_nbr,
                                         state_variables, system_variables, bmcps, dynamicws, model, preamble_eqs, 
                                         predetermined_variables, preamble_lagged_variables, forward_system_variables,
                                         forward_expressions_eqs, other_expressions_eqs, ids)                    
        end
    
        for ii1 in 1:aNumAdd
            @views state = aPoints1[:, ii1]
            @views pol .= pol_guess[:, ii1]
            (status, results, info) =     (status, results, info) = solve_path!(f1!, JA1!, JJ, lb, ub, pol)
            polInt[:, ii1] .= results
        end
    end
end

# Using PATHSolver
function solve_path!(f!, JA!, JJ::SparseMatrixCSC, lb, ub, initial_values; kwargs...)
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

        JA!(JJ, z)  # SparseMatrixCSC

        # Pouring JJ::SparseMatrixCSC to the format used in PATH
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

    nnz = SparseArrays.nnz(JJ)

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

function solve_path!(f!, JA!, JJ::Matrix, lb, ub, initial_values; kwargs...)
    n = length(initial_values)
    @assert n == length(lb) == length(ub)

    function function_callback(n::Cint, z::Vector{Cdouble}, F_val::Vector{Cdouble})
        @assert n == length(z) == length(F_val)
        f!(F_val, z)
        return Cint(0)
    end

    NNZ = 0
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

        JA!(JJ, z)  # Matrix

        i = 1
        for c in 1:n
            col_start[c] = i
            col_len[c] = n
            for r in 1:n
                row[i] = r
                i += 1
            end
        end
        copyto!(data, JJ)

        return Cint(0)
    end

    # var_name, F_name in RawIndex
    var_name = Vector{String}(undef, 0)
    F_name = Vector{String}(undef, 0)

    nnz = n*n

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

end #module
