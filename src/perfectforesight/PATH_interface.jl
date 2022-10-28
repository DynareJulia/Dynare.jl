using PATHSolver

function mcp_perfectforesight_core!(
    perfect_foresight_ws::PerfectForesightWs,
    context::Context,
    periods::Int64,
    guess_values::Matrix{Float64},
    initialvalues::Vector{Float64},
    dynamic_ws::DynamicWs,
)
    m = context.models[1]
#    ddf = context.dynarefunctions
    results = context.results.model_results[1]
    work = context.work
    residuals = zeros(periods * m.endogenous_nbr)
    dynamic_variables = dynamic_ws.dynamic_variables
    temp_vec = dynamic_ws.temporary_values
    steadystate = results.trends.endogenous_steady_state
    terminalvalues = view(guess_values, :, periods)
    params = work.params
    JJ = perfect_foresight_ws.J
    lb = perfect_foresight_ws.lb
    ub = perfect_foresight_ws.ub
    permutations = perfect_foresight_ws.permutations
    mcp!(lb, ub, permutations, m.mcps, context, periods)

    exogenous = perfect_foresight_ws.x



    ws_threaded = [
        Dynare.DynamicWs(
            m.endogenous_nbr,
            m.exogenous_nbr,
            length(dynamic_variables),
            length(temp_vec),
        ) for i = 1:Threads.nthreads()
    ]

    function f!(residuals, y)
        @debug "$(now()): start f!"
        get_residuals!(
            residuals,
            vec(y'),
            initialvalues,
            terminalvalues,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
 #           ddf,
            periods,
            temp_vec,
            permutations = permutations,
        )
        @debug "$(now()): end f!"
    end

    function JA!(
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
        y::AbstractVecOrMat{Float64},
    )
        @debug "$(now()): start J!"
        A = makeJacobian!(
            JJ,
            vec(y),
            initialvalues,
            terminalvalues,
            exogenous,
            context,
            periods,
            ws_threaded,
            permutations = permutations,
        )
        @debug count(!iszero, A) / prod(size(A))
        @debug "$(now()): end J!"
        return A
    end

    function fj!(residuals, JJ, y)
        f!(residuals, vec(y))
        JA!(JJ, vec(y))
    end

    @debug "$(now()): start makeJacobian"
    A0 = makeJacobian!(
        JJ,
        vec(guess_values),
        initialvalues,
        terminalvalues,
        exogenous,
        context,
        periods,
        ws_threaded,
    )
    @debug "$(now()): end makeJacobian"
    @debug "$(now()): start f!"
    f!(residuals, vec(guess_values))
    @debug "$(now()): end f!"
    @debug "$(now()): start J!"
    JA!(A0, guess_values)

    function J!(y::AbstractVecOrMat{Float64})
        JA!(A0, y)
        return JA!(A0, y)
    end

    @debug "$(now()): end J!"
    df = OnceDifferentiable(f!, J!, vec(guess_values), residuals, A0)
    @debug "$(now()): start nlsolve"

    rr = copy(residuals)
    F = lu(A0)
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

# Using PATHSolver
function solve_path!(F, J, lb, ub, initial_values; kwargs...)
    n = length(initial_values)
    @assert n == length(lb) == length(ub)

    function function_callback(n::Cint, z::Vector{Cdouble}, F_val::Vector{Cdouble})
        @assert n == length(z) == length(F_val)
        F(F_val, z)
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

        Jac = J(z)  # SparseMatrixCSC

        # Pouring Jac::SparseMatrixCSC to the format used in PATH
        for i = 1:n
            col_start[i] = Jac.colptr[i]
            col_len[i] = Jac.colptr[i+1] - Jac.colptr[i]
        end

        rv = rowvals(Jac)
        nz = nonzeros(Jac)
        num_nonzeros = SparseArrays.nnz(Jac)
        for i = 1:num_nonzeros
            row[i] = rv[i]
            data[i] = nz[i]
        end

        return Cint(0)
    end

    # var_name, F_name in RawIndex
    var_name = Vector{String}(undef, 0)
    F_name = Vector{String}(undef, 0)

    # overestimating number of nonzeros in Jacobian
    nnz = max(SparseArrays.nnz(J(lb)), SparseArrays.nnz(J(ub)))
    nnz = max(nnz, SparseArrays.nnz(J(initial_values)))
    n = length(lb)
    for i = 1:2
        z_rand = max.(lb, min.(ub, rand(Float64, n)))
        nnz_rand = SparseArrays.nnz(J(z_rand))
        nnz = max(nnz, nnz_rand)
    end
    nnz = min(2 * nnz, n^2)


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

    # This function has changed the content of m already.
    return Int(status), z, info
end


function ramsey_constraints!(context, field)
    for c in field["ramsey_model_constraints"]
        constraint = split(c["constraint"], limit = 3)
        eq = context.symboltable[constraint[1]].orderintype
        push!(context.models[1].mcps, (eq, constraint...))
    end
end

function get_mcps!(mcps, model)
    for (i, eq) in enumerate(model)
        tags = get(eq, "tags", "")
        if !isempty(tags)
            tag = get(eq["tags"], "mcp", "")
            if !isempty(tag)
                push!(mcps, (i, split(tag, limit = 3)...))
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
        ivar = context.symboltable[var].orderintype
        boundary = Dynare.dynare_parse_eval(String(expr), context)
        if op[1] == '<'
            for i = ivar:n:n*periods
                ub[i] = boundary
            end
        elseif op[1] == '>'
            for i = ivar:n:n*periods
                lb[i] = boundary
            end
        else
            error("MCP operator must be <, <=, > or >=")
        end
        if ivar != eqn
            push!(permutations, (ivar, eqn))
        end
    end
end
