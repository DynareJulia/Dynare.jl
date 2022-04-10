using PATHSolver

function mcp_perfectforesight_core!(perfect_foresight_ws::PerfectForesightWs,
                                    context::Context,
                                    periods::Int64,
                                    y0::Matrix{Float64},
                                    dynamic_ws::DynamicWs)
    m = context.models[1]
    ddf = context.dynarefunctions
    results = context.results.model_results[1]
    work = context.work
    residuals = zeros(periods*m.endogenous_nbr)
    dynamic_variables = dynamic_ws.dynamic_variables
    temp_vec = dynamic_ws.temporary_values
    steadystate = results.trends.endogenous_steady_state
    initialvalues = steadystate
    terminalvalues = view(y0, :, periods)
    params = work.params
    JJ = perfect_foresight_ws.J
    
    exogenous = perfect_foresight_ws.x

   

    ws_threaded = [Dynare.DynamicWs(m.endogenous_nbr,
                                    m.exogenous_nbr,
                                    length(dynamic_variables),
                                    length(temp_vec))
                   for i = 1:Threads.nthreads()]

    function f!(residuals, y)
        @debug "$(now()): start f!"
        get_residuals!(residuals,
                       vec(y),
                       initialvalues,
                       terminalvalues,
                       exogenous,
                       dynamic_variables,
                       steadystate,
                       params,
                       m,
                       ddf,
                       periods,
                       temp_vec,
                       )
        @debug "$(now()): end f!"
    end
    
    function J!(A::SparseArrays.SparseMatrixCSC{Float64, Int64}, y::AbstractVecOrMat{Float64})
        @debug "$(now()): start J!"
        A = makeJacobian!(JJ, vec(y), initialvalues, terminalvalues, exogenous, context, periods, ws_threaded)
        @debug count(!iszero, A)/prod(size(A))
        @debug "$(now()): end J!"
    end
    
    function fj!(residuals, JJ, y)
        f!(residuals, vec(y))
        J!(JJ, vec(y))
    end

    @debug "$(now()): start makeJacobian"
    A0 = makeJacobian!(JJ, vec(y0), initialvalues, terminalvalues, exogenous, context, periods, ws_threaded)
    @debug "$(now()): end makeJacobian"
    @debug "$(now()): start f!"
    f!(residuals, vec(y0))
    @debug "$(now()): end f!"
    @debug "$(now()): start J!"
    J!(A0, y0)
    @debug "$(now()): end J!"
    df = OnceDifferentiable(f!, J!, vec(y0), residuals, A0)
    @debug "$(now()): start nlsolve"

    rr = copy(residuals)
    F = lu(A0)
    res = nlsolve(df, vec(y0), method=:robust_trust_region, show_trace=true)
    @debug "$(now()): end nlsolve"
    endogenous_names = get_endogenous_longname(context.symboltable)
    push!(context.results.model_results[1].simulations,
          Simulation("Sim1", "", TimeDataFrame(DataFrame(transpose(reshape(res.zero, m.endogenous_nbr, periods)),
                                                         endogenous_names),
                                               UndatedDate(1))))
end
end


# Using PATHSolver
function solve_path!(F, J, lb, ub, initial_values; kwargs...)

    function function_callback(n::Cint, z::Vector{Cdouble}, F_val::Vector{Cdouble})
        @assert n == length(z) == length(F_val)
        F_val .= F(z)
         return Cint(0)
    end

    function jacobian_callback(        
        n::Cint,
        nnz::Cint,
        z::Vector{Cdouble},
        col_start::Vector{Cint},
        col_len::Vector{Cint},
        row::Vector{Cint},
        data::Vector{Cdouble}
    )
        @assert n == length(z) == length(col_start) == length(col_len)
        @assert nnz == length(row) == length(data)    

        Jac = J(z)  # SparseMatrixCSC

        # Pouring Jac::SparseMatrixCSC to the format used in PATH
        for i in 1:n
            col_start[i] = Jac.colptr[i]
            col_len[i] = Jac.colptr[i + 1] - Jac.colptr[i]
        end

        rv = rowvals(Jac)
        nz = nonzeros(Jac)
        num_nonzeros = SparseArrays.nnz(Jac)
        for i in 1:num_nonzeros
            row[i] = rv[i]
            data[i] = nz[i]
        end

        return Cint(0)
    end

    # var_name, F_name in RawIndex
    var_name = Vector{String}(undef, 0)
    F_name = Vector{String}(undef, 0)

    # overestimating number of nonzeros in Jacobian
    nnz = max( SparseArrays.nnz(J(lb)), SparseArrays.nnz(J(ub)) )
    nnz = max( nnz, SparseArrays.nnz(J(initial_values)) )
    for i in 1:2
        z_rand = max.(lb, min.(ub, rand(Float64, n)))
        nnz_rand = SparseArrays.nnz(J(z_rand))
        nnz = max( nnz, nnz_rand )
    end
    nnz = min( 2 * nnz, n^2 )


    F_val = zeros(3)
    @show function_callback(Int32(3), initial_values, F_val)
    @show F_val
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
        kwargs...
    )

    # This function has changed the content of m already.
    return Int(status), z, info
end


function ramsey_contraints(context, field)
    for c in field["ramsey_model_constraints"]
        push!(context.work.mcps, c["constraint"])
    end
end
    
function get_tags!(mcps, modeljson)
    for (i, eq) in enumerate(modeljson["model"])
        tag = get(eq["tags"], "mcp", "")
        if !isempty(tag)
            push!(mcps, (i, tag))
        end
    end
end

function mcp!(lb, ub, permutations, mcps, context)
    for m in mcps
        (m1, m2) = m 
        (var, op, expr) = split(m2)
        
        if op[1] == '<'
            ub[m1] = Dynare.dynare_parse_eval(String(expr), context)
        elseif op[1] == '>'
            lb[m1] = Dynare.dynare_parse_eval(String(expr), context)
        else
            error("MCP operator must be <, <=, > or >=")
        end
        iva = context.symboltable[var].orderintype
        if iva != m1
            push!(permutations, (iva, m1))
        end
    end
end

mcps = []
get_tags!(mcps, modeljson)
lb = repeat([-Inf], 3)
ub = repeat([Inf], 3)
permutations = []
mcp!(lb, ub, permutations, mcps, context)

function reorder!(x, permutations, offset)
    for p in permutations
        p1 = p[1] + offset
        p2 = p[2] + offset
        x[p2], x[p1] = x[p1], x[p2]
    end
end

x = collect(1:10)
reorder!(x, [(1, 3), (2, 6)], 0)
@show x
reorder!(x, [(1, 3), (2, 6)], 0)
@show x
x = collect(1:20)
reorder!(x, [(1, 3), (2, 6)], 10)
@show x
