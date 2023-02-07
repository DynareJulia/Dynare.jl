include("makeA.jl")

"""
PerfectForesightSetupOptions type 
    periods::Int64 - number of periods in the simulation [required]
    datafile::String - optional data filename with guess values for the simulation and initial and terminal values
"""
struct PerfectForesightSetupOptions
    periods::Int64
    datafile::String
    function PerfectForesightSetupOptions(options::Dict{String,Any})
        periods = 0
        datafile = ""
        for (k, v) in pairs(options)
            if k == "periods"
                periods = v::Int64
            elseif k == "datafile"
                datafile = v::String
            end
        end
        if periods == 0
            throw(DomainError(periods, "periods must be set to a number greater than zero"))
        end
        new(periods, datafile)
    end
end


function perfect_foresight_setup!(context, field)
    options = PerfectForesightSetupOptions(get(field, "options", Dict{String,Any}()))
    context.work.perfect_foresight_setup["periods"] = options.periods
    context.work.perfect_foresight_setup["datafile"] = options.datafile
end

@enum PerfectForesightAlgo trustregionA
@enum LinearSolveAlgo ilu pardiso
@enum InitializationAlgo initvalfile steadystate firstorder linearinterpolation

struct PerfectForesightSolverOptions
    algo::PerfectForesightAlgo
    display::Bool
    homotopy::Bool
    linear_solve_algo::LinearSolveAlgo
    maxit::Int64
    tolf::Float64
    tolx::Float64
    function PerfectForesightSolverOptions(context, field)
        algo = trustregionA
        display = true
        homotopy = false
        linear_solve_algo = ILU
        maxit = 50
        tolf = 1e-5
        tolx = 1e-5
        for (k, v) in pairs(options)
            if k == "stack_solve_algo"
                algo = v::Int64
            elseif k == "noprint"
                display = false
            elseif k == "print"
                display = true
            elseif k == "homotopy"
                homotopy = true
            elseif k == "no_homotopy"
                homotopy = false
            elseif k == "solve_algo"
                linear_solve_algo = v
            elseif k == "maxit"
                maxit = v
            elseif k == "tolf"
                tolf = v
            elseif k == "tolx"
                tolf = v
            end
        end
        new(algo, display, homotopy, linear_solve_algo, maxit, tolf, tolx)
    end
end

struct PerfectForesightWs
    y::Vector{Float64}
    x::Vector{Float64}
    shocks::Matrix{Float64}
    J::SparseMatrixCSC
    lb::Vector{Float64}
    ub::Vector{Float64}
    permutationsR::Vector{Tuple{Int64,Int64}}
    permutationsJ::Vector{Tuple{Int64,Int64}}
    function PerfectForesightWs(context, periods)
        m = context.models[1]
        modfileinfo = context.modfileinfo
        trends = context.results.model_results[1].trends
        y = Vector{Float64}(undef, m.endogenous_nbr)
        if m.exogenous_nbr > 0
            if modfileinfo.has_endval
                exogenous_steady_state = trends.exogenous_terminal_steady_state
            else
                exogenous_steady_state = trends.exogenous_steady_state
            end
            x = repeat(exogenous_steady_state, periods)
        else
            x = Vector{Float64}(undef, 0, 0)
        end
        if length(context.work.shocks) > 0
            shocks_tmp = context.work.shocks
            pmax = Int64(length(shocks_tmp) / m.exogenous_nbr)
            shocks = Matrix{Float64}(undef, pmax, m.exogenous_nbr)
            shocks .= reshape(shocks_tmp, (pmax, m.exogenous_nbr))
            # adding shocks to exogenous variables
            view(x, 1:pmax*m.exogenous_nbr) .= vec(transpose(shocks))
        else
            shocks = Matrix{Float64}(undef, 0, 0)
        end
        permutationsR = [(p[1], p[2]) for p in m.mcps]
        colptr = m.dynamic_g1_sparse_colptr
        rowval = m.dynamic_g1_sparse_rowval
        (J, permutationsJ) = makeJacobian(colptr,
                                          rowval,
                                          m.endogenous_nbr,
                                          periods,
                                          m.mcps)
        lb = Float64[]
        ub = Float64[]
        new(y, x, shocks, J, lb, ub, permutationsR, permutationsJ)
    end
end

function mcp_parse(mcps, context)
    mcp1 = Tuple{Int64, Int64, String, Float64}[]
    for m in mcps
        m1 = (m[1], context.symboltable[m[2]].orderintype, m[3], dynare_parse_eval(m[4], context))
        push!(mcp1, m1)
    end
    return mcp1
end

function perfect_foresight_solver!(context, field)
    periods = context.work.perfect_foresight_setup["periods"]
    datafile = context.work.perfect_foresight_setup["datafile"]
    m = context.models[1]
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr =  m.dynamic_tmp_nbr::Vector{Int64}
    dynamic_ws = DynamicWs(m.endogenous_nbr,
                           m.exogenous_nbr,
                           sum(tmp_nbr[1:2]),
                           m.dynamic_g1_sparse_colptr,
                           m.dynamic_g1_sparse_rowval)

    perfect_foresight_ws = PerfectForesightWs(context, periods)
    X = perfect_foresight_ws.shocks
    initial_values = get_dynamic_initialvalues(context)
    terminal_values = get_dynamic_terminalvalues(context, periods)
    guess_values = perfect_foresight_initialization!(
        context,
        periods,
        datafile,
        X,
        perfect_foresight_ws,
        steadystate,
        dynamic_ws,
    )
    if haskey(field, "options") && get(field["options"], "lmmcp.status", false)
        mcp_perfectforesight_core!(
            perfect_foresight_ws,
            context,
            periods,
            guess_values,
            initial_values,
            terminal_values,
            dynamic_ws,
        )
    else
        perfectforesight_core!(
            perfect_foresight_ws,
            context,
            periods,
            guess_values,
            initial_values,
            terminal_values,
            dynamic_ws,
        )
    end
end

function get_dynamic_initialvalues(context::Context)
    endo_nbr = context.models[1].endogenous_nbr 
    work = context.work
    modfileinfo = context.modfileinfo
    y0 = zeros(context.models[1].endogenous_nbr)
    if modfileinfo.has_histval
        @views for i in eachindex(skipmissing(work.histval[lastindex(work.histval, 1), 1:endo_nbr]))
            y0[i] = work.histval[end, i]
        end
        return y0
    elseif modfileinfo.has_initval_file
        @views for i in eachindex(skipmissing(work.initval[lastindex(work.initval, 1), 1:endo_nbr]))
            y0[i] = work.initval[end, i]
        end
        return y0
    else
        trends = context.results.model_results[1].trends
        if isempty(trends.endogenous_steady_state)
            compute_steady_state!(context, Dict{String,Any}())
        end
        return trends.endogenous_steady_state
    end
end

function get_dynamic_terminalvalues(context::Context, periods)
    work = context.work
    modfileinfo = context.modfileinfo
    yT = zeros(context.models[1].endogenous_nbr)
    if modfileinfo.has_initval_file
        @views for i in eachindex(skipmissing(view(work.initval, size(work.initval, 1), :)))
            yT[i] = work.initval[end, i]
        end
        return yT
    else
        trends = context.results.model_results[1].trends
        if isempty(trends.endogenous_steady_state)
            compute_steady_state!(context, Dict{String,Any}())
        end
        if modfileinfo.has_endval
            return trends.endogenous_terminal_steady_state
        else
            return trends.endogenous_steady_state
        end
    end
end

function perfect_foresight_initialization!(
    context,
    periods,
    datafile,
    exogenous,
    perfect_foresight_ws,
    algo::InitializationAlgo,
    dynamic_ws::DynamicWs,
)
    modfileinfo = context.modfileinfo
    trends = context.results.model_results[1].trends
    if algo == initvalfile
        initval = work.initval
        guess_values = view(initval, :, 2:periods - 1)
    elseif algo == linearinterpolation
    elseif algo == steadystate 
        if isempty(trends.endogenous_steady_state)
            compute_steady_state!(context, Dict{String,Any}())
        end
        if modfileinfo.has_endval
            x = trends.endogenous_terminal_steady_state
        else
            x = trends.endogenous_steady_state
        end
        guess_values = repeat(x, periods)
    elseif algo == firstorder
        guess_values = simul_first_order!(context, periods, exogenous, dynamic_ws)
    end
    return guess_values
end

function simul_first_order!(
    context::Context,
    periods::Int64,
    X::AbstractVector{Float64},
    dynamic_ws::DynamicWs,
)
    pre_options = Dict{String,Any}("periods" => periods)
    options = StochSimulOptions(pre_options)
    m = context.models[1]
    results = context.results.model_results[1]
    params = context.work.params
    compute_stoch_simul!(context, dynamic_ws, params, options)
    steadystate = results.trends.endogenous_steady_state
    linear_trend = results.trends.endogenous_linear_trend
    y0 = zeros(m.endogenous_nbr)
    simulresults = Matrix{Float64}(undef, m.endogenous_nbr, periods + 1)
    work = context.work
    histval = work.histval
    modfileinfo = context.modfileinfo
    if modfileinfo.has_histval
        for i in eachindex(skipmissing(view(work.histval, size(work.histval, 1), :)))
            y0[i] = work.histval[end, i]
        end
    else
        if work.model_has_trend[1]
            y0 .= steadystate - linear_trend
        else
            y0 .= steadystate
        end
    end
    A = zeros(m.endogenous_nbr, m.endogenous_nbr)
    B = zeros(m.endogenous_nbr, m.exogenous_nbr)
    make_A_B!(A, B, m, results)
    simul_first_order!(simulresults, y0, steadystate, A, B, X)
    return (view(simulresults, :, 1),
            view(simulresults, :, 2:periods - 1),
            view(simulresults, :, periods))
end


function perfectforesight_core!(
    perfect_foresight_ws::PerfectForesightWs,
    context::Context,
    periods::Int64,
    y0::AbstractVector{Float64},
    initialvalues::Vector{<:Real},
    terminalvalues::Vector{<:Real},
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

    exogenous = perfect_foresight_ws.x
    ws_threaded = [
        Dynare.DynamicWs(
            m.endogenous_nbr,
            m.exogenous_nbr,
            length(temp_vec),
            m.dynamic_g1_sparse_colptr,
            m.dynamic_g1_sparse_rowval
        ) for i = 1:Threads.nthreads()
    ]
    nzval = similar(m.dynamic_g1_sparse_rowval, Float64)
    nzval1 = similar(nzval)
    
    f! = make_pf_residuals(
            initialvalues,
            terminalvalues,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
            periods,
            temp_vec,
            perfect_foresight_ws.permutationsR,
    )

    J! = make_pf_jacobian(
            DFunctions.dynamic_derivatives!,
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

#    function fj!(residuals, JJ, y)
#        f!(residuals, vec(y))
#        J!(JJ, vec(y))
#    end
    J!(JJ, y0)
    df = OnceDifferentiable(f!, J!, y0, residuals, JJ)
    @debug "$(now()): start nlsolve"


    res = nlsolve(df, y0, method = :robust_trust_region, show_trace = false, ftol=cbrt(eps()))
    print_nlsolver_results(res)
    @debug "$(now()): end nlsolve"
    endogenous_names = get_endogenous_longname(context.symboltable)
    push!(
        context.results.model_results[1].simulations,
        Simulation(
            "Sim1",
            "",
            TimeDataFrame(
                DataFrame(
                    transpose(reshape(res.zero, m.endogenous_nbr, periods)),
                    endogenous_names,
                ),
                UndatedDate(1),
            ),
        ),
    )
end

function make_pf_residuals(
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
        ) where T <: Real
    function f!(residuals::AbstractVector{T}, y::AbstractVector{T})
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

function make_pf_jacobian(
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
    nzval1::AbstractVector{T}
) where {T <: Real, N <: Integer}
    function J!(
        A::SparseArrays.SparseMatrixCSC{Float64,Int64},
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

function get_residuals!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    m::Model,
    periods::Int64,
    temp_vec::AbstractVector{Float64};
    permutations::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[],
)
    n = m.endogenous_nbr

    rx = 1:m.exogenous_nbr
    @views get_residuals_1!(
        residuals,
        endogenous,
        initialvalues,
        exogenous[rx],
        dynamic_variables,
        steadystate,
        params,
        temp_vec,
        permutations = permutations,
    )
    for t = 2:periods - 1
        rx = rx .+ m.exogenous_nbr
        @views get_residuals_2!(
            residuals,
            endogenous,
            exogenous[rx],
            steadystate,
            params,
            temp_vec,
            t,
            permutations = permutations,
        )
    end
    rx = rx .+ m.exogenous_nbr
    @views get_residuals_3!(
        residuals,
        endogenous,
        terminalvalues,
        exogenous[rx],
        dynamic_variables,
        steadystate,
        params,
        temp_vec,
        periods,
        permutations = permutations,
    )
    return residuals
end

@inline function reorder!(x, permutations)
    for p in permutations
        p1 = p[1]
        p2 = p[2]
        x[p2], x[p1] = x[p1], x[p2]
    end
end

@inline function reorder!(x, permutations, offset)
    for p in permutations
        p1 = p[1] + offset
        p2 = p[2] + offset
        x[p2], x[p1] = x[p1], x[p2]
    end
end

@inline function reorder1!(x, permutations, n)
    isempty(permutations) && return
    reorder!(x, permutations, 0)
    reorder!(x, permutations, n)
end

@inline function reorder2!(x, permutations, t, n)
    isempty(permutations) && return
    reorder!(x, permutations, (t - 2) * n)
    reorder!(x, permutations, (t - 1) * n)
    reorder!(x, permutations, t * n)
end

@inline function reorder3!(x, permutations, t, n)
    isempty(permutations) && return
    reorder!(x, permutations, (t - 2) * n)
    reorder!(x, permutations, (t - 1) * n)
end


function get_residuals_1!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    temp_vec::AbstractVector{Float64};
    permutations::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[],
)
    n = length(steadystate)
    copyto!(dynamic_variables, initialvalues)
    copyto!(dynamic_variables, n + 1, endogenous, 1, 2*n)
    vr = view(residuals, 1:n)
    DFunctions.dynamic!(
        temp_vec,
        vr,
        dynamic_variables,
        exogenous,
        params,
        steadystate,
    )
    reorder!(vr, permutations)
end

function get_residuals_2!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    temp_vec::AbstractVector{Float64},
    t::Int64;
    permutations::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[],
)
    n = length(steadystate)
    k1 = (t - 1)*n  .+ (1:n)
    k2 = (t - 2)*n .+ (1:3*n)
    @views vr = residuals[k1]
    DFunctions.dynamic!(
        temp_vec,
        vr,
        endogenous[k2],
        exogenous,
        params,
        steadystate,
    )
    reorder!(vr, permutations)
end

function get_residuals_3!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    temp_vec::AbstractVector{Float64},
    periods::Int64;
    permutations::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[],
)
    n = length(steadystate)
    copyto!(dynamic_variables, 1, endogenous, (periods - 2)*n + 1, 2*n)
    copyto!(dynamic_variables, 2*n + 1, terminalvalues)
    @views vr = residuals[(periods - 1)*n .+ (1:n)]
    DFunctions.dynamic!(
        temp_vec,
        vr,
        dynamic_variables,
        exogenous,
        params,
        steadystate,
    )
    reorder!(vr, permutations)
end

function print_nlsolver_results(r)
    @printf "Results of Nonlinear Solver Algorithm\n"
    @printf " * Algorithm: %s\n" r.method
    @printf " * Inf-norm of residuals: %f\n" r.residual_norm
    @printf " * Iterations: %d\n" r.iterations
    @printf " * Convergence: %s\n" converged(r)
    @printf "   * |x - x'| < %.1e: %s\n" r.xtol r.x_converged
    @printf "   * |f(x)| < %.1e: %s\n" r.ftol r.f_converged
    @printf " * Function Calls (f): %d\n" r.f_calls
    return
end
include("PATH_interface.jl")
