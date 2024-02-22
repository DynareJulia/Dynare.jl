include("makeA.jl")

@enum PerfectForesightAlgo trustregionA
@enum LinearSolveAlgo ilu pardiso
@enum InitializationAlgo initvalfile steadystate firstorder linearinterpolation

struct DynarePerfectForesightSimulationFailed <: Exception end
function Base.showerror(io::IO, e::DynarePerfectForesightSimulationFailed)
    print("""
    Dynare couldn't compute the perfectf foresight simulation.
    Either there is no solution or the guess values
    are too far from the solution
     """)
end

abstract type LinearSolver end
struct IluLS <: LinearSolver end
struct PardisoLS <: LinearSolver end
    
function linear_solver!(::IluLS,
                        x::AbstractVector{Float64},
                        A::AbstractMatrix{Float64},
                        b::AbstractVector{Float64})
    x .= A\b
end

abstract type NonLinearSolver end
struct PathNLS <: NonLinearSolver end
struct DefaultNLS <: NonLinearSolver end

mutable struct PerfectForesightOptions
    algo::PerfectForesightAlgo
    datafile::String
    display::Bool
    homotopy::Bool
    initialization_algo::InitializationAlgo
    solve_algo
    linear_solve_algo::LinearSolveAlgo
    maxit::Int
    mcp::Bool
    method
    periods::Int
    show_trace::Bool
    tolf::Float64
    tolx::Float64
end

function PerfectForesightOptions(context::Context, field::Dict{String,Any})
    algo = trustregionA
    datafile = context.work.perfect_foresight_setup["datafile"]
    display = true
    homotopy = false
    initialization_algo = steadystate
    solve_algo = NewtonRaphson(linesearch = LineSearchesJL(method = NonlinearSolve.BackTracking()))
    linear_solve_algo = ilu
    maxit = 50
    mcp = false
    method = TrustRegion(radius_update_scheme=RadiusUpdateSchemes.NLsolve)
    periods =  context.work.perfect_foresight_setup["periods"]
    show_trace = false
    tolf = 1e-5
    tolx = 1e-5
    if haskey(field, "options")
        for (k, v) in pairs(field["options"])
            if k == "noprint"
                display = false
            elseif k == "print"
                display = true
            elseif k == "homotopy"
                homotopy = true
            elseif k == "lmmcp.status"
                mcp = true
            elseif k == "no_homotopy"
                homotopy = false
            elseif k == "maxit"
               maxit = v
            elseif k == "tolf"
                tolf = v
            elseif k == "tolx"
                tolf = v
            end
        end
    end 
    return PerfectForesightOptions(algo, datafile, display, homotopy,
                                   initialization_algo, solve_algo,
                                   linear_solve_algo, maxit, mcp,
                                   method, periods, show_trace, tolf, tolx)
end

struct FlipInformation
    flipped_variables::Vector{Bool}
    ix_period::SparseVector{Int, Int}
    ix_stack::SparseVector{Int, Int}
    iy_period::SparseVector{Int, Int}
    
    function FlipInformation(context, periods, infoperiod)
        m = context.models[1]
        symboltable = context.symboltable
        endogenous_nbr = m.endogenous_nbr
        exogenous_nbr = m.exogenous_nbr
        if isempty(context.work.scenario)
            flipped_variables = spzeros(Int, 0)
            ix_period = spzeros(Int, 0)
            ix_stack = spzeros(Int, 0)
            iy_period = spzeros(Int, 0)
        else
            N = periods*m.endogenous_nbr
            flipped_variables = repeat([false], N)
            ix_period = spzeros(Int, N)
            ix_stack = similar(ix_period)
            iy_period = similar(ix_period)
            s1 = context.work.scenario[infoperiod]
            for (period, f1) in s1
                p = length(infoperiod:period)
                for (s, f2) in f1
                    ss = string(s)
                    if is_endogenous(ss, symboltable)
                        iy1 = symboltable[ss].orderintype
                        iy2 = (p - 1)*endogenous_nbr + iy1
                        ix1 = symboltable[string(f2[2])].orderintype
                        ix2 = (p - 1)*exogenous_nbr + ix1
                        ix_period[iy2] = 3*endogenous_nbr + ix1
                        ix_stack[iy2] = ix2
                        iy_period[iy2] = iy1
                        flipped_variables[iy2] = true
                    end
                end
            end
        end
        new(flipped_variables, ix_period, ix_stack, iy_period)
    end
end

struct PerfectForesightWs
    y::Vector{Float64}
    x::Vector{Float64}
    shocks::Vector{Float64}
    J::SparseMatrixCSC{Float64, Int64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    permutationsR::Vector{Tuple{Int64,Int64}}  # permutations indices for residuals
    permutationsJ::Vector{Tuple{Int64,Int64}}  # permutations indices for Jacobian rows
    flipinfo::FlipInformation
    function PerfectForesightWs(context::Context, periods=Int; infoperiod = Undated(1))
        m = context.models[1]
        modfileinfo = context.modfileinfo
        trends = context.results.model_results[1].trends
        w = context.work
        y = Vector{Float64}(undef, m.endogenous_nbr)
        if m.exogenous_nbr > 0
            if modfileinfo.has_endval
                if isempty(trends.exogenous_terminal_steady_state)
                    if !isempty(w.endval_exogenous)
                        exogenous_steady_state = zeros(m.exogenous_nbr)
                    else
                        exogenous_steady_state = w.endval_exogenous
                    end
                else
                    exogenous_steady_state = trends.exogenous_terminal_steady_state
                end
            else
                if isempty(trends.exogenous_steady_state)
                    if isempty(w.initval_exogenous)
                        exogenous = zeros(m.exogenous_nbr)
                    else
                        exogenous_steady_state = w.initval_exogenous
                    end
                else
                    exogenous_steady_state = trends.exogenous_steady_state
                end
            end
            x = repeat(exogenous_steady_state, periods)
        else
            x = Vector{Float64}(undef, 0)
        end
        if length(context.work.shocks) > 0
            shocks = context.work.shocks
            pmax = Int64(length(shocks) / m.exogenous_nbr)
            # adding shocks to exogenous variables
            view(x, 1:pmax*m.exogenous_nbr) .= shocks
        else
            shocks = Vector{Float64}(undef, 0)
        end
        permutationsR = [(p[1], p[2]) for p in m.mcps]
        colptr = m.dynamic_g1_sparse_colptr
        rowval = m.dynamic_g1_sparse_rowval
        flipinfo = FlipInformation(context, periods, infoperiod)
        if isempty(context.work.scenario)
            (J, permutationsJ) = makeJacobian(colptr,
                                              rowval,
                                              m.endogenous_nbr,
                                              periods,
                                              m.mcps)
        else
            (J, permutationsJ) = makeJacobian(colptr,
                                              rowval,
                                              m.endogenous_nbr,
                                              periods,
                                              m.mcps,
                                              flipinfo)
        end            
            lb = Float64[]
        ub = Float64[]
        new(y, x, shocks, J, lb, ub, permutationsR, permutationsJ, flipinfo)
    end
end

include("conditional_simulation.jl")

function perfect_foresight_setup!(context::Context, field=Dict{String, Any})
    periods = 0
    datafile = ""
    for (k, v) in pairs(field["options"])
        if k == "periods"
            periods = v::Int64
        elseif k == "datafile"
            datafile = v::String
        end
    end
    if periods == 0
        throw(DomainError(periods, "periods must be set to a number greater than zero"))
    end
    context.work.perfect_foresight_setup["periods"] = periods
    context.work.perfect_foresight_setup["datafile"] = datafile
end

function perfect_foresight_solver!(context, field)
    options = PerfectForesightOptions(context, field)
    _perfect_foresight!(context, options)
end
"""
    perfect_foresight!(; periods, context = context, display = true,
                       linear_solve_algo=ilu, maxit = 50, mcp = false,
                       tolf = 1e-5, tolx = 1e-5)
 
 # Keyword arguments
 - `periods::Int`: number of periods in the simulation [required]
 - `context::Context=context`: context in which the simulation is computed
 - `display::Bool=true`: whether to display the results
 - `solve_algo = NewtonRaphson(linesearch = LineSearchesJL(method = NonlinearSolve.BackTracking()))`: nonlinear solver algorithm
 - `linear_solve_algo::LinearSolveAlgo=ilu`: algorithm used for the solution of the linear
   problem. Either `ilu` or `pardiso`. `ilu` is the sparse linear solver used by default in Julia.
   To use the Pardiso solver, write `using Pardiso` before running Dynare.
 - `maxit::Int=50` maximum number of iterations
 - `method`: solver method inside `solve_algo` package
 - `mcp::Bool=false`L whether to solve a mixed complementarity problem with
   occasionally binding constraints
 - `show_trace::Bool=false`: whether to show solver iterations
 - `tolf::Float64=1e-5`: tolerance for the norm of residualts
 - `tolx::Float64=1e-5`: tolerance for the norm of the change in the result
"""
function perfect_foresight!(;context::Context = context,
                            algo::PerfectForesightAlgo = trustregionA,
                            datafile::String = "",
                            display::Bool = true,
                            homotopy::Bool = false,
                            initialization_algo::InitializationAlgo = steadystate,
                            solve_algo = NewtonRaphson(linesearch = LineSearchesJL(method = NonlinearSolve.BackTracking())),
                            linear_solve_algo::LinearSolveAlgo = ilu,
                            maxit::Int = 50,
                            mcp::Bool = false,
                            method = TrustRegion(radius_update_scheme=RadiusUpdateSchemes.NLsolve),
                            periods::Int,
                            show_trace::Bool = false,
                            tolf::Float64 = 1e-5,
                            tolx::Float64 = 1e-5)

    options = PerfectForesightOptions(algo, datafile, display,
                                      homotopy, initialization_algo,
                                      solve_algo, linear_solve_algo, maxit,
                                      mcp, method, periods, show_trace, tolf, tolx)

    scenario = context.work.scenario
    check_scenario(scenario)
    if isempty(scenario) || length(scenario) == 1
        _perfect_foresight!(context, options)
    else
        recursive_perfect_foresight!(context, options)
    end
end

function check_scenario(scenario)
    # all infoperiods
    k = collect(keys(scenario))
    # for each infoperiod > 1 
    for (i, k1) in enumerate(sort(k))
        if k1 == 1
            continue
        end
        # period shocked for that infoperiod
        sk1 = keys(scenario[k1])
        # for all shock periods in the current infoperiod
        for p1 in sk1
            # check that future relevant infoperiods contain a comparable shock
            for k2 = k[i+1:end]
                # check only infoperiods <= this shock period
                if k2 > p1
                    break
                end
                p2 = collect(keys(scenario[k2]))
                if !(p1 in p2)
                    error("A shock must be explicitly confirmed or modified in all relevant subsequent infoperiods")
                else
                    for v in p1
                        if !(v in keys(scenario[k2][p1]))
                            error("A shock must be explicitly confirmed or modified in all relevant subsequent infoperiods")

                        end
                    end
                end
            end
        end
    end
end
                
function mcp_parse(mcps::Vector{Tuple{Int64,Int64,String,String}}, context::Context)
    mcp1 = Tuple{Int64, Int64, String, Float64}[]
    for m in mcps
        m1 = (m[1], m[2], m[3], dynare_parse_eval(m[4], context))
        push!(mcp1, m1)
    end
    return mcp1
end

function _perfect_foresight!(context::Context, options::PerfectForesightOptions)
    datafile = options.datafile
    periods = options.periods
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
        terminal_values,
        periods,
        datafile,
        X,
        perfect_foresight_ws,
        options.initialization_algo,
        dynamic_ws,
    )
    if options.mcp
        if isnothing(Base.get_extension(Dynare, :PathSolver))
            error("You must load PATH with 'using PATHSolver'")
        end
        if isempty(context.work.scenario)
            mcp_perfectforesight_core!(
                PathNLS(),
                perfect_foresight_ws,
                context,
                periods,
                guess_values,
                initial_values,
                terminal_values,
                dynamic_ws,
                maxit = options.maxit,
                tolf = options.tolf,
                tolx = options.tolx
            )
        else
            mcp_perfectforesight_core!(
                PathNLS(),
                perfect_foresight_ws,
                context,
                periods,
                guess_values,
                initial_values,
                terminal_values,
                dynamic_ws,
                perfect_foresight_ws.flipinfo,
                1,
                maxit = options.maxit,
                tolf = options.tolf,
                tolx = options.tolx
            )
        end            
    elseif !isempty(context.work.scenario)
        perfectforesight_core_conditional!(
            perfect_foresight_ws,
            context,
            periods,
            guess_values,
            initial_values,
            terminal_values,
            dynamic_ws,
            perfect_foresight_ws.flipinfo,
            1,
            linear_solve_algo = nothing,
            maxit = options.maxit,
            method = options.method,
            nonlinear_solve_algo = options.solve_algo,
            show_trace = options.show_trace,
            tolf = options.tolf,
            tolx = options.tolx
        )
    else
        perfectforesight_core!(
            perfect_foresight_ws,
            context,
            periods,
            guess_values,
            initial_values,
            terminal_values,
            options.solve_algo,
            options.linear_solve_algo,
            dynamic_ws,
            linear_solve_algo = options.linear_solve_algo,
            maxit = options.maxit,
            method = options.method,
            nonlinear_solve_algo = "NonlinearSolve",
            show_trace = false,
            tolf = options.tolf,
            tolx = options.tolx
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
            return work.initval_endogenous
        else
            return trends.endogenous_steady_state
        end
    end
end

function get_dynamic_terminalvalues(context::Context, periods::Int)
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
            return work.endval_endogenous
        end
        if modfileinfo.has_endval
            return trends.endogenous_terminal_steady_state
        else
            return trends.endogenous_steady_state
        end
    end
end

function perfect_foresight_initialization!(
    context::Context,
    terminal_values::Vector{Float64},
    periods::Int,
    datafile::String,
    exogenous::Vector{Float64},
    perfect_foresight_ws::PerfectForesightWs,
    algo::InitializationAlgo,
    dynamic_ws::DynamicWs,
    )
    work = context.work
    modfileinfo = context.modfileinfo
    trends = context.results.model_results[1].trends
    if algo == initvalfile
        initval = work.initval
        guess_values = view(initval, :, 2:periods - 1)
    elseif algo == linearinterpolation
    elseif algo == steadystate 
        if isempty(trends.endogenous_steady_state)
            compute_steady_state!(context)
        end
        if modfileinfo.has_endval
            x = trends.endogenous_terminal_steady_state
        else
            x = trends.endogenous_steady_state
        end
        guess_values = repeat(x, periods)
    elseif algo == firstorder
        if isempty(trends.endogenous_steady_state)
            compute_steady_state!(context)
        end
        y = simul_first_order!(context, periods, perfect_foresight_ws.x, dynamic_ws)
        guess_values = vec(y[1])
        terminal_values .= y[2]
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
    endogenous = results.trends.endogenous_steady_state
    endogenous3 = repeat(endogenous, 3)
    exogenous = results.trends.exogenous_steady_state
    model = context.models[1]
    compute_first_order_solution!(
        context,
        endogenous3,
        exogenous,
        endogenous,
        params,
        model,
        dynamic_ws,
        options
    )   
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

    exogenous = transpose(reshape(X, m.exogenous_nbr, periods))
    simul_first_order!(simulresults, y0, steadystate, A, B, exogenous)

    return (view(simulresults, :, 1:periods - 1),
            view(simulresults, :, periods))
end

using LinearSolve

function perfectforesight_core!(
    perfect_foresight_ws::PerfectForesightWs,
    context::Context,
    periods::Int64,
    y0::AbstractVector{Float64},
    initialvalues::Vector{<:Real},
    terminalvalues::Vector{<:Real},
    solve_algo,
    linear_solve_algo::LinearSolveAlgo,
    dynamic_ws::DynamicWs;
    maxit = 50,
    nonlinear_solve_algo = "NonlinearSolve",
    linear_solve_algo = nothing,
    method = TrustRegion(),
    show_trace = true,
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

    j! = make_pf_jacobian(
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
    function jj!(JJ, y, p)
        JJt = copy(transpose(JJ))
        j!(JJt, y, p)
        JJ .= transpose(JJt)
    end

    results = dynare_nonlinear_solvers(f!, j!, JJ, [], [], y0, exogenous, work.params,
                                       linear_solve_algo = linear_solve_algo,
                                       maxit = maxit,
                                       method = method,
                                       solve_algo = nonlinear_solve_algo,
                                       show_trace = show_trace,
                                       tolf = tolf,
                                       tolx = tolx)
    @debug "$(now()): end nlsolve"
    make_simulation_results!(context::Context, results, exogenous, terminalvalues, periods)
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
    function f!(residuals::AbstractVector{T}, y::AbstractVector{T}, params::AbstractVector{T})
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
        params
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
    @info @sprintf "Results of Nonlinear Solver Algorithm"
    @info @sprintf " * Algorithm: %s" r.alg
    @info @sprintf " * Inf-norm of residuals: %f" norm(r.resid)
    @info @sprintf " * Convergence: %s" r.retcode
    @info @sprintf " * Function Calls (f): %d" r.stats.nf
    return
end

function make_simulation_results!(context::Context, y, x,terminalvalues, periods)
    m = context.models[1]
    trends = context.results.model_results[1].trends
    work = context.work
    endogenous_names = get_endogenous(context.symboltable)
    exogenous_names = get_exogenous(context.symboltable)
    if !isempty(work.histval)
        initialvalues = work.histval
    else 
        if !isempty(work.initval_exogenous)
            initialvalues_x = work.initval_exogenous
        else
            initialvalues_x = trends.exogenous_steady_state
        end
        if !isempty(work.initval_endogenous)
            initialvalues = vcat(work.initval_endogenous,
                                 initialvalues_x)
        else
            initialvalues = vcat(trends.endogenous_steady_state,
                                 initialvalues_x)
        end
    end

    if size(initialvalues, 2) == 1
        initialvalues = reshape(initialvalues, 1, length(initialvalues))
    end 

    data = vcat(initialvalues,
                hcat(
                    transpose(reshape(y, m.endogenous_nbr, periods)),
                    transpose(reshape(x, m.exogenous_nbr, periods))
                ),
                transpose(vcat(terminalvalues,
                               trends.exogenous_steady_state))
            )
    push!(
        context.results.model_results[1].simulations,
        Simulation(
            Undated(1),
            Undated(periods),
            "Sim1",
            "",
            AxisArrayTable(
                data,
                Undated(1 - size(initialvalues, 1)):Undated(periods + 1),
                vcat(
                    [Symbol(s) for s in endogenous_names], 
                    [Symbol(s) for s in exogenous_names]
                    )
            )
        )
    )
end

function recursive_perfect_foresight!(context::Context, options::PerfectForesightOptions)
    periods = options.periods
    endogenous_names = get_endogenous(context.symboltable)
    endogenous_nbr = length(endogenous_names)
    exogenous_names = get_exogenous(context.symboltable)
    scenario = context.work.scenario
    if context.modfileinfo.has_histval
        initial_periods = lastindex(context.work.histval, 1)
    else
        initial_periods = 1
    end
    let simulation = Matrix{Float64}(undef, 0, 0)
        for (j, i) in enumerate(sort(collect(keys(scenario))))
            if i == 1
                initial_values = get_dynamic_initialvalues(context)
            else
                initial_values = Vector{Float64}(simulation[i + initial_periods - 1, 1:endogenous_nbr])
            end
            this_simulation = _recursive_perfect_foresight!(context, i, initial_values, options)
            data = Matrix(this_simulation[j].data)
            if i > 1
                simulation = vcat(simulation[1:(i + initial_periods - 1),:],  data[initial_periods + 1:end, :])
            else
                simulation = data
            end
        end
        push!(
            context.results.model_results[1].simulations,
            Simulation(
                Undated(1),
                Undated(periods),
                "Sim1",
                "",
                AxisArrayTable(
                    simulation,
                    Undated(1-initial_periods):Undated(periods + 1),
                    vcat(
                        [Symbol(s) for s in endogenous_names], 
                        [Symbol(s) for s in exogenous_names]
                    )
                )
            )
        )
    end
end

function _recursive_perfect_foresight!(context::Context, infoperiod, initial_values, options)
    datafile = options.datafile
    periods = options.periods - infoperiod + 1
    
    m = context.models[1]
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr =  m.dynamic_tmp_nbr::Vector{Int64}
    dynamic_ws = DynamicWs(m.endogenous_nbr,
                           m.exogenous_nbr,
                           sum(tmp_nbr[1:2]),
                           m.dynamic_g1_sparse_colptr,
                           m.dynamic_g1_sparse_rowval)

    perfect_foresight_ws = PerfectForesightWs(context, periods, infoperiod = infoperiod)
    X = perfect_foresight_ws.shocks
    terminal_values = get_dynamic_terminalvalues(context, periods)
    guess_values = perfect_foresight_initialization!(
        context,
        terminal_values,
        periods,
        datafile,
        X,
        perfect_foresight_ws,
        options.initialization_algo,
        dynamic_ws,
    )
    if options.mcp
        if isnothing(Base.get_extension(Dynare, :PathSolver))
            error("You must load PATH with 'using PATHSolver'")
        end
        if isempty(context.work.scenario)
            mcp_perfectforesight_core!(
                PathNLS(),
                perfect_foresight_ws,
                context,
                periods,
                guess_values,
                initial_values,
                terminal_values,
                dynamic_ws,
                maxit = options.maxit,
                tolf = options.tolf,
                tolx = options.tolx
            )
        else
            mcp_perfectforesight_core!(
                PathNLS(),
                perfect_foresight_ws,
                context,
                periods,
                guess_values,
                initial_values,
                terminal_values,
                dynamic_ws,
                perfect_foresight_ws.flipinfo,
                infoperiod,
                maxit = options.maxit,
                tolf = options.tolf,
                tolx = options.tolx
            )
        end
    elseif !isempty(context.work.scenario)
        perfectforesight_core_conditional!(
            perfect_foresight_ws,
            context,
            periods,
            guess_values,
            initial_values,
            terminal_values,
            dynamic_ws,
            perfect_foresight_ws.flipinfo,
            infoperiod,
            linear_solve_algo = options.linear_solve_algo,
            maxit = options.maxit,
            method = options.method,
            nonlinear_solve_algo = options.solve_algo,
            show_trace = options.show_trace,
            tolf = options.tolf,
            tolx = options.tolx
        )
    else
        perfectforesight_core!(
            perfect_foresight_ws,
            context,
            periods,
            guess_values,
            initial_values,
            terminal_values,
            options.solve_algo,
            options.linear_solve_algo,
            dynamic_ws,
            maxit = options.maxit,
            tolf = options.tolf,
            tolx = options.tolx
        )
    end
end


