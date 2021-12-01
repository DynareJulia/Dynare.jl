using NLsolve

include("makeA.jl")

"""
PerfectForesightSetupOptions type 
    periods::Int64 - number of periods in the simulation [required]
    datafile::String - optional data filename with guess values for the simulation
                       and initial and terminal values
"""
struct PerfectForesightSetupOptions
    periods::Int64
    datafile::String
    function PerfectForesightSetupOptions(options::Dict{String, Any})
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

@enum PerfectForesightAlgo trustregion
@enum LinearSolveAlgo ilu pardiso

struct PerfectForesightSolverOptions
    algo::PerfectForesightAlgo
    display::Bool
    homotopy::Bool
    linear_solve_algo::LinearSolveAlgo
    maxit::Int64
    tolf::Float64
    tolx::Float64
    function PerfectForesightSolverOptions(context, field)
        algo = trustregion
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
        new(algo, display, homotopy,  linear_solve_algo, maxit, tolf, tolx)
    end
end

struct PerfectForesightWs
    y::Vector{Float64}
    x::Vector{Float64}
    J::Jacobian
    function PerfectForesightWs(context, periods)
        m = context.models[1]
        y = Vector{Float64}(undef, (periods+2)*m.endogenous_nbr)
        x = Vector{Float64}(undef, (periods+2)*m.exogenous_nbr)
        J = Jacobian(context, periods)
        new(y, x, J)
    end
end

function perfect_foresight_solver!(context, field)
    periods = context.work.perfect_foresight_setup["periods"]
    datafile = context.work.perfect_foresight_setup["datafile"]
    perfect_foresight_ws = PerfectForesightWs(periods)
    perfect_foresight_initialization!(context, periods, datafile, perfect_foresight_ws)
end

function perfect_foresight_initialization!(context, periods, datafile, perfect_foresight_ws)
    y = perfect_foresight_ws.y
    x = perfect_foresight_ws.x
    simul_first_order!(context, periods)
end

function simul_first_order!(context::Context, periods::Int64)
    pre_options = Dict{String, Any}("periods" => periods)
    options = Dynare.StochSimulOptions(pre_options)
    m = context.models[1]
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr = m.dynamic!.tmp_nbr::Vector{Int64}
    ws = Dynare.DynamicWs(m.endogenous_nbr, m.exogenous_nbr, ncol, sum(tmp_nbr[1:2]))
    Dynare.stoch_simul_core!(context, ws, options)
end


function perfectforesight!(perfect_foresight_ws::PerfectForesightWs)
    f!(y) = get_residuals!(
            residuals,
            vec(y),
            initialvalues,
            terminalvalues,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            md,
            periods,
            temp_vec,
        )
    J!(y) = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)

    res = nlsolve(f!, J!, x0)
end

function get_residuals!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    m::Model,
    periods::Int64,
    temp_vec::AbstractVector{Float64},
)
    lli = m.lead_lag_incidence
    dynamic! = m.dynamic!.dynamic!

    get_residuals_1!(
        residuals,
        endogenous,
        initialvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        m,
        periods,
        temp_vec,
    )
    t1 = n + 1
    t2 = 2 * n
    for t = 2:periods-1
        get_residuals_2!(
            residuals,
            endogenous,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
            periods,
            temp_vec,
            t,
            t1,
            t2,
        )
        t1 += n
        t2 += n
    end
    get_residuals_3!(
        residuals,
        endogenous,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        m,
        periods,
        temp_vec,
        periods,
        t1,
        t2,
    )
    return residuals
end

function get_residuals_1!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    m::Model,
    periods::Int64,
    temp_vec::AbstractVector{Float64},
)
    lli = m.lead_lag_incidence
    dynamic! = m.dynamic!.dynamic!

    get_initial_dynamic_endogenous_variables!(
        dynamic_variables,
        endogenous,
        initialvalues,
        lli,
        2,
    )
    vr = view(residuals, 1:n)
    @inbounds Base.invokelatest(
        dynamic!,
        temp_vec,
        vr,
        dynamic_variables,
        exogenous,
        params,
        steadystate,
        2,
    )
end

function get_residuals_2!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    m::Model,
    periods::Int64,
    temp_vec::AbstractVector{Float64},
    t::Int64,
    t1::Int64,
    t2::Int64,
)
    lli = m.lead_lag_incidence
    dynamic! = m.dynamic!.dynamic!

    get_dynamic_endogenous_variables!(dynamic_variables, endogenous, lli, t)
    vr = view(residuals, t1:t2)
    @inbounds Base.invokelatest(
        dynamic!,
        temp_vec,
        vr,
        dynamic_variables,
        exogenous,
        params,
        steadystate,
        t,
    )
end

function get_residuals_3!(
    residuals::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    dynamic_variables::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    m::Model,
    periods::Int64,
    temp_vec::AbstractVector{Float64},
    t::Int64,
    t1::Int64,
    t2::Int64,
)
    lli = m.lead_lag_incidence
    dynamic! = m.dynamic!.dynamic!

    get_terminal_dynamic_endogenous_variables!(
        dynamic_variables,
        endogenous,
        terminalvalues,
        lli,
        t,
    )
    vr = view(residuals, t1:t2)
    @inbounds Base.invokelatest(
        dynamic!,
        temp_vec,
        vr,
        dynamic_variables,
        exogenous,
        params,
        steadystate,
        t,
    )
end


