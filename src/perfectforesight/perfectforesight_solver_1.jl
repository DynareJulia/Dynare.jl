#using Dynare
using IncompleteLU
using LinearAlgebra
using SuiteSparse
using SuiteSparse.UMFPACK

#include("perfectforesight_solvers.jl")
#include("gmres_solver.jl")
include("makeA.jl")


abstract type LinearSolver end
struct LuSolver <: LinearSolver end
struct GmresSolver <: LinearSolver end
abstract type NLStrategy end
struct FullNewtonStep <: NLStrategy end
struct LineSearch <: NLStrategy end
struct DogLeg <: NLStrategy end

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
    n = m.endogenous_nbr

    Dynare.get_initial_dynamic_endogenous_variables!(
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

    Dynare.get_dynamic_endogenous_variables!(dynamic_variables, endogenous, lli, t)
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

    Dynare.get_terminal_dynamic_endogenous_variables!(
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
    n = m.endogenous_nbr

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

function solve1(
    residuals,
    y,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
    JJ,
    context,
    ws_threaded,
    n,
)
    dy = zeros(n * periods)
    for i = 1:10
        get_residuals!(
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
        A = makeJacobian!(
            JJ,
            vec(y),
            initialvalues,
            terminalvalues,
            exogenous,
            context,
            periods,
            ws_threaded,
        )
        dy .= A \ vec(residuals)
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end
end

function solve2(
    residuals,
    y,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
    JJ,
    context,
    ws_threaded,
    n,
    ws,
)
    tmp = zeros(n * periods)
    dy = zeros(n * periods)
    for i = 1:10
        get_residuals!(
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
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        fill!(tmp, 0.0)
        Dynare.gmres!(tmp, A, residuals, Pr = ws.P, maxiter = 8)
        ldiv!(dy, ws.P, tmp)
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end
    0
end

function solve3(
    τ,
    residuals,
    y,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
    JJ,
    context,
    ws_threaded,
    n,
    verbose = false,
)
    tmp = zeros(n * periods)
    dy = zeros(n * periods)
    for i = 1:10
        get_residuals!(
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
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        lu = ilu(A, τ = τ)
        fill!(tmp, 0.0)
        Dynare.gmres!(dy, A, residuals, Pr = lu, verbose = verbose, maxiter = 4)
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end

end

function solve4(
    residuals,
    y,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
    JJ,
    context,
    ws_threaded,
    n,
)
    tmp = zeros(n * periods)
    dy = zeros(n * periods)
    local F::UmfpackLU
    for i = 1:10
        get_residuals!(
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
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        F = lu(A)
        ldiv!(dy, F, vec(residuals))
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end

end


function solve5(
    residuals,
    y,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
    JJ,
    context,
    ws_threaded,
    n,
    ps,
)
    tmp = zeros(n * periods)
    dy = zeros(n * periods)
    local F::UmfpackLU
    for i = 1:10
        get_residuals!(
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
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        solve!(ps, dy, A, vec(residuals))
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end

end
