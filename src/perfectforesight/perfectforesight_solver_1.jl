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


"""
geometric_series!(Y, A, B, X, T, tmp1, tmp2) -> Y

computes ``y_t = ∑_{i=0}^{T-1}A^iBx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X``.
"""
function geometric_series!(
    Y::StridedVecOrMat{Float64},
    A::StridedVecOrMat{Float64},
    B::StridedVecOrMat{Float64},
    X::StridedVecOrMat{Float64},
)
    n, periods = size(Y)
    mul!(Y, B, transpose(X))
    for t = periods-1:-1:1
        vt1 = view(Y, :, t + 1)
        vt2 = view(Y, :, t)
        mul!(vt2, A, vt1, 1.0, 1.0)
    end
    return Y
end

"""
simul_first_order_1!(Y, y0, G, H, X, T) -> Y

computes the solution of a perfect foresight linear model:
``y_t = Gy_{t-1} + ∑_{i=0}^{T-1}G^iHx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X`` containing shocks known at the beginning of the simulation.
vector ``y0`` contains the initial values of the simulation and ``y_0 = y0``.
Matrix ``Y`` doesn't contain the initial values.
 ``G`` and ``H`` matrices must have been computed by a linear rational expectation model solver.
"""
function simul_first_order_1!(
    Y::StridedMatrix{Float64},
    y0::AbstractVector{Float64},
    G::StridedMatrix{Float64},
    H::StridedMatrix{Float64},
    X::StridedMatrix{Float64},
)
    window = size(X, 1)
    fill!(Y, 0.0)
    vy = view(Y, :, 1:window)
    geometric_series!(vy, G, H, X)
    vy_1 = view(Y, :, 1)
    vy_1 = y0
    for t = 1:size(Y, 2)
        vy = view(Y, :, t)
        mul!(vy, G, vy_1, 1.0, 1.0)
        vy_1 = vy
    end
end

"""
simul_first_order_1!(Y, y0, c, G, H, X, T1, T, tmp1, tmp2) -> Y

computes the solution of a perfect foresight linear model:
``y_t -c = G(y_{t-1} - c) + ∑_{i=0}^{T-1}G^iHx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X`` containing shocks known at the beginning of the simulation.
 ``c`` is a constant vector.
vector ``y0`` contains the initial values of the simulation and ``y_0 = y0``.
Matrix ``Y`` doesn't contain the initial values.
 ``G`` and ``H`` matrices must have been computed by a linear rational expectation model solver.
"""
function simul_first_order_1!(
    Y::StridedMatrix{Float64},
    y0::AbstractVector{Float64},
    c::AbstractVector{Float64},
    G::StridedMatrix{Float64},
    H::StridedMatrix{Float64},
    X::StridedMatrix{Float64},
)
    y0 .-= c
    simul_first_order_1!(Y, y0, G, H, X)
    y0 .+= c
    Y .+= c
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
        A = makeJacobian!(JJ, vec(y), initialvalues, terminalvalues, exogenous, context, periods, ws_threaded)
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
