#using Dynare
using IncompleteLU
using LinearAlgebra
using SuiteSparse
using SuiteSparse.UMFPACK

include("perfectforesight_solvers.jl")
include("gmres_solver.jl")
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
function geometric_series!(Y::StridedVecOrMat{Float64},
                           A::StridedVecOrMat{Float64},
                           B::StridedVecOrMat{Float64},
                           X::StridedVecOrMat{Float64},
                           T::Integer,
                           tmp1::StridedVecOrMat{Float64},
                           tmp2::StridedVecOrMat{Float64})
    n, periods = size(Y)
    vr = view(Y, :, 1:T)
    vx = view(X, 1:T, :)
    mul!(vr, B, transpose(vx))
    vt1 = view(tmp1, :, 2:T)
    copyto!(tmp1, 1, Y, n+1, (periods-1)*n)
    for t = 2:T
        last = periods - t + 1
        vt1 = view(tmp1, :, 1:last)
        vt2 = view(tmp2, :, 1:last )
        mul!(vt2, A, vt1, 1.0, 0.0)
        vr = view(Y, :, 1:last)
        vr .+= vt2
        if t < T
            copyto!(tmp1, 1, tmp2, n + 1, (last - 1)*n)
        end
    end
end

"""
simul_first_order_1!(Y, y0, G, H, X, T1, T, tmp1, tmp2) -> Y

computes the solution of a perfect foresight linear model:
``y_t = Gy_{t-1} + ∑_{i=0}^{T-1}G^iHx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X`` containing shocks known at the beginning of the simulation.
``y0`` are the initial values of the simulation and ``y_0 = y0``.
 ``G`` and ``H`` matrices must have been computed by a linear rational expectation model solver.
"""
function simul_first_order_1!(Y::StridedMatrix{Float64},
                              y0::AbstractVector{Float64},
                              G::StridedMatrix{Float64},
                              H::StridedMatrix{Float64},
                              X::StridedMatrix{Float64},
                              T1::Int64,
                              T::Int64,
                              temp1::StridedMatrix{Float64},
                              temp2::StridedMatrix{Float64})
    vy = view(Y, :, 2:T1 + 1) 
    geometric_series!(vy, G, H, X, T1, tmp1, tmp2)
    vy_1 = view(Y, :, 1)
    vy_1 .= y0
    for t = 2:T1 + 1
        vy = view(Y, :, t)
        mul!(vy, G, vy_1, 1.0, 1.0)
        vy_1 = vy
    end
    for t = T1 + 2:periods
        vy = view(Y, :, t)
        mul!(vy, G, vy_1)
        vy_1 = vy
    end
end

"""
simul_first_order_1!(Y, y0, c, G, H, X, T1, T, tmp1, tmp2) -> Y

computes the solution of a perfect foresight linear model:
``y_t -c = G(y_{t-1} - c) + ∑_{i=0}^{T-1}G^iHx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X`` containing shocks known since period 1. ``c`` is a constant vector.
 ``G`` and ``H`` matrices must have been computed by a linear rational expectation model solver.
"""
function simul_first_order_1!(Y::StridedMatrix{Float64},
                              y0::AbstractVector{Float64},
                              c::AbstractVector{Float64},
                              G::StridedMatrix{Float64},
                              H::StridedMatrix{Float64},
                              X::StridedMatrix{Float64},
                              T1::Int64,
                              T::Int64,
                              temp1::StridedMatrix{Float64},
                              temp2::StridedMatrix{Float64})

    Y .-= c
    y0 .-= c
    simul_first_order_1!(Y, y0, G, H, X, T1, T, temp1, temp2)
    y0 .+= c
    Y .+= c
end

function get_residuals!(residuals::AbstractMatrix{Float64},
                        endogenous::Vector{Float64},
                        exogenous::AbstractMatrix{Float64},
                        dynamic_variables::AbstractVector{Float64},
                        steadystate::AbstractVector{Float64},
                        params::AbstractVector{Float64},
                        m::Model,
                        periods::Int64,
                        temp_vec)
    lli = m.lead_lag_incidence
    dynamic! = m.dynamic!.dynamic!
    for t = 2:periods
        get_dynamic_endogenous_variables!(dynamic_variables,
                                          endogenous, lli, m, t)
        vr = view(residuals, :, t - 1)
        @inbounds Base.invokelatest(dynamic!,
                                    temp_vec,
                                    vr,
                                    dynamic_variables,
                                    exogenous,
                                    params,
                                    steadystate,
                                    t)
    end
end

function solve1(residuals, y, exogenous, dynamic_variables,
                steadystate, params, md, periods, temp_vec,
                JJ, context, ws_threaded, n)
    dy = zeros(n*periods)
    for i = 1:10
        get_residuals!(residuals, vec(y), exogenous, dynamic_variables, steadystate, params, md, periods, temp_vec)
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        dy .= A\vec(residuals)
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end
end

function solve2(residuals, y, exogenous, dynamic_variables,
                steadystate, params, md, periods, temp_vec,
                JJ, context, ws_threaded, n, ws)
    tmp = zeros(n*periods)
    dy = zeros(n*periods)
    for i = 1:10
        get_residuals!(residuals, vec(y), exogenous, dynamic_variables, steadystate, params, md, periods, temp_vec)
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        fill!(tmp, 0.0)
        Dynare.gmres!(tmp, A, residuals, Pr=ws.P, maxiter=8)
        ldiv!(dy, ws.P, tmp)
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end
    0
end

function solve3(τ, residuals, y, exogenous, dynamic_variables,
                steadystate, params, md, periods, temp_vec,
                JJ, context, ws_threaded, n, verbose=false)
    tmp = zeros(n*periods)
    dy = zeros(n*periods)
    for i = 1:10
        get_residuals!(residuals, vec(y), exogenous, dynamic_variables, steadystate, params, md, periods, temp_vec)
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        lu = ilu(A, τ = τ)
        fill!(tmp, 0.0)
        Dynare.gmres!(dy, A, residuals, Pr=lu, verbose = verbose, maxiter=4)
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end
    
end

function solve4(residuals, y, exogenous, dynamic_variables,
                steadystate, params, md, periods, temp_vec,
                JJ, context, ws_threaded, n)
    tmp = zeros(n*periods)
    dy = zeros(n*periods)
    local F::UmfpackLU
    for i = 1:10
        get_residuals!(residuals, vec(y), exogenous, dynamic_variables, steadystate, params, md, periods, temp_vec)
        A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
        F = lu(A)
        ldiv!(dy, F, vec(residuals))
        view(y, n+1:(periods+1)*n) .-= dy
        if norm(residuals) < 1e-10
            break
        end
    end
    
end



