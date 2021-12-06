"""
geometric_series!(Y, A, B, X) -> Y

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
simul_first_order!(Y, y0, G, H, X) -> Y

computes the solution of a perfect foresight linear model:
``y_t = Gy_{t-1} + ∑_{i=0}^{T-1}G^iHx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X`` containing shocks known at the beginning of the simulation.
vector ``y0`` contains the initial values of the simulation and ``y_0 = y0``.
Matrix ``Y`` doesn't contain the initial values.
 ``G`` and ``H`` matrices must have been computed by a linear rational expectation model solver.
"""
function simul_first_order!(
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
simul_first_order!(Y, y0, c, G, H, X) -> Y

computes the solution of a perfect foresight linear model:
``y_t -c = G(y_{t-1} - c) + ∑_{i=0}^{T-1}G^iHx_{t-1}`` for ``t=1,…,T``
where ``y_t`` are the columns of matrix ``Y`` and 
    ``x_t`` the columns of matrix ``X`` containing shocks known at the beginning of the simulation.
 ``c`` is a constant vector.
vector ``y0`` contains the initial values of the simulation and ``y_0 = y0``.
Matrix ``Y`` doesn't contain the initial values.
 ``G`` and ``H`` matrices must have been computed by a linear rational expectation model solver.
"""
function simul_first_order!(
    Y::StridedMatrix{Float64},
    y0::AbstractVector{Float64},
    c::AbstractVector{Float64},
    G::StridedMatrix{Float64},
    H::StridedMatrix{Float64},
    X::StridedMatrix{Float64},
)
    y0 .-= c
    simul_first_order!(Y, y0, G, H, X)
    y0 .+= c
    Y .+= c
end

"""
    function simul_first_order!(Y, y0, x, c, A, B, T)

simulates the linear model
```
    y_t - c = A (y_{t-1} - c) + Bx_t
```
for t = 1, ...,T
""" 
function simul_first_order!(
    Y::AbstractMatrix{Float64},
    y0::AbstractVector{Float64},
    x::AbstractVecOrMat{Float64},
    c::AbstractVector{Float64},
    A::StridedVecOrMat{Float64},
    B::StridedVecOrMat{Float64},
    T::Int64,
)
    mul!(view(Y, 2:T+1, :), view(x, 2:T+1, :), transpose(B))
    r_1 = view(Y, 1, :)
    r_1 .= y0 .- c
    for t = 2:T+1
        r = view(Y, t, :)
        mul!(r, A, r_1, 1.0, 1.0)
        r_1 .+= c
        r_1 = r
    end
    r_1 .+= c
end
