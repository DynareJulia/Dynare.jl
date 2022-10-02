using BenchmarkTools
using ForwardDiff

struct DiffCache{T<:AbstractArray,S<:AbstractArray}
    du::T
    dual_du::S
end

function DiffCache(u::AbstractArray{T}, siz, chunk_sizes) where {T}
    x = zeros(T, prod(chunk_sizes .+ 1) * prod(siz))
    DiffCache(u, x)
end

"""

`dualcache(u::AbstractArray, N::Int = ForwardDiff.pickchunksize(length(u)); levels::Int = 1)`
`dualcache(u::AbstractArray; N::AbstractArray{<:Int})`

Builds a `DualCache` object that stores both a version of the cache for `u`
and for the `Dual` version of `u`, allowing use of pre-cached vectors with
forward-mode automatic differentiation. Supports nested AD via keyword `levels`
or specifying an array of chunk_sizes.

"""
function dualcache(
    u::AbstractArray,
    N::Int = ForwardDiff.pickchunksize(length(u));
    levels::Int = 1,
)
    DiffCache(u, size(u), N * ones(Int, levels))
end
dualcache(u::AbstractArray, N::AbstractArray{<:Int}) = DiffCache(u, size(u), N)
function dualcache(u::AbstractArray, ::Type{Val{N}}; levels::Int = 1) where {N}
    dualcache(u, N; levels)
end
dualcache(u::AbstractArray, ::Val{N}; levels::Int = 1) where {N} = dualcache(u, N; levels)

"""

`get_tmp(dc::DiffCache, u)`

Returns the `Dual` or normal cache array stored in `dc` based on the type of `u`.

"""
function get_tmp(dc::DiffCache, u::T) where {T<:ForwardDiff.Dual}
    reinterpret(T, dc.dual_du)
end

function get_tmp(dc::DiffCache, u::AbstractArray{T}) where {T<:ForwardDiff.Dual}
    reinterpret(T, dc.dual_du)
end

get_tmp(dc::DiffCache, u::Number) = dc.du
get_tmp(dc::DiffCache, u::AbstractArray) = dc.du


function f(r, y, x)
    r = get_tmp(r, y)
    r[1] = y[1] * y[2]
    r[2] = y[1] + y[2] * x[1]
    return r
end

y = rand(2)
r = similar(y)
rd = dualcache(r)

@btime ForwardDiff.jacobian(y -> f(rd, y, [3]), ones(2))
