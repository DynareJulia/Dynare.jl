using BenchmarkTools
using LinearAlgebra

@show Threads.nthreads()

n = 10
T = 100000

A = randn(n, n)

y = randn(T*n)
z = similar(y)

function f1(y, z, A, n, T)
    @inbounds for i=1:T
        vy = view(y, (i-1)*n+1:i*n)
        vz = view(z, (i-1)*n+1:i*n)
        mul!(vz, A, vy)
    end
end

f1(y, z, A, n, T)
@btime f1(y, z, A, n, T)

function f2(y, z, A, n, T)
    @inbounds Threads.@threads for i=1:T
        vy = view(y, (i-1)*n+1:i*n)
        vz = view(z, (i-1)*n+1:i*n)
        mul!(vz, A, vy)
    end
end

f2(y, z, A, n, T)
@btime f2(y, z, A, n, T)

