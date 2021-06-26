using BenchmarkTools
import LinearAlgebra
import LinearAlgebra.BLAS: @blasfunc, BlasInt

include("../src/perfectforesight/gmres_solver.jl")

function preconditioner1!(rout::AbstractVector{Float64},
                      rin::AbstractVector{Float64},
                      g::AbstractMatrix{Float64},
                      hh::AbstractMatrix{Float64},
                      preconditioner_window::Int64,
                      periods::Int64)

    m = size(hh, 1)
    mk = m*preconditioner_window
    @inbounds mul!(rout, 1, hh, 1, m, mk, rin, 1)
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir = (i - 1)*m + 1
        mul!(rout, ir, hh, 1, m, mk, rin, ir)
    end
    ir = m + 1
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
    mk -= m
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir = (i - 1)*m + 1
        ir_m = ir - m
        mul!(rout, ir, hh, 1, m, mk, rin, ir)
        mk -= m
    end
    ir = (periods - preconditioner_window + 1)*m + 1 
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
end

# flexible mul!
function preconditioner2!(rout::AbstractVector{Float64},
                      rin::AbstractVector{Float64},
                      g::AbstractMatrix{Float64},
                      hh::AbstractMatrix{Float64},
                      preconditioner_window::Int64,
                      periods::Int64)

    m = size(hh, 1)
    mk = m*preconditioner_window
    @inbounds mul!(rout, 1, hh, 1, m, mk, rin, 1)
    @inbounds @Threads.threads for i = 2:(periods - preconditioner_window + 1)
        ir1 = (i - 1)*m + 1
        mul!(rout, ir1, hh, 1, m, mk, rin, ir1)
    end
    ir = m + 1
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
    @inbounds @Threads.threads for i = periods - preconditioner_window + 2:periods
        ir1 = (i - 1)*m + 1
        mk1 = m*(periods - i + 1)
        mul!(rout, ir1, hh, 1, m, mk1, rin, ir1)
    end
    ir = (periods - preconditioner_window + 1)*m + 1 
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
end

# views
function preconditioner3!(rout::AbstractVector{Float64},
                      rin::AbstractVector{Float64},
                      g::AbstractMatrix{Float64},
                      hh::AbstractMatrix{Float64},
                      preconditioner_window::Int64,
                      periods::Int64)

    m = size(hh, 1)
    mk = m*preconditioner_window
    @inbounds mul!(rout, 1, hh, 1, m, mk, rin, 1)
    @inbounds @Threads.threads for i = 2:(periods - preconditioner_window + 1)
        vrout = view(rout, (i - 1)*m + 1:i*m)
        vrin = view(rin, (i - 1)*m + 1:(i + preconditioner_window - 1)*m)
        mul!(vrout, hh, vrin)
    end
    ir = m + 1
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
    @inbounds @Threads.threads for i = periods - preconditioner_window + 2:periods
        vrout = view(rout, (i - 1)*m + 1: i*m)
        vrin = view(rin, (i - 1)*m + 1: periods*m)
        vhh = view(hh, :, 1:(periods - i + 1)*m)
        mul!(vrout, vhh, vrin)
    end
    ir = (periods - preconditioner_window + 1)*m + 1 
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
end

n = 50
g = randn(n, n)
periods = 100
rin = rand(periods*n)
rout = similar(rin)
preconditioner_window = 3
hh = randn(n,preconditioner_window*n)

@btime preconditioner1!(rout, rin, g, hh, preconditioner_window, periods)
if periods == 6
    target = similar(rout)
    target[1:6] = hh*rin[1:18]
    target[7:12] = g*target[1:6]
    target[7:12] += hh*rin[7:24]
    target[13:18] = g*target[7:12]
    target[13:18] += hh*rin[13:30]
    target[19:24] = g*target[13:18]
    target[19:24] += hh*rin[19:36]
    target[25:30] = g*target[19:24]
    target[25:30] += hh[:,1:12]*rin[25:36]
    target[31:36] = g*target[25:30]
    target[31:36] += hh[:,1:6]*rin[31:36]
    @test rout[1:12] ≈ target[1:12]
    @test rout ≈ target
end

BLAS.set_num_threads(1)

@btime preconditioner2!(rout, rin, g, hh, preconditioner_window, periods)
if periods == 6
    target = similar(rout)
    target[1:6] = hh*rin[1:18]
    target[7:12] = g*target[1:6]
    target[7:12] += hh*rin[7:24]
    target[13:18] = g*target[7:12]
    target[13:18] += hh*rin[13:30]
    target[19:24] = g*target[13:18]
    target[19:24] += hh*rin[19:36]
    target[25:30] = g*target[19:24]
    target[25:30] += hh[:,1:12]*rin[25:36]
    target[31:36] = g*target[25:30]
    target[31:36] += hh[:,1:6]*rin[31:36]
    @test rout[1:12] ≈ target[1:12]
    @test rout ≈ target
end

@btime preconditioner3!(rout, rin, g, hh, preconditioner_window, periods)
if periods == 6
    target = similar(rout)
    target[1:6] = hh*rin[1:18]
    target[7:12] = g*target[1:6]
    target[7:12] += hh*rin[7:24]
    target[13:18] = g*target[7:12]
    target[13:18] += hh*rin[13:30]
    target[19:24] = g*target[13:18]
    target[19:24] += hh*rin[19:36]
    target[25:30] = g*target[19:24]
    target[25:30] += hh[:,1:12]*rin[25:36]
    target[31:36] = g*target[25:30]
    target[31:36] += hh[:,1:6]*rin[31:36]
    @test rout[1:12] ≈ target[1:12]
    @test rout ≈ target
end
