using Dynare
using FastLapackInterface
using FastLapackInterface.LinSolveAlgo
using IterativeSolvers
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearMaps
import Base.\
import LinearAlgebra.mul!
import LinearAlgebra.ldiv!
import LinearAlgebra.BLAS: @blasfunc, BlasInt, BlasFloat, libblas
using LinearRationalExpectations
using SparseArrays
using Test

function mul!(y::StridedVector{Float64}, offset_y::Int64,
                  a::StridedMatrix{Float64}, offset_a::Int64,
                  ma::Int64, na::Int64, x::StridedVector{Float64},
                  offset_x::Int64, α::Real, β::Real)
    ccall((@blasfunc(dgemv_), libblas), Cvoid,
          (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'N', ma, na,
          α, Ref(a, offset_a), stride(a, 2),
          Ref(x, offset_x), stride(x, 1), β, Ref(y, offset_y),
          stride(y, 1))
end

function mul!(y::StridedVector{Float64}, offset_y::Int64,
                  a::Transpose{Float64, <: StridedMatrix}, offset_a::Int64,
                  ma::Int64, na::Int64, x::StridedVector{Float64},
                  offset_x::Int64, α::Real, β::Real)
    ccall((@blasfunc(dgemv_), libblas), Cvoid,
          (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
           Ref{BlasInt}),
          'T', ma, na,
          α, Ref(a, offset_a), stride(a, 2),
          Ref(x, offset_x), stride(x, 1), β, Ref(y, offset_y),
          stride(y, 1))
end

mul!(y::StridedVector{Float64}, offset_y::Int64,
     a::StridedMatrix{Float64}, offset_a::Int64, ma::Int64,
     na::Int64, x::StridedVector{Float64}, offset_x::Int64) =
         mul!(y, offset_y, a, offset_a, ma, na, x, offset_x,
              1.0, 0.0)

mul!(y::StridedVector{Float64}, offset_y::Int64,
     a::Transpose{Float64, <: StridedMatrix}, offset_a::Int64, ma::Int64,
     na::Int64, x::StridedVector{Float64}, offset_x::Int64) =
         mul!(y, offset_y, a, offset_a, ma, na, x, offset_x,
              1.0, 0.0)


function get_dynamic_endogenous_variables!(y::AbstractVector{Float64},
                                           data::AbstractVector{Float64},
                                           lli::Matrix{Int64},
                                           m::Model, period::Int64)
    m, n = size(lli)
    p = (period - 2)*n
    @inbounds for i = 1:m
        @simd for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p + j]
            end
        end
        p += n
    end
end

function get_dynamic_endogenous_variables_0!(y::AbstractVector{Float64},
                                             initial_values::AbstractVector{Float64},
                                             data::AbstractVector{Float64},
                                             lli::Matrix{Int64},
                                             m::Model, period::Int64)
    m, n = size(lli)
    for j = 1:n
        k = lli[1, j]
        if k > 0
            y[k] = initial_values[j]
        end
    end
    p = (period - 2)*n
    for i = 2:m
        for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p + j]
            end
        end
        p += 1
    end
end

function get_dynamic_endogenous_variables_T!(y::AbstractVector{Float64},
                                             terminal_values::AbstractVector{Float64},
                                             data::AbstractVector{Float64},
                                             lli::Matrix{Int64},
                                             m::Model, period::Int64)
    m, n = size(lli)
    p = (period - 2)*n
    @inbounds for i = 1:m-1
        @simd for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p + j]
            end
        end
        p += 1
    end
    @inbounds for j = 1:n
        k = lli[3, j]
        if k > 0
            y[k] = terminal_values[j]
        end
    end
end

function get_jacobian!(work::Work,
                       endogenous::AbstractVector{Float64},
                       exogenous::Matrix{Float64},
                       steadystate::Vector{Float64},
                       m::Model,
                       period::Int64)
    dynamic_variables = work.dynamic_variables
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables!(dynamic_variables, endogenous,
                                      lli, m, period)
    compute_jacobian(work, dynamic_variables, exogenous,
                     steadystate, m, period)
end

function get_jacobian_0!(work::Work,
                         initial_values::AbstractVector{Float64},
                         endogenous::AbstractVector{Float64},
                         exogenous::Matrix{Float64},
                         steadystate::Vector{Float64},
                         m::Model,
                         period::Int64)
    dynamic_variables = work.dynamic_variables
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables_0!(dynamic_variables,
                                        initial_values,
                                        endogenous,
                                        lli, m, period)
    dynamic! = m.dynamic!.dynamic!
    temporary_values = work.temporary_values
    residuals = work.residuals
    jacobian = work.jacobian
    params = work.params
    compute_jacobian(work, dynamic_variables, exogenous,
                     steadystate, m, period)
end

function get_jacobian_T!(work::Work,
                         terminalal_values::AbstractVector{Float64},
                         endogenous::AbstractVector{Float64},
                         exogenous::Matrix{Float64},
                         steadystate::Vector{Float64},
                         m::Model,
                         period::Int64)
    dynamic_variables = work.dynamic_variables
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables_T!(dynamic_variables,
                                        terminal_values,
                                        endogenous,
                                        lli, m, period)
    compute_jacobian(work, dynamic_variables, exogenous,
                          steady_state, m, period)
end

function compute_jacobian(work::Work,
                          dynamic_variables::AbstractVector{Float64},
                          exogenous::AbstractMatrix{Float64},
                          steadystate::AbstractVector{Float64},
                          m::Model,
                          period::Int64)
    dynamic! = m.dynamic!.dynamic!
    temporary_values = work.temporary_values
    residuals = work.residuals
    jacobian = work.jacobian
    params = work.params
    @inbounds Base.invokelatest(dynamic!,
                                temporary_values,
                                residuals,
                                jacobian,
                                dynamic_variables,
                                exogenous,
                                params,
                                steadystate,
                                period)  
end

function get_direction!(y::AbstractVector{Float64},
                        work::Dynare.Work,
                        residuals::AbstractVector{Float64},
                        endogenous::AbstractVector{Float64},
                        exogenous::AbstractMatrix{Float64},
                        steadystate::AbstractVector{Float64},
                        presiduals::AbstractVector{Float64},
                        m::Dynare.Model,
                        n::Int64)
    dynamic_variables = work.dynamic_variables
    ndyn = length(dynamic_variables)
    npred = m.n_bkwrd + m.n_both
    nfwrd = m.n_fwrd + m.n_both
    nendo = m.endogenous_nbr
    lli = m.lead_lag_incidence
    offset_y = 1
    @inbounds for period = 1:n
        get_jacobian!(work, endogenous, exogenous, steadystate, m,
                      period + 1)
        get_dynamic_endogenous_variables!(presiduals, residuals,
                                          lli, m, period + 1)
        mul!(y, offset_y, work.jacobian, 1, nendo, ndyn,
             presiduals, 1)
        offset_y += nendo
    end
end

function fan_columns!(y::Matrix{Float64}, x::Matrix{Float64}, columns::Vector{Int64},
                      offset_x::Int64)
    my, ny = size(y)
    mx, nx = size(x)
    source = offset_x*mx + 1
    @inbounds for i = 1:length(columns)
        k = columns[i]
        destination = (k - 1)*my + 1
        copyto!(y, destination, x, source, mx)
        source += mx
    end
end

function get_abc!(a::Matrix{Float64}, b::Matrix, c, jacobian, m)
    i_bkwrd = m.i_bkwrd_b
    i_current = m.i_current
    i_fwrd = m.i_fwrd_b
    fan_columns!(a, jacobian, m.i_bkwrd_b, 0)
    offset = m.n_bkwrd + m.n_both
    fan_columns!(b, jacobian, m.i_current, offset)
    offset += m.n_current
    fan_columns!(c, jacobian, m.i_fwrd_b, offset)
end

"""
    function h0!(h0::Matrix{Float64}, a::Matrix{Float64}, b::Matrix{Float64}, c::Matrix{Float64})

computes h0 = inv(a+b*c)
"""
function h0!(h0::AbstractMatrix{Float64}, a::AbstractMatrix{Float64},
             b::AbstractMatrix{Float64}, c::AbstractMatrix{Float64},
             work::AbstractMatrix{Float64}, linsolve_ws)
    copy!(work, a)
    mul!(work, b, c, 1.0, 1.0)
    fill!(h0, 0.0)
    n = size(h0, 1)
    m = 1
    @inbounds for i = 1:n
        h0[m] = 1.0
        m += n + 1
    end
    linsolve_core!(work, h0, linsolve_ws)
end

"""
function hh!(hh::AbstractMatrix{Float64}, h::AbstractMatrix{Float64}, 
             f::AbstractMatrix{Float64}, hf::AbstractMatrix{Float64},
             n::Int64, work1::AbstractMatrix{Float64},
             work2::AbstractMatrix{Float64})

computes hh = [h, h(1), h(2), ..., h(n-1)] and hh(i) = -h*f*h(-1)
"""
function hh!(hh::AbstractMatrix{Float64}, h::AbstractMatrix{Float64},
             f::AbstractMatrix{Float64}, hf::AbstractMatrix{Float64},
             n::Int64, work1::AbstractMatrix{Float64},
             work2::AbstractMatrix{Float64})
    m = size(h, 1)
    m2 = m*m
    @inbounds mul!(hf, h, f, -1, 0)
    @inbounds copy!(work1, h)
    k = 1
    @inbounds for i = 1:n
        copyto!(hh, k, work1, 1, m2)
        if i < n
            mul!(work2, hf, work1)
            copy!(work1, work2)
            k += m2
        end
    end
end

function preconditioner!(rout::AbstractVector{Float64},
                      rin::AbstractVector{Float64},
                      g::AbstractMatrix{Float64},
                      hh::AbstractMatrix{Float64},
                      k::Int64,
                      n::Int64)

    m = size(hh, 1)
    mk = m*k
    @inbounds mul!(rout, 1, hh, 1, m, mk, rin, 1)
    ir = m + 1
    m2 = m*m
    ihh = m2 + 1
    @inbounds for i = 2:(n - k + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m)
        mul!(rout, ir, hh, 1, m, mk, rin, ir, 1, 1)
        ir += m
        
    end
    mk -= m
    @inbounds for i = n - k + 2:n
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m)
        mul!(rout, ir, hh, 1, m, mk, rin, ir, 1, 1)
        ir += m
        mk -= m
    end
end

function makeA(jacobian::Matrix{Float64},
               n::Int64)
    i, j, v = findnz(jacobian)
    nvar = size(jacobian, 1)
    m = length(i)
    nm = n*m - count(j .<= nvar) - count(j .> 2*nvar)
    i1 = zeros(Int64, nm)
    j1 = zeros(Int64, nm)
    v1 = zeros(nm)
    r = 1
    for el = 1:m
        if j[el] > nvar
            i1[r] = i[el]
            j1[r] = j[el] - nvar
            v1[r] = v[el]
            r += 1
        end
    end        
    offset = nvar
    for k = 2:(n - 1)
        for el = 1:m
            i1[r] = i[el] + offset
            j1[r] = j[el] + offset - nvar
            v1[r] = v[el]
            r += 1
        end
        offset += nvar
    end
    for el = 1:m
        if j[el] <= 2*nvar
            i1[r] = i[el] + offset
            j1[r] = j[el] + offset - nvar
            v1[r] = v[el]
            r += 1
        end
    end
    return sparse(i1, j1, v1)
end

struct LREprecond
    y::AbstractVector{Float64}
    g::AbstractMatrix{Float64}
    h::AbstractMatrix{Float64}
    hf::AbstractMatrix{Float64}
    hh::AbstractMatrix{Float64}
    work1::AbstractMatrix{Float64}
    work2::AbstractMatrix{Float64}
    k::Int64
    n::Int64
    linsolve_ws::LinSolveWs
    function LREprecond(k::Int64,
                     n::Int64,
                     a::AbstractMatrix{Float64},
                     b::AbstractMatrix{Float64},
                     c::AbstractMatrix{Float64},
                     g::AbstractMatrix{Float64})
        m = size(g, 1)
        mn = m*n
        h = Matrix{Float64}(undef, m, m)
        hf = similar(h)
        hh = Matrix{Float64}(undef, m, mn)
        work1 = similar(h)
        work2 = similar(h)
        linsolve_ws = LinSolveWs(n)
        copy!(work2, b)
        h0!(h, work2, c, g, work1, linsolve_ws)
        hh!(hh, h, c, hf, n, work1, work2)
        y = Vector{Float64}(undef, mn)
        new(y, g, h, hf, hh, work1, work2, k, n, linsolve_ws)
    end
end

function ldiv!(y::AbstractVector{Float64}, P::LREprecond, x::AbstractVector{Float64})
    preconditioner!(y, x, P.g, P.hh, P.k, P.n)
end

function ldiv!(P::LREprecond, x::AbstractVector{Float64})
    preconditioner!(P.y, x, P.g, P.hh, P.k, P.n)
    x .= P.y
end

(\)(P::LREprecond, x::AbstractVector{Float64}) = ldiv!(P, copy(x))

function gmres_solver!(endogenous, exogenous, steadystate,
                       preconditioner_window, periods, work)

    LRE = LinearMap(k*n) do C, B
        copyto!(work.residuals, n + 1, B, 1, k*n) 
        get_direction!(C, work, residuals, endogenous, exogenous,
                       steadystate, work.presiduals, md, k)
    end

    x, h = gmres!(rout, LRE, res, log=true, verbose=true, Pl=P)
end
