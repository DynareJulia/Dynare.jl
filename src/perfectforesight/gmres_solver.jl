using BenchmarkTools
using Dynare
using FastLapackInterface
using FastLapackInterface.LinSolveAlgo
using IterativeSolvers
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearMaps
using LinearRationalExpectations
import Base.\
import LinearAlgebra.mul!
import LinearAlgebra.ldiv!
import LinearAlgebra.BLAS: @blasfunc, BlasInt, BlasFloat, libblas
using LinearRationalExpectations
using SparseArrays

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
        for j = 1:n
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
    @inbounds for j = 1:n
        k = lli[1, j]
        if k > 0
            y[k] = initial_values[j]
        end
    end
    p = (period - 2)*n
    @inbounds for i = 2:m
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
        for j = 1:n
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
                         initial_values,
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
    fill!(jacobian, 0.0)
    Base.invokelatest(dynamic!,
                      temporary_values,
                      residuals,
                      jacobian,
                      dynamic_variables,
                      exogenous,
                      params,
                      steadystate,
                      period)
end

function jacobian_time_vec!(y::AbstractVector{Float64},
                            dynamic_variables::AbstractVector{Float64},
                            residuals::AbstractVector{Float64},
                            endogenous::AbstractVector{Float64},
                            exogenous::AbstractMatrix{Float64},
                            steadystate::AbstractVector{Float64},
                            presiduals::AbstractVector{Float64},
                            g::AbstractMatrix{Float64},
                            temp_vec::AbstractVector{Float64},
                            work::Work,
                            m::Model,
                            n::Int64)
    #=
    y .= NaN
    if any(isnan.(residuals))
        throw(ArgumentError("contains NaN in residuals"))
    end
    =#
    npred = m.n_bkwrd + m.n_both
    nfwrd = m.n_fwrd + m.n_both
    n_current = m.n_current
    nendo = m.endogenous_nbr
    lli = m.lead_lag_incidence
    offset_y = 1
    @inbounds for period = 1:n
        get_jacobian!(work, endogenous, exogenous, steadystate, m,
                      period + 1)
#=
        if any(isnan.(work.residuals))
            @show "Periods $period"
            @show endogenous
            @show exogenous
            @show work.jacobian
            @show work.residuals
            throw(ArgumentError("contains NaN in residuals"))
        end
=# 
        get_dynamic_endogenous_variables!(presiduals, residuals,
                                          lli, m, period + 1)

        ndyn = length(presiduals)
        @inbounds mul!(y, offset_y, work.jacobian, 1, nendo, ndyn,
                       presiduals, 1)
#=        
        if any(isnan.(y[1:period*nendo]))
            @show period
            @show y[(period-1)*nendo + 1:period*nendo]
            display(work.jacobian)
            println(' ')
            @show presiduals
            @show endogenous[1:18]
            @show m.params
            throw(ArgumentError("contains NaN in y"))
        end
=#
        offset_y += nendo
    end
    # setting terminal period according to linear approximation
    offset_y -= nendo
    #select forward looking variables
    k = 1 
    @inbounds for i = 1:nendo
        if lli[3, i] > 0
            mul!(temp_vec, k, g, i, 1, nendo, presiduals, npred + 1)
            k += 1
        end
    end
    mul!(y, offset_y, work.jacobian, (npred+n_current)*nendo + 1, nendo, nfwrd,
         temp_vec, 1, 1.0, 1.0)
    #=
    if any(isnan.(y))
    throw(ArgumentError("contains NaN in whole y"))
    end
    =#
    return 0
end

function makeA0(jacobian::AbstractMatrix{Float64},
               g::AbstractMatrix{Float64},
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
    A = sparse(i1, j1, v1, n*nvar, n*nvar)
    # terminal condition
    k = (n-1)*nvar + 1:n*nvar
    A[k,k] .+= jacobian[:,2*nvar+1:end]*g 
    return A
end

function makeA(jacobian::AbstractMatrix{Float64},
               g::AbstractMatrix{Float64},
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
    A = sparse(i1, j1, v1, n*nvar, n*nvar)
    # terminal condition
    k = (n-1)*nvar + 1:n*nvar
    A[k,k] .+= jacobian[:,2*nvar+1:end]*g 
    return A
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

function get_abc!(a::Matrix{Float64}, b::Matrix, c::Matrix{Float64}, jacobian::Matrix{Float64}, m::Model)
    fill!(a, 0.0)
    fill!(b, 0.0)
    fill!(c, 0.0)
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
    @inbounds copy!(work, a)
    @inbounds mul!(work, b, c, 1.0, 1.0)
    fill!(h0, 0.0)
    n = size(h0, 1)
    m = 1
    @inbounds for i = 1:n
        h0[m] = 1.0
        m += n + 1
    end
    @inbounds linsolve_core!(work, h0, linsolve_ws)
    if any(isnan.(h0))
        if any(isnan.(a))
            throw(ArgumentError("NaN in a"))
        end
        if any(isnan.(b))
            throw(ArgumentError("NaN in b"))
        end
        if any(isnan.(c))
            throw(ArgumentError("NaN in c"))
        end
        throw(ArgumentError("NaN in h0"))
    end
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
             preconditioner_window::Int64, work1::AbstractMatrix{Float64},
             work2::AbstractMatrix{Float64})
    if any(isnan.(h))
        throw(ArgumentError("NaN in h"))
    end
    if any(isnan.(f))
        throw(ArgumentError("NaN in f"))
    end
    m = size(h, 1)
    m2 = m*m
    @inbounds mul!(hf, h, f, -1, 0)
    @inbounds copy!(work1, h)
    k = 1
    @inbounds for i = 1:preconditioner_window
        copyto!(hh, k, work1, 1, m2)
        if i < preconditioner_window
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
                      preconditioner_window::Int64,
                      periods::Int64)

    m = size(hh, 1)
    mk = m*preconditioner_window
    @inbounds mul!(rout, 1, hh, 1, m, mk, rin, 1)
    ir = m + 1
    m2 = m*m
    ihh = m2 + 1
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m)
        mul!(rout, ir, hh, 1, m, mk, rin, ir, 1, 1)
        ir += m
    end
    mk -= m
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m)
        mul!(rout, ir, hh, 1, m, mk, rin, ir, 1, 1)
        ir += m
        mk -= m
    end
end

struct LREprecond
    y::AbstractVector{Float64}
    g::AbstractMatrix{Float64}
    h::AbstractMatrix{Float64}
    hf::AbstractMatrix{Float64}
    hh::AbstractMatrix{Float64}
    work1::AbstractMatrix{Float64}
    work2::AbstractMatrix{Float64}
    preconditioner_window::Int64
    periods::Int64
    linsolve_ws::LinSolveWs
    function LREprecond(periods::Int64,
                        preconditioner_window::Int64,
                        a::AbstractMatrix{Float64},
                        b::AbstractMatrix{Float64},
                        c::AbstractMatrix{Float64},
                        g::AbstractMatrix{Float64})
        m = size(g, 1)
        mp = m*preconditioner_window
        h = Matrix{Float64}(undef, m, m)
        hf = similar(h)
        hh = Matrix{Float64}(undef, m, mp)
        work1 = similar(h)
        work2 = similar(h)
        linsolve_ws = LinSolveWs(m)
        @inbounds copy!(work2, b)
        h0!(h, work2, c, g, work1, linsolve_ws)
        hh!(hh, h, c, hf, preconditioner_window, work1, work2)
        y = zeros(m*periods)
        new(y, g, h, hf, hh, work1, work2, preconditioner_window, periods, linsolve_ws)
    end
end

function ldiv!(y::AbstractVector{Float64}, P::LREprecond, x::AbstractVector{Float64})
    preconditioner!(y, x, P.g, P.hh, P.preconditioner_window, P.periods)
end

function ldiv!(P::LREprecond, x::AbstractVector{Float64})
    preconditioner!(P.y, x, P.g, P.hh, P.preconditioner_window, P.periods)
    return P.y
end

(\)(P::LREprecond, x::AbstractVector{Float64}) = ldiv!(P, copy(x))

struct GmresWs
    a::Matrix{Float64}
    b::Matrix{Float64}
    c::Matrix{Float64}
    g::Matrix{Float64}
    dynamic_variables::Vector{Float64}
    jacobian::Matrix{Float64}
    residuals::Vector{Float64}
    presiduals::Vector{Float64}
    periods::Int64
    endogenous::Vector{Float64}
    exogenous::Matrix{Float64}
    steadystate::Vector{Float64}
    temp_vec::Vector{Float64}
    P::LREprecond
    LREMap::LinearMap
    function GmresWs(periods::Int64, preconditioner_window::Int64,
                     context::Context, algo::String)
        m = context.models[1]
        n = m.endogenous_nbr
        a = Matrix{Float64}(undef, n, n)
        b = Matrix{Float64}(undef, n, n)
        c = Matrix{Float64}(undef, n, n)
        g = zeros(n, n)
        steadystate = Vector{Float64}(undef, n)
        steadystate .= context.results.model_results[1].trends.endogenous_steady_state
        endogenous = repeat(steadystate, periods + 2)
        exogenous = zeros(periods + 2, m.exogenous_nbr)
        lli = m.lead_lag_incidence
        dynamic_variables = zeros(nnz(sparse(lli)))
        steadystate = context.results.model_results[1].trends.endogenous_steady_state
        temp_vec=Vector{Float64}(undef, n)
        work = context.work
        LREWs = LinearRationalExpectationsWs(
            algo, n, m.exogenous_nbr, m.exogenous_deterministic_nbr,
            m.i_fwrd_b, m.i_current, m.i_bkwrd_b, m.i_both, m.i_static)
        LREresults = LinearRationalExpectationsResults(n,
                                                       m.exogenous_nbr,
                                                       LREWs.backward_nbr)
        options = context.options["stoch_simul"]
        if algo == "GS"
            options["generalized_schur"]["criterium"] = 1 + 1e-6
        else
            options["cyclic_reduction"] = Dict(["tol" => 1e-8])
        end

        first_order_solver!(LREresults,
                            algo,
                            work.jacobian,
                            options,
                            LREWs)
        @inbounds for i = 1:LREWs.backward_nbr
            k = LREWs.backward_indices[i]
            for j = 1:n
                g[j, k] = LREresults.g1_1[j, i]
            end
        end
        if any(isnan.(g))
            throw(ArgumentError("NaN in g"))
        end
        residuals = zeros((periods+2)*n)
        presiduals = zeros(m.n_bkwrd + m.n_current + m.n_fwrd + 2*m.n_both)
        get_jacobian!(work, endogenous, exogenous, steadystate, m, 2)
        jacobian = Matrix{Float64}(undef, size(work.jacobian))
        copy!(jacobian, work.jacobian)
        n = size(g, 1)
        if any(isnan.(jacobian))
            throw(ArgumentError("NaN in jacobian"))
        end
        @inbounds get_abc!(a, b, c, jacobian, m)
        P = LREprecond(periods, preconditioner_window, a, b, c, g)
        LREMap = LinearMap(periods*n) do C, B
            @inbounds copyto!(residuals, n + 1, B, 1, periods*n)
            jacobian_time_vec!(C, dynamic_variables, residuals,
                               endogenous, exogenous,
                               steadystate, presiduals, g,
                               temp_vec, work, m, periods)
        end
        new(a, b, c, g, dynamic_variables, jacobian, residuals, presiduals,
            periods, endogenous, exogenous, steadystate, temp_vec, P, LREMap)
    end
end

function gmres_solver!(rout::Vector{Float64}, res::Vector{Float64},
                       periods::Int64, preconditioner_window::Int64,
                       model::Model, work::Work, ws::GmresWs;
                       log=false, verbose=false)
    @inbounds x, h = gmres!(rout, LREMap, res, log=log,
                            verbose=verbose, Pr=P)
    @show x
    @show h
    return 0
end

function gmres_solver_test!(rout::Vector{Float64}, res::Vector{Float64},
                            periods::Int64, preconditioner_window::Int64,
                            model::Model, work::Work, ws::GmresWs;
                            log=false, verbose=false)
    n = size(ws.g, 1)
    @inbounds get_abc!(ws.a, ws.b, ws.c, ws.jacobian, model)
    
    P = LREprecond(periods, preconditioner_window, ws.a, ws.b, ws.c, ws.g)
    ldiv!(P, res)
    LREMap = LinearMap(periods*n) do C, B
        @inbounds copyto!(ws.residuals, n + 1, B, 1, periods*n) 
        jacobian_time_vec!(C, ws.dynamic_variables, ws.residuals, ws.endogenous, ws.exogenous,
                           ws.steadystate, ws.presiduals, work, model, periods)
    end
    @time x, h = gmres!(rout, LREMap, res, log=log, verbose=verbose, Pl=P)
    @time x, h = gmres!(rout, LREMap, res, log=log, verbose=verbose, Pl=P)
    @show norm(LREMap*rout - res)
    A = SparseArrays.sparse(LREMap)
    @time A = SparseArrays.sparse(LREMap)
    A\res
    @time A\res
    return 0
end

function get_first_order_solution!(context::Context)
    results = context.results.model_results[1]
    model = context.models[1]
    work = context.work
    options = context.options
    endogenous = results.trends.endogenous_steady_state
    exogenous = results.trends.exogenous_steady_state
    compute_first_order_solution!(results.linearrationalexpectations,
                                  endogenous, exogenous, endogenous,
                                  model, work, options)
end
