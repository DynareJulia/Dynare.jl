using LinearAlgebra
using FastLapackInterface
using FastLapackInterface.LinSolveAlgo

import LinearAlgebra.BLAS: @blasfunc, BlasInt, BlasFloat, libblas

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

struct JacTimesVec
    jacobian::Matrix{Float64}
    dynamic_variables::Vector{Float64}
    temp_vec::Vector{Float64}
    residuals::Vector{Float64}
    params::Vector{Float64}
    function JacTimesVec(context)
        md = context.models[1]
        work = context.work
        jacobian = similar(work.jacobian)
        dynamic_variables =
            zeros(md.n_bkwrd + md.n_current + md.n_fwrd + 2*md.n_both)
        temp_vec = Vector{Float64}(undef, sum(md.dynamic!.tmp_nbr[1:2]))
        residuals = Vector{Float64}(undef, md.endogenous_nbr)
        params = work.params
        new(jacobian, dynamic_variables,
            temp_vec, residuals, params)
    end
end

function jacobian_time_vec_period!(y::AbstractVector{Float64},
                                   ws::JacTimesVec,
                                   residuals::AbstractVector{Float64},
                                   endogenous::AbstractVector{Float64},
                                   exogenous::AbstractMatrix{Float64},
                                   steadystate::AbstractVector{Float64},
                                   lli::Matrix{Int64},
                                   ndyn::Int64,
                                   md::Model,
                                   period::Int64)
    nvar = md.endogenous_nbr
    get_jacobian!(ws, endogenous, exogenous, steadystate, md,
                  period + 1)
    get_dynamic_endogenous_variables!(ws.dynamic_variables ,
                                      residuals, lli, md,
                                      period + 1)
    offset_y = (period - 1)*nvar + 1
    @inbounds mul!(y, offset_y, ws.jacobian, 1, nvar, ndyn,
                   ws.dynamic_variables, 1)
end

function jacobian_time_vec!(y::AbstractVector{Float64},
                            residuals::AbstractVector{Float64},
                            endogenous::AbstractVector{Float64},
                            exogenous::AbstractMatrix{Float64},
                            steadystate::AbstractVector{Float64},
                            g::AbstractMatrix{Float64},
                            m::Model,
                            periods::Int64,
                            ws::Vector{JacTimesVec})
    npred = m.n_bkwrd + m.n_both
    nfwrd = m.n_fwrd + m.n_both
    n_current = m.n_current
    nvar = m.endogenous_nbr
    lli = m.lead_lag_incidence
    ndyn = m.n_bkwrd + m.n_current + m.n_fwrd + 2*m.n_both
    oldthreadnbr = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    @inbounds @Threads.threads for period = 1:periods
        k = Threads.threadid()
        jacobian_time_vec_period!(y, ws[k], residuals, endogenous,
                                  exogenous, steadystate, lli, ndyn,
                                  m, period)
    end
    BLAS.set_num_threads(oldthreadnbr)
    # setting terminal period according to linear approximation
    offset_y = (periods - 1)*nvar + 1
    #select forward looking variables
    get_dynamic_endogenous_variables!(ws[1].dynamic_variables ,
                                      residuals, lli, m,
                                      periods + 1)
    k = 1 
    @inbounds for i = 1:nvar
        if lli[3, i] > 0
            mul!(ws[1].temp_vec, k, g, i, 1, nvar, ws[1].dynamic_variables,
                 npred + 1)
            k += 1
        end
    end
    @inbounds mul!(y, offset_y, ws[1].jacobian, (npred+n_current)*nvar + 1, nvar,
         nfwrd, ws[1].temp_vec, 1, 1.0, 1.0)
    return 0
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
    oldthreadnbr = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    @inbounds @Threads.threads for i = 2:(periods - preconditioner_window + 1)
        ir1 = (i - 1)*m + 1
        mul!(rout, ir1, hh, 1, m, mk, rin, ir1)
    end
    BLAS.set_num_threads(oldthreadnbr)
    ir = m + 1 
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
    BLAS.set_num_threads(1)
    @inbounds @Threads.threads for i = periods - preconditioner_window + 2:periods
        ir2 = (i - 1)*m + 1
        mk1 = m*(periods - i + 1)
        mul!(rout, ir2, hh, 1, m, mk1, rin, ir2)
    end
    BLAS.set_num_threads(oldthreadnbr)
    ir = (periods - preconditioner_window + 1)*m + 1 
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
    return rout
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
    residuals::Vector{Float64}
    presiduals::Vector{Float64}
    periods::Int64
    endogenous::Vector{Float64}
    exogenous::Matrix{Float64}
    steadystate::Vector{Float64}
    temp_vec::Vector{Float64}
    J::Jacobian
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
        steadystate_exo = Vector{Float64}(undef, m.exogenous_nbr)
        steadystate .= context.results.model_results[1].trends.endogenous_steady_state
        steadystate_exo .= context.results.model_results[1].trends.exogenous_steady_state
        endogenous = repeat(steadystate, periods + 2)
        exogenous = repeat(steadystate_exo',periods + 2)
        lli = m.lead_lag_incidence
        dynamic_variables = zeros(nnz(sparse(lli)))
        temp_vec=Vector{Float64}(undef, n)
        work = context.work
        LREWs = LinearRationalExpectationsWs(
            algo, n, m.exogenous_nbr, m.exogenous_deterministic_nbr,
            m.i_fwrd_b, m.i_current, m.i_bkwrd_b, m.i_both, m.i_static)
        LREresults = LinearRationalExpectationsResults(n,
                                                       m.exogenous_nbr,
                                                       LREWs.backward_nbr)
        options = LinearRationalExpectationsOptions()
        first_order_solver!(LREresults,
                            algo,
                            work.jacobian,
                            options,
                            LREWs)
        g1_1 = LREresults.g1_1
        @inbounds for i = 1:LREWs.backward_nbr
            vg = view(g, :, LREWs.backward_indices[i])
            vg1 = view(g1_1, :, i)
            copy!(vg, vg1)
        end
        if any(isnan.(g))
            throw(ArgumentError("NaN in g"))
        end
        residuals = zeros((periods+2)*n)
        presiduals = zeros(m.n_bkwrd + m.n_current + m.n_fwrd + 2*m.n_both)
        ws_threaded = [JacTimesVec(context) for i=1:Threads.nthreads()]
        get_jacobian!(ws_threaded[1], endogenous, exogenous, steadystate, m, 2)
        n = size(g, 1)
        if any(isnan.(ws_threaded[1].jacobian))
            throw(ArgumentError("NaN in jacobian"))
        end
        @inbounds get_abc!(a, b, c, ws_threaded[1].jacobian, m)
        P = LREprecond(periods, preconditioner_window, a, b, c, g)
        LREMap = LinearMap(periods*n) do C, B
            @inbounds copyto!(residuals, n + 1, B, 1, periods*n)
            jacobian_time_vec!(C, residuals, endogenous, exogenous,
                               steadystate, g, m, periods,
                               ws_threaded)
        end
        J = Jacobian(context, periods)
        new(a, b, c, g, dynamic_variables, residuals, presiduals,
            periods, endogenous, exogenous, steadystate, temp_vec, J, P, LREMap)
    end
end

function gmres_solver!(rout::Vector{Float64}, res::Vector{Float64},
                       periods::Int64, preconditioner_window::Int64,
                       model::Model, work::Work, ws::GmresWs;
                       log=false, verbose=false)
    @inbounds x, h = gmres!(rout, LREMap, res, log=log,
                            verbose=verbose, Pr=ws.P)
    return 0
end

function get_first_order_solution!(context::Context)
    results = context.results.model_results[1]
    model = context.models[1]
    work = context.work
    options = StochSimulOptions(Dict{String, Any}())
    endogenous = results.trends.endogenous_steady_state
    exogenous = results.trends.exogenous_steady_state
    compute_first_order_solution!(results.linearrationalexpectations,
                                  endogenous, exogenous, endogenous,
                                  model, work, options)
end
