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



function jacobian_time_vec!(y::AbstractVector{Float64},
                            residuals::AbstractVector{Float64},
                            endogenous::AbstractVector{Float64},
                            exogenous::AbstractMatrix{Float64},
                            steadystate::AbstractVector{Float64},
                            g::AbstractMatrix{Float64},
                            m::Model,
                            n::Int64,
                            ws::Vector{JacTimesVec})
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
    ndyn = m.n_bkwrd + m.n_current + m.n_fwrd + 2*m.n_both
    oldthreadnbr = BLAS.get_thread_nbr()
    BLAS.set_num_threads(1)
    @inbounds @Threads.threads for period = 1:n
        k = Threads.threadid()
        get_jacobian!(ws[k], endogenous, exogenous, steadystate, m,
                      period + 1)
        get_dynamic_endogenous_variables!(ws[k].dynamic_variables ,
                                          residuals, lli, m,
                                          period + 1)
        offset_y = (period - 1)*nendo + 1
        @inbounds mul!(y, offset_y, ws[k].jacobian, 1, nendo, ndyn,
                       ws[k].dynamic_variables, 1)
    end
    # setting terminal period according to linear approximation
    offset_y = (n - 1)*nendo + 1
    #select forward looking variables
    k = 1 
    @inbounds for i = 1:nendo
        if lli[3, i] > 0
            mul!(ws[1].temp_vec, k, g, i, 1, nendo, ws[1].dynamic_variables,
                 npred + 1)
            k += 1
        end
    end
    mul!(y, offset_y, ws[1].jacobian, (npred+n_current)*nendo + 1, nendo,
         nfwrd, ws[1].temp_vec, 1, 1.0, 1.0)
    BLAS.set_num_threads(oldthreadnbr)
    return 0
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
    oldthreadnbr = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    @inbounds @Threads.threads for i = 2:(periods - preconditioner_window + 1)
        ir = (i - 1)*m + 1
        mul!(rout, ir, hh, 1, m, mk, rin, ir)
    end
    BLAS.set_num_threads(oldthreadnbr)
    ir = (periods - preconditioner_window + 1)*m + 1 
    @inbounds for i = 2:(periods - preconditioner_window + 1)
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
    end
    BLAS.set_num_threads(1)
    @inbounds @Threads.threads for i = periods - preconditioner_window + 2:periods
        ir1 = (i - 1)*m + 1
        mk1 = m*(periods - i + 1)
        mul!(rout, ir, hh, 1, m, mk1, rin, ir)
    end
    BLAS.set_num_threads(oldthreadnbr)
    ir = (periods - preconditioner_window + 1)*m + 1 
    @inbounds for i = periods - preconditioner_window + 2:periods
        ir_m = ir - m
        mul!(rout, ir, g, 1, m, m, rout, ir_m, 1, 1)
        ir += m
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
    options = context.options
    endogenous = results.trends.endogenous_steady_state
    exogenous = results.trends.exogenous_steady_state
    compute_first_order_solution!(results.linearrationalexpectations,
                                  endogenous, exogenous, endogenous,
                                  model, work, options)
end
