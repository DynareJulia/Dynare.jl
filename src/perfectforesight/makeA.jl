using SparseArrays
import SparseArrays

struct SparseStorage{Tv,Ti<:Integer}
    I::Vector{Ti}
    J::Vector{Ti}
    V::Vector{Tv}
    klasttouch::Vector{Ti}
    csrrowptr::Vector{Ti}
    csrcolval::Vector{Ti}
    csrnzval::Vector{Tv}
    csccolptr::Vector{Ti}
    function SparseStorage{Tv,Ti}(
        m::Integer,
        n::Integer,
        nz::Integer,
    ) where {Tv,Ti<:Integer}
        I = Vector{Ti}(undef, nz)
        J = Vector{Ti}(undef, nz)
        V = Vector{Tv}(undef, nz)
        klasttouch = Vector{Ti}(undef, n)
        csrrowptr = Vector{Ti}(undef, m + 1)
        csrcolval = Vector{Ti}(undef, nz)
        csrnzval = Vector{Tv}(undef, nz)
        csccolptr = Vector{Ti}(undef, n + 1)
        new(I, J, V, klasttouch, csrrowptr, csrcolval, csrnzval, csccolptr)
    end
end

SparseStorage(m::Integer, n::Integer, nz::Integer) = SparseStorage{Float64,Int64}(m, n, nz)

function SparseStorage(
    m::Integer,
    n::Integer,
    nz::Integer,
    I::Vector{Ti},
    J::Vector{Ti},
    V::Vector{Tv},
) where {Tv,Ti}
    nz = length(I)
    klasttouch = Vector{Ti}(undef, n)
    csrrowptr = Vector{Ti}(undef, m + 1)
    csrcolval = Vector{Ti}(undef, nz)
    csrnzval = Vector{Tv}(undef, nz)
    csccolptr = Vector{Ti}(undef, n + 1)
    SparseStorage(I, J, V, klasstouch, csrrowptr, csrcolval, csrnzval, csccolptr)
end

SparseStorage(
    m::Integer,
    n::Integer,
    nz::Integer,
    I::Vector{Int64},
    J::Vector{Int64},
    V::Vector{Float64},
) = SparseStorage{Float64,Int64}(
    m::Integer,
    n::Integer,
    nz::Integer,
    I::Vector{Int64},
    J::Vector{Int64},
    V::Vector{Float64},
)

sparse!(m::Integer, n::Integer, SS::SparseStorage) = sparse!(
    SS.I,
    SS.J,
    SS.V,
    m,
    n,
    (x, y) -> x + y,
    SS.klasttouch,
    SS.csrrowptr,
    SS.csrcolval,
    SS.csrnzval,
    SS.csccolptr,
    SS.J,
    SS.V,
)

struct Jacobian
    nrow::Int64
    ss::SparseStorage
    nzval::Vector{Float64}
    maxcol::Int64
    steadystate::Vector{Float64}
    tmp_nvar_npred::Matrix{Float64}
    tmp_nvar_nfwrd::Matrix{Float64}
    tmp_nfwrd_npred::Matrix{Float64}
    function Jacobian(context, periods)
        md = context.models[1]
        nvar = md.endogenous_nbr
        steadystate = context.results.model_results[1].trends.endogenous_steady_state
        steadystate_exo = context.results.model_results[1].trends.exogenous_steady_state
        work = context.work
        endogenous = repeat(steadystate, 3)
        exogenous = repeat(steadystate_exo', 2)
        maxcol = md.n_bkwrd + md.n_current + md.n_fwrd + 2 * md.n_both
        nz = periods * md.NNZDerivatives[1]
        nrow = periods * md.endogenous_nbr
        ss = SparseStorage(nrow, nrow, nz)
        @inbounds fill!(ss.I, 1)
        @inbounds fill!(ss.J, 1)
        @inbounds fill!(ss.V, 0.0)
        nzval = Vector{Float64}(undef, nz)
        npred = md.n_bkwrd + md.n_both
        nfwrd = md.n_fwrd + md.n_both
        tmp_nvar_npred = zeros(nvar, npred)
        tmp_nvar_nfwrd = zeros(nvar, nfwrd)
        tmp_nfwrd_npred = zeros(nfwrd, npred)
        new(
            nrow,
            ss,
            nzval,
            maxcol,
            steadystate,
            tmp_nvar_npred,
            tmp_nvar_nfwrd,
            tmp_nfwrd_npred,
        )
    end
end

struct PeriodJacobianWs
    jacobian::Matrix{Float64}
    dynamic_variables::Vector{Float64}
    temp_vec::Vector{Float64}
    residuals::Vector{Float64}
    params::Vector{Float64}
    function PeriodJacobianWs(context)
        md = context.models[1]
        work = context.work
        jacobian = similar(work.jacobian)
        dynamic_variables = zeros(md.n_bkwrd + md.n_current + md.n_fwrd + 2 * md.n_both)
        temp_vec = Vector{Float64}(undef, sum(md.dynamic!.tmp_nbr[1:2]))
        residuals = Vector{Float64}(undef, md.endogenous_nbr)
        params = work.params
        new(jacobian, dynamic_variables, temp_vec, residuals, params)
    end
end


#=
function get_dynamic_endogenous_variables!(y::AbstractVector{Float64},
                                           data::AbstractVector{Float64},
                                           lli::Matrix{Int64},
                                           m::Model,
                                           period::Int64)
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


function get_jacobian!(ws::JacTimesVec,
                       endogenous::AbstractVector{Float64},
                       exogenous::Matrix{Float64},
                       steadystate::Vector{Float64},
                       m::Model,
                       period::Int64)
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables!(ws.dynamic_variables,
                                      endogenous, lli, m, period)
    compute_jacobian(ws, exogenous,
                     steadystate, m, period)
end


function compute_jacobian(ws::JacTimesVec,
                          exogenous::AbstractMatrix{Float64},
                          steadystate::AbstractVector{Float64},
                          m::Model,
                          period::Int64)
    dynamic! = m.dynamic!.dynamic!
    fill!(ws.jacobian, 0.0)
    @inbounds Base.invokelatest(dynamic!,
                      ws.temp_vec,
                      ws.residuals,
                      ws.jacobian,
                      ws.dynamic_variables,
                      exogenous,
                      ws.params,
                      steadystate,
                      period)
end
=#

function permute_row(ir, permutations)
    for p in permutations
        ir == p[1] && return p[2]
        ir == p[2] && return p[1]
    end
    return ir
end

#=
function make_one_period!(
    JA::Jacobian,
    params::AbstractVector{Float64},
    endogenous::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    periods::Int64,
    md::Model,
    jacobian_columns::Vector{Int64},
    nnz_period::Int64,
    maxcol::Int64,
    ws::DynamicWs,
    t::Int64;
    permutations::Tuple{Int64, Int64} = Tuple{Int64, Int64}[]
)
    nvar = md.endogenous_nbr
    oc = (t - 1) * md.endogenous_nbr
    jacobian = get_dynamic_jacobian!(ws, params, endogenous, exogenous, JA.steadystate, md, t)
    i, j, v = findnz(sparse(jacobian))
    r1 = (t - 1) * nnz_period + 1
    @inbounds for el = 1:length(i)
        if j[el] <= maxcol
            kjel = jacobian_columns[j[el]]
            JA.ss.I[r1] = permute_row(i[el], permutations) + oc
            JA.ss.J[r1] = kjel + oc - nvar
            JA.ss.V[r1] = v[el]
            r1 += 1
        end
    end
end
=#


                            
    
function makeJacobian!(
    JA::Jacobian,
    endogenous::AbstractVector{Float64},
    initialvalues::AbstractVector{Float64},
    terminalvalues::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    context::Context,
    periods::Int64,
    ws::Vector{DynamicWs};
    permutations::Vector{Tuple{Int64,Int64}} = Tuple{Int64,Int64}[],
)
    I = JA.ss.I
    J = JA.ss.J
    V = JA.ss.V
    klasttouch = JA.ss.klasttouch
    csrcolval = JA.ss.csrcolval
    csrnzval = JA.ss.csrnzval
    csrrowptr = JA.ss.csrrowptr
    csccolptr = JA.ss.csccolptr
    md = context.models[1]
    df = context.dynarefunctions
    steadystate = JA.steadystate
    maxcol = JA.maxcol
    nvar = md.endogenous_nbr
    npred = md.n_bkwrd + md.n_both
    nnz_period = md.NNZDerivatives[1]
    lengthI = periods * nnz_period
    resize!(I, lengthI)
    resize!(J, lengthI)
    resize!(V, lengthI)
    resize!(csrcolval, lengthI)
    resize!(csrnzval, lengthI)
    params = context.work.params
    @debug "any(isnan.(params))=$(any(isnan.(params)))"
    @debug "any(isnan.(endogenous))=$(any(isnan.(endogenous)))"
    @debug "any(isnan.(initialvalues))=$(any(isnan.(initialvalues)))"
    @debug "any(isnan.(exogenous))=$(any(isnan.(exogenous)))"
    @debug "any(isnan.(steadystate))=$(any(isnan.(steadystate)))"
    function make_one_period!(r::Int64, t::Int64, tid::Int64)
        jacobian = get_dynamic_jacobian!(
            ws[tid],
            params,
            endogenous,
            exogenous,
            JA.steadystate,
            md,
            df,
            t,
        )
        i, j, v = findnz(sparse(jacobian))
        oc = (t - 1) * nvar
        @inbounds for el = 1:length(i)
            if j[el] <= maxcol
                kjel = jacobian_columns[j[el]]
                I[r] = permute_row(i[el], permutations) + oc
                J[r] = kjel + oc - nvar
                V[r] = v[el]
                r += 1
            end
        end
        return r
    end
    jacobian = get_initial_jacobian!(
        ws[1],
        params,
        endogenous,
        initialvalues,
        exogenous,
        steadystate,
        md,
        df,
        2,
    )
    r = 1
    @debug "any(isnan.(ws[1].jacobian))=$(any(isnan.(jacobian)))"
    jacobian_columns =
        [i for (i, x) in enumerate(transpose(md.lead_lag_incidence)) if x > 0]
    i, j, v = findnz(jacobian)
    @debug "period 2: isnan.(v) = $(findall(isnan.(v)))"
    fill!(V, 0.0)
    @inbounds for el = 1:length(i)
        if j[el] <= maxcol
            kjel = jacobian_columns[j[el]]
            if kjel > nvar
                I[r] = permute_row(i[el], permutations)
                J[r] = kjel - nvar
                V[r] = v[el]
                r += 1
            end
        end
    end
    #    @Threads.threads
    for t = 2:(periods-1)
        #        tid = Threads.threadid()
        tid = 1
        r = make_one_period!(r, t, tid)
    end
    jacobian = get_terminal_jacobian!(
        ws[1],
        params,
        endogenous,
        terminalvalues,
        exogenous,
        steadystate,
        md,
        df,
        periods,
    )
    i, j, v = findnz(sparse(jacobian))
    @debug "period $periods: isnan.(v) = $(findall(isnan.(v)))"
    oc = (periods - 1) * nvar
    @inbounds for el = 1:length(i)
        if j[el] <= maxcol
            kjel = jacobian_columns[j[el]]
            if kjel <= 2 * nvar
                I[r] = permute_row(i[el], permutations) + oc
                J[r] = kjel + oc - nvar
                V[r] = v[el]
                r += 1
            end
        end
    end
    #=
    # terminal condition
    ic = md.n_bkwrd + md.n_both + md.n_current .+ (1:md.n_fwrd+md.n_both)
    g = context.results.model_results[1].linearrationalexpectations.g1_1
    @inbounds CmultG!(
        JA.tmp_nvar_npred,
        JA.tmp_nvar_nfwrd,
        JA.tmp_nfwrd_npred,
        ws[1].jacobian,
        g,
        ic,
        md.i_fwrd_b,
    )
    (i, j, v) = findnz(sparse(JA.tmp_nvar_npred))
    needed_space = r + length(i) - 1 
    extra_space =  needed_space - length(I) 
    @show extra_space
    if extra_space > 0
        resize!(I, needed_space)
        resize!(J, needed_space)
        resize!(V, needed_space)
        resize!(csrcolval, needed_space)
        resize!(csrnzval, needed_space)
    end
    @inbounds for el = 1:length(i)
        I[r] = i[el] + oc
        J[r] = jacobian_columns[j[el]] + oc
        V[r] = v[el]
        r += 1
    end
    =#
    resize!(I, r - 1)
    resize!(J, r - 1)
    resize!(V, r - 1)
    resize!(csrcolval, r - 1)
    resize!(csrnzval, r - 1)
    n = periods * nvar
    #length(colptr) == n + 1 && colptr[end] - 1 == length(rowval) == length(nzval)

    A = SparseArrays.sparse!(
        I,
        J,
        V,
        n,
        n,
        (x, y) -> x + y,
        klasttouch,
        csrrowptr,
        csrcolval,
        csrnzval,
        csccolptr,
        J,
        V,
    )
    return A
end

function CmultG!(
    y::AbstractMatrix{Float64},
    tmp_nvar_nfwrd::AbstractMatrix{Float64},
    tmp_nfwrd_npred::AbstractMatrix{Float64},
    jacobian::AbstractMatrix{Float64},
    g::AbstractMatrix{Float64},
    ic::AbstractVector{Int64},
    ir::AbstractVector{Int64},
)
    fill!(tmp_nvar_nfwrd, 0.0)
    nvar, npred = size(y)
    nfwrd = length(ir)
    for i = 1:nfwrd
        for j = 1:nvar
            tmp_nvar_nfwrd[j, i] = jacobian[j, ic[i]]
        end
    end
    fill!(tmp_nfwrd_npred, 0.0)
    for i = 1:npred
        for j = 1:nfwrd
            tmp_nfwrd_npred[j, i] = g[ir[j], i]
        end
    end
    mul!(y, tmp_nvar_nfwrd, tmp_nfwrd_npred)
end
