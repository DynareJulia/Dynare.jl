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
        dynamic_variables = zeros(md.n_bkwrd + md.n_current + md.n_fwrd + 2 * md.n_both)
        temp_vec = Vector{Float64}(undef, sum(md.dynamic!.tmp_nbr[1:2]))
        residuals = Vector{Float64}(undef, md.endogenous_nbr)
        params = work.params
        new(jacobian, dynamic_variables, temp_vec, residuals, params)
    end
end

struct Jacobian
    nrow::Int64
    I::Vector{Int64}
    J::Vector{Int64}
    V::Vector{Float64}
    klasstouch::Vector{Int64}
    colptr::Vector{Int64}
    rowptr::Vector{Int64}
    colval::Vector{Int64}
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
        vj = view(work.jacobian, :, 1:maxcol)
        nz = periods * md.NNZDerivatives[1]
        nrow = periods * md.endogenous_nbr
        I = Vector{Int64}(undef, nz)
        @inbounds fill!(I, 1)
        J = Vector{Int64}(undef, nz)
        @inbounds fill!(J, 1)
        V = Vector{Float64}(undef, nz)
        @inbounds fill!(V, 0.0)
        klasstouch = Vector{Int64}(undef, nz)
        colptr = Vector{Int64}(undef, nrow + 1)
        rowptr = Vector{Int64}(undef, nrow + 1)
        colval = Vector{Int64}(undef, nz)
        nzval = Vector{Float64}(undef, nz)
        npred = md.n_bkwrd + md.n_both
        nfwrd = md.n_fwrd + md.n_both
        tmp_nvar_npred = zeros(nvar, npred)
        tmp_nvar_nfwrd = zeros(nvar, nfwrd)
        tmp_nfwrd_npred = zeros(nfwrd, npred)
        new(
            nrow,
            I,
            J,
            V,
            klasstouch,
            colptr,
            rowptr,
            colval,
            nzval,
            maxcol,
            steadystate,
            tmp_nvar_npred,
            tmp_nvar_nfwrd,
            tmp_nfwrd_npred,
        )
    end
end

function get_dynamic_endogenous_variables!(
    y::AbstractVector{Float64},
    data::AbstractVector{Float64},
    lli::Matrix{Int64},
    m::Model,
    period::Int64,
)
    m, n = size(lli)
    p = (period - 2) * n
    @inbounds for i = 1:m
        for j = 1:n
            k = lli[i, j]
            if k > 0
                y[k] = data[p+j]
            end
        end
        p += n
    end
end


function get_jacobian!(
    ws::JacTimesVec,
    endogenous::AbstractVector{Float64},
    exogenous::Matrix{Float64},
    steadystate::Vector{Float64},
    m::Model,
    period::Int64,
)
    lli = m.lead_lag_incidence
    get_dynamic_endogenous_variables!(ws.dynamic_variables, endogenous, lli, m, period)
    compute_jacobian(ws, exogenous, steadystate, m, period)
end


function compute_jacobian(
    ws::JacTimesVec,
    exogenous::AbstractMatrix{Float64},
    steadystate::AbstractVector{Float64},
    m::Model,
    period::Int64,
)
    dynamic! = m.dynamic!.dynamic!
    fill!(ws.jacobian, 0.0)
    @inbounds Base.invokelatest(
        dynamic!,
        ws.temp_vec,
        ws.residuals,
        ws.jacobian,
        ws.dynamic_variables,
        exogenous,
        ws.params,
        steadystate,
        period,
    )
end

include("makeA.jl")
include("gmres_solver.jl")
