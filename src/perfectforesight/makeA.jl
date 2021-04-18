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
        get_jacobian!(work, endogenous, exogenous, steadystate, md, 2)
        maxcol = size(work.jacobian, 2) - md.exogenous_nbr
        vj = view(work.jacobian, :, 1:maxcol)
        nz = periods * nnz(sparse(vj))
        nrow = periods * md.endogenous_nbr
        I = Vector{Int64}(undef, nz)
        J = Vector{Int64}(undef, nz)
        V = Vector{Float64}(undef, nz)
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

function makeJacobian!(
    JA::Jacobian,
    endogenous::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    context::Context,
    periods::Int64,
)
    md = context.models[1]
    work = context.work
    steadystate = JA.steadystate
    maxcol = JA.maxcol
    nvar = md.endogenous_nbr
    npred = md.n_bkwrd + md.n_both
    get_jacobian!(work, endogenous, exogenous, steadystate, md, 2)
    r = 1
    k = [i for (i, x) in enumerate(transpose(md.lead_lag_incidence)) if x > 0]
    i, j, v = findnz(sparse(work.jacobian))
    for el = 1:length(i)
        if j[el] <= maxcol
            kjel = k[j[el]]
            if kjel > nvar
                JA.I[r] = i[el]
                JA.J[r] = k[j[el]] - nvar
                JA.V[r] = v[el]
                r += 1
            end
        end
    end
    offset = nvar
    for t = 2:(periods-1)
        get_jacobian!(work, endogenous, exogenous, steadystate, md, 2)
        i, j, v = findnz(sparse(work.jacobian))
        for el = 1:length(i)
            if j[el] <= maxcol
                kjel = k[j[el]]
                JA.I[r] = i[el] + offset
                JA.J[r] = kjel + offset - nvar
                JA.V[r] = v[el]
                r += 1
            end
        end
        offset += nvar
    end
    get_jacobian!(work, endogenous, exogenous, steadystate, md, 2)
    i, j, v = findnz(sparse(work.jacobian))
    for el = 1:length(i)
        if j[el] <= maxcol
            kjel = k[j[el]]
            if kjel <= 2*nvar
                JA.I[r] = i[el] + offset
                JA.J[r] = kjel + offset - nvar
                JA.V[r] = v[el]
                r += 1
            end
        end
    end
    # terminal condition
     ic = md.n_bkwrd + md.n_both + md.n_current .+ (1:md.n_fwrd+md.n_both)
    g = context.results.model_results[1].linearrationalexpectations.g1_1
    CmultG!(
        JA.tmp_nvar_npred,
        JA.tmp_nvar_nfwrd,
        JA.tmp_nfwrd_npred,
        work.jacobian,
        g,
        ic,
        md.i_fwrd_b,
    )
    (i, j, v) = findnz(sparse(JA.tmp_nvar_npred))
    needed_space = r + length(i) - 1 
    extra_space =  needed_space - length(JA.I) 
    if extra_space > 0
        resize!(JA.I, needed_space)
        resize!(JA.J, needed_space)
        resize!(JA.V, needed_space)
        resize!(JA.klasstouch, needed_space)
        resize!(JA.colval, needed_space)
        resize!(JA.nzval, needed_space)
    end
    for el = 1:length(i)
        JA.I[r] = i[el] + offset
        JA.J[r] = k[j[el]] + offset
        JA.V[r] = v[el]
        r += 1
    end
    resize!(JA.I, r - 1)
    resize!(JA.J, r - 1)
    resize!(JA.V, r - 1)
    n = periods * nvar
    return SparseArrays.sparse!(
        JA.I,
        JA.J,
        JA.V,
        n,
        n,
        (x, y) -> x + y,
        JA.klasstouch,
        JA.rowptr,
        JA.colval,
        JA.nzval,
        JA.colptr,
        JA.J,
        JA.V,
    )
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

