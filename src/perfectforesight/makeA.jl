function make_one_period!(JA::Jacobian,
                          endogenous::AbstractVector{Float64},
                          exogenous::AbstractMatrix{Float64},
                          periods::Int64,
                          md::Model,
                          jacobian_columns::Vector{Int64},
                          nnz_period::Int64,
                          maxcol::Int64,
                          ws::JacTimesVec,
                          t::Int64
                          )
    nvar = md.endogenous_nbr
    oc = (t - 1)*md.endogenous_nbr
    get_jacobian!(ws, endogenous, exogenous, JA.steadystate, md, t)
    i, j, v = findnz(sparse(ws.jacobian))
    r1 = (t - 1)*nnz_period + 1
    for el = 1:length(i)
        if j[el] <= maxcol
            kjel = jacobian_columns[j[el]]
            JA.I[r1] = i[el] + oc
            JA.J[r1] = kjel + oc - nvar
            JA.V[r1] = v[el]
            r1 += 1
        end
    end
end

function makeJacobian!(
    JA::Jacobian,
    endogenous::AbstractVector{Float64},
    exogenous::AbstractMatrix{Float64},
    context::Context,
    periods::Int64,
    ws::Vector{JacTimesVec}
)
    md = context.models[1]
    steadystate = JA.steadystate
    maxcol = JA.maxcol
    nvar = md.endogenous_nbr
    npred = md.n_bkwrd + md.n_both
    nnz_period = md.NNZDerivatives[1]
    get_jacobian!(ws[1], endogenous, exogenous, steadystate, md, 2)
    r = 1
    jacobian_columns = [i for (i, x) in enumerate(transpose(md.lead_lag_incidence)) if x > 0]
    i, j, v = findnz(sparse(ws[1].jacobian))
    for el = 1:length(i)
        if j[el] <= maxcol
            kjel = jacobian_columns[j[el]]
            if kjel > nvar
                JA.I[r] = i[el]
                JA.J[r] = kjel - nvar
                JA.V[r] = v[el]
                r += 1
            end
        end
    end
    @Threads.threads for t = 2:(periods-1)
        tid = Threads.threadid()
        make_one_period!(JA, endogenous, exogenous, periods, md,
                         jacobian_columns, nnz_period, maxcol,
                         ws[tid], t)
    end
    get_jacobian!(ws[1], endogenous, exogenous, steadystate, md, periods)
    i, j, v = findnz(sparse(ws[1].jacobian))
    oc = (periods - 1)*nvar
    r = (periods - 1)*nnz_period + 1
    for el = 1:length(i)
        if j[el] <= maxcol
            kjel = jacobian_columns[j[el]]
            if kjel <= 2*nvar
                JA.I[r] = i[el] + oc
                JA.J[r] = kjel + oc - nvar
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
        ws[1].jacobian,
        g,
        ic,
        md.i_fwrd_b,
    )
    (i, j, v) = findnz(sparse(JA.tmp_nvar_npred))
    needed_space = periods*nnz_period + length(i) - 1 
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
        JA.I[r] = i[el] + oc
        JA.J[r] = jacobian_columns[j[el]] + oc
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

