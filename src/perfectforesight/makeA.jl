using SparseArrays

function bigindex!(bigrowval, r, rowval, c1, c2, offset)
    c1 > length(rowval) && return r
    for c in c1:c2
        bigrowval[r] = rowval[c] + offset
        r += 1
    end
    return r
end

function permutation(permutations, colptr, rowval)
    p0 = Tuple{Int64, Int64}[]
    rowval1 = copy(rowval)
    rowval2 = copy(rowval)
    k = 1
    for i in 1:length(colptr) - 1
        if colptr[i + 1] > colptr[i]
            vr1 = view(rowval1, colptr[i]:colptr[i+1]-1)
            vr2 = view(rowval2, colptr[i]:colptr[i+1]-1)
            for p in permutations
                p1, p2 = p
                p1 > p2 && ((p2, p1) = (p1, p2))
                for (j, r) in enumerate(vr1)
                    if r == p1
                        vr1[j] = p2
                    elseif r == p2
                        vr1[j] = p1
                        break
                    end
                end
            end
            copy!(vr2, vr1)
            if !issorted(vr1)
                isort = sortperm(vr1)
                for (j, k) in enumerate(isort)
                    if j != k
                        push!(p0, (j, k) .+ colptr[i] .- 1)
                        vr2[j] = vr1[k]
                    end
                end
            end
        end
    end
    return(p0, rowval2)
end

function makeJacobian(colptr, rowval, endogenous_nbr, periods, permutations)
    startv = colptr[endogenous_nbr + 1]
    endv = colptr[2*endogenous_nbr + 1]

    rowval_length = colptr[3*endogenous_nbr + 1] - 1
    nnz = (periods - 1)*rowval_length - startv + endv
    bigcolptr = Vector{Int64}(undef, periods*endogenous_nbr + 1)
    bigrowval = Vector{Int64}(undef, nnz)
    bignzval = similar(bigrowval, Float64)
    permutations1 = []

    # compute permutations for one period Jacobian
    if !isempty(permutations)
        (permutations1, rowval) = permutation(permutations, colptr, rowval)
    end
    
    r = 1
    c = 1
    # first periods
    bigcolptr[1] = 1
    for i in 1:endogenous_nbr
        c += 1
        offset = 0
        bigcolptr[c] = (bigcolptr[c - 1]
                        + colptr[endogenous_nbr + i + 1]
                        - colptr[endogenous_nbr + i]
                        + colptr[i + 1]
                        - colptr[i])
        r1 = colptr[endogenous_nbr + i]
        r= bigindex!(bigrowval,
                     r,
                     rowval,
                     colptr[endogenous_nbr + i],
                     colptr[endogenous_nbr + i + 1] - 1,
                     offset)
        r= bigindex!(bigrowval,
                     r,
                     rowval,
                     colptr[i],
                     colptr[i + 1] - 1,
                     offset + endogenous_nbr)
    end

    # intermediary periods
    for p in 2:periods - 1
        for i in 1:endogenous_nbr
            c += 1
            offset = (p - 2)*endogenous_nbr
            bigcolptr[c] = (bigcolptr[c - 1]
                            + colptr[2*endogenous_nbr + i + 1]
                            - colptr[2*endogenous_nbr + i]
                            + colptr[endogenous_nbr + i + 1]
                            - colptr[endogenous_nbr + i]
                            + colptr[i + 1]
                            - colptr[i])
            r1 = colptr[endogenous_nbr + i]
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[2*endogenous_nbr + i],
                         colptr[2*endogenous_nbr + i + 1] - 1,
                         offset)
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[endogenous_nbr + i],
                         colptr[endogenous_nbr + i + 1] - 1,
                         offset + endogenous_nbr)
            r= bigindex!(bigrowval,
                         r,
                         rowval,
                         colptr[i],
                         colptr[i + 1] - 1,
                         offset + 2*endogenous_nbr)
        end
        offset += endogenous_nbr
    end
    #terminal period
    for i in 1:endogenous_nbr
        c += 1
        offset = (periods - 2)*endogenous_nbr
        bigcolptr[c] = (bigcolptr[c - 1]
                        + colptr[2*endogenous_nbr + i + 1]
                        - colptr[2*endogenous_nbr + i]
                        + colptr[endogenous_nbr + i + 1]
                        - colptr[endogenous_nbr + i])
        r1 = colptr[endogenous_nbr + i]
        r= bigindex!(bigrowval,
                     r,
                     rowval,
                     colptr[2*endogenous_nbr + i],
                     colptr[2*endogenous_nbr + i + 1] - 1,
                     offset)
        r= bigindex!(bigrowval,
                     r,
                     rowval,
                     colptr[endogenous_nbr + i],
                     colptr[endogenous_nbr + i + 1] - 1,
                     offset + endogenous_nbr)

    end

    b1 = colptr[endogenous_nbr + 1]
    b2 = (periods - 1)*rowval_length + colptr[2*endogenous_nbr + 1] - 1
    nperiods = periods*endogenous_nbr
    return (SparseMatrixCSC(nperiods, nperiods, bigcolptr, bigrowval, bignzval), permutations1)
end

function reorder_derivatives!(nzval, permutations, ws)
    copy!(ws, nzval)
    for p in permutations
        p1, p2 = p
        nzval[p1] = ws[p2]
    end
end

function updateJacobian!(J::SparseMatrixCSC,
                         G1!,
                         endogenous::AbstractVector{<: Real},
                         initialvalues::AbstractVector{<: Real},
                         terminalvalues::AbstractVector{<: Real},
                         dynamic_variables::AbstractVector{<: Real},
                         exogenous::AbstractVector{<: Real},
                         periods,
                         temporary_var::AbstractVector{<: Real},
                         params::AbstractVector{<: Real},
                         steady_state::AbstractVector{<: Real},
                         colptr::AbstractVector{Int},
                         nzval::AbstractVector{<: Real},
                         endogenous_nbr,
                         exogenous_nbr,
                         permutations,
                         ws::AbstractVector{<: Real})
    bigcolptr = J.colptr
    offset = 1
    rx = 1:exogenous_nbr
    @views begin
        copyto!(dynamic_variables, initialvalues)
        copyto!(dynamic_variables, endogenous_nbr + 1, endogenous, 1, 2*endogenous_nbr)
        G1!(temporary_var, nzval, dynamic_variables, exogenous[rx], params, steady_state)
        !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
        oy = endogenous_nbr
        ox = exogenous_nbr
        for c in 1:2*endogenous_nbr
            k = bigcolptr[c]
            n = colptr[endogenous_nbr + c+1] - colptr[endogenous_nbr + c]
            copyto!(J.nzval, k, nzval, colptr[endogenous_nbr + c], n)
        end
        ry = 1:3*endogenous_nbr
        rx = rx .+ ox
        G1!(temporary_var, nzval, endogenous[ry], exogenous[rx], params, steady_state)
        !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
        for c in 1:2*endogenous_nbr
            k = bigcolptr[c] + colptr[endogenous_nbr + c + 1] - colptr[endogenous_nbr + c]
            n = colptr[c+1] - colptr[c]
            copyto!(J.nzval, k, nzval, colptr[c], n)
        end
        for c in 2*endogenous_nbr + 1:3*endogenous_nbr
            k = bigcolptr[c]
            n = colptr[c+1] - colptr[c]
            copyto!(J.nzval, k, nzval, colptr[c], n)
        end

        for p in 1:periods - 3
            ry = ry .+ oy
            rx = rx .+ ox
            G1!(temporary_var, nzval, endogenous[ry], exogenous[rx], params, steady_state)
            !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
            for c in 1:endogenous_nbr
                k = (bigcolptr[c + p*endogenous_nbr]
                     + colptr[c + 2*endogenous_nbr + 1] - colptr[c + 2*endogenous_nbr]
                     + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
            for c in endogenous_nbr + 1:2*endogenous_nbr
                k = (bigcolptr[c + p*endogenous_nbr]
                     + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
            for c in 2*endogenous_nbr + 1: 3*endogenous_nbr
                k = bigcolptr[c + p*endogenous_nbr]
                n = colptr[c+1] - colptr[c]
                copyto!(J.nzval, k, nzval, colptr[c], n)
            end
        end
        
        rx = rx .+ ox
        copyto!(dynamic_variables, 1, endogenous, (periods - 2)*endogenous_nbr + 1, 2*endogenous_nbr)
        copyto!(dynamic_variables, 2*endogenous_nbr + 1, terminalvalues)
        G1!(temporary_var, nzval, dynamic_variables, exogenous[rx], params, steady_state)
        !isempty(permutations) && reorder_derivatives!(nzval, permutations, ws)
        for c in 1:endogenous_nbr
            k = (bigcolptr[c + (periods - 2)*endogenous_nbr]
                 + colptr[c + 2*endogenous_nbr + 1] - colptr[c + 2*endogenous_nbr]
                 + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
            n = colptr[c+1] - colptr[c]
            copyto!(J.nzval, k, nzval, colptr[c], n)
        end
        for c in endogenous_nbr + 1:2*endogenous_nbr
            k = (bigcolptr[c + (periods - 2)*endogenous_nbr]
                 + colptr[c + endogenous_nbr + 1] - colptr[c + endogenous_nbr])
            n = colptr[c+1] - colptr[c]
            copyto!(J.nzval, k, nzval, colptr[c], n)
        end
    end
    return J
end

function CmultG!(
    y::AbstractMatrix{Float64},
    tmp_nvar_nfwrd::AbstractMatrix{Float64},
    tmp_nfwrd_npred::AbstractMatrix{Float64},
    jacobian::AbstractMatrix{Float64},
    g::AbstractMatrix{Float64},
    ic::AbstractVector{Int},
    ir::AbstractVector{Int},
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
