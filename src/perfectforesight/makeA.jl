using SparseArrays

function permute_row(ir, permutations)
    for p in permutations
        ir == p[1] && return p[2]
        ir == p[2] && return p[1]
    end
    return ir
end

function bigindex!(bigrowval, r, rowval, c1, c2, offset)
    c1 > length(rowval) && return r
    @show c1, c2
    offset1 = offset - rowval[c1] + 1
    for c in c1:c2
        @show r, rowval[c] + offset
        bigrowval[r] = rowval[c] + offset
        r += 1
    end
    return r
end

        
function makeJacobian(colptr, rowval, endogenous_nbr, periods, permutations)
    startv = colptr[endogenous_nbr + 1]
    endv = colptr[2*endogenous_nbr + 1]

    rowval_length = colptr[3*endogenous_nbr + 1] - 1
    nnz = (periods - 1)*rowval_length - startv + endv
    @show nnz
    bigcolptr = Vector{Int64}(undef, periods*endogenous_nbr + 1)
    bigrowval = Vector{Int64}(undef, nnz)
    bignzval = similar(bigrowval, Float64)
    
    for (i, r) in enumerate(rowval)
        rowval[i] = permute_row(r, permutations)
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
        @show bigcolptr
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
        @show p
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
            @show bigcolptr
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
        @show bigcolptr
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
    @show size(bigcolptr)
    @show size(bigrowval)
    @show size(bignzval)
    return SparseMatrixCSC(nperiods, nperiods, bigcolptr, bigrowval, bignzval)
end

function updateJacobian!(J::SparseMatrixCSC, G1!, endogenous, exogenous, periods, temporary_var, params, steady_state, colptr, nzval, endogenous_nbr, exogenous_nbr)
    bigcolptr = J.colptr
    offset = 1
    ry = 1:3*endogenous_nbr
    rx = 1:exogenous_nbr
    @views begin
        G1!(temporary_var, nzval, endogenous[ry], exogenous[rx], params, steady_state, true)
        oy = endogenous_nbr
        ox = exogenous_nbr
        for c in 1:2*endogenous_nbr
            k = bigcolptr[c]
            n = colptr[endogenous_nbr + c+1] - colptr[endogenous_nbr + c]
            copyto!(J.nzval, k, nzval, colptr[endogenous_nbr + c], n)
        end
        ry1 = ry .+ oy
        rx1 = rx .+ ox

        G1!(temporary_var, nzval, endogenous[ry1], exogenous[rx1], params, steady_state, true)
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
            ry1 = ry .+ oy
            rx1 = rx .+ ox
            G1!(temporary_var, nzval, endogenous[ry1], exogenous[rx1], params, steady_state, true)
            oy += endogenous_nbr
            ox += exogenous_nbr
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
        
        ry1 = ry .+ oy
        rx1 = rx .+ ox
        G1!(temporary_var, nzval, endogenous[ry1], exogenous[rx1], params, steady_state, true)
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
