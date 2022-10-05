using Dynare
using LinearAlgebra
using SparseArrays
using Test

function bigindex!(bigindex, r, c1, c2)
    for c in c1:c2
        bigrowval[r] = j
        r += 1
    end
    return r
end

function reorder_blocks!(y, x, k, periods, offset)
    y[k1] .= 
    for p in 1:periods
        
function makeJacobian(colptr, rowval, endogenous_nbr, periods)
    startv = colptr[endogenous_nbr + 1]
    endv = colptr[2*endogenous_nbr + 1]

    rowval_length = colptr[3*endogenous_nbr + 1] - 1
    nnz = (periods - 1)*rowval_length - startv + endv
    bigcolptr = Vector{Int64}(undef, periods*endogenous_nbr + 1)
    bigrowval = Vector{Int64}(undef, nnz)
    bignzval = similar(bigrowval, Float64)
    index1 = Vector{Int64}(undef, rowval_length - colptr[endogenous_nbr + 1] + 1)
    index2 = Vector{Int64}(undef, rowval_length)
    index3 = Vector{Int64}(undef, colptr[2*endogenous_nbr + 1] - 1)
    
                           
    r = 1
    c = 1
    # first periods
    bigcolptr[1] = 1
    for i in 1:endogenous_nbr
        c += 1
        bigcolptr[c] = (bigcolptr[c - 1]
                        + colptr[endogenous_nbr + i + 1]
                        - colptr[endogenous_nbr + i]
                        + colptr[i + 1]
                        - colptr[i])
        r1 = colptr[endogenous_nbr + i]
        r= bigindex!(index1,
                     r,
                     colptr[endogenous_nbr + i],
                     colptr[endogenous_nbr + i + 1] - 1)
        r= bigindex!(index1,
                     r,
                     colptr[2*endogenous_nbr + i],
                     colptr[2*endogenous_nbr + i + 1] - 1)
        r= bigindex!(index2,
                     r,
                     colptr[i],
                     colptr[i + 1] - 1)
        r= bigindex!(index2,
                     r,
                     colptr[endogenous_nbr + i],
                     colptr[endogenous_nbr + i + 1] - 1)
        r= bigindex!(index2,
                     r,
                     colptr[2*endogenous_nbr + i],
                     colptr[2*endogenous_nbr + i + 1] - 1)
        r= bigindex!(index3,
                     r,
                     colptr[i],
                     colptr[i + 1] - 1)
        r= bigindex!(index3,
                     r,
                     colptr[endogenous_nbr + i],
                     colptr[endogenous_nbr + i + 1] - 1)
    end

    # intermediary periods
    offset = 0
    for p in 2:periods - 1
        for i in 1:endogenous_nbr
            c += 1
            bigcolptr[c] = (bigcolptr[c - 1]
                            + colptr[2*endogenous_nbr + i + 1]
                            - colptr[2*endogenous_nbr + i]
                            + colptr[endogenous_nbr + i + 1]
                            - colptr[endogenous_nbr + i]
                            + colptr[i + 1]
                            - colptr[i])
            r1 = colptr[endogenous_nbr + i]
            r= rowval_column!(bigrowval,
                              r,
                              rowval,
                              colptr[2*endogenous_nbr + i],
                              colptr[2*endogenous_nbr + i + 1] - 1,
                              offset)
            r= rowval_column!(bigrowval,
                              r,
                              rowval,
                              colptr[endogenous_nbr + i],
                              colptr[endogenous_nbr + i + 1] - 1,
                              offset + endogenous_nbr)
            r= rowval_column!(bigrowval,
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
        bigcolptr[c] = (bigcolptr[c - 1]
                        + colptr[2*endogenous_nbr + i + 1]
                        - colptr[2*endogenous_nbr + i]
                        + colptr[endogenous_nbr + i + 1]
                        - colptr[endogenous_nbr + i])
        r1 = colptr[endogenous_nbr + i]
        r= rowval_column!(bigrowval,
                          r,
                          rowval,
                          colptr[2*endogenous_nbr + i],
                          colptr[2*endogenous_nbr + i + 1] - 1,
                          offset)
        r= rowval_column!(bigrowval,
                          r,
                          rowval,
                          colptr[endogenous_nbr + i],
                          colptr[endogenous_nbr + i + 1] - 1,
                          offset + endogenous_nbr)

    end

    b1 = colptr[endogenous_nbr + 1]
    b2 = (periods - 1)*rowval_length + colptr[2*endogenous_nbr + 1] - 1
    nperiods = periods*endogenous_nbr
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


A  = sprand(5, 15, 0.1)
AA = vcat(hcat(A[:, 6:15], zeros(5, 10)),
          hcat(A, zeros(5, 5)),
          hcat(zeros(5, 5), A),
          hcat(zeros(5, 10), A[:, 1:10]))


J = makeJacobian(A.colptr, A.rowval, 5, 4)

@test J.colptr == AA.colptr
@test J.rowval == AA.rowval


context = @dynare "test/models/example1pf/example1pf_sparse.mod"

m = context.models[1]

periods = 4

results = context.results.model_results[1]
steady_state = results.trends.endogenous_steady_state
endogenous = repeat(steady_state, periods + 2)
exogenous = repeat(results.exogenous_steady_state, periods + 2)
temporary_var = Vector{Float64}(undef, sum(m.dynamic_tmp_nbr[1:2]))
params = context.work.params
rowval = m.dynamic_g1_sparse_rowval
colptr = m.dynamic_g1_sparse_colptr

J = makeJacobian(colptr,
                 rowval,
                 m.endogenous_nbr,
                 periods)

r1 = colptr[m.endogenous_nbr + 1]
n = colptr[3*m.endogenous_nbr + 1] - 1
n1 = colptr[3*m.endogenous_nbr  + 1] - colptr[m.endogenous_nbr + 1] + 1
n2 = colptr[2*m.endogenous_nbr + 1]
df = Dynare.DFunctions
nzval = Vector{Float64}(undef, colptr[end] - 1)
df.SparseDynamicG1!(temporary_var, nzval, endogenous[1:18], exogenous[1:2], params, steady_state, true)
A = SparseMatrixCSC(m.endogenous_nbr, 3*m.endogenous_nbr + m.exogenous_nbr, colptr, rowval, nzval)
AA = vcat(hcat(A[:, 7:18], zeros(6, 12)),
          hcat(A[:,1:18], zeros(6, 6)),
          hcat(zeros(6, 6), A[:,1:18]),
          hcat(zeros(6, 12), A[:, 1:12]))

@test J.colptr == AA.colptr
@test J.rowval == AA.rowval
updateJacobian!(J, df.SparseDynamicG1!, endogenous, exogenous, periods, temporary_var, params, steady_state, colptr, nzval, m.endogenous_nbr, m.exogenous_nbr)

@test J == AA
