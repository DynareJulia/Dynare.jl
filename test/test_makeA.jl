using Dynare
using LinearAlgebra
using SparseArrays
using Test

function rowval_column!(bigrowval, r, rowval, c1, c2, offset)
    for j in c1:c2
        bigrowval[r] = rowval[j] + offset 
        r += 1
    end
    return r
end

function makeJacobian(colptr, rowval, endogenous_nbr, periods)
    startv = colptr[endogenous_nbr + 1]
    endv = colptr[2*endogenous_nbr + 1]

    rowval_length = colptr[3*endogenous_nbr + 1] - 1
    nnz = (periods - 1)*rowval_length - startv + endv
    bigcolptr = Vector{Int64}(undef, periods*endogenous_nbr + 1)
    bigrowval = Vector{Int64}(undef, nnz)
    bignzval = Vector{Float64}(undef, nnz)
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
        r= rowval_column!(bigrowval,
                          r,
                          rowval,
                          colptr[endogenous_nbr + i],
                          colptr[endogenous_nbr + i + 1] - 1,
                          0)
        r= rowval_column!(bigrowval,
                          r,
                          rowval,
                          colptr[i],
                          colptr[i + 1] - 1,
                          endogenous_nbr)

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

function updateJacobian!(J::SparseMatrixCSC, G1!, endogenous, exogenous, periods, temporary_var, params, steady_state, nzval, r1, n, n1, n2, endogenous_nbr, exogenous_nbr)
    offset = 1
    ry = 1:3*endogenous_nbr
    rx = 1:exogenous_nbr
    oy = 0
    ox = 0
    @views begin
        @show endogenous[ry]
        @show exogenous[rx]
        G1!(temporary_var, nzval, endogenous[ry], exogenous[rx], params, steady_state, true)
        copyto!(J.nzval, offset, nzval, r1, n1)
        display(Matrix(SparseMatrixCSC(endogenous_nbr, 3*endogenous_nbr + exogenous_nbr, m.dynamic_g1_sparse_colptr, m.dynamic_g1_sparse_rowval, nzval)))
        offset += n1
        @show offset
        oy += endogenous_nbr
        ox += exogenous_nbr
        
        for p in 2:periods - 1
            ry1 = ry .+ oy
            rx1 = rx .+ ox
            G1!(temporary_var, nzval, endogenous[ry1], exogenous[rx1], params, steady_state, true)
            copyto!(J.nzval, offset, nzval, 1, n)
            offset += n
            @show offset
            oy += endogenous_nbr
            ox += exogenous_nbr
        end
        
        ry1 = ry .+ oy
        rx1 = rx .+ ox
        G1!(temporary_var, nzval, endogenous[ry1], exogenous[rx1], params, steady_state, true)
        @show offset, n2
        copyto!(J.nzval, offset, nzval, 1, n2)
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
@show params
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
updateJacobian!(J, df.SparseDynamicG1!, endogenous, exogenous, periods, temporary_var, params, steady_state, nzval, r1, n, n1, n2, m.endogenous_nbr, m.exogenous_nbr)                  
