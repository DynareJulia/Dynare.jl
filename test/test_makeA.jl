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

    rowval_length = length(rowval)
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
            @show colptr[endogenous_nbr + i]:colptr[endogenous_nbr + i + 1] - 1
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
        @show colptr[endogenous_nbr + i]:colptr[endogenous_nbr + i + 1] - 1
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
    @show length(bignzval)
    @show length(bigrowval)
    @show b2
    nperiods = periods*endogenous_nbr
    return SparseMatrixCSC(nperiods, nperiods, bigcolptr, bigrowval, bignzval)

end

function updateJacobian!(J::SparseMatrixCSC, G1!, endogenous, exogenous, periods, temporary_var, params, steady_state, rowval, r1, r2, n, n1, n2)
    offset = 1
    ry = 1:3*endogenous_nbr
    rx = 1:exogenous_nbr
    @views begin
        G1!(temporary_var, rowval, endogenous[ry], exogenous[rx], params, steady_state, true)
        copyto!(J.nzval, offset, rowval, c1, n1)
        offset += n1
        ry .+= 1
        rx .+= 1
        
        for p in 2:periods - 1
            G1!(temporary_var, rowval, endogenous[ry], exogenous[rx], params, steady_state, true)
            copyto!(J.nzval, offset, rowval, 1, n)
            offset += n
            ry .+= 1
            rx .+= 1
        end
        
        G1!(temporary_var, rowval, endogenous[ry], exogenous[rx], params, steady_state, true)
        copyto!(J.nzval, 1, rowval, c2, n2)
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


context = @dynare "test/models/example1/example1"

m = context.models[1]

periods = 4
J = makeJacobian(m.dynamic_g1_sparse_colptr,
                 m.dynamic_g1_sparse_rowval,
                 m.endoegnous_nbr,
                 periods)

results = context.results.model_results[1]
steady_state = results.endogenous_steady_state
endogenous = repeat(steady_state, periods)
exogenous = repeat(results.exogenous_steady_state, periods)
temporary_var = Vector{Float64}(undef, sum(m.dynamic_tmp_nbr[1:2]))
params = context.work.params
rowval = m.dynare_g1_sparse_rowval
colptr = m.dynare_g1_sparse_colptr
r1 = colptr[m.endogenous_nbr + 1]
r2 = colptr[2*m.endogenous_nbr + 1] - 1
n = colptr[3*m.endogenous_nbr + 1] - 1
n1 = colptr[3*m.endogenous_nbr  + 1] - colptr[m.endogenous_nbr + 1]
n2 = r2
updateJacobian!(J, df.SparseDynamicG1!, endogenous, exogenous, periods, temporary_var, params, steady_state, rowval, r1, r2, n, n1, n2)                  
