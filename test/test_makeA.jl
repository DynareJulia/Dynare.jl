using Dynare
using LinearAlgebra
using SparseArrays
using Test


A  = sprand(5, 15, 0.1)
AA = vcat(hcat(A[:, 6:15], zeros(5, 10)),
          hcat(A, zeros(5, 5)),
          hcat(zeros(5, 5), A),
          hcat(zeros(5, 10), A[:, 1:10]))
display(AA[1:10,1:10])

(J, permutations1) = Dynare.makeJacobian(A.colptr, A.rowval, 5, 4, [])

@test J.colptr == AA.colptr
@test J.rowval == AA.rowval

permutations =[(1,3), (2, 5)]
sort!(permutations, by=x->x[1])

(permutations1, rowval1) = Dynare.permutation(permutations, A.colptr, A.rowval)
A1 = A[[3, 5, 1, 4, 2], :]
@test rowval1 == A1.rowval
@show permutations1

AP = copy(A)
AP.rowval .= rowval1
for p in permutations1
    p1, p2 = p
    @show p1, p2
    @show AP.nzval[p1], A.nzval[p2]
    AP.nzval[p1] = A.nzval[p2]
end
@test AP == A1

@show A.colptr


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

@show colptr
(J, permutations1) = Dynare.makeJacobian(colptr,
                                         rowval,
                                         m.endogenous_nbr,
                                         periods,
                                         [])

r1 = colptr[m.endogenous_nbr + 1]
n = colptr[3*m.endogenous_nbr + 1] - 1
n1 = colptr[3*m.endogenous_nbr  + 1] - colptr[m.endogenous_nbr + 1] + 1
n2 = colptr[2*m.endogenous_nbr + 1]
df = Dynare.DFunctions
nzval = Vector{Float64}(undef, colptr[end] - 1)
@show steady_state
residual = zeros(m.endogenous_nbr)

#df.SparseDynamicG1!(temporary_var, nzval, endogenous[1:18], exogenous[1:2], params, steady_state, false)
A = SparseMatrixCSC(m.endogenous_nbr, 3*m.endogenous_nbr + m.exogenous_nbr, colptr, rowval, nzval)
df.dynamic!(temporary_var, residual, A, endogenous[1:18], exogenous[1:2], params, steady_state) 
AA = vcat(hcat(A[:, 7:18], zeros(6, 12)),
          hcat(A[:,1:18], zeros(6, 6)),
          hcat(zeros(6, 6), A[:,1:18]),
          hcat(zeros(6, 12), A[:, 1:12]))

@test J.colptr == AA.colptr
@test J.rowval == AA.rowval
@show nzval
Dynare.updateJacobian!(J, df.SparseDynamicG1!, endogenous, exogenous, periods, temporary_var, params, steady_state, colptr, nzval, m.endogenous_nbr, m.exogenous_nbr, [], [])
@show nzval
@test J == AA

permutations =[(1,3), (2, 5)]
sort!(permutations, by=x->x[1])

(J1, permutations1) = Dynare.makeJacobian(A.colptr,
                                          A.rowval,
                                          m.endogenous_nbr,
                                          periods,
                                          permutations)
ws = similar(A.nzval)
Dynare.updateJacobian!(J1, df.SparseDynamicG1!, endogenous, exogenous, periods, temporary_var, params, steady_state, colptr, A.nzval, m.endogenous_nbr, m.exogenous_nbr, permutations1, ws)
let kk = []
    for i = 1:periods
        kk = vcat(kk, [3, 5, 1, 4, 2, 6] .+ (i - 1)*m.endogenous_nbr)
    end
    @show kk
    global J2 = J[kk, :]
end

@test J1.colptr == J2.colptr
@test J1.rowval == J2.rowval
@test J1.nzval == J2.nzval
@test J1 == J2 

