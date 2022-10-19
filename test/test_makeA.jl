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

J = Dynare.makeJacobian(A.colptr, A.rowval, 5, 4, [])

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

@show colptr
J = Dynare.makeJacobian(colptr,
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
Dynare.updateJacobian!(J, df.SparseDynamicG1!, endogenous, exogenous, periods, temporary_var, params, steady_state, colptr, nzval, m.endogenous_nbr, m.exogenous_nbr)
@show nzval
@test J == AA

permutations =[(1,3), (2, 5)]
sort!(permutations, by=x->x[1])

rowval1 = copy(rowval)
permutations1 = Tuple{Int64, Int64}[]
k = 1
for i in 1:length(colptr) - 1
    if colptr[i + 1] > colptr[i]
        vr = view(rowval1, colptr[i]:colptr[i+1]-1)
        @show vr
        for p in permutations
            p1, p2 = p
            p0 = Tuple{Int64, Int64}[]
            for (j,r) in enumerate(vr)
                if r == p1
                    vr[j] = p2
                    push!(p0, p)
                elseif r == p2
                    vr[j] = p1
                    push!(p0, p)
                    break
                end
            end
            if !issorted(vr)
                sort!(vr)
                np0 = length(p0)
                if np0 > 0
                    push!(permutations1, colptr[i] - 1 .+ p0[1])
                end
                if np0 > 1
                    push!(permutations1, colptr[i] - 1 .+ (p0[2][2], p0[2][1]))
                end
            end
        end
    end
end


nzval1 = copy(J.nzval)
