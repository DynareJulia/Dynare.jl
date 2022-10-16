using Dynare
using SparseArrays

Dynare.dynare_preprocess("models/example1/example1_sparse.mod", [])

DF= Dynare.DFunctions.load_model_functions("models/example1/example1")

context = @dynare "models/example1/example1_sparse";

m = context.models[1]
df = context.dynarefunctions
endo_nbr = m.endogenous_nbr
exo_nbr = m.exogenous_nbr
ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
tmp_nbr = sum(m.dynamic_tmp_nbr[1:2])
nzval_nbr = length(m.dynamic_g1_sparse_rowval)

endogenous = context.results.model_results[1].trends.endogenous_steady_state
exogenous = zeros(2,1)
steadystate = endogenous




ws = Dynare.DynamicWs(endo_nbr, exo_nbr, tmp_nbr, m.dynamic_g1_sparse_colptr, m.dynamic_g1_sparse_rowval)
T = zeros(tmp_nbr)
resid = zeros(endo_nbr)
y = repeat(steadystate, 3)
x = zeros(exo_nbr)

nzval = similar(m.dynamic_g1_sparse_rowval, Float64)
g1 = SparseMatrixCSC(endo_nbr, 3*endo_nbr + exo_nbr, m.dynamic_g1_sparse_colptr, m.dynamic_g1_sparse_rowval, nzval)
Dynare.DFunctions.dynamic!(T, resid, y, x, context.work.params, steadystate)
Dynare.DFunctions.dynamic!(T, resid, g1, y, x, context.work.params, steadystate)
