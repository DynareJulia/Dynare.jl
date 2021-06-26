using Dynare
using SparseArrays
using Test

include("../src/perfectforesight/gmres_solver.jl")

nvar = 20
nfwrd = 5
npred = 6
y = zeros(nvar, npred)
a = randn(nvar, 15)
b = randn(9, npred)
ic = 11:15
ir = 1:2:9
tmp1 = zeros(nvar, nfwrd)
tmp2 = zeros(nfwrd, npred)
CmultG!(y, tmp1, tmp2, a, b, collect(ic), collect(ir))
@test y ≈ a[:, ic] * b[ir, :]

context = @dynare "test/models/example1/example1.mod"

md = context.models[1]
nvar = md.endogenous_nbr
periods = 4
work = context.work
steadystate = context.results.model_results[1].trends.endogenous_steady_state
endogenous = repeat(steadystate, periods)
exogenous = zeros(periods, md.exogenous_nbr)

J = Jacobian(context, periods)
A = makeJacobian!(J, endogenous, exogenous, context, periods)

a = zeros(6, 6)
b = zeros(6, 6)
c = zeros(6, 6)
a[:, [3, 4, 6]] .= work.jacobian[:, 1:3]
b .= work.jacobian[:, 4:9]
c[:, [1, 2, 6]] = work.jacobian[:, 10:12]

target =
    vcat([b c zeros(6, 12)], [a b c zeros(6, 6)], [zeros(6, 6) a b c], [zeros(6, 12) a b])
g = context.results.model_results[1].linearrationalexpectations.g1_1
target[19:24, 18 .+ [3, 4, 6]] .+= c * g
@test A ≈ target
