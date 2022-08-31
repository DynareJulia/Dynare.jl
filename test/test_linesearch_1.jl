using Dynare
using Test
include("../src/perfectforesight/perfectforesight_solver_1.jl")
include("../src/perfectforesight/linesearch.jl")

context = @dynare "test/models/cgg/cgg.mod"

md = context.models[1]
results = context.results.model_results[1]
n = md.endogenous_nbr
p = md.exogenous_nbr
A = zeros(n, n)
A[:, md.i_bkwrd_b] .= results.linearrationalexpectations.g1_1
B = copy(results.linearrationalexpectations.g1_2)
periods = 102
exogenous = zeros(periods, p)
exogenous[2, :] = [0.02, 0.02]
window = 1
y = zeros(n, periods)
tmp1 = similar(y)
tmp2 = similar(y)
c = context.results.model_results[1].trends.endogenous_steady_state
simul_first_order_1!(y, zeros(n), A, B, exogenous, window, periods, tmp1, tmp2)

steadystate = context.results.model_results[1].trends.endogenous_steady_state
y .+= steadystate
residuals = zeros(n * (periods - 2))
dynamic_variables = zeros(7)
temp_vec = context.work.temporary_values

initialvalues = y[1:n]
terminalvalues = y[(periods-1)*n+1:periods*n]
periods = 100
tmp = zeros(n * periods)
dy = zeros(n * periods)
params = context.work.params
x = view(vec(y), n+1:(periods+1)*n)
get_residuals!(
    residuals,
    x,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
)
@show residuals[1:12]
@show 0.5 * dot(residuals, residuals)
#=
JJ = Jacobian(context, periods)
ws_threaded = [JacTimesVec(context) for i=1:Threads.nthreads()]
A = makeJacobian!(JJ, vec(y), exogenous, context, periods, ws_threaded)
τ = 0.01
my_lu = ilu(A, τ = τ)
fill!(tmp, 0.0)
Dynare.gmres!(dy, A, residuals, Pr=my_lu, verbose = true, maxiter=4)

stpmx = 100
stpmax = stpmx*max(sqrt(dot(y, y)), n*periods)
tolx = 1e-5
g = zeros(n*periods)
mul!(g, transpose(A), vec(residuals))
x0 = view(vec(y), n+1:(periods+1)*n)
x = zeros(n*periods)
params = context.work.params

linesearch!(x, residuals, x0, g, -dy, stpmax, get_residuals!, 1e-8,
            initialvalues, terminalvalues, exogenous,
            dynamic_variables, steadystate, params, md, periods,
            temp_vec)

=#
