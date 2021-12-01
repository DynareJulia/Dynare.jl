using BenchmarkTools
using Test
include("perfectforesight_solver_1.jl")

#tests
A = randn(4, 4)
B = randn(4, 3)
x = randn(5, 3)
results = randn(4, 5)
tmp1 = similar(results)
tmp2 = similar(results)
results_orig = copy(results)
z = results
geometric_series!(z, A, B, x)

@test z[:, 5] ≈ B * x[5, :]
@test z[:, 4] ≈ B * x[4, :] + A * B * x[5, :]
@test z[:, 3] ≈ B * x[3, :] + A * B * x[4, :] + A * A * B * x[5, :]
@test z[:, 2] ≈
      B * x[2, :] + A * B * x[3, :] + A * A * B * x[4, :] + A * A * A * B * x[5, :]
@test z[:, 1] ≈
      B * x[1, :] +
      A * B * x[2, :] +
      A * A * B * x[3, :] +
      A * A * A * B * x[4, :] +
      A * A * A * A * B * x[5, :]

initial_values = ones(4)
c = randn(4)
periods = 6
results = rand(4, periods)
tmp1 = similar(results)
tmp2 = similar(results)
simul_first_order_1!(results, initial_values, c, A, B, x)

z = zeros(4, 5)
geometric_series!(z, A, B, x)
for i = 2:5
    @test results[:, i] ≈ c .+ A * (results[:, i-1] - c) .+ z[:, i]
end
for i = 6:6
    @test results[:, i] ≈ c .+ A * (results[:, i-1] - c)
end

md = context.models[1]
results = context.results.model_results[1]
n = md.endogenous_nbr
p = md.exogenous_nbr
A = zeros(n, n)
A[:, md.i_bkwrd_b] .= results.linearrationalexpectations.g1_1
B = copy(results.linearrationalexpectations.g1_2)
periods = 102
exogenous = zeros(periods, p)
#exogenous[2, :] = [0.02, 0.02]
window = 2
y = zeros(n, periods)
tmp1 = similar(y)
tmp2 = similar(y)
c = context.results.model_results[1].trends.endogenous_steady_state
simul_first_order_1!(y, zeros(6), A, B, exogenous)
steadystate = context.results.model_results[1].trends.endogenous_steady_state
y .+= steadystate
initialvalues = steadystate
terminalvalues = y[:, periods]
residuals = zeros(n, periods - 2)
dynamic_variables = zeros(12)
temp_vec = context.work.temporary_values

preconditioner_window = 10
algo = "CR"
rout = zeros((periods - 2) * md.endogenous_nbr)
tmp = similar(rout)
JJ = Jacobian(context, periods - 2)
endo_nbr = 6
exo_nbr = 2
dynamic_nbr = 12
tmp_nbr = length(temp_vec)

ws_threaded = [Dynare.DynamicWs(endo_nbr,
                                exo_nbr,
                                dynamic_nbr,
                                tmp_nbr)
               for i = 1:Threads.nthreads()]

params = context.work.params
@btime begin
    simul_first_order_1!(y, zeros(6), A, B, exogenous)
    y .= steadystate
    solve1(
        vec(residuals),
        y,
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        md,
        periods - 2,
        temp_vec,
        JJ,
        context,
        ws_threaded,
        n,
    )
end

@btime begin
    ws = GmresWs(periods - 2, preconditioner_window, context, algo)
    y .+= steadystate
    solve2(
        residuals,
        y,
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        md,
        periods - 2,
        temp_vec,
        JJ,
        context,
        ws_threaded,
        n,
        ws,
    )
end

@btime begin
    simul_first_order_1!(y, zeros(6), A, B, exogenous)
    y .= steadystate
    solve3(
        0.01,
        residuals,
        y,
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        md,
        periods - 2,
        temp_vec,
        JJ,
        context,
        ws_threaded,
        n,
    )
end

@btime begin
    simul_first_order_1!(y, zeros(6), A, B, exogenous)
    y .= steadystate
    solve4(
        residuals,
        y,
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        md,
        periods - 2,
        temp_vec,
        JJ,
        context,
        ws_threaded,
        n,
    )
end

using Pardiso
ps = PardisoSolver()

@btime begin
    simul_first_order_1!(y, zeros(6), A, B, exogenous)
    y .= steadystate
    solve5(
        residuals,
        y,
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        md,
        periods - 2,
        temp_vec,
        JJ,
        context,
        ws_threaded,
        n,
        ps,
    )
end
