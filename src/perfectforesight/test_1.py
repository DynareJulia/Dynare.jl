using Test

#tests
A = randn(4, 4)
B = randn(4, 3)
x = randn(5, 3)
window = size(x, 1)
results = randn(4, 5)
tmp1 = similar(results)
tmp2 = similar(results)
results_orig = copy(results)
z = results
geometric_series!(z, A, B, x, window, tmp1, tmp2)

@test z[:, 5] ≈ B*x[5, :]
@test z[:, 4] ≈ B*x[4, :] + A*B*x[5, :]
@test z[:, 3] ≈ B*x[3, :] + A*B*x[4, :] + A*A*B*x[5, :]
@test z[:, 2] ≈ B*x[2, :] + A*B*x[3, :] + A*A*B*x[4, :] + A*A*A*B*x[5, :]
@test z[:, 1] ≈ B*x[1, :] + A*B*x[2, :] + A*A*B*x[3, :] + A*A*A*B*x[4, :] + A*A*A*A*B*x[5, :]

initial_values = ones(4)
c = randn(4)
periods = 8
results = rand(4, periods)
tmp1 = similar(results)
tmp2 = similar(results)
simul_first_order_1!(results, initial_values, x, c, A, B, window, periods, tmp1, tmp2)

z = zeros(4, 5)
geometric_series!(z, A, B, x, window, tmp1, tmp2)
@test results[:, 1] == initial_values .- c
for i = 2:6
    @test results[:, i] ≈ c .+ A*results[:, i-1] .+ z[:, i-1]
end
for i = 7:8
    @test results[:, i] ≈ c .+ A*results[:, i-1]
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
exogenous[2, :] = [0.02, 0.02]
window = 2
y = zeros(n, periods)
tmp1 = similar(y)
tmp2 = similar(y)
c = context.results.model_results[1].trends.endogenous_steady_state
#simul_first_order_1!(y, zeros(6), exogenous, zeros(6), A, B, window, periods, tmp1, tmp2)
steadystate = context.results.model_results[1].trends.endogenous_steady_state
y .+= steadystate
residuals = zeros(n, periods - 2)
dynamic_variables = zeros(12)
temp_vec = context.work.temporary_values
