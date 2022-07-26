using Dynare
using Test

A = randn(4, 4)
B = randn(4, 3)
x = randn(5, 3)
results = randn(4, 5)
tmp1 = similar(results)
tmp2 = similar(results)
results_orig = copy(results)
z = results
Dynare.geometric_series!(z, A, B, x)

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
Dynare.simul_first_order!(results, initial_values, c, A, B, x)

z = zeros(4, 5)
Dynare.geometric_series!(z, A, B, x)
for i = 2:5
    @test results[:, i] ≈ c .+ A * (results[:, i-1] - c) .+ z[:, i]
end
for i = 6:6
    @test results[:, i] ≈ c .+ A * (results[:, i-1] - c)
end

x = [[0, 0.1] [-0.1, 0] [0, 0]]
Dynare.simul_first_order!(results, initial_values, c, A, B, x)
@test results[:, 1] ≈ c + A * (ones(4) - c) + B * x[1, :] + A * B * x[2, :]
@test results[:, 2] ≈ c + A * (results[:, 1] - c) + B * x[2, :]
@test results[:, 3] ≈ c + A * (results[:, 2] - c)
@test results[:, 5] ≈ c + A * A * A * (results[:, 2] - c)
