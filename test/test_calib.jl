using CSV
using Dynare
using Test

context = @dynare "models/stochastic_trend_drift.mod/calib1";
results = context.results.model_results[1]
data = CSV.File("data.csv")
Yorig = [data.y data.infl]'
Yorig = Yorig[:, 2:end]
Y = Matrix{Float64}(undef, size(Yorig))
varobs_ids = [1, 4]
Dynare.remove_linear_trend!(
    Y,
    Yorig,
    results.trends.endogenous_steady_state[varobs_ids],
    results.trends.endogenous_linear_trend[varobs_ids],
)

@test Y[:, 1:3] ≈ context.results.model_results[1].smoother["alphah"][[1, 4], 1:3]
