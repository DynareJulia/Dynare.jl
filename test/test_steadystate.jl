using Dynare
using Test

context1 = @dynare "models/example3ss/example3ss_analytical.mod"
Dynare.compute_steady_state!(context1, Dict{String,Any}())
target = context1.results.model_results[1].trends.endogenous_steady_state

context2 = @dynare "models/example3ss/example3ss.mod"
@test isapprox(context2.results.model_results[1].trends.endogenous_steady_state, target)

context3 = @dynare "models/example3ss/example3ss_partial.mod"
@test isapprox(context3.results.model_results[1].trends.endogenous_steady_state, target)
