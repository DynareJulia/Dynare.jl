using Dynare
using Test

context1 = @dynare "models/example3/example3.mod"
Dynare.compute_steady_state!(context1)
target = context1.results.model_results[1].trends.endogenous_steady_state

context2 = @dynare "models/example3ss/example3ss.mod"
x0 = Float64.(vec(view(context2.work.initval_endogenous, 1, :)))
Dynare.solve_steady_state!(context2, x0)
@test context1.results.model_results[1].trends.endogenous_steady_state ≈ target


