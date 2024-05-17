using Dynare
using Test

# CSV.write fails with AxisArrays
#include("test_initialization.jl")
#include("test_nonlinear_initialization.jl")
include("test_steadystate.jl")
include("test_scenario.jl")
include("test_deterministic_trends.jl")

context = @dynare "models/example1/example1.mod"
irfs = irf()
@test length(irfs) == 2
irfs_e = irf(:e)
@test size(irfs_e) == (40, 6)
irfs_u = irf(:u)
@test size(irfs_u) == (40, 6)
irfs_e_y = irf(:e, :y)
#@test Dynare.periods(irfs_e) ==
#      range(Dynare.UndatedDate(1), stop = Dynare.UndatedDate(40), step = Dynare.Undated(1))
#@test Dynare.dataframe(irfs_e_y) ==
#      Dynare.dataframe(context.results.model_results[1].irfs[:e])[!, [:y]]
#irfs_u_b = irf(:u, :b)
#@test Dynare.dataframe(irfs_u_b) ==
#      Dynare.dataframe(context.results.model_results[1].irfs[:u])[!, [:b]]
#irfs_u_b = irf("u", "b")
#@test Dynare.dataframe(irfs_u_b) ==
#      Dynare.dataframe(context.results.model_results[1].irfs[:u])[!, [:b]]
@dynare "models/example2/example2.mod"
# CSV.write fails with AxisArrays
#@dynare "models/example3/example3.mod"
@dynare "models/example3ss/example3ss.mod"
@dynare "models/example3ss/example3ss_analytical.mod"
@dynare "models/example3ss/example3ss_partial.mod"
@dynare "models/example3report/example3report.mod"
@dynare "models/cgg/cgg_ramsey.mod"
#@dynare "models/stochastic_trend_drift/trend1.mod"

context = @dynare "models/example1pf/example1pf"
sim = Dynare.simulation()
@test length(sim) == 1
@test size(sim[1].data) == (1002, 8)
@test Dynare.AxisArrayTables.AxisArrays.axes(getfield(sim[1].data, :data), 1).val ==
      range(Dynare.Undated(0), stop = Dynare.Undated(1001), step = Dynare.Undated(1))
sim_a = Dynare.simulation(:a)
@test sim_a ==
      context.results.model_results[1].simulations[1].data.a[2:end-1]
sim_a = Dynare.simulation("a")
@test sim_a ==
      context.results.model_results[1].simulations[1].data.a[2:end-1]
trends = context.results.model_results[1].trends
@test Dynare.simulation(1)[1000, :] ≈ transpose(vcat(trends.endogenous_steady_state, trends.exogenous_steady_state))

context = @dynare "models/example1pf/example1pf_endval"
@test Dynare.simulation(1)[300, 1:6] ≈ transpose(context.results.model_results[1].trends.endogenous_terminal_steady_state)
context = @dynare "models/initialization/neoclassical1"
context = @dynare "models/initialization/neoclassical5"
using PATHSolver
context = @dynare "models/irreversible/irbc2a"
@test all(Matrix(simulation("i"))[10:52] .≈ 0)
