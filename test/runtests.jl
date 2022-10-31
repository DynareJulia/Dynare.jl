using Dynare
using Test

include("test_initialization.jl")
include("test_nonlinear_initialization.jl")
include("test_steadystate.jl")

context = @dynare "models/example1/example1.mod"
irfs = irf()
@test length(irfs) == 2
irfs_e = irf(:e)
@test size(irfs_e) == (40, 6)
irfs_u = irf(:u)
@test size(irfs_u) == (40, 6)
irfs_e_y = irf(:e, :y)
@test Dynare.periods(irfs_e) ==
      range(Dynare.UndatedDate(1), stop = Dynare.UndatedDate(40), step = Dynare.Undated(1))
@test Dynare.dataframe(irfs_e_y) ==
      Dynare.dataframe(context.results.model_results[1].irfs[:e])[!, [:y]]
irfs_u_b = irf(:u, :b)
@test Dynare.dataframe(irfs_u_b) ==
      Dynare.dataframe(context.results.model_results[1].irfs[:u])[!, [:b]]
irfs_u_b = irf("u", "b")
@test Dynare.dataframe(irfs_u_b) ==
      Dynare.dataframe(context.results.model_results[1].irfs[:u])[!, [:b]]
@dynare "models/example2/example2.mod"
@dynare "models/example3/example3.mod"
@dynare "models/example3ss/example3ss.mod"
@dynare "models/example3ss/example3ss_analytical.mod"
@dynare "models/example3ss/example3ss_partial.mod"
@dynare "models/example3report/example3report.mod"
@dynare "models/cgg/cgg_ramsey.mod"
#@dynare "models/stochastic_trend_drift/trend1.mod"
context = @dynare "models/example1pf/example1pf"
sim = Dynare.simulation()
@test length(sim) == 1
@test size(Dynare.dataframe(sim[1].data)) == (300, 6)
@test Dynare.periods(sim[1].data) ==
      range(Dynare.UndatedDate(1), stop = Dynare.UndatedDate(300), step = Dynare.Undated(1))
sim_a = Dynare.simulation(:a)
@test Dynare.dataframe(sim_a) ==
      Dynare.dataframe(context.results.model_results[1].simulations[1].data)[!, [:a]]
sim_a = Dynare.simulation("a")
@test Dynare.dataframe(sim_a) ==
      Dynare.dataframe(context.results.model_results[1].simulations[1].data)[!, [:a]]
