using Dynare

include("test_initialization.jl")
include("test_nonlinear_initialization.jl")

@dynare "models/example1/example1.mod"
@dynare "models/example2/example2.mod"
@dynare "models/example3/example3.mod"
@dynare "models/example3report/example3report.mod"
@dynare "models/example3ss/example3ss.mod"
@dynare "models/cgg/cgg_ramsey.mod"
@dynare "models/stochastic_trend_drift/trend1.mod"
@dynare "models/example1pf/example1pf"
