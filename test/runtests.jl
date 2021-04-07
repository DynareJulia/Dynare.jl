using Dynare

@dynare "models/example1/example1.mod"
@dynare "models/example2/example2.mod"
@dynare "models/example3/example3.mod"

@dynare "models/cgg/cgg_ramsey.mod"

# The following needs support for “deterministic_trends” in the preprocessor.
# On 2021-04-07, this feature is only present on Michel’s private preprocessor branch.
#@dynare "models/stochastic_trend_drift/trend1.mod"
