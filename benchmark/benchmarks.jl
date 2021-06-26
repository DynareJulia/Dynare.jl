using BenchmarkTools
using Dynare

SUITE = BenchmarkGroup()

SUITE["example1"] = BenchmarkGroup(["example1", "order1", "cycle reduction"]) 
SUITE["example2"] = BenchmarkGroup(["example2", "order1", "schur", "auxiliary variable"]) 
SUITE["example3"] = BenchmarkGroup(["example3", "order1", "schur", "simulation", "smoother", "csv"]) 
SUITE["cgg_ramsey"] = BenchmarkGroup(["cgg", "ramsey"]) 
SUITE["trend1"] = BenchmarkGroup(["trend1", "nonstationary"]) 

SUITE["example1"]["run"] = @benchmarkable @dynare "test/models/example1/example1";
SUITE["example2"]["run"] = @benchmarkable @dynare "test/models/example2/example2";
SUITE["example3"]["run"] = @benchmarkable @dynare "test/models/example3/example3";
SUITE["cgg_ramsey"]["run"] = @benchmarkable @dynare "test/models/cgg/cgg_ramsey.mod"
SUITE["trend1"]["run"] = @benchmarkable @dynare "test/models/stochastic_trend_drift/trend1.mod"



