using BenchmarkTools
using Dynare

SUITE = BenchmarkGroup()

SUITE["irbc2"] = BenchmarkGroup(["irbc", "N2"])
SUITE["irbc20"] = BenchmarkGroup(["irbc", "N20"])
SUITE["irbc100"] = BenchmarkGroup(["irbc", "N100"])

SUITE["irbc2"]["run"] = @benchmarkable @dynare "../test/models/irbc/irbc2";
SUITE["irbc20"]["run"] = @benchmarkable @dynare "../test/models/irbc/irbc20";
SUITE["irbc100"]["run"] = @benchmarkable @dynare "../test/models/irbc/irbc100";
