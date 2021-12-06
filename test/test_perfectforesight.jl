using Dynare
include("../src/perfectforesight/perfectforesight.jl")

context = @dynare "dynare.jl/test/models/example1pf/example1pf"

periods = 10
datafile = ""
perfect_foresight_ws = PerfectForesightWs(context, periods)

perfect_foresight_initialization!(context, periods, datafile, perfect_foresight_ws)
