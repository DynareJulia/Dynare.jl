using PkgBenchmark

import Dynare

myscript = "/home/michel/projects/julia-packages/dynare.jl/benchmark/benchmark1.jl"
benchmarkpkg(Dynare; script=myscript)
benchmarkpkg(Dynare, "29d7476"; script=myscript)
benchmarkpkg(Dynare, "5625c09"; script=myscript)
#benchmarkpkg(Dynare, "1e788c8"; script=myscript)
