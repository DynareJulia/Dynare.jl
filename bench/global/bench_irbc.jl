using CpuId
using Dynare
using Statistics
using Tasmanian

function bench1(ncountries, depth)
    context = dynare("irbc1", "-DN=$ncountries");
    (grid, state_variables, policy_variables) = sparsegridapproximation(scaleCorrExclude=["lambda"], gridDepth = depth, maxRef = 0, tol_ti=0.001);
    # skipping JIT compilation time of first iteration
    return mean([t.value for (i, t) in enumerate(context.timings["sparsegrids"]) if i > 1]), grid
end

function bench()
    depth = 4
    open("bench_irbc.md", "w+") do outfile
        buffer = IOBuffer()
        println(buffer, "# Dynare sparsegrid benchmark")
        println(buffer, cpubrand(), cpucores())
        println(buffer, "Number of threads: $(Threads.nthreads())")

        try
            ompthreads = ENV["OMP_NUM_THREADS"]
            println(buffer, "Number of OMP threads: $ompthreads")
        catch
        end

        for N in [2, 3 ] #4, 8, 12, 20]
            m, grid = bench1(N, depth)
            println(buffer, "$N countries, depth $depth, $(getNumPoints(grid)) grid points, mean iteration time: $m milliseconds")
        end
        write(outfile, String(take!(buffer)))
    end
end

bench()
