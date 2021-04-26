module Dynare
using Plots
using GR
using Revise

gr()

include("dynare_containers.jl")
include("model.jl")
export get_abc, get_de
include("symboltable.jl")
include("initialization.jl")
include("deterministic_trends.jl")
include("DynareParser.jl")
export parser, get_jacobian_at_steadystate!
include("DynarePreprocessor.jl")
export dynare_preprocess
include("steady_state/SteadyState.jl")
export steady_state!
include("dynare_table.jl")
include("reporting/report.jl")
include("graphics.jl")
include("data.jl")
include("filters/kalman/kalman.jl")
include("optimal_policy.jl")
include("perturbations.jl")

export @dynare

mutable struct CommandLineOptions
    compilemodule::Bool
    function CommandLineOptions()
        compilemodule = true
        new(compilemodule)
    end
end
        
macro dynare(modfilename, args...)
    modfilename = string(modfilename)
    if  occursin(r"\.mod$", modfilename)
        modname = modfilename[1:length(modfilename)-4]
    else
        modname = modfilename
    end
    options = CommandLineOptions()
    for (i, a) in enumerate(args)
        if a == "nocompile"
            options.compilemodule = false
            deleteat!(args, i)
        end
    end
    dynare_preprocess(modname*".mod", args)
    context = parser(modname, options)
    return context
end

end # module
