module Dynare
using Revise


@Base.kwdef struct CommandLineOptions
    compilemodule::Bool = true
end

include("utils.jl")
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
#include("perfectforesight/perfectforesight_solvers.jl")
export @dynare

macro dynare(modfile_arg::String, args...)
    modname = get_modname(modfile_arg)
    arglist = []
    compilemodule = true
    for (i, a) in enumerate(args)
        if a == "nocompile"
            compilemodule = false
        else
            push!(arglist, a)
        end
    end
    options = CommandLineOptions(compilemodule)
    modfilename = modname*".mod"
    dynare_preprocess(modfilename, arglist)
    context = parser(modname, options)
    return context
end

function get_modname(s::String)
    modfilename = string(s)
    mdoname = split(modfilename, ".")[1]
    if  occursin(r"\.mod$", modfilename)
        modname::String = modfilename[1:length(modfilename)-4]
    else
        modname = modfilename
    end
    return modname
end

include("precompile_Dynare.jl")
_precompile_()
end # module
