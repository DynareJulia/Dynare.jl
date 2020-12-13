module Dynare
using Plots
using GR

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
include("Utils.jl")
export get_power_deriv
include("dynare_table.jl")
include("reporting/report.jl")
include("graphics.jl")
include("data.jl")
include("filters/kalman/kalman.jl")
include("optimal_policy.jl")
include("perturbations.jl")

export @dynare

macro dynare(modfilename, args...)
    modfilename = string(modfilename)
    if  occursin(r"\.mod$", modfilename)
        modname = modfilename[1:length(modfilename)-4]
    else
        modname = modfilename
    end
    dynare_preprocess(modname*".mod", args)
    context = parser(modname)
    return context
end

end # module
