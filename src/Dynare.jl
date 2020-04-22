module Dynare

include("model.jl")
export get_abc, get_de
include("parser/DynareParser.jl")
export parser, get_jacobian_at_steadystate!
include("parser/DynarePreprocessor.jl")
export dynare_preprocess
include("steady_state/SteadyState.jl")
export steady_state!
include("Utils.jl")
export get_power_deriv

export @dynare

macro dynare(modfilename, args...)
    if  occursin(r"\.mod$", modfilename)
        modname = modfilename[1:length(modfilename)-4]
    else
        modname = modfilename
    end
    dynare_preprocess(modname*".mod", args)
    context = Dynare.Context()
    parser(modname, context)
    return context
end

end # module
