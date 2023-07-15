module Dynare

using ExtendedDates
using Logging
using Printf

Base.@kwdef struct CommandLineOptions
    compilemodule::Bool = true
end

using LinearRationalExpectations
const LRE = LinearRationalExpectations

using KalmanFilterTools
using KOrderPerturbations

include("utils.jl")
include("dynare_functions.jl")
include("dynare_containers.jl")
include("accessors.jl")
export irf, simulation
include("model.jl")
export get_abc, get_de
include("symboltable.jl")
include("data.jl")
include("distributions/distribution_parameters.jl")
include("estimation/parse_prior.jl")
include("initialization.jl")
include("deterministic_trends.jl")
include("DynareParser.jl")
export parser
include("DynarePreprocessor.jl")
export dynare_preprocess
include("steady_state/SteadyState.jl")
export steadystate!
include("dynare_table.jl")
export round
include("reporting/report.jl")
include("graphics.jl")
include("filters/kalman/kalman.jl")
include("optimal_policy.jl")
include("perturbations.jl")
include("perfectforesight/perfectforesight.jl")
include("estimation/priorprediction.jl")
include("simulations.jl")
include("nonlinear/NLsolve.jl")
using .NLsolve
include("estimation/estimation.jl")
export mh_estimation

export @dynare, dynare

macro dynare(modfile_arg::String, args...)
    dynare(modfile_arg, args...)
end

function dynare(modfile_arg::String, args...)
    @info "Dynare version: $(module_version(Dynare))"
    modname = get_modname(modfile_arg)
    @info "$(now()): Starting @dynare $modfile_arg"
    arglist = []
    compilemodule = true
    preprocessing = true
    for (i, a) in enumerate(args)
        if a == "nocompile"
            compilemodule = false
        elseif a == "nopreprocessing"
            preprocessing = false
        else
            push!(arglist, a)
        end
    end
    if preprocessing
        modfilename = modname * ".mod"
        dynare_preprocess(modfilename, arglist)
    end
    @info "$(now()): End of preprocessing"

    # onlymodel option performs only preprocessing
    "onlymodel" in arglist && return nothing
    
    options = CommandLineOptions(compilemodule)
    context = parser(modname, options)
    return context
end

function get_modname(modfilename::String)
    if occursin(r"\.mod$", modfilename)
        modname::String = modfilename[1:length(modfilename)-4]
    else
        modname = modfilename
    end
    return modname
end

if !isdefined(Base, :get_extension)
  include("../ext/PardisoSolver.jl")
  include("../ext/PathSolver.jl")
end

using PrecompileTools
@compile_workload begin
    # redirect output to avoid startling a user during precompilation 
    redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        modelpath = dirname(@__DIR__) * "/test/models/example1/example1"
        dynare(modelpath)
    end
    end
end

#include("precompile_Dynare.jl")
#_precompile_()
end # module
