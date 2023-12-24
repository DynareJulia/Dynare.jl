module Dynare

using Reexport
@reexport using ExtendedDates

using LoggingExtras
using NLsolve
using NonlinearSolve
using Printf

Base.@kwdef struct CommandLineOptions
    compilemodule::Bool = true
    stoponerror::Bool = false
end

using LinearRationalExpectations
const LRE = LinearRationalExpectations

using KalmanFilterTools
using KOrderPerturbations

include("logging.jl")
include("utils.jl")
include("dynare_functions.jl")
include("dynare_containers.jl")
include("accessors.jl")
export forecast, irf, simulation, smoother 
include("model.jl")
export get_abc, get_de
include("symboltable.jl")
include("data.jl")
include("distributions/distribution_parameters.jl")
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
export Report, add_graph!, add_model!, add_page!, add_paragraph!, add_table! 
include("graphics.jl")
export plot_forecast, plot_recursive_forecast
include("filters/kalman/kalman.jl")
export smoother
include("optimal_policy.jl")
include("perturbations.jl")
export localapproximation!
include("perfectforesight/perfectforesight.jl")
export perfect_foresight!, scenario!
include("estimation/priorprediction.jl")
export priorprediction
include("simulations.jl")
include("estimation/estimation.jl")
export covariance, mode_compute!, output_MCMCChains, plot_MCMCChains, plot_priors 
export plot_prior_posterior, prior!, rwmh_compute!, sms_compute!
export calibsmoother!
include("forecast.jl") 
export forecasting!, recursive_forecasting!
export @dynare, dynare
include("macros.jl")
export limits!, @limits
include("global/sparsegrids.jl")
export sparsegridapproximation

macro dynare(modfile_arg::String, args...)
    dynare(modfile_arg, args...)
end

function dynare(modfile_arg::String, args...)
    modname = get_modname(modfile_arg)
    set_logging(modname)
    @info "Dynare version: $(module_version(Dynare))"
    @info "$(now()): Starting @dynare $modfile_arg"
    arglist = []
    compilemodule = true
    preprocessing = true
    stoponerror = false
    for (i, a) in enumerate(args)
        if a == "nocompile"
            compilemodule = false
        elseif a == "nopreprocessing"
            preprocessing = false
        elseif a == "stoponerror"
            stoponerror = true
        else
            push!(arglist, a)
        end
    end
    if preprocessing
        modfilename = modname * ".mod"
        try
            dynare_preprocess(modfilename, arglist)
        catch
            return []
        end
    end
    @info "$(now()): End of preprocessing"

    # onlymodel option performs only preprocessing
    "onlymodel" in arglist && return nothing
    
    options = CommandLineOptions(compilemodule, stoponerror)
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
