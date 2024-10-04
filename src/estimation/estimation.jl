include("smc.jl")

using AdvancedMH
using AbstractMCMC
using Distributions
using FillArrays
using FiniteDiff: finite_difference_gradient, finite_difference_hessian
using .Iterators
using KernelDensity
using LogDensityProblems
using MCMCChains
using MCMCDiagnosticTools
using Metaheuristics
using Optim
using Plots
using PolynomialMatrixEquations: UndeterminateSystemException, UnstableSystemException
using Random
using Serialization
using StatsPlots
using TransformVariables
using TransformedLogDensities

import LogDensityProblems: dimension, logdensity, logdensity_and_gradient, capabilities

include("data.jl")
include("priors.jl")

Base.@kwdef struct EstimationOptions
    config_sig::Float64 = 0.8
    data::AxisArrayTable = AxisArrayTable(AxisArrayTables.AxisArray(Matrix{Float64}(undef, 0, 0), PeriodsSinceEpoch[], Symbol[]))
    datafile::String = ""
    diffuse_filter::Bool = false
    display::Bool = true
    fast_kalman_filter::Bool = true
    first_obs::PeriodsSinceEpoch = Undated(typemin(Int))
    last_obs::PeriodsSinceEpoch = Undated(typemin(Int))
    mcmc_chains::Int = 1
    mcmc_init_scale::Float64 = 0
    mcmc_jscale::Float64 = 0.2
    mcmc_replic::Int =  0
    mode_check::Bool = false
    mode_compute::Bool = true
    mode_file::String = ""
    nobs::Int = 0
    order::Int = 1
    plot_priors::Bool = true
    presample::Int = 0
end

function translate_estimation_options(options)
    new_options = copy(options)
    for (k, v) in options
        if k == "first_obs"
            new_options["first_obs"] = v[1]
        elseif k == "last_obs"
            new_options["last_obs"] = v[1]
        elseif k == "mh_jscale"
            new_options["mcmc_jscale"] = v
            delete!(new_options, "mh_jscale")
        elseif k == "mh_nblck"
            new_options["mcmc_chains"] = v
            delete!(new_options, "mh_nblck")
        elseif k == "mh_replic"
            new_options["mcmc_replic"] = v
            delete!(new_options, "mh_replic")
        elseif k == "plot_priors"
            v == 0 && (options["plot_priors"] = false)
        elseif k == "nobs"
            new_options["nobs"] = v[1]
        end
    end
    return new_options 
end

function estimation!(context, field::Dict{String, Any})
    opt = translate_estimation_options(field["options"])
    ff = NamedTuple{Tuple(Symbol.(keys(opt)))}(values(opt))
    options = EstimationOptions(; ff...)
    symboltable = context.symboltable
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    estimation_results = results.estimation
    varobs = context.work.observed_variables
    observations  = get_transposed_data!(context, options.datafile, options.data, varobs, options.first_obs, options.last_obs, options.nobs)
    nobs = size(observations, 2)
    estimated_parameters = context.work.estimated_parameters
    initial_parameter_values = get_initial_value_or_mean()
    set_estimated_parameters!(context, initial_parameter_values)
    if options.plot_priors
        plot_priors()
    end 
    if options.mode_compute
        (res, mode, tstdh, mode_covariance) = posterior_mode!(context,  initial_parameter_values, observations, varobs, transformed_parameters = true)
        options.display && log_result(res)
        estimation_result_table(
            get_parameter_names(estimated_parameters),
            mode,
            tstdh,
            "Posterior mode"
        )
    end

    if !isempty(options.mode_file)
        mode, mode_covariance = get_mode_file(options.mode_file)
    end

    if options.mcmc_replic > 0
        chain = mh_estimation(context, observations, varobs, mode, 
                              covariance = options.mcmc_jscale*mode_covariance,
                              mcmc_replic = options.mcmc_replic,
                              mcmc_chains = options.mcmc_chains)
        output_MCMCChains(context, chain, options.display, options.display)
        options.display && plot_prior_posterior(chain)
    end 
       
    return nothing
end

"""
     mode_compute!(; 
                 context=context,
                 data = AxisArrayTable(AxisArrayTables.AxisArray(Matrix(undef, 0, 0))),
                 datafile = "",
                 diffuse_filter::Bool = false,
                 display::Bool = false,
                 fast_kalman_filter::Bool = true,
                 first_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
                 initial_values = get_initial_value_or_mean(),
                 last_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
                 mode_check::Bool = false,
                 nobs::Int = 0,
                 order::Int = 1,
                 presample::Int = 0,
                 transformed_parameters = true)

computes the posterior mode.

# Keyword arguments
- `context::Context=context`: context of the computation
- `data::AxisArrayTable`: AxisArrayTable containing observed variables
- `datafile::String:  data filename (can't be used as the same time as `dataset`)
- `first_obs::PeriodsSinceEpoch`: first observation (default: first observation in the dataset)
- `initial_values`: initival parameter values for optimization algorithm (default: `estimated_params_init` block if present or prior mean)
- `last_obs::PeriodsSinceEpoch`: last period (default: last period of the dataset)
- `nobs::Int = 0`: number of observations (default: entire dataset)
- `transformed_parameters = true`: whether to transform estimated parameter so as they take their value on R

Either `data` or `datafile` must be specified.
"""
function mode_compute!(; algorithm = BFGS,
                       context=context,
                       data = AxisArrayTable(AxisArrayTables.AxisArray(Matrix(undef, 0, 0))),
                       datafile = "",
                       diffuse_filter::Bool = false,
                       display::Bool = false,
                       fast_kalman_filter::Bool = true,
                       first_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
                       initial_values = get_initial_value_or_mean(),
                       last_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
                       mode_check::Bool = false,
                       nobs::Int = 0,
                       order::Int = 1,
                       presample::Int = 0,
                       show_trace = false,
                       transformed_parameters = true
)
    symboltable = context.symboltable
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    estimation_results = results.estimation

    varobs = context.work.observed_variables
    observations = get_transposed_data!(context, datafile, data, varobs, first_obs, last_obs, nobs)
    (res, mode, tstdh, mode_covariance) = posterior_mode!(context, initial_values, observations, varobs, algorithm = algorithm, show_trace = show_trace, transformed_parameters = transformed_parameters)
    if display
        log_result(res)
        estimation_result_table(
            get_parameter_names(context.work.estimated_parameters),
            mode,
            tstdh,
            "Posterior mode"
        )
    end
end

"""
    rwmh_compute!(;context::Context=context,
             back_transformation::Bool = true,
             datafile::String = "",
             data::AxisArrayTable = AxisArrayTable(AxisArrayTables.AxisArray(Matrix(undef, 0, 0))),
             diffuse_filter::Bool = false,
             display::Bool = true,
             fast_kalman_filter::Bool = true,
             first_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
             initial_values = prior_mean(context.work.estimated_parameters),
             covariance::AbstractMatrix{Float64} = Matrix(prior_variance(context.work.estimated_parameters)),
             transformed_covariance::Matrix{Float64} = Matrix{Float64}(undef, 0,0),
             last_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
             mcmc_chains::Int = 1,
             mcmc_init_scale::Float64 = 0.0,
             mcmc_jscale::Float64 = 0.0,
             mcmc_replic::Int =  0,
             mode_compute::Bool = true,
             nobs::Int = 0,
             order::Int = 1,
             plot_chain::Bool = true,
             plot_posterior_density::Bool = false, 
             presample::Int = 0,
             transformed_parameters::Bool = true
)

runs random walk Monte Carlo simulations of the posterior

# Keyword arguments
- `context::Context=context`: context of the computation
- `covariance::AbstractMatrix{Float64}`: 
- `data::AxisArrayTable`: AxisArrayTable containing observed variables
- `datafile::String`:  data filename
- `first_obs::PeriodsSinceEpoch`: first observation (default: first observation in the dataset)
- `initial_values`: initival parameter values for optimization algorithm (default: `estimated_params_init` block if present or prior mean)
- `last_obs::PeriodsSinceEpoch`: last period (default: last observation in the dataset)
- `mcmc_chains::Int` number of MCMC chains (default: 1)
- `mcmc_jscale::Float64`: scale factor of proposal
- `mcmc_replic::Int`: =  0,
- `nobs::Int = 0`: number of observations (default: entire dataset)
- `plot_chain::Bool`: whether to display standard MCMC chain output (default:true)
- `plot_posterior_density::Bool`: wether to display plots with prior and posterior densities (default: false)
- `transformed_covariance::Matrix{Float64}`: covariance of transformed parameters (default: empty)
- `transformed_parameters = true`: whether to transform estimated parameter so as they take their value on R

Either `data` or `datafile` must be specified.
"""
function rwmh_compute!(;context::Context=context,
             back_transformation::Bool = true,
             datafile::String = "",
             data::AxisArrayTable = AxisArrayTable(AxisArrayTables.AxisArray(Matrix(undef, 0, 0))),
             diffuse_filter::Bool = false,
             display::Bool = true,
             fast_kalman_filter::Bool = true,
             first_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
             initial_values = prior_mean(context.work.estimated_parameters),
             covariance::AbstractMatrix{Float64} = Matrix(prior_variance(context.work.estimated_parameters)),
             transformed_covariance::Matrix{Float64} = Matrix{Float64}(undef, 0,0),
             last_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
             mcmc_chains::Int = 1,
             mcmc_init_scale::Float64 = 0.0,
             mcmc_jscale::Float64 = 0.0,
             mcmc_replic::Int =  0,
             mode_compute::Bool = true,
             nobs::Int = 0,
             order::Int = 1,
             plot_chain::Bool = true,
             plot_posterior_density::Bool = false, 
             presample::Int = 0,
             transformed_parameters::Bool = true
)
    symboltable = context.symboltable
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    estimation_results = results.estimation
    varobs = context.work.observed_variables
    observations = get_transposed_data!(context, datafile, data, varobs, first_obs, last_obs, nobs)
    chain = mh_estimation(context, 
        observations,
        varobs, 
        initial_values, 
        covariance = mcmc_jscale*covariance, 
        transformed_covariance = mcmc_jscale*transformed_covariance,
        back_transformation=back_transformation, 
        first_obs=first_obs, 
        last_obs=last_obs, 
        mcmc_chains=mcmc_chains, 
        mcmc_replic=mcmc_replic, 
        transformed_parameters=transformed_parameters
    )
    output_MCMCChains(context, chain, display, plot_chain)
    plot_posterior_density && plot_prior_posterior(chain)
    return chain
end

function output_MCMCChains(context, chain, display, plot_chain)
    estimation_results = context.results.model_results[1].estimation
    n = estimation_results.posterior_mcmc_chains_nbr += 1
    path = mkpath(joinpath(context.modfileinfo.modfilepath, "output"))
    serialize("$path/mcmc_chain_$n.jls", chain)
    display && log_result(chain)
    plot_chain && plot_MCMCChains(chain, n, "$path/graphs", display)
end

function plot_MCMCChains(chain, n, path, display)
    mkpath(path)
    filename = "$(path)/MCMC_chains_$n"
    pl = StatsPlots.plot(chain)
    display && graph_display(pl)
    savefig(filename)
end 

function smc_compute!(;context=context,
                      datafile = "",
                      first_obs = 1,
                      last_obs = 0,
                      nobs = 0
                      )
    symboltable = context.symboltable
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    estimation_results = results.estimation
    trends = results.trends
    has_trends = context.modfileinfo.has_trends

    ep = context.work.estimated_parameters
    np = length(ep.prior)
    observations = get_observations(context, datafile, data, first_obs, last_obs, nobs)
    ssws = SSWs(context, size(observations, 1), context.work.observed_variables)    
    
    pdraw!(θ) = prior_draw!(θ, ep)
    logprior(θ) = logpriordensity(θ, ep)
    function llikelihood(θ)
        try 
            return loglikelihood(θ, context, observations, ssws)
        catch
            return -Inf
        end 
    end 

    smc(pdraw!, logprior, llikelihood, ep.name, length(ep.prior))
end 

function prior_draw!(pdraw, ep)
    for (i, p) in enumerate(ep.prior)
        pdraw[i] = rand(p, 1)[1]
    end
end

prior_mean(ep::EstimatedParameters) = [mean(p) for p in ep.prior]

function prior_variance(ep::EstimatedParameters)
    d = Vector{Float64}(undef, length(ep.prior))
    for (i, p) in enumerate(ep.prior)
        v = var(p)
        d[i] = v > 1 ? 1 : v
    end 
    return Diagonal(d)
end 

"""
 covariance(chain:Chains)

computes the covariance matrix of MCMC `chain`
"""
function covariance(chain::Chains)
    c = copy(chain.value.data[:,1:end-1,1])
    m = mean(c, dims=1)
    c .-= m
    return Symmetric(c'*c/size(c,1))
end

struct SSWs{D<:AbstractFloat,I<:Integer}
    a0::Vector{D}
    dynamicws::Dynare.DynamicWs
    cssws::ComputeStochSimulWs
    stochsimulws::Dynare.StochSimulWs
    H::Matrix{D}
    lagged_state_ids::Vector{I}
    P::Matrix{D}
    Q::Matrix{D}
    R::Matrix{D}
    state_ids::Vector{I}
    stoch_simul_options::Dynare.StochSimulOptions
    T::Matrix{D}
    obs_idx::Vector{I}
    obs_idx_state::Vector{I}
    Y::Matrix{Union{D,Missing}}
    Z::Matrix{D}
    kalman_ws::KalmanLikelihoodWs{D, I}
    function SSWs(context, observations, varobs)
        model = context.models[1]
        D = eltype(context.work.params)
        symboltable = context.symboltable
        dynamicws = Dynare.DynamicWs(context)
        cssws = Dynare.ComputeStochSimulWs(context)
        stochsimulws = Dynare.StochSimulWs(context, 1)
        stoch_simul_options = Dynare.StochSimulOptions(Dict{String,Any}())
        obs_idx = [
            symboltable[String(v)].orderintype for
            v in varobs]
        state_ids = sort!(union(obs_idx, model.i_bkwrd_b))
        obs_idx_state = [Base.findfirst(isequal(i), 
        state_ids) for i in obs_idx]
        lagged_state_ids = findall(in(model.i_bkwrd_b), state_ids)
        np = model.exogenous_nbr
        ns = length(state_ids)
        ny = length(varobs)
        nobs = size(observations, 2)
        Y = similar(observations)
        Z = zeros(D, (ny, ns))
        H = zeros(D, (ny, ny))
        T = zeros(D, (ns, ns))
        R = zeros(D, (ns, np))
        Q = zeros(D, (np, np))
        a0 = zeros(D, ns)
        P = zeros(D, (ns, ns))
        T = zeros(D, (ns, ns))
        new{D,Int64}(
            a0,
            dynamicws,
            cssws,
            stochsimulws,
            H,
            lagged_state_ids,
            P,
            Q,
            R,
            state_ids,
            stoch_simul_options,
            T,
            obs_idx,
            obs_idx_state,
            Y,
            Z,
            KalmanLikelihoodWs(
                ny,
                ns,
                np,
                nobs
            )
        )
    end
end

## Estimation problems
struct DSGELogLikelihood{F<:Function, UT}
    dimension::Int
    context::Context
    observations::Matrix{UT}
    ssws::SSWs
end

struct DSGELogPosteriorDensity{UT}
    dimension::Integer
    context::Context
    observations::Matrix{Union{Missing, UT}}
    ssws::SSWs
end

function DSGELogPosteriorDensity(context, observations, obsvarnames)
    n = length(context.work.estimated_parameters)
    ssws = SSWs(context, observations, obsvarnames)
    DSGELogPosteriorDensity{eltype(observations)}(n, context, observations, ssws)
end

function logpriordensity(x, estimated_parameters)::Float64
    lpd = 0.0
    k = 1
    for (k, p) in enumerate(estimated_parameters.prior)
        lpd += logpdf(p, x[k])
    end
    return lpd
end

## Parameters transformation on ℝ
function DSGETransformation(ep::EstimatedParameters)
    tvec = []
    for p in ep.prior
        if typeof(p) <: Distributions.Beta
            push!(tvec, as𝕀)
        elseif typeof(p) <: Distributions.Gamma
            push!(tvec, asℝ₊)
        elseif typeof(p) <: Distributions.InverseGamma
            push!(tvec, asℝ₊)
        elseif typeof(p) <: InverseGamma1
            push!(tvec, asℝ₊)
        elseif typeof(p) <: Distributions.Normal
            push!(tvec, asℝ)
        elseif typeof(p) <: Distributions.Uniform
            push!(tvec, as(Real, p.a, p.b))
        elseif typeof(p) <: Distributions.Weibull
            push!(tvec, asℝ₊)
        else
            error("Unknown prior distribution")
        end 
    end
    return as((tvec...,))
end

function (problem::DSGELogLikelihood)(θ)
    @debug θ
    lpd = -Inf
    context = problem.context
    try
      lpd = logpriordensity(θ, context.work.estimated_parameters)
      lpd += loglikelihood(θ, context, problem.observations, problem.ssws)
    catch e
        @debug e
        lpd = -Inf
    end
    return lpd
end    

function (problem::DSGELogPosteriorDensity)(θ)
    @debug θ
    lpd = -Inf
    context = problem.context
    try
      lpd = logpriordensity(θ, context.work.estimated_parameters)
        lpd += loglikelihood(θ, context, problem.observations, problem.ssws)
    catch e
        @debug e
        lpd = -Inf
    end
    return lpd
end    

function set_estimated_parameters!(context::Context, value::AbstractVector{T}) where {T<:Union{Missing,Real}}
    ep = context.work.estimated_parameters
    for (k, p) in enumerate(zip(ep.index, ep.parametertype))
        set_estimated_parameters!(context, p[1], value[k], Val(p[2]))
    end
end

function set_estimated_parameters!(
    context::Context,
    index::N,
    value::T,
    ::Val{EstParameter},
) where {T<:Real, N<:Integer}
    context.work.params[index] = value
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{EstSDMeasurement},
) where {T<:Real, N<:Integer}
    context.work.Sigma_m[index[1], index[2]] = value*value
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{EstSDShock},
) where {T<:Real, N<:Integer}
    context.models[1].Sigma_e[index[1], index[2]] = value*value
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{EstVarMeasurement},
) where {T<:Real, N<:Integer}
    context.work.Sigma_m[index[1], index[2]] = value
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{EstVarShock},
) where {T<:Real, N<:Integer}
    context.models[1].Sigma_e[index[1], index[2]] = value
    return nothing
end

function set_covariance!(corr, covariance, index1, index2)
    V1 = covariance[index1, index1]
    V2 = covariance[index2, index2]
    covariance[index1, index2] = corr*sqrt(V1*V2)
    covariance[index2, index1] = corr*sqrt(V1*V2)
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{EstCorrMeasurement},
) where {T<:Real, N<:Integer}
    set_covariance!(value, context.work.Sigma_m,
    index[1], index[2])
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{EstCorrShock},
) where {T<:Real, N<:Integer}
    set_covariance!(value, context.models[1].Sigma_e,
    index[1], index[2])
    return nothing
end

function loglikelihood(
    θ,
    context::Context,
    observations::AbstractMatrix{D},
    ssws::SSWs{T,I},
)::T where {T<:Real,D<:Union{Missing,<:Real},I<:Integer}
    @debug θ
    parameters = collect(θ)
    model = context.models[1]
    results = context.results.model_results[1]
    work = context.work
    estimated_parameters = work.estimated_parameters
    model_parameters = work.params
    set_estimated_parameters!(context, parameters)
    fill!(context.results.model_results[1].trends.exogenous_steady_state, 0.0)
    #compute steady state and first order solution
    compute_steady_state!(context)
    compute_stoch_simul!(
        context,
        ssws.dynamicws,
        ssws.cssws,
        ssws.stochsimulws,
        model_parameters,
        ssws.stoch_simul_options;
        variance_decomposition = false,
    )
    LRE = LinearRationalExpectations
    LRE_results = results.linearrationalexpectations
    compute_variance!(LRE_results, model.Sigma_e, workspace(LRE.VarianceWs, context))
    # build state space representation
    ssws.Y .= observations
    detrend_data!(
        ssws.Y,
        context,
        context.work.observed_variables,
        dim = 1
    )
    ns = length(ssws.state_ids)
    np = model.exogenous_nbr
    ny, nobs = size(ssws.Y)
    obs_idx_state = ssws.obs_idx_state
    for i = 1:ny
        ssws.Z[i, obs_idx_state[i]] = T(1)
    end
    ssws.H .= work.Sigma_m
    vg1 = view(LRE_results.g1_1, ssws.state_ids, :)
    view(ssws.T, :, ssws.lagged_state_ids) .= vg1
    vg2 = view(LRE_results.g1_2, ssws.state_ids, :)
    ssws.R .= vg2
    ssws.Q .= model.Sigma_e
    fill!(ssws.a0, T(0))
    ssws.P .= view(
        LRE_results.endogenous_variance,
        ssws.state_ids,
        ssws.state_ids,
    )
    start = 1
    last = nobs
    presample = 0
    # workspace for kalman_likelihood
    klws = ssws.kalman_ws
    Y = ssws.Y
    if any(ismissing.(Y))
        # indices of non-missing observations are in data_pattern
        data_pattern = Vector{Vector{Int}}(undef, 0)
        for i = 1:nobs
            push!(data_pattern, findall(.!ismissing.(view(Y[:, i]))))
        end
        return Dynare.kalman_likelihood(
            Y,
            ssws.Z,
            ssws.H,
            ssws.T,
            ssws.R,
            ssws.Q,
            ssws.a0,
            ssws.P,
            start,
            last,
            presample,
            klws,
            data_pattern,
        )
    else
        return Dynare.kalman_likelihood(
            Y,
            ssws.Z,
            ssws.H,
            ssws.T,
            ssws.R,
            ssws.Q,
            ssws.a0,
            ssws.P,
            start,
            last,
            presample,
            klws,
        )
    end
end

## Smooth penalty functions
"""
function penalty(e::UndeterminateSystemException, eigenvalues, backward_nbr)

returns the sum of the modulus of the forward_nbr largest eigenvalues smaller or equal to 1.0 in modulus 
"""      
function smooth_penalty(e::UndeterminateSystemException,  
                 eigenvalues::AbstractVector{<:Union{T,Complex{T}}},
                 backward_nbr::Integer,
                 ) where {T<:Real}
    n = length(eigenvalues)
    sort!(eigenvalues, by=abs)
    first_unstable_root = findfirst(x -> x > 1.0, eigenvalues)
    return sum(eigenvalues, backward_nbr + 1, first_unstable_root)
end

"""
function penalty(e::UnstableSystemException, eigenvalues, backward_nbr)

returns the sum of the modulus of the backward_nbr smallest eigenvalues larger than 1.0 in modulus 
"""      
function smooth_penalty(e::UnstableSystemException,  
                 eigenvalues::AbstractVector{<:Union{T,Complex{T}}},
                 backward_nbr::Integer,
                 ) where {T<:Real}
    n = length(eigenvalues)
    sort!(eigenvalues, by=abs)
    first_unstable_root = findfirst(x -> x > 1.0, eigenvalues)
    return sum(view(eigenvalues, first_unstable_root, backward_nbr))
end


"""
function penalty(e::DomainError, eigenvalues, backward_nbr)

returns the absolute value of the argument triggering a DomainError exception
"""      
function smooth_penalty(e::DomainError, 
                 eigenvalues::AbstractVector{<:Union{T,Complex{T}}},
                 backard_nbr::Integer,
                 ) where {T<:Real}
    msg = sprint(showerror, e, catch_backtrace())
    x = parse(Float64, rsplit(rsplit(msg, ":")[1], " ")[3])
    return abs(x)
end

"""
function penalty(e::LinearAlgebra.LAPACKException, eigenvalues, backward_nbr)

catches unexpected exceptions and rethrow
"""      
function smooth_penalty(e::LinearAlgebra.LAPACKException, 
                 eigenvalues::AbstractVector{<:Union{T,Complex{T}}},
                 backard_nbr::Integer,
                 ) where {T<:Real}
    return -Inf
end

"""
function penalty(e::Exception, eigenvalues, backward_nbr)

catches unexpected exceptions and rethrow
"""      
function smooth_penalty(e::Exception, 
                 eigenvalues::AbstractVector{<:Union{T,Complex{T}}},
                 backard_nbr::Integer,
                 ) where {T<:Real}
    return -Inf
end

function get_symbol(symboltable, indx)
    for k in symboltable
        if k[2].symboltype == SymbolType(3) && k[2].orderintype == indx
            return k[2].longname
        end
    end
end

function dynare_optimize(objective, initialvalue, algorithm; algoptions=(), optimoptions=() )
    let res, minimizer
        if algorithm in [ABC]
            n = length(initialvalue)
            bounds = BoxConstrainedSpace(lb = -10*ones(n), ub=10*ones(n))
            res = Metaheuristics.optimize(objective, bounds, algorithm(;algoptions...))
            minimizer = Metaheuristics.minimizer(res)
        elseif algorithm in [BFGS, LBFGS, NelderMead, Newton, NewtonTrustRegion]
            res = Optim.optimize(objective, initialvalue, algorithm(;algoptions...), Optim.Options(;optimoptions...) )
            minimizer = res.minimizer
        else
            error("Unknown algorithm $algorithm")
        end
        return res, minimizer
    end 
end 

## Estimation
function ml_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    kwargs...,
)
    ep = context.work.estimated_parameters
    problem = DSGELogLikelihood(context, datafile, first_obs, last_obs)
    problem(θ) = -Loglikelihood(θ)
    transformation = DSGETransformation(ep)
    ℓ = TransformedLogDensity(transformation, problem)
    objective(θ) = -LogDensityProblems.logdensity(ℓ, θ)
    p0 = ep.initialvalue
    ip0 = collect(TransformVariables.inverse(transformation, tuple(p0...)))

    (res, minimizer) = dynare_optimize(
        tansformed_problem,
        p0,
        algorithm,
        optimoptions = (f_tol = 1e-5, show_trace = true)
    )

    hess = finite_difference_hessian(problem.f, minimizer)
    inv_hess = inv(hess)
    hsd = sqrt.(diag(hess))
    invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
    stdh = sqrt.(diag(invhess))
#=
    maximum_likelihood_result_table(
        context.symboltable,
        params_indices,
        res.minimizer,
        stdh,
    )
    results = maximum_likelihood(
        params,
        params_indices,
        shock_variance_indices,
        measurement_variance_indices,
        observed_variables,
        observations,
        context,
        ssws,
    )
=#
end

function posterior_mode!(
    context,
    initialvalue,
    observations,
    obsvarnames;
    algorithm = BFGS,
    iterations = 1000,
    show_trace = false,
    transformed_parameters = true
)
    ep = context.work.estimated_parameters
    results = context.results.model_results[1].estimation
    problem = DSGELogPosteriorDensity(context, observations, obsvarnames)
    if transformed_parameters
        transformation = DSGETransformation(ep)
        ℓ = TransformedLogDensity(transformation, problem)
        objective(θ) = -LogDensityProblems.logdensity(ℓ, θ)
        p0 = collect(TransformVariables.inverse(transformation, tuple(initialvalue...)))
        @debug objective(p0)
        (res, minimizer) = dynare_optimize(
            objective,
            p0,
            algorithm,
            optimoptions = (show_trace = show_trace, f_tol = 1e-5, iterations = iterations),
        )
        hess = finite_difference_hessian(objective, minimizer)
        hsd = sqrt.(diag(hess))
        invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
        stdh = sqrt.(diag(invhess))
        results.transformed_posterior_mode = copy(minimizer)
        results.transformed_posterior_mode_std = copy(stdh)
        results.transformed_posterior_mode_covariance = copy(invhess)
        results.posterior_mode = copy(collect(TransformVariables.transform(transformation, minimizer)))
        objective2(θ) = -problem(θ)
        hess = finite_difference_hessian(objective2, results.posterior_mode)
        hsd = sqrt.(diag(hess))
        invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
        d = diag(invhess)
        if any(d .< 0)
            k = findall(d .< 0)
            d[k] .= NaN
        end
        results.posterior_mode_std = copy(sqrt.(d))
        results.posterior_mode_covariance = copy(invhess)
    else
        objective1(θ) = -problem(θ)
        p0 = initialvalue
        @debug objective1(p0)
        (res, minimizer) = dynare_optimize(
            objective1,
            p0,
            algorithm,
            optimoptions = (show_trace = show_trace, f_tol = 1e-5, iterations = iterations),
        )
        hess = finite_difference_hessian(objective1, minimizer)
        if !isposdef((hess + hess')/2)
            @show minimizer
            error("Hessian isn't positive definite")
        end 
        hsd = sqrt.(diag(hess))
        invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
        results.posterior_mode = copy(minimizer)
        results.posterior_mode_std = copy(sqrt.(diag(invhess)))
        results.posterior_mode_covariance = copy(invhess)
    end
    
    return(res,
           results.posterior_mode,
           results.posterior_mode_std,
           results.posterior_mode_covariance
           )
end

function mh_estimation(
    context,
    observations,
    obsvarnames,
    initial_values;
    covariance,
    back_transformation = true,
    first_obs = 1,
    last_obs = 0,
    mcmc_chains = 1,
    mcmc_replic = 100000,
    transformed_covariance = Matrix(undef, 0, 0),
    transformed_parameters = true,
    kwargs...,
)
    ep = context.work.estimated_parameters
    param_names = ep.name
    problem = DSGELogPosteriorDensity(context, observations, obsvarnames)

    if transformed_parameters
        transformation = DSGETransformation(context.work.estimated_parameters)
 #       transformed_density(θ) = problem(collect(TransformVariables.transform(transformation, θ)))
        transformed_problem = TransformedLogDensity(transformation, problem)
        if back_transformation
            initial_values = Vector(collect(TransformVariables.inverse(transformation, tuple(initial_values...))))
        end
        let proposal_covariance
            if !isempty(transformed_covariance)
                proposal_covariance = transformed_covariance
            else
                proposal_covariance = inverse_transform_variance(transformation, initial_values, Matrix(covariance))    
            end
            chain = run_mcmc(transformed_problem, initial_values, proposal_covariance, param_names, mcmc_replic, mcmc_chains)
        end 
        if back_transformation
            transform_chains!(chain, transformation, problem)
        end 
    else
        chain = run_mcmc(problem, initial_values, Matrix(covariance), param_names, mcmc_replic, mcmc_chains)
    end 
    return chain
end 

function run_mcmc(posterior_density, initial_values, proposal_covariance, param_names, mcmc_replic, mcmc_chains)    
    model = DensityModel(posterior_density)
    proposal_covariance = (proposal_covariance + proposal_covariance')/2
    spl = RWMH(MvNormal(zeros(length(initial_values)), proposal_covariance))
    if mcmc_chains == 1
        AbstractMCMC.setprogress!(false, silent=true)
        chain = MCMCChains.sample(posterior_density, spl, mcmc_replic,
                   initial_params = Float64.(initial_values),
                   param_names = get_parameter_names(context.work.estimated_parameters),
                   chain_type = Chains)
    else    
        old_blas_threads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        if ndims(initial_values) == 1 || size(initial_values, 2) == 1
            initial_values_ = FillArrays.Fill(initial_values, mcmc_chains)
        end
        chain = MCMCChains.sample(posterior_density, spl, MCMCThreads(), mcmc_replic,
                   mcmc_chains,
                   initial_params = initial_values_,
                   param_names = get_parameter_names(context.work.estimated_parameters),
                   chain_type = Chains)
        BLAS.set_num_threads(old_blas_threads)
    end 
    display_acceptance_rate(chain)
    return chain
end 

## Results
function estimation_result_table(param_names, estimated_value, stdh, title)
    table = Matrix{Any}(undef, length(param_names) + 1, 4)
    table[1, 1] = "Parameter"
    table[1, 2] = "Estimated value"
    table[1, 3] = "Standard error"
    table[1, 4] = "80% confidence interval"
    for (i, n) in enumerate(param_names)
        table[i+1, 1] = n
        table[i+1, 2] = estimated_value[i]
        table[i+1, 3] = stdh[i]
        table[i+1, 4] = estimated_value[i] ± (1.28 * stdh[i])
    end
    dynare_table(table, title)
end

function maximum_likelihood_result_table(symboltable, param_indices, estimated_params, stdh)
    table = Matrix{Any}(undef, length(param_indices) + 1, 4)
    table[1, 1] = "Parameter"
    table[1, 2] = "Estimated value"
    table[1, 3] = "Standard error"
    table[1, 4] = "80% confidence interval"
    for (i, k) in enumerate(param_indices)
        table[i+1, 1] = get_symbol(symboltable, k)
        table[i+1, 2] = estimated_params[i]
        table[i+1, 3] = stdh[i]
        table[i+1, 4] = estimated_params[i] ± (1.28 * stdh[i])
    end
    dynare_table(table, "Results from maximum likelihood estimation")
end

function display_acceptance_rate(chains)
    ar = acceptance_rate(chains)
    println("")
    for (i, s) in enumerate(ar) 
        println("Acceptance rate chain $i: $(ar[i])")
    end
    println("")
end 


function acceptance_rate(chains::Chains)
    data = chains.value.data
    n1, n2, n3 = size(data)
    rr = zeros(n3)
    k = drop(axes(data, 1), 1)
    @inbounds for j in 1:n3      
        for i in 2:n1
            # 2 identical draws indicate a rejected proposal
            if data[i, n2, j] == data[i-1, n2, j]
                rr[j] += 1
            end
        end
    end
    return 1 .- rr./n1
end    

# Utilities
function get_eigenvalues(context)
    eigenvalues = context.results.model_results[1].linearrationalexpectations.eigenvalues
    return eigenvalues
end

function get_initial_value_or_mean(; context = context, 
                                   initialvalues = context.work.estimated_parameters.initialvalue)
    epprior = context.work.estimated_parameters.prior
    if isempty(initialvalues)
        return [mean(p) for p in epprior]
    else
        @assert length(initialvalues) == length(epprior)
        return [ismissing(initialvalue) ? mean(prior) : initialvalue for (initialvalue, prior) in zip(initialvalues, epprior)]
    end
end

function get_mode_file(modfilepath)
    modfilename = strip(basename(modfilepath), ".mod")
    context = deserialize(joinpath(modfilepath, "output", modfilename))["context"]
    results = context.results.model_results[1]
    mode = results.posterior_mode
    mode_covariance = results.posterior_mode_covariance
    return (mode, mode_covariance)
end

function get_parameter_names(estimated_parameters::EstimatedParameters)
    parametertype = estimated_parameters.parametertype
    name = estimated_parameters.name
    return [get_parameter_name(Val(t), n) for (t, n) in zip(parametertype, name)]
end

get_parameter_name(::Val{EstParameter}, name) = name
get_parameter_name(::Val{EstSDShock}, name) = "Std($(name))"
get_parameter_name(::Val{EstSDMeasurement}, name) = "Std($(name))"
get_parameter_name(::Val{EstVarShock}, name) = "Var($(name))"
get_parameter_name(::Val{EstVarMeasurement}, name) = "Var($(name))"
get_parameter_name(::Val{EstCorrShock}, name) = "Corr($(name[1]), $(name[2]))"
get_parameter_name(::Val{EstCorrMeasurement}, name) = "Corr($(name[1]), $(name[2]))"



function get_observations(context, datafile, data, first_obs, last_obs, nobs)
    results = context.results.model_results[1]
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    if datafile != "" || length(data) > 0
        Yorig = get_data!(context, datafile, data, varobs, first_obs, last_obs, nobs)
    else
        error("The procedure needs a data file or a AxisArrayTable!")
    end
    return transpose(Yorig), column_labels(Yorig)
end

# Pretty printing name of estimated parameters

# transformation of standard errors using delta approach
function transform_std(T::TransformVariables.TransformTuple, x, std)
    ts = []
    for (i, t) in enumerate(T.transformations)
        push!(ts, transform_std(t, x[i], std[i]))
    end
    return ts
end

function transform_variance(T::TransformVariables.TransformTuple, m, var)
    var1 = copy(var)
    for (i, t) in enumerate(T.transformations)
        transform_std(t, m, view(var1, i, :))
    end
    for (i, t) in enumerate(T.transformations)
        transform_std(t, m, view(var1, i, :))
    end
    return var1
end

transform_std(::TransformVariables.Identity, x, std) = std
transform_std(::TransformVariables.ShiftedExp{true, Int64}, x, std) = exp.(x).*std
transform_std(::TransformVariables.ScaledShiftedLogistic{<:Real}, x, std) = logistic.(x).*(1-logistic.(x)).*std

# inverse transformation of standard errors using delta approach
function inverse_transform_std(T::TransformVariables.TransformTuple, x, std)
    ts = []
    for (i, t) in enumerate(T.transformations)
        push!(ts, inverse_transform_std(t, x[i], std[i]))
    end
    return ts
end

function inverse_transform_variance(T::TransformVariables.TransformTuple, m, var)
    var1 = copy(var)
    for (i, t) in enumerate(T.transformations)
        inverse_transform_std(t, m, view(var1, i, :))
    end
    for (i, t) in enumerate(T.transformations)
        inverse_transform_std(t, m, view(var1, i, :))
    end
    return var1
end

inverse_transform_std(::TransformVariables.Identity, x, std) = std
inverse_transform_std(::TransformVariables.ShiftedExp{true, Int64}, x, std) = std./x
inverse_transform_std(::TransformVariables.ScaledShiftedLogistic{<:Real}, x, std) = std./((1 .- x).*x)

function find(names::Vector, s)
    for (i, n) in enumerate(names)
        if n == s
            return i
        end
    end
end

has_posterior_mode(context) = length(context.work.estimated_parameters.posterior_mode)>0

function mcmc_diagnostics(chains, context, names)
    indices = [find(context.work.estimated_parameters.name, e) for e in names]
    iterations, column_nbr = size(chains.value.data)
    prior_pdfs = []
    n_plots = length(indices)
    posterior_pdfs = []
    posterior_x_axes = []
    prior_x_axes = []
    transformation = Dynare.DSGETransformation(context.work.estimated_parameters)
    transformed_data = Matrix{Float64}(undef, iterations, column_nbr - 1)
    
    for it in 1:iterations
        @views transformed_data[it, :] .= TransformVariables.transform(transformation, chains.value.data[it, 1:end-1])
    end
    
    for i in indices
        U = kde(transformed_data[:, i])
        prior = context.work.estimated_parameters.prior[i]
        m, v = mean(prior), var(prior)
        prior_x_axis = LinRange(m-15*v, m+15*v, length(U.x))
        prior_pdf = [pdf(prior, e) for e in prior_x_axis]
        push!(posterior_pdfs, U.density*(U.x[2] - U.x[1]))
        push!(posterior_x_axes, U.x)
        push!(prior_pdfs, prior_pdf/sum(prior_pdf))
        push!(prior_x_axes, prior_x_axis)
    end

    posterior_pdfs = hcat(posterior_pdfs...)
    posterior_x_axes = hcat(posterior_x_axes...)
    prior_pdfs = hcat(prior_pdfs...)
    prior_x_axes = hcat(prior_x_axes...)
    names = hcat(names...)

    f = Plots.plot(posterior_x_axes, posterior_pdfs, layout=n_plots, linecolor=:darkblue, legend=false, linewidth=3)
    plot!(f, prior_x_axes, prior_pdfs, layout=n_plots, title=names, linecolor=:darkgrey, legend=false, linewidth=3)
    if has_posterior_mode(context)
        posterior_modes = [context.work.estimated_parameters.posterior_mode[i] for i in indices]
        Plots.plot!(f, hcat(posterior_modes...), layout=n_plots, seriestype=:vline, line=(:dot,3), linecolor=:green, legend=false, linewidth=3)
    end
    f
end

function transform_chains!(chains, t, posterior_density)
    y = chains.value.data
    nparams = size(y, 2) - 1
    tmp = Vector{Float64}(undef, size(y, 1))
    for i in axes(y, 3)
        for j in axes(y, 1)
            @views vj1 = y[j, 1:nparams, i]  
            @views vj2 = y[j, nparams + 1, i]  
            tparam, dp = TransformVariables.transform_and_logjac(t, vj1)
            @views vj1 .= tparam
            y[j, nparams + 1, i] -= dp
        end          
    end
end 

function log_result(result)
    io = IOBuffer()
    show(io, "text/plain", result)
    @info String(take!(io))
end
