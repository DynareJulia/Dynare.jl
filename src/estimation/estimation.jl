using Distributions
using DynamicHMC
using FiniteDiff: finite_difference_gradient, finite_difference_hessian
using KernelDensity
using LogDensityProblems
using MCMCChains
using MCMCDiagnosticTools
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
include("estimated_parameters.jl")

Base.@kwdef struct EstimationOptions
    config_sig::Float64 = 0.8
    data::AxisArrayTable = AxisArrayTable(AxisArrayTables.AxisArray([;;]))
    datafile::String = ""
    diffuse_filter::Bool = false
    display::Bool = false
    fast_kalman_filter::Bool = true
    first_obs::PeriodsSinceEpoch = Undated(1)
    last_obs::PeriodsSinceEpoch = Undated(0)
    mcmc_chains::Int = 1
    mcmc_init_scale::Float64 = 0
    mcmc_jscale::Float64 = 0
    mcmc_replic::Int =  0
    mode_check::Bool = false
    mode_compute::Bool = true
    nobs::Int = 0
    order::Int = 1
    plot_prior::Bool = false
    presample::Int = 0
end

function translate_estimation_options(options)
    new_options = copy(options)
    for (k, v) in options
        if k == "mh_jscale"
            new_options["mcmc_jscale"] = v
            delete!(new_options, "mh_jscale")
        elseif k == "mh_nblck"
            new_options["mcmc_chains"] = v
            delete!(new_options, "mh_nblck")
        elseif k == "mh_replic"
            new_options["mcmc_replic"] = v
            delete!(new_options, "mh_replic")
        end
    end
    return new_options 
end

function get_observables(datafile, varobs, first_obs, last_obs, symboltable, has_trends, steady_state, linear_trend)
    if (filename = datafile) != ""
        Yorig = get_data(filename, varobs, start = first_obs, last = last_obs)
    else
        error("estimation needs a data file or an AxisArrayTable!")
    end
    #=
    varobs_ids =
        [symboltable[v].orderintype for v in varobs if is_endogenous(v, symboltable)]
    Y = Matrix{Union{Float64,Missing}}(undef, size(Yorig
    if has_trends
        remove_linear_trend!(
            Y,
            Yorig,
            steady_state[varobs_ids],
            linear_trend[varobs_ids],
        )
    else
        Y .= Yorig .- steady_state[varobs_ids]
    end
    =#
    return Yorig
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
    trends = results.trends
    has_trends = context.modfileinfo.has_trends
        
    observations = get_observables(options.datafile, varobs, options.first_obs, options.last_obs, symboltable, has_trends, trends.endogenous_steady_state, trends.endogenous_linear_trend)
    nobs = size(observations, 2)
    ssws = SSWs(context, nobs, varobs)
    estimated_parameters = context.work.estimated_parameters
    initial_values = get_initial_value_or_mean(estimated_parameters)
    
    set_estimated_parameters!(context, initial_values)
    
    if options.mode_compute
        (res, mode, tstdh, mode_covariance) = posterior_mode(context,  initial_values, observations)
        @show res
    end

    if options.mcmc_replic > 0
        chain = mh_estimation(context, observations, mode, 
        #options.mcmc_jscale*Matrix(prior_variance(context.work.estimated_parameters),
        options.mcmc_jscale*mode_covariance,
        mcmc_replic=options.mcmc_replic)
        display(chain)
        StatsPlots.plot(chain)
    end 
       
    return nothing
end

function mode_compute(; context=context,
                 datafile = "",
                 data = AxisArrayTable(AxisArrayTables.AxisArray([;;])),
                 diffuse_filter::Bool = false,
                 display::Bool = false,
                 fast_kalman_filter::Bool = true,
                 first_obs::PeriodsSinceEpoch = Undated(1),
                 initial_values = prior_mean(context.work.estimated_parameters),
                 last_obs::PeriodsSinceEpoch = Undated(0),
                 mode_check::Bool = false,
                 nobs::Int = 0,
                 order::Int = 1,
                 presample::Int = 0
)
    symboltable = context.symboltable
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    estimation_results = results.estimation
    varobs = context.work.observed_variables
    trends = results.trends
    has_trends = context.modfileinfo.has_trends
    
    observations = get_observables(datafile, varobs, first_obs, last_obs, symboltable, has_trends, trends.endogenous_steady_state, trends.endogenous_linear_trend)
    (res, mode, tstdh, mode_covariance) = posterior_mode(context, initial_values, observations)
    @show res
end

function rwmh_compute(;context=context,
             datafile = "",
             back_transformation = true,
             data = AxisArrayTable(AxisArrayTables.AxisArray([;;])),
             diffuse_filter::Bool = false,
             display::Bool = true,
             fast_kalman_filter::Bool = true,
             first_obs::PeriodsSinceEpoch = Undated(1),
             initial_values = prior_mean(context.work.estimated_parameters),
             covariance = Matrix(prior_variance(context.work.estimated_parameters)),
             last_obs::PeriodsSinceEpoch = Undated(0),
             mcmc_chains::Int = 1,
             mcmc_init_scale::Float64 = 0.0,
             mcmc_jscale::Float64 = 0.0,
             mcmc_replic::Int =  0,
             mode_compute::Bool = true,
             nobs::Int = 0,
             order::Int = 1,
             plot_chain::Bool = true,
             plot_posterior_density = false, 
             presample::Int = 0,
             transformed_parameters = true
)
    symboltable = context.symboltable
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    estimation_results = results.estimation
    varobs = context.work.observed_variables
    trends = results.trends
    has_trends = context.modfileinfo.has_trends

    observations = get_observables(datafile, varobs, first_obs, last_obs, symboltable, has_trends, trends.endogenous_steady_state, trends.endogenous_linear_trend)
    (chain, back_transformed_chain) = mh_estimation(context, observations, initial_values, mcmc_jscale*covariance, back_transformation=back_transformation, first_obs=first_obs, last_obs=last_obs, mcmc_chains=mcmc_chains, mcmc_replic=mcmc_replic, transformed_parameters=transformed_parameters)
    output_mcmc_chain!(context, chain, display, plot_chain)
    plot_posterior_density && plot_prior_posterior(context, back_transformed_chain)
    plot_chain && StatsPlots.plot(chain)
    return chain
end

function output_mcmc_chain!(context, chain, display, plot_chain)
    estimation_results = context.results.model_results[1].estimation
    n = estimation_results.posterior_mcmc_chains_nbr += 1
    path = mkpath(joinpath(context.modfileinfo.modfilepath, "output"))
    serialize("$path/mcmc_chain_$n.jls",
    chain) 
    display && Base.display(chain)
    path 
    plot_chain && plot_MCMCChains(chain, n, "$path/graphs", display)
end

function plot_MCMCChains(chain, n, path, display)
    mkpath(path)
    filename = "$(path)/MCMC_chains_$n"
    pl = StatsPlots.plot(chain)
    display && graph_display(pl)
    savefig(filename)
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

function covariance(chain::Chains)
    c = copy(chain.value.data[:,1:end-1,1])
    m = mean(c)
    c .-= m
    return Symmetric(c'*c/length(c))
end

struct SSWs{D<:AbstractFloat,I<:Integer}
    a0::Vector{D}
    dynamicws::Dynare.DynamicWs
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
    function SSWs(context, nobs, varobs)
        model = context.models[1]
        D = eltype(context.work.params)
        symboltable = context.symboltable
        dynamicws = Dynare.DynamicWs(context)
        stoch_simul_options = Dynare.StochSimulOptions(Dict{String,Any}())
        obs_idx = [
            symboltable[v].orderintype for
            v in varobs]
        state_ids = sort!(union(obs_idx, model.i_bkwrd_b))
        obs_idx_state = [Base.findfirst(isequal(i), 
        state_ids) for i in obs_idx]
        lagged_state_ids = findall(in(model.i_bkwrd_b), state_ids)
        np = model.exogenous_nbr
        ns = length(state_ids)
        ny = length(varobs)
        Y = Matrix{Union{D,Missing}}(undef, ny, nobs)
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
struct DSGENegativeLogLikelihood{F<:Function, UT}
    f::F
    dimension::Integer
    context::Context
    observations::Matrix{UT}
    ssws::SSWs
    function DSGENegativeLogLikelihood(context, datafile, first_obs, last_obs)
        n = length(context.work.estimated_parameters)
        observations = get_observations(context, datafile, first_obs, last_obs)
        ssws = SSWs(context, size(observations, 2), context.work.observed_variables)
        f = make_negativeloglikelihood(context, observations, first_obs, last_obs, ssws)
        new{typeof(f),eltype(observations)}(f, n, context, observations, ssws)
    end
end

struct DSGELogPosteriorDensity{F <: Function, UT}
    f::F
    dimension::Integer
    context::Context
    observations::Matrix{Union{Missing,UT}}
    ssws::SSWs
    function DSGELogPosteriorDensity(context, observations, first_obs, last_obs)
        n = length(context.work.estimated_parameters)
        ssws = SSWs(context, size(observations, 2), context.work.observed_variables)
        f = make_logposteriordensity(context, observations, first_obs, last_obs, ssws)
        new{typeof(f), eltype(observations)}(f, n, context, observations, ssws)
    end
end

struct DSGENegativeLogPosteriorDensity{F <:Function, UT}
    f::F
    dimension::Integer
    context::Context
    observations::Matrix{Union{Missing,UT}}
    ssws::SSWs
    function DSGENegativeLogPosteriorDensity(context, datafile, first_obs, last_obs)
        n = length(context.work.estimated_parameters)
        observations = get_observations(context, datafile, first_obs, last_obs)
        ssws = SSWs(context, size(observations, 2), context.work.observed_variables)
        f = make_negativelogposteriordensity(context, observations, first_obs, last_obs, ssws)
        new{typeof(f), eltype(observations)}(f, n, context, observations, ssws)
    end
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

function make_negativeloglikelihood(context, observations, first_obs, last_obs, ssws)
    previous_value = 0.0
    function lognegativelikelihood(x)
        lll = Inf
        try
            lll = loglikelihood(x, context, observations, ssws)
            previous_value  = lll
        catch e
            eigenvalues = get_eigenvalues(context)
            model = context.models[1]
            backward_nbr = model.n_bkwrd + model.n_both
            #            lll = previous_value + penalty(e, eigenvalues, backward_nbr)
            @show e
            lll = Inf
        end
        return -lll
    end
    return lognegativelikelihood
end

function make_logposteriordensity(context, observations, first_obs, last_obs, ssws)
    function logposteriordensity(x)::Float64
        lpd = logpriordensity(x, context.work.estimated_parameters)
        if abs(lpd) == Inf
            return lpd
        end
        try
            lpd += loglikelihood(x, context, observations, ssws)
        catch e
            @debug e
            lpd = -Inf
        end
        return lpd
    end
    return logposteriordensity
end

## methods necessary for the interface, 
dimension(ld::DSGELogPosteriorDensity) = ld.dimension
logdensity(ld::DSGELogPosteriorDensity, x) = ld.f(x)
logdensity_and_gradient(ld::DSGELogPosteriorDensity, x) =
    ld.f(x), finite_difference_gradient(ld.f, x)
logdensity(tld::TransformedLogDensity, x) = tld.log_density_function(collect(TransformVariables.transform(tld.transformation, x)))

function logdensity_and_gradient(tld::TransformedLogDensity, x)
    tx = collect(TransformVariables.transform(tld.transformation, x))
    return (tld.log_density_function(tx), finite_difference_gradient(tld.log_density_function, tx))
end

capabilities(ld::DSGELogPosteriorDensity) = LogDensityProblems.LogDensityOrder{1}() ## we provide only first order derivative
capabilities(ld::TransformedLogDensity) = LogDensityProblems.LogDensityOrder{1}() ## we provide only first order derivative

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

function get_initial_value_or_mean(ep::EstimatedParameters)
    return [ismissing(initialvalue) ? mean(prior) : initialvalue for (initialvalue, prior) in zip(ep.initialvalue, ep.prior)]
end

function set_estimated_parameters!(context::Context, value::AbstractVector{T}) where {T<:Real}
    ep = context.work.estimated_parameters
    for (k, p) in enumerate(zip(ep.index, ep.parametertype))
        set_estimated_parameters!(context, p[1], value[k], Val(p[2]))
    end
end

function set_estimated_parameters!(
    context::Context,
    index::Integer,
    value::T,
    ::Val{Parameter},
) where {T<:Real}
    context.work.params[index[1]] = value
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{Endogenous},
) where {T<:Real,N<:Integer}
    context.work.Sigma_m[index[1], index[2]] = value*value
    return nothing
end

function set_estimated_parameters!(
    context::Context,
    index::Pair{N,N},
    value::T,
    ::Val{Exogenous},
) where {T<:Real,N<:Integer}
    context.models[1].Sigma_e[index[1], index[2]] = value*value
    return nothing
end

function loglikelihood(
    θ,
    context::Context,
    observations::Matrix{D},
    ssws::SSWs{T,I},
)::T where {T<:Real,D<:Union{Missing,<:Real},I<:Integer}
    parameters = collect(θ)
    model = context.models[1]
    results = context.results.model_results[1]
    work = context.work
    estimated_parameters = work.estimated_parameters
    model_parameters = work.params

    set_estimated_parameters!(context, parameters)

    fill!(context.results.model_results[1].trends.exogenous_steady_state, 0.0)
    #compute steady state and first order solution
    Dynare.compute_stoch_simul!(
        context,
        ssws.dynamicws,
        model_parameters,
        ssws.stoch_simul_options;
        variance_decomposition = false,
    )
    LRE = LinearRationalExpectations
    LRE_results = results.linearrationalexpectations
    compute_variance!(LRE_results, model.Sigma_e, workspace(LRE.VarianceWs, context))
    # build state space representation
    steady_state = results.trends.endogenous_steady_state[ssws.obs_idx]
    n = size(ssws.Y, 2)
    row = 1
    Dynare.remove_linear_trend!(
        ssws.Y,
        observations,
        results.trends.endogenous_steady_state[ssws.obs_idx],
        results.trends.endogenous_linear_trend[ssws.obs_idx],
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
        data_pattern = Vector{Vector{I}}(undef, I(0))
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

## Estimation
function ml_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    kwargs...,
)
    problem = DSGENegativeLogLikelihood(context, datafile, first_obs, last_obs)
    transformation = DSGETransformation(context.work.estimated_parameters)
    transformed_problem = TransformedLogDensity(transformation, problem)
    (p0, v0) = get_initial_values(context.work.estimated_parameters)
    ip0 = collect(TransformVariables.inverse(transformation, tuple(p0...)))

    res = optimize(
        problem.f,
        p0,
        NelderMead(),
        Optim.Options(f_tol = 1e-5, show_trace = true),
    )

    hess = finite_difference_hessian(problem.f, res.minimizer)
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

function posterior_mode(
    context,
    initialvalue,
    observations;
    first_obs = 1,
    last_obs = 0,
    iterations = 1000,
    show_trace = false
)
    ep = context.work.estimated_parameters
    results = context.results.model_results[1].estimation
    problem = DSGELogPosteriorDensity(context, observations, first_obs, last_obs)
    transformation = DSGETransformation(ep)
    transformed_density(θ) = -problem.f(collect(Dynare.TransformVariables.transform(transformation, θ)))
    transformed_density_gradient!(g, θ) = (g = finite_difference_gradient(transformed_density, θ))
    p0 = initialvalue
    ip0 = collect(TransformVariables.inverse(transformation, tuple(p0...)))
    res = optimize(
        transformed_density,
        ip0,
        LBFGS(),
        Optim.Options(show_trace = show_trace, f_tol = 1e-5, iterations = iterations),
    )
    @show res
    @show res.minimizer
    hess = finite_difference_hessian(transformed_density, res.minimizer)
    inv_hess = inv(hess)
    @show diag(hess)
    @show TransformVariables.transform(transformation, diag(hess))
    @show diag(inv_hess)
    hsd = sqrt.(diag(hess))
    invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
    stdh = sqrt.(diag(invhess))
    std = transform_std(transformation, res.minimizer, stdh)
    mode = collect(TransformVariables.transform(transformation, res.minimizer))
    
    results.posterior_mode = mode
    results.posterior_mode_std = std
    results.posterior_mode_covariance = invhess
    
    estimation_result_table(
        get_parameter_names(ep),
        mode,
        std,
        "Posterior mode"
    )
    return(res, mode, std, inv_hess)
end

function mh_estimation(
    context,
    observations,
    initial_values,
    covariance;
    back_transformation = true,
    first_obs = 1,
    last_obs = 0,
    mcmc_chains = 1,
    mcmc_replic = 100000,
    transformed_parameters = true,
    kwargs...,
)
    ep = context.work.estimated_parameters
    param_names = ep.name
    problem = DSGELogPosteriorDensity(context, observations, first_obs, last_obs)
    #(p0, v0) = get_initial_values(context.work.estimated_parameters)

    if transformed_parameters
        transformation = DSGETransformation(context.work.estimated_parameters)
        transformed_density(θ) = problem.f(collect(TransformVariables.transform(transformation, θ)))
        transformed_density_gradient!(g, θ) = (g = finite_difference_gradient(transformed_density, θ))
        initial_values = collect(TransformVariables.inverse(transformation, tuple(initial_values...)))
        posterior_density(θ) = transformed_density(θ)
        proposal_covariance = inverse_transform_variance(transformation, initial_values, Matrix(covariance))
        chain = run_mcmc(transformed_density, initial_values, proposal_covariance, param_names, mcmc_replic, mcmc_chains)
        if back_transformation
            transform_chains(chain, transformation, problem.f)
            # Recompute posterior density
            nparams = length(ep.prior)
            n1 = nparams + 1
            back_transformed_chain = sample(chain, 1000)
            y = back_transformed_chain.value.data
            for k in axes(y, 1)
                @views begin
                θ1 = y[k, 1:nparams, 1]
                y[k, n1, 1] = posterior_density(θ1)
                end
            end 
        end 
    else
        chain = run_mcmc(problem.f, initial_values, Matrix(covariance), param_names, mcmc_replic, mcmc_chains)
        back_transformed_chain  = chain
    end 
    imode1 = argmax(chain.value.data[:,end,1])
    mode1 = chain.value.data[imode1, 1:end-1, 1]
    return chain, back_transformed_chain
end 

function run_mcmc(posterior_density, initial_values, proposal_covariance, param_names, mcmc_replic, mcmc_chains)    
    model = DensityModel(posterior_density)
    spl = RWMH(MvNormal(zeros(length(initial_values)), proposal_covariance))
    if mcmc_chains == 1
        chain = AbstractMCMC.sample(model, spl, mcmc_replic,
                   init_params = initial_values,
                   param_names = context.work.estimated_parameters.name,
                   chain_type = Chains)
                   display_acceptance_ratio(spl, mcmc_replic)
    else    
        old_blas_threads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        chain = AbstractMCMC.sample(model, spl, MCMCThreads(), mcmc_replic,
                   mcmc_chains,
                   init_params = Iterators.repeated(initial_values),
                   param_names = context.work.estimated_parameters.name,
                   chain_type = Chains)
        BLAS.set_num_threads(old_blas_threads)
    end 
    return chain
end 

function hmc_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    iterations = 1000,
    initial_values = Vector{Float64}(undef, 0),
    initial_energy = Vector{Float64}(undef, 0),
    kwargs...,
)
    problem = DSGELogPosteriorDensity(context, datafile, first_obs, last_obs)
    transformation = DSGETransformation(context.work.estimated_parameters)
    transformed_problem = TransformedLogDensity(transformation, problem.f)
    transformed_logdensity(θ) = TransformVariables.transform_logdensity(transformed_problem.transformation,
    transformed_problem.log_density_function,
    θ)
    function logdensity_and_gradient(ℓ, x)
        return(transformed_logdensity(x),
               FiniteDiff.finite_difference_gradient(transformed_logdensity, x))
    end

    if ismissing(initial_values)
        (p0, v0) = get_initial_values(context.work.estimated_parameters)
        ip0 = collect(TransformVariables.inverse(transformation, tuple(p0...)))
    else
        p0, v0 = collect(initial_values), collect(initial_energy)
    end
    @show v0
    results = DynamicHMC.mcmc_keep_warmup(
        Random.GLOBAL_RNG,
        #        transformed_problem,
        transformed_problem,
        30;
        initialization = (q = p0, κ = GaussianKineticEnergy(diagm(0 => Vector{Float64}(v0)))),
        warmup_stages = default_warmup_stages(),
#        reporter = ProgressMeterReport(),
    )
    @show results
    parameter_names = get_parameter_names(context.work.estimated_parameters)
    estimation_result_table(
        parameter_names,
        mean(results.inference.chain),
        sqrt.(diag(results.κ.M⁻¹)),
        "Results from Bayesian estimation",
    )
    return results
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

function hmcmc_result_table(parameter_names)
    mean(results.chain)
    sqrt.(diag(results.κ.M⁻¹))
    "Results from Bayesian estimation"
end

function display_acceptance_ratio(spl::MetropolisHastings,replic)
    println("")
    println("Acceptance ratio MCMC chain: $(spl.n_acceptances/replic)")
    println("")
end 

function display_acceptance_ratio(spl::Vector{MetropolisHastings},replic)
    println("")
    for (i, s) in enumerate(spl) 
        println("Acceptance ratio chain $i: $(s.n_acceptances/replic)")
    end
    println("")
end 

# Utilities
function get_eigenvalues(context)
    eigenvalues = context.results.model_results[1].linearrationalexpectations.eigenvalues
    return eigenvalues
end

function get_initial_values(
    estimated_parameters::EstimatedParameters,
)::Tuple{Vector{Float64},Vector{Float64}}
    return (mean.([p for p in estimated_parameters.prior]),
            var.([p for p in estimated_parameters.prior])
            )
end

function get_parameter_names(estimated_parameters::EstimatedParameters)
    return [get_parameter_name(n) for n in estimated_parameters.name]
end

get_parameter_name(name::String) = name
get_parameter_name(name::Pair{String, String}) = name[1]

function get_observations(context, datafile, first_obs, last_obs)
    results = context.results.model_results[1]
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    obs_idx =
        [symboltable[v].orderintype for v in varobs if is_endogenous(v, symboltable)]

    if datafile != ""
        varnames = [v for v in varobs if is_endogenous(v, symboltable)]
        Yorig = get_data(datafile, varnames, start = first_obs, last = last_obs)
    else
        error("calib_smoother needs a data file or a TimeDataFrame!")
    end
    Y = Matrix{Union{Float64,Missing}}(undef, size(Yorig))

    remove_linear_trend!(
        Y,
        Yorig,
        results.trends.endogenous_steady_state[obs_idx],
        results.trends.endogenous_linear_trend[obs_idx],
    )
    return Y
end

# Pretty printing name of estimated parameters
print_parameter_name(::Val{Parameter}, name) = name
print_parameter_name(::Val{Endogenous}, name) = "Std($(name[1]))"
print_parameter_name(::Val{Exogenous}, name) = "Std($(name[1]))"

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


function plot_priors(context, names, n_points = 100)
    indices = [find(context.work.estimated_parameters.name, e) for e in names]
    prior_pdfs = []
    n_plots = length(indices)
    prior_x_axes = []
    for i in indices
        prior = context.work.estimated_parameters.prior[i]
        m, v = mean(prior), var(prior)
        prior_x_axis = LinRange(m-15*v, m+15*v, n_points)
        prior_pdf = [pdf(prior, e) for e in prior_x_axis]
        push!(prior_pdfs, prior_pdf/sum(prior_pdf))
        push!(prior_x_axes, prior_x_axis)
    end
    prior_pdfs = hcat(prior_pdfs...)
    prior_x_axes = hcat(prior_x_axes...)
    names = hcat(names...)
    f = Plots.plot(prior_x_axes, prior_pdfs, layout=n_plots, title=names, linecolor=:darkgrey, labels=false, linewidth=3)
    f
end

function transform_chains(chains, t, posterior_density)
    y = chains.value.data
    nparams = size(y, 2) - 1
    tmp = Vector{Float64}(undef, size(y, 1))
    for i in axes(y, 3)
        for j in 1:nparams
            vy = view(y, :, j, i)   
            tmp .= map(TransformVariables.transform(t.transformations[j]),
                     vy)
            vy .= tmp
        end          
    end
end 

# alias for InverseGamma distribution
InverseGamma2 = InverseGamma

struct stdev
    s::Symbol
    function stdev(s_arg::Symbol)
        s = s_arg
        new(s)
    end
end

struct corr
    s1::Symbol
    s2::Symbol
    function corr(s1_arg, s2_arg)
        s1 = s1_arg
        s2 = s2_arg
        new(s1, s2)
    end 
end 

function prior(s::Symbol; shape, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    name = string(s)
    push!(ep.index, symboltable[name].orderintype)
    push!(ep.name, name)
    push!(ep.initialvalue, initialvalue)
    push!(ep.parametertype, Parameter)
    prior_(s, shape, mean, stdev, domain, variance, ep)
end
    
function prior(s::stdev; shape, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    name = string(s.s)
    push!(ep.index, symboltable[name].orderintype)
    push!(ep.name, name)
    push!(ep.initialvalue, initialvalue)
    push!(ep.parametertype, symboltable[name].symboltype)
    prior_(s, shape, mean, stdev, domain, variance, ep)
end

function prior(s::corr; shape, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    name1 = string(s.s1)
    name2 = string(s.s2)
    index1 = symboltable[name1].orderintype
    index2 = symboltable[name2].orderintype
    push!(ep.index, (index1 =>index2))
    push!(ep.name, (name1 => name2))
    push!(ep.initialvalue, initialvalue)
    push!(ep.parametertype, symboltable[name1].symboltype)
    prior_(s, shape, mean, stdev, domain, variance, ep)
end

function prior_(s, shape, mean, stdev, domain, variance, ep)
    if ismissing(variance) && !ismissing(stdev)
        variance = stdev*stdev
    end 
    if shape == Beta
        α, β = beta_specification(mean, variance)
        push!(ep.prior, Beta(α, β))
    elseif shape == Gamma
        α, β = gamma_specification(mean, variance) 
        push!(ep.prior, Gamma(α, β))
    elseif shape == InverseGamma1
        α, β = inverse_gamma_1_specification(mean, variance) 
        push!(ep.prior, InverseGamma1(α, β))
    elseif shape == InverseGamma2
        α, β = inverse_gamma_2_specification(mean, variance)
        push!(ep.prior, InverseGamma(α, β))
    elseif shape == Normal
        push!(ep.prior, Normal(mean, stdev))
    elseif shape == Uniform
        if isempty(domain)
            a, b = uniform_specification(mean, variance)
        else
            a, b = domain
        end 
        push!(ep.prior, Uniform(a, b))
    elseif shape == Weibull
        α, β = weibull_specification(mean, stdev)
        push!(ep.prior, Weibull(α, β))
    else
        erro("Unknown prior shape")
    end 
end         

function get_index_name(s::Symbol, symboltable::SymbolTable)
    name = String(s)
    index = symboltable[name].orderintype
    return (index, name)
end 

function plot_prior_posterior(context::Context, chains)
    ep = context.work.estimated_parameters
    mode = context.results.model_results[1].estimation.posterior_mode
    @assert length(ep.prior) > 0 "There is no defined prior"
    @assert length(chains) > 0 "There is no MCMC chain"

    path = "$(context.modfileinfo.modfilepath)/graphs/"
    mkpath(path)
    filename = "$(path)/PriorPosterior"
    nprior = length(ep.prior)
    (nbplt, nr, nc, lr, lc, nstar) = pltorg(nprior)
    ivars = collect(1:nr*nc)
    for p = 1:nbplt-1
        filename1 = "$(filename)_$(p).png"
        plot_panel_prior_posterior(
            ep.prior,
            chains,
            mode,
            "Priors",
            ep.name,
            ivars,
            nr,
            nc,
            nr * nc,
            filename1
        )
        ivars .+= nr * nc
    end
    ivars = ivars[1:nstar]
    filename = "$(filename)_$(nbplt).png"
    plot_panel_prior_posterior(
        ep.prior,
        chains,
        mode,
        "Priors",
        ep.name,
        ivars,
        lr,
        lc,
        nstar,
        filename
    )
end

function plot_panel_prior_posterior(
    prior,
    chains,
    mode,
    title,
    ylabels,
    ivars,
    nr,
    nc,
    nstar,
    filename;
    kwargs...
)
    sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]
    for (i, j) in enumerate(ivars)
        posterior_density = kde(vec(get(chains, Symbol(ylabels[j]))[1].data))
        
        sp[i] = Plots.plot(prior[j], title = ylabels[j], labels = "Prior", kwargs...)
        plot!(posterior_density, labels = "Posterior")
    end

    pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title = "Prior and posterior distributions")
    graph_display(pl)
    savefig(filename)
end

