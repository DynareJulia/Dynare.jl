using AdvancedMH
using Distributions
using DynamicHMC
using FiniteDiff: finite_difference_gradient, finite_difference_hessian
using KernelDensity
using LogDensityProblems
using MCMCChains
using Optim
using Plots
using PolynomialMatrixEquations: UndeterminateSystemException, UnstableSystemException
using Random
using TransformVariables
using TransformedLogDensities

import LogDensityProblems: dimension, logdensity, logdensity_and_gradient, capabilities

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
    varobs_ids::Vector{I}
    Y::Matrix{Union{D,Missing}}
    Z::Matrix{D}
    kalman_ws::KalmanLikelihoodWs{D, I}
    function SSWs(context, nobs, varobs)
        model = context.models[1]
        D = eltype(context.work.params)
        symboltable = context.symboltable
        tmp_nbr = sum(model.dynamic_tmp_nbr[1:2])
        ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both
        dynamicws = Dynare.DynamicWs(
            model.endogenous_nbr,
            model.exogenous_nbr,
            tmp_nbr,
            model.dynamic_g1_sparse_colptr,
            model.dynamic_g1_sparse_rowval
        )
        stoch_simul_options = Dynare.StochSimulOptions(Dict{String,Any}())
        varobs_ids0 = [
            symboltable[v].orderintype for
            v in varobs]
        state_ids = sort!(union(varobs_ids0, model.i_bkwrd_b))
        varobs_ids = [findfirst(isequal(symboltable[v].orderintype), state_ids) for v in varobs]
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
            varobs_ids,
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
    function DSGELogPosteriorDensity(context, datafile, first_obs, last_obs)
        n = length(context.work.estimated_parameters)
        observations = get_observations(context, datafile, first_obs, last_obs)
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
        push_prior!(tvec, Val(typeof(p)))
    end
    return as((tvec...,))
end

push_prior!(tvec, ::Val{Distributions.Beta{Float64}}) = push!(tvec, as𝕀) 
push_prior!(tvec, ::Val{Distributions.Gamma{Float64}}) = push!(tvec, asℝ₊) 
push_prior!(tvec, ::Val{Dynare.InverseGamma1{Float64}}) = push!(tvec, asℝ₊) 
push_prior!(tvec, ::Val{Distributions.InverseGamma{Float64}}) = push!(tvec, asℝ₊) 
push_prior!(tvec, ::Val{Distributions.Normal{Float64}}) = push!(tvec, asℝ) 
push_prior!(tvec, ::Val{Distributions.Uniform{Float64}}) = push!(tvec, as(Real, p.domain...)) 
push_prior!(tvec, ::Val{Distributions.Weibull{Float64}}) = push!(tvec, asℝ₊) 

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
            # eigenvalues = get_eigenvalues(context)
            # model = context.models[1]
            # backward_nbr = model.n_bkwrd + model.n_both
            # @debug e
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
        @show e
        @debug e
        lpd = -Inf
    end
    return lpd
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
    varobs = work.observed_variables

    set_estimated_parameters!(context, parameters)
    #compute steady state and first order solution
    Dynare.compute_stoch_simul!(
        context,
        ssws.dynamicws,
        model_parameters,
        ssws.stoch_simul_options;
        variance_decomposition = false,
    )
    # build state space representation
    steady_state = results.trends.endogenous_steady_state[ssws.varobs_ids]
    linear_trend_coeffs = results.trends.endogenous_linear_trend[ssws.varobs_ids]
    n = size(ssws.Y, 2)
    row = 1
    linear_trend = collect(row - 1 .+ (1:n))
    Dynare.remove_linear_trend!(
        ssws.Y,
        observations,
        results.trends.endogenous_steady_state[ssws.varobs_ids],
        results.trends.endogenous_linear_trend[ssws.varobs_ids],
    )
    ns = length(ssws.state_ids)
    np = model.exogenous_nbr
    ny, nobs = size(ssws.Y)
    varobs_ids = ssws.varobs_ids
    for i = 1:ny
        ssws.Z[i, ssws.varobs_ids[i]] = T(1)
    end
    vg1 = view(results.linearrationalexpectations.g1_1, ssws.state_ids, :)
    view(ssws.T, :, ssws.lagged_state_ids) .= vg1
    vg2 = view(results.linearrationalexpectations.g1_2, ssws.state_ids, :)
    ssws.R .= vg2
    ssws.Q .= model.Sigma_e
    fill!(ssws.a0, T(0))
    ssws.P .= view(
        context.results.model_results[1].linearrationalexpectations.endogenous_variance,
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
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    iterations = 1000,
    show_trace = false
)
    @debug show_trace = true
    problem = DSGELogPosteriorDensity(context, datafile, first_obs, last_obs)
    transformation = DSGETransformation(context.work.estimated_parameters)
    transformed_problem = TransformedLogDensity(transformation, problem)
    transformed_density(θ) = -problem.f(collect(Dynare.TransformVariables.transform(transformation, θ)))
    transformed_density_gradient!(g, θ) = (g = finite_difference_gradient(transformed_density, θ))
    (p0, v0) = get_initial_values(context.work.estimated_parameters)
    ip0 = collect(TransformVariables.inverse(transformation, tuple(p0...)))
    res = optimize(
        transformed_density,
        ip0,
        LBFGS(),
#        Optim.Options(f_tol = 1e-5, show_trace = show_trace, iterations = iterations),
        Optim.Options(show_trace = show_trace, iterations = iterations),
    )
    @show res
    hess = finite_difference_hessian(transformed_density, res.minimizer)
    inv_hess = inv(hess)
    hsd = sqrt.(diag(hess))
    invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
    stdh = sqrt.(diag(invhess))
    tstdh = transform_std(transformation, res.minimizer, stdh)
    mode = collect(TransformVariables.transform(transformation, res.minimizer))
    
    if length(context.work.estimated_parameters.posterior_mode) == 0
        resize!(context.work.estimated_parameters.posterior_mode, length(mode))
    end
    context.work.estimated_parameters.posterior_mode .= mode
    
    estimation_result_table(
        context.work.estimated_parameters.name,
        context.work.estimated_parameters.parametertype,
        mode,
        tstdh,
        "Posterior mode"
    )
    return(res, mode, tstdh, inv_hess)
end

function mh_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    iterations = 100000,
    initial_values = Vector{Float64}(undef, 0),
    covariance = Matrix{Float64}(undef, 0, 0),
    kwargs...,
)

    problem = DSGELogPosteriorDensity(context, datafile, first_obs, last_obs)
    transformation = DSGETransformation(context.work.estimated_parameters)
    transformed_problem = TransformedLogDensity(transformation, problem)
    transformed_density(θ) = problem.f(collect(Dynare.TransformVariables.transform(transformation, θ)))
    transformed_density_gradient!(g, θ) = (g = finite_difference_gradient(transformed_density, θ))
    #    (p0, v0) = get_initial_values(context.work.estimated_parameters)
    #    ip0 = collect(TransformVariables.inverse(transformation, tuple(initial_values...)))

    model = DensityModel(transformed_density)

    spl = RWMH(MvNormal(zeros(length(initial_values)), Matrix(covariance)))

    # Sample from the posterior.
    chain = sample(model, spl, iterations,
                   init_params = initial_values,
                   param_names = context.work.estimated_parameters.name,
                   chain_type = Chains)
    @show spl.n_acceptances
    display(my_chains)
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
function estimation_result_table(param_names, param_type, estimated_value, stdh, title)
    table = Matrix{Any}(undef, length(param_names) + 1, 4)
    table[1, 1] = "Parameter"
    table[1, 2] = "Estimated value"
    table[1, 3] = "Standard error"
    table[1, 4] = "80% confidence interval"
    for (i, k) in enumerate(zip(param_type, param_names))
        table[i+1, 1] = print_parameter_name(Val(k[1]), k[2])
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
    varobs_ids =
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
        results.trends.endogenous_steady_state[varobs_ids],
        results.trends.endogenous_linear_trend[varobs_ids],
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

#transform_std(::asℝ, x, std) = std
transform_std(::TransformVariables.ShiftedExp{true, Int64}, x, std) = exp(x)*std
transform_std(::TransformVariables.ScaledShiftedLogistic{Int64}, x, std) = logistic(x)*(1-logistic(x))*std

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
    f = Plots.plot(prior_x_axes, prior_pdfs, layout=n_plots, title=names, linecolor=:darkgrey, legend=false, linewidth=3)
    f
end

