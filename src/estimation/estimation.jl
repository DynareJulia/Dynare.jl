using AdvancedMH
using Distributions
using DynamicHMC
using FiniteDiff: finite_difference_gradient, finite_difference_hessian
using LogDensityProblems
using MCMCChains
using Optim
using PolynomialMatrixEquations: UndeterminateSystemException, UnstableSystemException
using Random
using TransformVariables

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
    function SSWs(context, nobs, varobs)
        model = context.models[1]
        D = eltype(context.work.params)
        symboltable = context.symboltable
        tmp_nbr = context.dynarefunctions.dynamic!.tmp_nbr
        ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both
        dynamicws = Dynare.DynamicWs(
            model.endogenous_nbr,
            model.exogenous_nbr,
            ncol,
            sum(tmp_nbr[1:2]),
        )
        stoch_simul_options = Dynare.StochSimulOptions(Dict{String,Any}())
        varobs_ids = [
            symboltable[v].orderintype for
            v in varobs if Dynare.is_endogenous(v, symboltable)
        ]
        state_ids = union(varobs_ids, model.i_bkwrd_b)
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
        )
    end
end

struct DSGELogPosteriorDensity
    f::Function
    dimension::Integer
    context::Context
    observations::Matrix{Union{Missing,<:Real}}
    ssws::SSWs
    function DSGELogPosteriorDensity(context, datafile, first_obs, last_obs)
        n = length(context.work.estimated_parameters)
        observations = get_observations(context, datafile, first_obs, last_obs)
        ssws = SSWs(context, size(observations, 2), context.work.observed_variables)
        f = make_logposteriordensity(context, observations, first_obs, last_obs, ssws)
        new(f, n, context, observations, ssws)
    end
end

function logpriordensity(x, estimated_parameters)
    lpd = 0.0
    k = 1
    for p in estimated_parameters.prior_R
        lpd += logpdf(p.prior, x[k])
        k += 1
    end
    for p in estimated_parameters.prior_Rplus
        lpd += logpdf(p.prior, x[k])
        k += 1
    end
    for p in estimated_parameters.prior_01
        lpd += logpdf(p.prior, x[k])
        k += 1
    end
    for p in estimated_parameters.prior_AB
        lpd += logpdf(p.prior, x[k])
        k += 1
    end
    return lpd
end

function DSGEtransformation(ep::EstimatedParameters)
    tvec = []
    for p in ep.prior_R
        push!(tvec, asℝ)
    end
    for p in ep.prior_Rplus
        push!(tvec, asℝ₊)
    end
    for p in ep.prior_01
        push!(tvec, as𝕀)
    end
    for p in ep.prior_AB
        push!(tvec, as(Real, p.domain...))
    end
    return as((tvec...,))
end

function make_logposteriordensity(context, observations, first_obs, last_obs, ssws)
    function logposteriordensity(x)
        lpd = -Inf
        try
            lpd = logpriordensity(x, context.work.estimated_parameters)
            lpd += loglikelihood(x, context, observations, ssws)
        catch e
            @show e
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
    ld.f(collect(x)), finite_difference_gradient(ld.f, collect(x))
capabilities(ld::DSGELogPosteriorDensity) = LogDensityProblems.LogDensityOrder{1}() ## we provide only first order derivative

(problem::DSGELogPosteriorDensity)(θ) =
    logposteriordensity(collect(θ), problem.context, problem.observations, problem.ssws)

function set_estimated_parameters!(context::Context, value::Vector{T}) where {T<:Real}
    ep = context.work.estimated_parameters
    k = 1
    for p in ep.prior_R
        set_estimated_parameters!(context, p.index, value[k], Val(p.parametertype))
        k += 1
    end
    for p in ep.prior_Rplus
        set_estimated_parameters!(context, p.index, value[k], Val(p.parametertype))
        k += 1
    end
    for p in ep.prior_01
        set_estimated_parameters!(context, p.index, value[k], Val(p.parametertype))
        k += 1
    end
    for p in ep.prior_AB
        set_estimated_parameters!(context, p.index, value[k], Val(p.parametertype))
        k += 1
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
    context.models[1].Sigma_e[index[1], index[2]] = value
    return nothing
end

function loglikelihood(
    parameters::Vector{T},
    context::Context,
    observations::Matrix{D},
    ssws::SSWs{T,I},
) where {T<:Real,D<:Union{Missing,<:Real},I<:Integer}

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
    kalman_statevar_ids = collect(1:model.endogenous_nbr)
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
    klws = Dynare.KalmanLikelihoodWs(ny, ns, np, nobs)
    Y = ssws.Y
    if any(ismissing.(Y))
        # indices of non-missing observations are in data_pattern
        data_pattern = Vector{Vector{I}}(undef, I(0))
        for i = 1:nobs
            push!(data_pattern, findall(.!ismissing.(view(Y[:, i]))))
        end
        return Dynare.kalman_likelihood(
            Y,
            swws.Z,
            ssws.H,
            ssws.T,
            ssws.R,
            ssws.Q,
            ssws.a0,
            swws.P,
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


pervious_optimum = Inf
function minimas_unstable!(m, x)
    fill!(m, Inf)
    for xx in x
        xxx = abs(xx)
        if xxx > 1
            x1 = xxx - 1
            for mm in m
                x1 = xxx - 1
                if x1 < mm
                    mm = x1
                end
            end
        end
    end
    return m
end

function maximas_stable!(m, x)
    fill!(m, 0.0)
    for xx in x
        xxx = abs(x)
        if xxx < 1
            x1 = 1 - xxx
            for mm in m
                if x1 > mm
                    mm = 1 - x1
                end
            end
        end
    end
    return m
end

function penalty(
    eigenvalues::AbstractVector{<:Union{T,Complex{T}}},
    forward_nbr::Integer,
) where {T<:Real}
    n = length(eigenvalues)
    unstable_nbr = count(abs.(eigenvalues) .> 1.0)
    excess_unstable_nbr = unstable_nbr - forward_nbr
    if excess_unstable_nbr > 0
        return sum(
            minimas_unstable!(view(optimum_work, 1:excess_unstable_nbr), eigenvalues),
        )
    else
        return sum(
            minimas_unstable!(view(optimum_work, 1:-excess_unstable_nbr), eigenvalues),
        )
    end
end

function get_symbol(symboltable, indx)
    for k in symboltable
        if k[2].symboltype == SymbolType(3) && k[2].orderintype == indx
            return k[2].longname
        end
    end
end


function maximum_likelihood(
    params::Vector{T},
    params_indices,
    shock_variance_indices,
    measurement_variance_indices,
    varobs,
    observations,
    context,
    ssws,
) where {T<:Real}

    optimum_work = Vector{Float64}(undef, context.models[1].endogenous_nbr)
    history = 0.0

    # objective function
    function negative_loglikelihood(p)
        try
            f =
                -logposteriordensity(
                    p,
                    params_indices,
                    shock_variance_indices,
                    measurement_variance_indices,
                    varobs,
                    observations,
                    context,
                    ssws,
                )
            history = f
            return f
        catch e
            if e isa Union{UndeterminateSystemException,UnstableSystemException}
                @debug e, p
                model = context.models[1]
                forward_nbr = model.n_fwrd + model.n_both
                return penalty(
                    eigvals(
                        context.results.model_results[1].linearrationalexpectations.gs1,
                    ),
                    forward_nbr,
                )
            elseif e isa DomainError
                @debug DomainError, p
                rethrow(e)
                msg = sprint(showerror, e, catch_backtrace())
                x = parse(T, rsplit(rsplit(msg, ":")[1], " ")[3])
                #println(x)
                return abs(x) + history
            else
                rethrow(e)
            end
        end
    end

    initial_values = get_initial_values(context.work.estimated_parameters)
    f(p) = negative_loglikelihood(p)
    res = optimize(
        f,
        initial_values,
        NelderMead(),
        Optim.Options(f_tol = 1e-5, show_trace = true),
    )

    hess = finite_difference_hessian(negative_loglikelihood, res.minimizer)
    inv_hess = inv(hess)
    hsd = sqrt.(diag(hess))
    invhess = inv(hess ./ (hsd * hsd')) ./ (hsd * hsd')
    stdh = sqrt.(diag(invhess))

    maximum_likelihood_result_table(
        context.symboltable,
        params_indices,
        res.minimizer,
        stdh,
    )
end


function metropolis_hastings(
    params,
    params_indices,
    shock_variance_indices,
    measurement_variance_indices,
    varobs,
    observations,
    mh_iter,
    invhess,
    context,
    ssws,
)
    function loglikelihood_density(params)
        try
            f = loglikelihood(
                params,
                params_indices,
                shock_variance_indices,
                measurement_variance_indices,
                varobs,
                observations,
                context,
                ssws,
            )
            return f
        catch
            return -Inf
        end
    end


    model = DensityModel(loglikelihood_density)

    spl = RWMH(MvNormal(zeros(length(params_indices)), 0.5 * (invhess + invhess')))

    # Sample from the posterior.
    chain = sample(model, spl, 50000; param_names = parameter_names, chain_type = Chains)
end

function hamiltonian_mcmc(problem, p0) end

function estimation(context)
    varobs = context.varobs
    observations = copy(transpose(Matrix(Dynare.simulation(varobs))))
    nobs = size(observations, 2)
    ssws = SSWs(context, nobs, varobs)
    # estimated parameters: rho, alpha, theta, tau
    params_indices = [2, 3, 5, 6]
    # no shock_variance estimated
    shock_variance_indices = Vector{Int}(undef, 0)
    # no measurement errors in the model
    measurement_variance_indices = Vector{Int}(undef, 0)
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

function get_indices(estimated_parameters, symboltable)
    parameters_indices = Vector{Int64}(undef, 0)
    shock_variance_indices = Vector{Int64}(undef, 0)
    measurement_variance_indices = Vector{Int64}(undef, 0)

    pp(i, ::Val{Endogenous}) = push!(measurement_variance_indices, i)
    pp(i, ::Val{Exogenous}) = push!(shock_variance_indices, i)
    pp(i, ::Val{Parameter}) = push!(parameters_indices, i)

    for name in estimated_parameters.name
        p = symboltable[name]
        pp(p.orderintype, Val(p.symboltype))
    end

    return (parameters_indices, shock_variance_indices, measurement_variance_indices)
end

function ml_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    optim_algo = "optim",
    kwargs...,
)
    problem = DSGELogPosteriorDensity(logposteriordensity, datafile, first_obs, last_obs)
    transformation = DSGETransformation(context.work.estimated_parameters)
    transformed_problem = TransformedLogDensity(transformation, problem)

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
end

function get_initial_values(
    estimated_parameters::EstimatedParameters,
)::Tuple{Vector{Float64},Vector{Float64}}
    return (
        [
            mean.([p.prior for p in estimated_parameters.prior_R])
            mean.([p.prior for p in estimated_parameters.prior_Rplus])
            mean.([p.prior for p in estimated_parameters.prior_01])
            mean.([p.prior for p in estimated_parameters.prior_AB])
        ],
        [
            var.([p.prior for p in estimated_parameters.prior_R])
            var.([p.prior for p in estimated_parameters.prior_Rplus])
            var.([p.prior for p in estimated_parameters.prior_01])
            var.([p.prior for p in estimated_parameters.prior_AB])
        ],
    )
end

function get_parameter_names(estimated_parameters::EstimatedParameters)
    return [
        [p.name for p in estimated_parameters.prior_R]
        [p.name for p in estimated_parameters.prior_Rplus]
        [p.name for p in estimated_parameters.prior_01]
        [p.name for p in estimated_parameters.prior_AB]
    ]
end

function mh_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    optim_algo = "optim",
    kwargs...,
)

    observations = get_observations(context, datafile, first_obs, last_obs)
    observed_variables = context.work.observed_variables
    estimated_parameters = context.work.estimated_parameters
    params = get_initial_values(estimated_parameters)
    (params_indices, shock_variance_indices, measurement_variance_indices) =
        get_indices(estimated_parameters, context.symboltable)
    ssws = SSWs(context, size(observations, 2), observed_variables)
    mh_iter = 10000

    metropolis_hastings(
        params,
        params_indices,
        shock_variance_indices,
        measurement_variance_indices,
        observed_variables,
        observations,
        mh_iter,
        invhess,
        context,
        ssws,
    )
end

function hmc_estimation(
    context;
    datafile = "",
    first_obs = 1,
    last_obs = 0,
    optim_algo = "optim",
    kwargs...,
)
    problem = DSGELogPosteriorDensity(context, datafile, first_obs, last_obs)
    transformation = DSGEtransformation(context.work.estimated_parameters)
    transformed_problem = TransformedLogDensity(transformation, problem)
    (p0, v0) = get_initial_values(context.work.estimated_parameters)
    results = mcmc_with_warmup(
        Random.GLOBAL_RNG,
        problem,
        1000;
        initialization = (q = p0, κ = GaussianKineticEnergy(I(problem.dimension))),
        reporter = ProgressMeterReport(),
    )
    parameter_names = get_parameter_names(context.work.estimated_parameters)
    estimation_result_table(
        parameter_names,
        mean(results.chain),
        sqrt.(diag(results.κ.M⁻¹)),
        "Results from Bayesian estimation",
    )
end

function estimation_result_table(param_names, estimated_params, stdh, title)
    table = Matrix{Any}(undef, length(param_names) + 1, 4)
    table[1, 1] = "Parameter"
    table[1, 2] = "Estimated value"
    table[1, 3] = "Standard error"
    table[1, 4] = "80% confidence interval"
    for (i, k) in enumerate(param_names)
        table[i+1, 1] = k
        table[i+1, 2] = estimated_params[i]
        table[i+1, 3] = stdh[i]
        table[i+1, 4] = estimated_params[i] ± (1.28 * stdh[i])
    end
    dynare_table(table, title)
end
