using Dynare
using PolynomialMatrixEquations: UndeterminateSystemException, UnstableSystemException
using Optim
using FiniteDiff: finite_difference_hessian
using LinearAlgebra: inv, diag, eigvals
#using ForwardDiff: hessian

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

function estimated_parameters!(
    context::Context,
    estimated_params::Vector{D},
    param_indices::Vector{I},
    shock_variance_indices::Vector{I},
    measurement_variance_indices::Vector{I},
)  where {D<:AbstractFloat,I<:Integer}
    k = 1
    for j in param_indices
        context.work.params[j] = estimated_params[k]
        k += 1
    end
    for j in shock_variance_indices
        # !!! linear indices to variance matrix
        context.model.Sigma_e[j] = estimated_params[k]
        k += 1
    end
    for j in measurement_variance_indices
        # not implemented
        k += 1
    end
end

function loglikelihood(
    estimated_params::Vector{D},
    params_indices::Vector{I},
    shock_variance_indices::Vector{I},
    measurement_variance_indices::Vector{I},
    varobs::Tuple{String,String},
    observations::Matrix{D},
    context::Dynare.Context,
    ssws::SSWs{D,I},
) where {D<:AbstractFloat,I<:Integer}
    estimated_parameters!(context,
                          estimated_params,
                          params_indices,
                          shock_variance_indices,
                          measurement_variance_indices)
    model = context.models[1]
    params = context.work.params
    results = context.results.model_results[1]

    #compute steady state and first order solution
    Dynare.compute_stoch_simul!(
        context,
        ssws.dynamicws,
        params,
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
        ssws.Z[i, ssws.varobs_ids[i]] = D(1)
    end
    vg1 = view(results.linearrationalexpectations.g1_1, ssws.state_ids, :)
    view(ssws.T, :, ssws.lagged_state_ids) .= vg1
    vg2 = view(results.linearrationalexpectations.g1_2, ssws.state_ids, :)
    ssws.R .= vg2
    ssws.Q .= model.Sigma_e
    fill!(ssws.a0, D(0))
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
        Dynare.kalman_likelihood(
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
        Dynare.kalman_likelihood(
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

# generate artificial data with model example5
# provide model parsing
context = @dynare "test/models/example5/example5_est_a.mod";

# simulation is an accessor to simulated series
# we assume that we observe output, y, and consumption, c
varobs = ("y", "c")
observations = copy(transpose(Matrix(Dynare.simulation(varobs))))
nobs = size(observations, 2)
ssws = SSWs(context, nobs, varobs)
# estimated parameters: rho, alpha, theta, tau
params_indices = [2, 3, 5, 6]
# no shock_variance estimated
shock_variance_indices = Vector{Int}(undef, 0)
# no measurement errors in the model
measurement_variance_indices = Vector{Int}(undef, 0)

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

optimum_work = Vector{Float64}(undef, context.models[1].endogenous_nbr)    
function penalty(eigenvalues::AbstractVector,
                 forward_nbr::Integer)
    n = length(eigenvalues)
    unstable_nbr = count(abs.(eigenvalues) .> 1.0)
    excess_unstable_nbr = unstable_nbr - forward_nbr
    if excess_unstable_nbr > 0
        return sum(minimas_unstable!(view(optimum_work, 1:excess_unstable_nbr), eigenvalues))
    else
        return sum(minimas_unstable!(view(optimum_work, 1:-excess_unstable_nbr), eigenvalues))
    end
end

function negative_loglikelihood(
    params::Vector{T}
) where {T<:AbstractFloat}
    try
        return -loglikelihood(params,
                              params_indices,
                              shock_variance_indices,
                              measurement_variance_indices,
                              varobs,
                              observations,
                              context,
                              ssws)
    catch e
        if e isa Union{UndeterminateSystemException,UnstableSystemException}
            model = context.models[1]
            forward_nbr = model.n_fwrd + model.n_both
            return penalty(eigvals(context.results.model_results[1].linearrationalexpectations.gs1),
                           forward_nbr)
        else
            rethrow(e)
        end
    end
end

# objective function
# estimated parameters: rho, alpha, theta,tau
init_guess = [0.95, 0.36, 2.95, 0.025] 
f(p) = negative_loglikelihood(p)
res = optimize(f, init_guess, NelderMead())

#hess = hessian(negative_loglikelihood, res.minimizer)
hess = finite_difference_hessian(negative_loglikelihood, res.minimizer)
#println(hess)
inv_hess = inv(hess)
println(diag(hess))
hsd = sqrt.(diag(hess))
invhess = inv(hess./(hsd*hsd'))./(hsd*hsd')
println(diag(invhess))
stdh = sqrt.(diag(invhess))
println("variance: ", stdh)

