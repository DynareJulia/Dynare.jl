using Dynare

struct SSWs
    a0::Vector{Float64}
    dynamicws::Dynare.DynamicWs
    H::Matrix{Float64}
    lagged_state_ids::Vector{Int64}
    P::Matrix{Float64}
    Q::Matrix{Float64}
    R::Matrix{Float64}
    state_ids::Vector{Int64}
    stoch_simul_options::Dynare.StochSimulOptions
    T::Matrix{Float64}
    varobs_ids::Vector{Int64}
    Y::Matrix{Union{Float64, Missing}}
    Z::Matrix{Float64}
    function SSWs(context, nobs, varobs )
        model = context.models[1]
        symboltable = context.symboltable
        tmp_nbr = context.dynarefunctions.dynamic!.tmp_nbr
        ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both
        dynamicws =  Dynare.DynamicWs(model.endogenous_nbr, model.exogenous_nbr, ncol, sum(tmp_nbr[1:2]))
        stoch_simul_options = Dynare.StochSimulOptions(Dict{String, Any}())
        varobs_ids =
            [symboltable[v].orderintype for v in varobs if Dynare.is_endogenous(v, symboltable)]
        state_ids = union(varobs_ids, model.i_bkwrd_b)
        lagged_state_ids = findall(in(model.i_bkwrd_b), state_ids)
        np = model.exogenous_nbr
        ns = length(state_ids)
        ny = length(varobs)
        Y = Matrix{Union{Float64, Missing}}(undef, ny, nobs)
        Z = zeros(ny, ns)
        H = zeros(ny, ny)
        T = zeros(ns, ns)
        R = zeros(ns, np)
        Q = zeros(np, np)
        a0 = zeros(ns)
        P = zeros(ns, ns)
        T = zeros(ns, ns)
        new(a0, dynamicws, H, lagged_state_ids, P, Q, R, state_ids, stoch_simul_options, T, varobs_ids, Y, Z)
    end
end

function loglikelihood(params::Vector{Float64}, varobs, observations::Matrix{Float64}, context::Dynare.Context, ssws::SSWs)
    model = context.models[1]
    results = context.results.model_results[1]

    #compute steady state and first order solution
    Dynare.compute_stoch_simul!(context, ssws.dynamicws, params, ssws.stoch_simul_options)

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
        ssws.Z[i, ssws.varobs_ids[i]] = 1.0
    end
    vg1 = view(
        results.linearrationalexpectations.g1_1,
        ssws.state_ids,
        :,
    )
    view(ssws.T, :, ssws.lagged_state_ids) .= vg1
    vg2 = view(
        results.linearrationalexpectations.g1_2,
        ssws.state_ids,
        :,
    )
    ssws.R .= vg2
    ssws.Q .= model.Sigma_e
    fill!(ssws.a0, 0.0)
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
        data_pattern = Vector{Vector{Int64}}(undef, 0)
        for i = 1:nobs
            push!(data_pattern, findall(.!ismissing.(view(Y[:, i]))))
        end
        Dynare.kalman_likelihood(Y, swws.Z, ssws.H, ssws.T, ssws.R, ssws.Q, ssws.a0, swws.P, start,
                                 last, presample, klws, data_pattern)
    else
        Dynare.kalman_likelihood(Y, ssws.Z, ssws.H, ssws.T, ssws.R, ssws.Q, ssws.a0, ssws.P, start,
                                 last, presample, klws)
    end
end

# generate artificial data with model example5
# provide model parsing
context = @dynare "models/example5/example5.mod";

# simulation is an accessor to simulated series
# we assume that we observe output, y, and consumption, c
varobs = ("y", "c")
observations = copy(transpose(Matrix(Dynare.simulation(varobs))))
nobs = size(observations, 2)
ssws = SSWs(context, nobs, varobs)
loglikelihood(params) = loglikelihood(params, varobs, observations, context, ssws)

# objective function
# context.work.params contains the parameter values of the DGP
loglikelihood(context.work.params)
