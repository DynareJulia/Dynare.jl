struct StochSimulOptions
    display::Bool
    dr_algo::String
    first_period::Int64
    irf::Int64
    LRE_options::LinearRationalExpectationsOptions
    nar::Int64
    nonstationary::Bool
    order::Int64
    periods::Int64
end

function StochSimulOptions(options::Dict{String,Any})
    display = true
    dr_algo = "GS"
    first_period = 1
    irf = 40
    LRE_options = LinearRationalExpectationsOptions()
    nar = 5
    nonstationary = false
    order = 1
    periods = 0
    print_results = true
    for (k, v) in pairs(options)
        if k == "noprint"
            display = false
        elseif k == "dr_cycle_reduction" && v::Bool
            dr_algo = "CR"
        elseif k == "first_period"
            first_period = v::Int64
        elseif k == "irf"
            irf = v::Int64
        elseif k == "nar"
            nar = v::Int64
        elseif k == "nonstationary"
            nonstationary = v::Bool
        elseif k == "order"
            order = v::Int64
        elseif k == "periods"
            periods = v::Int64
        end
    end
    StochSimulOptions(display, dr_algo, first_period, irf, LRE_options, nar, 
        nonstationary, order, periods)
end

function display_stoch_simul(context::Context, options::StochSimulOptions)
    m = context.models[1]
    endogenous_names = get_endogenous_longname(context.symboltable)
    exogenous_names = get_exogenous_longname(context.symboltable)
    results = context.results.model_results[1]
    LRE_results = results.linearrationalexpectations
    stationary_variables = LRE_results.stationary_variables
    display_solution_function(LRE_results.g1, endogenous_names, exogenous_names, m)
    display_mean_sd_variance(
        results.trends.endogenous_steady_state,
        diag(LRE_results.endogenous_variance),
        endogenous_names,
        m,
    )
    display_variance_decomposition(
        LRE_results,
        endogenous_names,
        exogenous_names,
        m,
        options,
    )
    display_correlation(LRE_results, endogenous_names, m)
    display_autocorrelation(LRE_results, endogenous_names, m, options)
end

function display_stoch_simul2(context::Context, options::StochSimulOptions)
    m = context.models[1]
    endogenous_names = get_endogenous_longname(context.symboltable)
    exogenous_names = get_exogenous_longname(context.symboltable)
    results = context.results.model_results[1]
    LRE_results = results.linearrationalexpectations
    stationary_variables = LRE_results.stationary_variables
    display_solution_function2(results.solution_derivatives, endogenous_names, exogenous_names, m)
    simulation_results = long_second_order_simulation(context)
    display_mean_sd_variance2(simulation_results, endogenous_names, m)
    display_correlation2(simulation_results, endogenous_names, m)
    display_autocorrelation2(simulation_results, endogenous_names, m, options)
end

function long_second_order_simulation(context; periods = 100_000, burning = 100)
    model = context.models[1]
    results = context.results.model_results[1]
    GD = results.solution_derivatives
    
    exo_nbr = model.exogenous_nbr
    original_endo_nbr = model.original_endogenous_nbr
    state_index = model.i_bkwrd_b

    gy1 = GD[1][:, 1]
    y0 = zeros(size(gy1, 1))

    active_exogenous = findall(diag(model.Sigma_e) .> 0)
    Sigma_e = view(model.Sigma_e, active_exogenous, active_exogenous)
    C = transpose(cholesky(Sigma_e).U)
    n_active = length(active_exogenous)

    u_shock = [zeros(exo_nbr) for _ in 1:periods]
    for i in 1:periods
        u_shock[i][active_exogenous] = C*randn(n_active) 
    end

    simWs = SimulateWs(GD, length(y0), state_index, model.exogenous_nbr)    
    simulation_result_vec = simulate(GD, y0, u_shock, periods, simWs)
    simulation_result_mat = stack(simulation_result_vec)'
    
    # ignore burning periods and useless endogenous variables
    simulation_result_clean = simulation_result_mat[burning:end, 1:original_endo_nbr]
    return simulation_result_clean
end

function display_solution_function(
    g1::AbstractMatrix{Float64},
    endogenous_names::AbstractVector{String},
    exogenous_names::AbstractVector{String},
    m::Model,
)
    title = "Coefficients of approximate solution function (reduced form)"
    data = Matrix{Any}(undef, m.n_states + m.exogenous_nbr + 1, m.endogenous_nbr + 1)
    data[1, 1] = ""
    # row headers
    for i = 1:m.n_states
        data[i+1, 1] = "ϕ($(endogenous_names[m.i_bkwrd_b[i]]))"
    end
    offset = m.n_states + 1
    for i = 1:m.exogenous_nbr
        data[i+offset, 1] = "$(exogenous_names[i])_t"
    end
    # columns
    for j = 1:m.endogenous_nbr
        # header
        data[1, j+1] = "$(endogenous_names[j])_t"
        # data
        for i = 1:m.n_states+m.exogenous_nbr
            data[i+1, j+1] = g1[j, i]
        end
    end
    # Note: ϕ(x) = x_{t-1} - \bar x
    #    note = string("Note: ϕ(x) = x\U0209C\U0208B\U02081 - ", "\U00305", "x")
    note = string("Note: ϕ(x) = x_{t-1} - steady_state(x)")
    println("\n")
    dynare_table(data, title, note = note)
end

function display_solution_function2(
    G::Vector{Matrix{Float64}},
    endogenous_names::AbstractVector{String},
    exogenous_names::AbstractVector{String},
    m::Model,
)
    title = "Coefficients of approximate solution function (reduced form, order 2)"
    n1 = m.n_states + m.exogenous_nbr
    n11 = n1 + 1
    nrows = 2 + n1 + Int(n1*(n1 + 1)/2)
    data = Matrix{Any}(undef, nrows, m.endogenous_nbr + 1)
    data[1, 1] = ""
    data[2, 1] = "Δσ"
    row = 3
    # row headers
    for i = 1:m.n_states
        data[row, 1] = "ϕ($(endogenous_names[m.i_bkwrd_b[i]]))"
        row += 1
    end
    offset = m.n_states + 2
    for i = 1:m.exogenous_nbr
        data[row, 1] = "$(exogenous_names[i])_t"
        row += 1
    end
    offset += m.exogenous_nbr
    for i = 1:m.n_states
        for j = i:m.n_states
            data[row, 1] = "ϕ($(endogenous_names[m.i_bkwrd_b[i]]))*ϕ($(endogenous_names[m.i_bkwrd_b[j]]))"
            row += 1
        end
        for j = 1:m.exogenous_nbr
            data[row, 1] = "ϕ($(endogenous_names[m.i_bkwrd_b[i]]))*$(exogenous_names[j])_t"
            row += 1
        end
    end
    offset += m.n_states*n1
    for i = 1:m.exogenous_nbr
        for j = i:m.exogenous_nbr
            data[row, 1] = "$(exogenous_names[i])_t*$(exogenous_names[j])_t"
            row += 1
        end
    end
    
    # columns
    for j = 1:m.endogenous_nbr
        # header
        data[1, j + 1] = "$(endogenous_names[j])_t"
        data[2, j + 1] = 0.5*G[2][j, n11*n11]
        row = 3
        # data
        for i = 1:n1
            data[row, j + 1] = G[1][j, i]
            row += 1
        end
        for i = 1:m.n_states
            for r = i:n1
                k = (i-1)*n11 + r
                    s = i == r ? 0.5 : 1  
                data[row, j + 1] = s*G[2][j, k]
                row += 1
            end
        end
        for i = 1:m.exogenous_nbr
            for r = i:m.exogenous_nbr
                k = (m.n_states + i-1)*n11 + m.n_states + r
                    s = i == r ? 0.5 : 1  
                data[row, j + 1] = s*G[2][j, k]
                row += 1
            end
        end
        
    end

    # Note: ϕ(x) = x_{t-1} - \bar x
    #    note = string("Note: ϕ(x) = x\U0209C\U0208B\U02081 - ", "\U00305", "x")
    note = string("Note: ϕ(x) = x_{t-1} - steady_state(x)")
    println("\n")
    dynare_table(data, title, note = note)
end

robustsqrt(x) = sqrt(x + eps())

function display_mean_sd_variance(
    steadystate::AbstractVector{Float64},
    variance::AbstractVector{Float64},
    endogenous_names::AbstractVector{String},
    m::Model,
)
    title = "THEORETICAL MOMENTS"

    std = robustsqrt.(variance)
    data = Matrix{Any}(undef, m.original_endogenous_nbr + 1, 4)
    data[1, 1] = "VARIABLE"
    # row headers
    for i = 1:m.original_endogenous_nbr
        data[i+1, 1] = "$(endogenous_names[i])"
    end
    # column headers
    data[1, 2] = "MEAN"
    data[1, 3] = "STD. DEV."
    data[1, 4] = "VARIANCE"
    # data
    for i = 1:m.original_endogenous_nbr
        data[i+1, 2] = isnan(std[i]) ? NaN : steadystate[i]
        data[i+1, 3] = std[i]
        data[i+1, 4] = variance[i]
    end
    println("\n")
    dynare_table(data, title)
end

function display_mean_sd_variance2(sim_results, endogenous_names, model::Model)
    original_endo_nbr = model.original_endogenous_nbr    
    
    title = "SIMULATED MOMENTS"
    data = Matrix{Any}(undef, original_endo_nbr + 1, 4)
    data[1, :] = ["VARIABLE", "MEAN", "STD. DEV.", "VARIANCE"]
    data[2:end, 1] = endogenous_names[1:original_endo_nbr]
    data[2:end, 2] .= map(mean, eachcol(sim_results))
    data[2:end, 3] .= map(std, eachcol(sim_results))
    data[2:end, 4] .= map(var, eachcol(sim_results))
    println("\n")
    dynare_table(data, title)
end

function display_variance_decomposition(
    LREresults::LinearRationalExpectationsResults,
    endogenous_names::AbstractVector{String},
    exogenous_names::AbstractVector{String},
    model::Model,
    options::StochSimulOptions,
)
    n = model.original_endogenous_nbr
    LREws = LinearRationalExpectations.LinearRationalExpectationsWs(
        options.dr_algo,
        model.exogenous_nbr,
        model.i_fwrd_b,
        model.i_current,
        model.i_bkwrd_b,
        model.i_static,
    )
    VD = variance_decomposition(
        LREresults,
        LREws,
        model.Sigma_e,
        model.endogenous_nbr,
        model.exogenous_nbr,
        model.n_bkwrd + model.n_both,
    )
    title = "VARIANCE DECOMPOSITION (in percent)"
    stationary_variables = LREresults.stationary_variables
    stationary_nbr = count(stationary_variables[1:model.original_endogenous_nbr])
    data = Matrix{Any}(undef, stationary_nbr + 1, model.exogenous_nbr + 1)
    data[1, 1] = "VARIABLE"
    # row headers
    k = 2
    for i = 1:model.original_endogenous_nbr
        if stationary_variables[i]
            data[k, 1] = "$(endogenous_names[i])"
            k += 1
        end
    end
    # columns
    for j = 1:model.exogenous_nbr
        # header
        data[1, j+1] = "$(exogenous_names[j])"
        # data
        k = 1
        for i = 1:model.original_endogenous_nbr
            if stationary_variables[i]
                data[k+1, j+1] = VD[i, j]
                k += 1
            end
        end
    end
    println("\n")
    dynare_table(data, title)
end

function display_correlation(
    LREresults::LinearRationalExpectationsResults,
    endogenous_names::Vector{String},
    model::Model,
)
    title = "CORRELATION MATRIX"
    corr = correlation(LREresults.endogenous_variance)
    n = model.original_endogenous_nbr
    stationary_variables = LREresults.stationary_variables
    stationary_nbr = count(stationary_variables[1:n])
    data = Matrix{Any}(undef, stationary_nbr + 1, stationary_nbr + 1)
    data[1, 1] = ""
    # row headers
    k = 2
    for i = 1:n
        if stationary_variables[i]
            data[k, 1] = "$(endogenous_names[i])"
            k += 1
        end
    end
    # columns
    k1 = 2
    for j = 1:n
        if stationary_variables[j]
            # header
            data[1, k1] = "$(endogenous_names[j])"
            # data
            k2 = 2
            for i = 1:n
                if stationary_variables[i]
                    data[k2, k1] = corr[i, j]
                    k2 += 1
                end
            end
            k1 += 1
        end
    end
    println("\n")
    dynare_table(data, title)
end

function display_correlation2(sim_results, endogenous_names, model::Model)
    original_endo_nbr = model.original_endogenous_nbr    

    title = "SIMULATED CORRELATION MATRIX"
    data = Matrix{Any}(undef, original_endo_nbr + 1, original_endo_nbr + 1)
    data[1, 1] = ""
    data[1, 2:end] = endogenous_names[1:original_endo_nbr] |> permutedims 
    data[2:end, 1] = endogenous_names[1:original_endo_nbr] 
    data[2:end, 2:end] = cor(sim_results)
    println("\n")
    dynare_table(data, title)
end

function display_autocorrelation(
    LREresults::LinearRationalExpectationsResults,
    endogenous_names::AbstractVector{String},
    model::Model,
    options::StochSimulOptions,
)
    title = "AUTOCORRELATION COEFFICIENTS"
    stationary_variables = LREresults.stationary_variables
    endogenous_nbr = model.endogenous_nbr
    n1 = count(stationary_variables)
    ar = [zeros(n1) for i = 1:options.nar]
    S1a = zeros(n1, endogenous_nbr)
    S1b = similar(S1a)
    S2 = zeros(endogenous_nbr - model.n_bkwrd - model.n_both, endogenous_nbr)
    ar = autocorrelation!(
        ar,
        LREresults,
        S1a,
        S1b,
        S2,
        model.i_bkwrd_b,
        stationary_variables,
    )
    n = model.original_endogenous_nbr
    stationary_nbr = count(stationary_variables[1:n])
    data = Matrix{Any}(undef, stationary_nbr + 1, options.nar + 1)
    data[1, 1] = ""
    # row headers
    k = 2
    for i = 1:n
        if stationary_variables[i]
            data[k, 1] = "$(endogenous_names[i])"
            k += 1
        end
    end
    # columns
    for j = 1:options.nar
        # header
        data[1, j+1] = j
        # data
        k = 1
        for i = 1:n
            if stationary_variables[i]
                data[k+1, j+1] = ar[j][k]
                k += 1
            end
        end
    end
    println("\n")
    dynare_table(data, title)
end

function display_autocorrelation2(sim_results, endogenous_names, model, options::StochSimulOptions)
    original_endo_nbr = model.original_endogenous_nbr    
    lags = [i for i in 1:options.nar]

    title = "SIMULATED AUTOCORRELATION COEFFICIENTS"
    data = Matrix{Any}(undef,  original_endo_nbr + 1, options.nar + 1)
    data[1, 1] = ""
    data[2:end, 1] = endogenous_names[1:original_endo_nbr]
    data[1, 2:end] = lags
    data[2:end, 2:end] = autocor(sim_results, lags)'
    println("\n")
    dynare_table(data, title)
end

function make_A_B!(
    A::Matrix{Float64},
    B::Matrix{Float64},
    model::Model,
    results::ModelResults,
)
    vA = view(A, :, model.i_bkwrd_b)
    vA .= results.linearrationalexpectations.g1_1
    B .= results.linearrationalexpectations.g1_2
end

function stoch_simul!(context::Context, field::Dict{String,Any})
    options = StochSimulOptions(field["options"])
    m = context.models[1]
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr = m.dynamic_tmp_nbr
    ws = DynamicWs(context, order=options.order)
    stoch_simul_core!(context, ws, options)
end

function localapproximation!(; context::Context=context,
                             display = true,
                             dr_algo = "GS",
                             first_period = 1,
                             irf = 40,
                             LRE_options = LinearRationalExpectationsOptions(),
                             nar = 5,
                             nonstationary = false,
                             order = 1,
                             periods = 0,
                             print_results = true
                             )
    options = StochSimulOptions(display, dr_algo, first_period, irf, LRE_options, nar, 
        nonstationary, order, periods)
    m = context.models[1]
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr = m.dynamic_tmp_nbr
    ws = DynamicWs(context, order=options.order)
    stoch_simul_core!(context, ws, options)
end

function stoch_simul_core!(context::Context, ws::DynamicWs, options::StochSimulOptions)
    model = context.models[1]
    modfileinfo = context.modfileinfo
    results = context.results.model_results[1]
    work = context.work
    #check_parameters(work.params, context.symboltable)
    #check_endogenous(results.trends.endogenous_steady_state)
    compute_stoch_simul!(context, ws, work.params, options; variance_decomposition = true)
    if options.display
        if options.order == 1
            display_stoch_simul(context, options)
        elseif options.order == 2
            display_stoch_simul2(context, options)
        end
    end
    if options.irf > 0
        path = "$(context.modfileinfo.modfilepath)/graphs/"
        mkpath(path)
        filename = "$(path)/irfs"
        
        if options.order == 1
            irfs!(context, options.irf)
        elseif options.order == 2
            irfs2!(context, options.irf)
        end

        plot_irfs(
            context.results.model_results[1].irfs,
            model,
            context.symboltable,
            filename,
        )
        display(context.results.model_results[1].irfs[:e])
    end


    if (periods = options.periods) > 0
        steadystate = results.trends.endogenous_steady_state
        linear_trend = results.trends.endogenous_linear_trend
        y0 = zeros(model.endogenous_nbr)
        simulresults = Matrix{Float64}(undef, periods + 1, model.endogenous_nbr)
        histval = work.histval
        if modfileinfo.has_histval
            for i in eachindex(skipmissing(view(work.histval, size(work.histval, 1), :)))
                y0[i] = work.histval[end, i]
            end
        else
            if work.model_has_trend[1]
                y0 = steadystate - linear_trend
            else
                y0 = steadystate
            end
        end
        C = cholesky(model.Sigma_e)
        x = vcat(zeros(1, model.exogenous_nbr), randn(periods, model.exogenous_nbr) * C.U)
        A = zeros(model.endogenous_nbr, model.endogenous_nbr)
        B = zeros(model.endogenous_nbr, model.exogenous_nbr)
        make_A_B!(A, B, model, results)
        simul_first_order!(simulresults, y0, x, steadystate, A, B, periods)
        if work.model_has_trend[1]
            simulresults .+= collect(0:periods) * transpose(linear_trend)
        end
        first_period = ExtendedDates.Undated(options.first_period)
        endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
        tdf = AxisArrayTable(simulresults, first_period .+ (0:options.periods), endogenous_names)
        push!(results.simulations, Simulation(first_period, first_period + periods - 1, "", "stoch_simul", tdf))
    end
end

function check!(context::Context, field::Dict{String,Any})
    Nothing
end

function compute_stoch_simul!(
    context::Context,
    ws::DynamicWs,
    params::Vector{Float64},
    options::StochSimulOptions;
    kwargs...
        )
    model = context.models[1]
    order = options.order
    results = context.results.model_results[1]
    # don't check the steady state if the model is nonstationary
    fill!(results.trends.exogenous_steady_state, 0.0)
    compute_steady_state!(context, nocheck = options.nonstationary)
    endogenous = results.trends.endogenous_steady_state
    endogenous3 = repeat(endogenous, 3)
    exogenous = results.trends.exogenous_steady_state
    compute_first_order_solution!(
        context,
        endogenous3,
        exogenous,
        endogenous,
        params,
        model,
        ws,
        options; kwargs...
            )
    if order == 2
        steadystate = context.results.model_results[1].trends.endogenous_steady_state
        values = zeros(size(model.dynamic_g2_sparse_indices, 1))
        get_dynamic_derivatives2!(
            ws,
            params,
            endogenous3,
            exogenous,
            steadystate,
            values,
            model,
        )
        m = model
        G = Vector{Matrix{Float64}}(undef, 0)
        push!(G, context.results.model_results[1].linearrationalexpectations.g1)
        push!(G, zeros(model.endogenous_nbr,
                       (model.n_states + model.exogenous_nbr + 1)^2))
        solverWs = KOrderPerturbations.KOrderWs(model.endogenous_nbr,
                                          model.n_fwrd + model.n_both,
                                          model.n_states,
                                          model.n_current,
                                          model.exogenous_nbr,
                                          model.i_fwrd_b,
                                          model.i_bkwrd_b,
                                          model.i_current,
                                          1:model.n_states,
                                          order)
        moments = [0, vec(model.Sigma_e)]
        KOrderPerturbations.k_order_solution!(G,
                                              ws.derivatives,
                                              moments,
                                              order,
                                              solverWs)
        copy!(results.solution_derivatives, G)
    end
end

function irfs2!(context, periods; burning = 100)
    model = context.models[1]
    results = context.results.model_results[1]
    GD = results.solution_derivatives
    
    exo_nbr = model.exogenous_nbr
    state_index = model.i_bkwrd_b
    
    model.Sigma_e
    active_exogenous = findall(diag(model.Sigma_e) .> 0)
    nx = length(active_exogenous)
    endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
    exogenous_names = [Symbol(n) for n in get_exogenous_longname(context.symboltable)]
    exogenous_names = view(exogenous_names, active_exogenous)
    Sigma_e = view(model.Sigma_e, active_exogenous, active_exogenous)
    C = transpose(cholesky(Sigma_e).U)
    gy1 = GD[1][:, 1]
    n = size(gy1)[1]
    y0 = zeros(n)

    simWs = SimulateWs(GD, n, state_index, model.exogenous_nbr)

    for (i, exo_var) in enumerate(active_exogenous) 
        u_shock = [zeros(exo_nbr) for _ in 1:periods + burning]
        det_shock = C[:, i]

        mean_sims = [zeros(model.endogenous_nbr) for _ in 1:periods]
        replic = 1000
        for n_sim in 1:replic
            for i in 1:periods + burning
                u_shock[i][active_exogenous] = C*randn(nx) 
            end
            rand_sim = simulate(GD, y0, u_shock, periods + burning, simWs)
            u_shock_det = deepcopy(u_shock)
            @views u_shock_det[burning + 1] .+= det_shock  
            rand_det_sim = simulate(GD, y0, u_shock_det, periods + burning, simWs)
            @views mean_sims .+= (rand_det_sim[burning + 1:end] .- rand_sim[burning + 1:end])/replic
        end
        
        mean_irf_matrix = reduce(vcat, transpose.(mean_sims))

        tdf = AxisArrayTable(mean_irf_matrix, Undated(1):Undated(periods), endogenous_names) 
        results.irfs[exogenous_names[exo_var]] = tdf
    end

end


function compute_first_order_solution!(
    context::Context,
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    model::Model,
    ws::DynamicWs,
    options::StochSimulOptions;
    variance_decomposition::Bool = true,
)

    results = context.results.model_results[1]
    LRE_results = results.linearrationalexpectations
    jacobian = get_dynamic_jacobian!(
        ws,
        params,
        endogenous,
        exogenous,
        steadystate,
        model,
        2,
    )
    lli = model.lead_lag_incidence
    @views jacobian = hcat(Matrix(jacobian[:, findall(lli[1, :] .> 0)]),
                           Matrix(jacobian[:, model.endogenous_nbr .+ findall(lli[2, :] .> 0)]),
                           Matrix(jacobian[:, 2*model.endogenous_nbr .+ findall(lli[3, :] .> 0)]),
                           Matrix(jacobian[:, 3*model.endogenous_nbr .+ collect(1:model.exogenous_nbr)]))
    #    LRE.remove_static!(jacobian, wsLRE)
    LRE.first_order_solver!(LRE_results, jacobian, options.LRE_options,
                            workspace(LRE.LinearRationalExpectationsWs, context, algo=options.dr_algo))
    if variance_decomposition
        compute_variance!(LRE_results, model.Sigma_e, workspace(LRE.VarianceWs, context, algo=options.dr_algo))
    end
end

function irfs!(context, periods)
    model = context.models[1]
    model.Sigma_e
    results = context.results.model_results[1]
    active_exogenous = findall(diag(model.Sigma_e) .> 0)
    endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
    exogenous_names = [Symbol(n) for n in get_exogenous_longname(context.symboltable)]
    exogenous_names = view(exogenous_names, active_exogenous)
    Sigma_e = view(model.Sigma_e, active_exogenous, active_exogenous)
    C = transpose(cholesky(model.Sigma_e).U)
    x = zeros(length(active_exogenous))
    A = zeros(model.endogenous_nbr, model.endogenous_nbr)
    B = zeros(model.endogenous_nbr, model.exogenous_nbr)
    make_A_B!(A, B, model, results)
    for (i, j) in enumerate(active_exogenous)
        fill!(x, 0)
        @views x .= C[:, j]
        yy = Matrix{Float64}(undef, size(A, 1), periods)
        @views mul!(view(yy, :, 1), B[:, active_exogenous], x)
        for j = 2:periods
            mul!(view(yy, :, j), A, view(yy, :, j - 1))
        end
        tdf = AxisArrayTable(transpose(yy), Undated(1):Undated(periods), endogenous_names) 
        results.irfs[exogenous_names[i]] = tdf
    end
end
