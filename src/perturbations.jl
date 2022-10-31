struct StochSimulOptions
    display::Bool
    dr_algo::String
    first_period::Int64
    irf::Int64
    LRE_options::LinearRationalExpectationsOptions
    nar::Int64
    order::Int64
    periods::Int64
    function StochSimulOptions(options::Dict{String,Any})
        display = true
        dr_algo = "GS"
        first_period = 1
        irf = 40
        LRE_options = LinearRationalExpectationsOptions()
        nar = 5
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
            elseif k == "order"
                order = v::Int64
            elseif k == "periods"
                periods = v::Int64
            end
        end
        new(display, dr_algo, first_period, irf, LRE_options, nar, order, periods)
    end
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
        data[i+1, 2] = steadystate[i]
        data[i+1, 3] = std[i]
        data[i+1, 4] = variance[i]
    end
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

function display_autocorrelation(
    LREresults::LinearRationalExpectationsResults,
    endogenous_names::AbstractVector{String},
    model::Model,
    options::StochSimulOptions,
)
    title = "AUTOCORRELATION COEFFICIENTS"
    ar = [zeros(model.endogenous_nbr) for i = 1:options.nar]
    stationary_variables = LREresults.stationary_variables
    # doesn't work yet for nonstationary models
    if model.endogenous_nbr != count(stationary_variables)
        return
    end
    n1 = count(stationary_variables)
    S1a = zeros(n1, model.endogenous_nbr)
    S1b = similar(S1a)
    S2 = zeros(model.endogenous_nbr - model.n_bkwrd - model.n_both, model.endogenous_nbr)
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
                data[k+1, j+1] = ar[j][i]
                k += 1
            end
        end
    end
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
    ws = DynamicWs(context)
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
        display_stoch_simul(context, options)
    end
    if options.irf > 0
        exogenous_names = get_exogenous_longname(context.symboltable)
        n = model.endogenous_nbr
        m = model.exogenous_nbr
        y = Dict{Symbol,TimeDataFrame}
        irfs!(context, options.irf)
        path = "$(context.modfileinfo.modfilepath)/graphs/"
        mkpath(path)
        filename = "$(path)/irfs"
        plot_irfs(
            context.results.model_results[1].irfs,
            model,
            context.symboltable,
            filename,
        )
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
        first_period = ExtendedDates.UndatedDate(options.first_period)
        endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
        tdf = TimeDataFrame(simulresults, first_period, endogenous_names)
        push!(results.simulations, Simulation("", "stoch_simul", tdf))
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
    variance_decomposition::Bool = false,
)
    model = context.models[1]
    results = context.results.model_results[1]
    compute_steady_state!(context, Dict{String, Any}())
    endogenous = results.trends.endogenous_steady_state
    endogenous3 = repeat(endogenous, 3)
    exogenous = results.trends.exogenous_steady_state
    fill!(exogenous, 0.0)
    compute_first_order_solution!(
        context,
        endogenous3,
        exogenous,
        endogenous,
        params,
        model,
        ws,
        options,
        variance_decomposition,
    )
end

function compute_first_order_solution!(
    context::Context,
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    params::AbstractVector{Float64},
    model::Model,
    ws::DynamicWs,
    options::StochSimulOptions,
    variance_decomposition::Bool,
)

    # abbreviations
    LRE = LinearRationalExpectations
    LREWs = LinearRationalExpectationsWs

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
    algo = options.dr_algo
    wsLRE = LinearRationalExpectationsWs(algo,
                                         model.exogenous_nbr,
                                         model.i_fwrd_b,
                                         model.i_current,
                                         model.i_bkwrd_b,
                                         model.i_static,
                                         )
    lli = model.lead_lag_incidence
    @views J = hcat(Matrix(jacobian[:, findall(lli[1, :] .> 0)]),
                    Matrix(jacobian[:, model.endogenous_nbr .+ findall(lli[2, :] .> 0)]),
                    Matrix(jacobian[:, 2*model.endogenous_nbr .+ findall(lli[3, :] .> 0)]),
                    Matrix(jacobian[:, 3*model.endogenous_nbr .+ collect(1:model.exogenous_nbr)]))
    LRE.remove_static!(J, wsLRE)
    LRE.first_order_solver!(LRE_results, J, options.LRE_options, wsLRE)
    lre_variance_ws = LRE.VarianceWs(
        model.endogenous_nbr,
        model.n_bkwrd + model.n_both,
        model.exogenous_nbr,
        wsLRE,
    )
    compute_variance!(LRE_results, model.Sigma_e, lre_variance_ws)
end


function irfs!(context, periods)
    model = context.models[1]
    results = context.results.model_results[1]
    endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
    exogenous_names = [Symbol(n) for n in get_exogenous_longname(context.symboltable)]
    C = cholesky(model.Sigma_e + 1e-14 * I)
    x = zeros(model.exogenous_nbr)
    A = zeros(model.endogenous_nbr, model.endogenous_nbr)
    B = zeros(model.endogenous_nbr, model.exogenous_nbr)
    make_A_B!(A, B, model, results)
    for i = 1:model.exogenous_nbr
        fill!(x, 0)
        x[i] = robustsqrt(model.Sigma_e[i, i])
        yy = Matrix{Float64}(undef, size(A, 1), periods)
        mul!(view(yy, :, 1), B, x)
        for j = 2:periods
            mul!(view(yy, :, j), A, view(yy, :, j - 1))
        end
        tdf = TimeDataFrame(transpose(yy), UndatedDate(1), endogenous_names)
        results.irfs[exogenous_names[i]] = tdf
    end
end
