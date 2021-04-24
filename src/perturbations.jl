function display_stoch_simul(x, title, context)
    endogenous_names = get_endogenous_longname(context.symboltable)
    emptyrow = ["" for _= 1:size(x,1)]
    column_header = []
    #    map(x -> push!(column_header, string(x, "\U0209C")), endogenous_names)
    map(x -> push!(column_header, "$(x)_t"), endogenous_names)
    row_header = [""]
    map(x -> push!(row_header, "ϕ($x)"), endogenous_names[context.models[1].i_bkwrd_b])
    map(x -> push!(row_header, "$(x)_t"), get_exogenous_longname(context.symboltable))
    data = hcat(row_header,
                vcat(reshape(column_header, 1, length(column_header)),
                     x))
    # Note: ϕ(x) = x_{t-1} - \bar x
    #    note = string("Note: ϕ(x) = x\U0209C\U0208B\U02081 - ", "\U00305", "x")
    note = string("Note: ϕ(x) = x_{t-1} - steady_state(x)")
    println("\n")
    dynare_table(data, title, column_header, row_header, note)
end

function make_A_B!(A, B, model, results)
    vA = view(A, :, model.i_bkwrd_b)
    vA .= results.linearrationalexpectations.g1_1
    B .= results.linearrationalexpectations.g1_2
end

function stoch_simul!(context, field)
    model = context.models[1]
    options = context.options
    results = context.results.model_results[1]
    work = context.work
    options["stoch_simul"] = Dict()
    copy!(options["stoch_simul"], field["options"])
    #check_parameters(work.params, context.symboltable)
    #check_endogenous(results.trends.endogenous_steady_state)
    compute_stoch_simul!(context)
    compute_variance!(context)
    x = results.linearrationalexpectations.g1
    vx = view(x, :, 1:size(x, 2) - 1)
    steadystate = results.trends.endogenous_steady_state
    linear_trend = results.trends.endogenous_linear_trend
    y0 = copy(steadystate)
    display_stoch_simul(vx', "Coefficients of approximate solution function", context)
    if (periods = get(options["stoch_simul"], "periods", 0)) > 0
        simulresults = Matrix{Float64}(undef, periods + 1, model.endogenous_nbr)
        histval = work.histval
        if size(histval, 1) == 0
            # no histval
            if work.model_has_trend
                y0 = steadystate - linear_trend
            else
                y0 = steadystate
            end
        else
            # histval present
            y0 = view(histval, 1, :)
        end
        C = cholesky(model.Sigma_e)
        x = vcat(zeros(1, model.exogenous_nbr),
                 randn(periods, model.exogenous_nbr)*C.U)
        A = zeros(model.endogenous_nbr, model.endogenous_nbr)
        B = zeros(model.endogenous_nbr, model.exogenous_nbr)
        make_A_B!(A, B, model, results)
        simul_first_order!(simulresults, y0, x, steadystate, A, B, periods)
        if work.model_has_trend
            simulresults .+= collect(0:periods) * transpose(linear_trend) 
        end
        first_period = get(options["stoch_simul"], "first_periods", 1)
        endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
        data = TimeDataFrame(simulresults, Periods.Undated, first_period, endogenous_names)
        push!(results.simulations, Simulation("", "stoch_simul", options["stoch_simul"], data))
    end
end

function check!(context, field)
end

function compute_stoch_simul!(context)
    model = context.models[1]
    results = context.results.model_results[1]
    options = context.options["stoch_simul"]
    work = context.work
    compute_steady_state!(context)
    endogenous = results.trends.endogenous_steady_state
    exogenous = results.trends.exogenous_steady_state
    fill!(exogenous, 0.0)
    compute_first_order_solution!(results.linearrationalexpectations,
                                  endogenous, exogenous, endogenous,
                                  model, work, options)
end

function compute_first_order_solution!(
    results::LinearRationalExpectationsResults,
    endogenous::AbstractVector{Float64},
    exogenous::AbstractVector{Float64},
    steadystate::AbstractVector{Float64},
    model::Model, work::Work, options)

    # abbreviations
    LRE = LinearRationalExpectations
    LREWs = LinearRationalExpectationsWs

    options["cyclic_reduction"] = Dict()
    options["generalized_schur"] = Dict()

    get_jacobian!(work, endogenous, exogenous, steadystate,
                  model, 2)
    if isnothing(get(options,"dr_cycle_reduction", nothing))
        algo = "GS"
    else
        algo = "CR"
    end
    ws = LREWs(algo,
               model.endogenous_nbr,
               model.exogenous_nbr,
               model.exogenous_deterministic_nbr,
               model.i_fwrd_b,
               model.i_current,
               model.i_bkwrd_b,
               model.i_both,
               model.i_static)
    LRE.remove_static!(work.jacobian, ws)
    if algo == "GS"
        LRE.get_de!(ws, work.jacobian)
        options["generalized_schur"]["criterium"] = 1 + 1e-6
    else
        LRE.get_abc!(ws, work.jacobian)
        options["cyclic_reduction"]["tol"] = 1e-8
    end
    LRE.first_order_solver!(results,
                            algo,
                            work.jacobian,
                            options,
                            ws)
end

function compute_variance!(context)
    m = context.models[1]
    results = context.results.model_results[1]
    Σy = results.endogenous_variance
    A = results.linearrationalexpectations.gs1
    B1 = zeros(m.n_states, m.exogenous_nbr)
    g1_1 = results.linearrationalexpectations.g1_1
    g1_2 = results.linearrationalexpectations.g1_2
    vr1 = view(g1_2, m.i_bkwrd_b, :)
    B1 .= vr1
    ws = LyapdWs(m.n_states)
    Σ = zeros(m.n_states, m.n_states)
    tmp = zeros(m.n_states, m.exogenous_nbr)
    mul!(tmp, B1, m.Sigma_e)
    B = zeros(m.n_states, m.n_states)
    mul!(B, tmp, B1')
    extended_lyapd!(Σ, A, B, ws)
    stationary_variables = results.stationary_variables
    fill!(stationary_variables, true)
    state_stationary_variables = view(stationary_variables, m.i_bkwrd_b)
    nonstate_stationary_variables = view(stationary_variables, m.i_non_states)
    if ws.stationary_model
        state_stationary_nbr = m.n_states 
        nonstate_stationary_nbr = m.endogenous_nbr - m.n_states
    else
        fill!(Σy, NaN)
        state_stationary_variables .= .!ws.nonstationary_variables
        state_stationary_nbr = count(state_stationary_variables)
        vr3 = view(results.linearrationalexpectations.g1_1, m.i_non_states, :)
        for i = 1:(m.endogenous_nbr - m.n_states)
            for j = 1:m.n_states
                if ws.nonstationary_variables[j] && abs(vr3[i, j]) > 1e-10
                    nonstate_stationary_variables[j] = false
                    break
                end
            end
        end
        nonstate_stationary_nbr = count(nonstate_stationary_variables)
    end
    # state / state
    stationary_nbr = state_stationary_nbr + nonstate_stationary_nbr
    A2 = zeros(nonstate_stationary_nbr, state_stationary_nbr)
    vtmp = view(results.linearrationalexpectations.g1_1, m.i_non_states, :)
    A2 .= view(vtmp, nonstate_stationary_variables, state_stationary_variables)
    vtmp = view(Σy, m.i_bkwrd_b, m.i_bkwrd_b)
    vΣy = view(vtmp, state_stationary_variables, state_stationary_variables)
    vΣ = view(Σ, state_stationary_variables, state_stationary_variables) 
    vΣy .= vΣ
    # endogenous / nonstate
    n_non_states = m.endogenous_nbr - m.n_states
    B2 = zeros(nonstate_stationary_nbr, m.exogenous_nbr)
    vtmp = view(results.linearrationalexpectations.g1_2, m.i_non_states, :)
    vr2 = view(vtmp, nonstate_stationary_variables, :)
    B2 .= vr2
    Σ = zeros(stationary_nbr, nonstate_stationary_nbr)
    vg1 = view(results.linearrationalexpectations.g1_1, stationary_variables, state_stationary_variables)
    tmp1 = zeros(stationary_nbr, state_stationary_nbr)
    mul!(tmp1, vg1, vΣy)
    vg1 = view(results.linearrationalexpectations.g1_1, m.i_non_states, :)
    vg11 = view(vg1, nonstate_stationary_variables, state_stationary_variables)
    mul!(Σ, tmp1, transpose(vg11))
    tmp2 = zeros(stationary_nbr, m.exogenous_nbr)
    vg2 = view(results.linearrationalexpectations.g1_2, stationary_variables, :)
    mul!(tmp2, vg2, m.Sigma_e)
    mul!(Σ, tmp2, transpose(B2), 1.0, 1.0)
    vtmp = view(Σy, :, m.i_non_states) 
    vΣy = view(vtmp, stationary_variables, nonstate_stationary_variables)
    vΣy .= Σ
    # nonstate / state
    vtmp1 = view(Σy, m.i_non_states, m.i_bkwrd_b)
    vΣy1 = view(vtmp1, nonstate_stationary_variables, state_stationary_variables)
    vtmp2 = view(Σy, m.i_bkwrd_b, m.i_non_states)
    vΣy2 = view(vtmp2, state_stationary_variables, nonstate_stationary_variables)
    vΣy1 .= transpose(vΣy2)
end
    
function simul_first_order!(results, initial_values, x, c, A, B, periods)
    r_1 = view(results, 1, :)
    r_1 .= initial_values .- c
    for t = 2:periods + 1
        r = view(results, t, :)
        e = view(x, t, :)
        mul!(r, B, e)
        mul!(r, A, r_1, 1.0, 1.0)
        r_1 .+=  c
        r_1 = r
    end
    r_1 .+= c
end
