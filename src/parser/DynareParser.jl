using CSV
using DataFrames
using GR
using JSON
using KalmanFilterTools
using LinearRationalExpectations
using Periods
using Plots
using TimeDataFrames

struct Simulation
    name::String
    statement::String
    options::Dict{String, Any}
    data::TimeDataFrame
end

struct ModelResults
    endogenous_steady_state::Vector{Float64}
    endogenous_trend::Matrix{Float64}
    endogenous_variance::Matrix{Float64}
    exogenous_steady_state::Vector{Float64}
    exogenous_deterministic_steady_state::Vector{Float64}
    linearrationalexpectations::LinearRationalExpectationsResults
    simulations::Vector{Simulation}
    smoother::Dict{String, Any}
end

struct Results
    model_results::Vector{ModelResults}
end

mutable struct Work
    params::Vector{Float64}
    residuals::Vector{Float64}
    temporary_values::Vector{Float64}
    dynamic_variables::Vector{Float64}
    jacobian::Matrix{Float64}
    qr_jacobian::Matrix{Float64}
    model_has_trend::Bool
    histval::Matrix{Float64}
end

mutable struct Context
    symboltable::Dict{String, DynareSymbol}
    models::Vector{Model}
    options::Dict
    results::Results
    work::Work
end

struct ModelInfo
    lead_lag_incidence::Array{Int64}
    nstatic
    nfwrd
    npred
    nboth
    nsfwrd
    nspred
    ndynamic
    maximum_endo_lag
    maximum_endo_lead
    maximum_exo_lag
    maximum_exo_lead
    maximum_exo_det_lag
    maximum_exo_det_lead
    maximum_lag
    maximum_lead
    orig_maximum_endo_lag
    orig_maximum_endo_lead
    orig_maximum_exo_lag
    orig_maximum_exo_lead
    orig_maximum_exo_det_lag
    orig_maximum_exo_det_lead
    orig_maximum_lag
    orig_maximum_lead
end

context = Context(Dict{String, DynareSymbol}(),
                  Vector{Model}(undef, 0),
                  Dict(),
                  Results(Vector{ModelResults}(undef, 0)),
                  Work(Vector{Float64}(undef, 0),
                       Vector{Float64}(undef, 0),
                       Vector{Float64}(undef, 0),
                       Vector{Float64}(undef, 0),
                       Matrix{Float64}(undef, 0, 0),
                       Matrix{Float64}(undef, 0, 0),
                       false,
                       Matrix{Float64}(undef, 0, 0)
                       )
                  )
                         
function parser(modfilename)
    modelstring = open(f -> read(f, String), modfilename*"/model/json/modfile.json")
    modeljson = JSON.parse(modelstring)

    symboltable = SymbolTable()
    endo_nbr = set_symbol_table!(symboltable, modeljson["endogenous"], Endogenous)
    exo_nbr = set_symbol_table!(symboltable, modeljson["exogenous"], Exogenous)
    exo_det_nbr = set_symbol_table!(symboltable,
                                    modeljson["exogenous_deterministic"],
                                    ExogenousDeterministic)
    param_nbr = set_symbol_table!(symboltable, modeljson["parameters"], Parameter)
    model_info = get_model_info(modeljson["model_info"])
    model = Model(modfilename,
                  endo_nbr,
                  model_info.lead_lag_incidence,
                  exo_nbr,
                  0,
                  exo_det_nbr,
                  param_nbr,
                  model_info.maximum_endo_lag,
                  model_info.maximum_endo_lead,
                  model_info.maximum_exo_lag,
                  model_info.maximum_exo_lead,
                  model_info.maximum_exo_det_lag,
                  model_info.maximum_exo_det_lead,
                  model_info.maximum_lag,
                  model_info.maximum_lead,
                  model_info.orig_maximum_endo_lag,
                  model_info.orig_maximum_endo_lead,
                  model_info.orig_maximum_exo_lag,
                  model_info.orig_maximum_exo_lead,
                  model_info.orig_maximum_exo_det_lag,
                  model_info.orig_maximum_exo_det_lead,
                  model_info.orig_maximum_lag,
                  model_info.orig_maximum_lead
                  )

    if "varobs" in keys(modeljson)
        varobs = modeljson["varobs"]
        varobs_ids = modeljson["varobs_ids"]
    else
        varobs = []
        varobs_ids = []
    end

    if "varexobs" in keys(modeljson)
        varexobs = modeljson["varexobs"]
        varexobs_ids = modeljson["varexobs_ids"]
    else
        varexobs = []
        varexobs_ids = []
    end

    order = 1
    modelresults = ModelResults(Vector{Float64}(undef, endo_nbr),
                                Matrix{Float64}(undef, endo_nbr, 2),
                                Matrix{Float64}(undef, endo_nbr, endo_nbr),
                                Vector{Float64}(undef, exo_nbr),
                                Vector{Float64}(undef, exo_det_nbr),
                                LinearRationalExpectationsResults(order,
                                                                  endo_nbr,
                                                                  exo_nbr,
                                                                  model.n_states),
                                Vector{Simulation}(undef, 0),
                                Dict{String, Any}())
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2*model.n_both
    ncol1 = ncol + model.exogenous_nbr
    work = Work(Vector{Float64}(undef, model.parameter_nbr),
                Vector{Float64}(undef, model.endogenous_nbr),
                Vector{Float64}(undef, sum(model.dynamic!.tmp_nbr[1:2])),
                Vector{Float64}(undef, ncol),
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                Matrix{Float64}(undef, model.endogenous_nbr, ncol1),
                false,
                Matrix{Float64}(undef, 0, 0))
    results = Results([modelresults])
    global context = Context(symboltable, [model], Dict(), results, work)
    for field in modeljson["statements"]
        if field["statementName"] == "calib_smoother"
            calib_smoother!(context, field, varobs, varobs_ids)
        elseif field["statementName"] == "check"
            check(field)
        elseif field["statementName"] == "deterministic_trends"
            deterministic_trends!(context, field)
        elseif field["statementName"] == "histval"
            histval!(context, field, symboltable)
        elseif field["statementName"] == "initval"
            initval!(field)
        elseif field["statementName"] == "native"
            @show field["string"]
            expr = Meta.parse(field["string"])
            eval(expr)
            native_statement(field)
        elseif field["statementName"] == "param_init"
            initialize_parameter!(work.params, field, symboltable)
        elseif field["statementName"] == "perfect_foresight_setup"
            perfect_foresight_setup!(options, field)
        elseif field["statementName"] == "perfect_foresight_solver"
            perfect_foresight_solver!(context, field)
        elseif field["statementName"] == "shocks"
            shocks!(model.Sigma_e, field, symboltable)
        elseif field["statementName"] == "stoch_simul"
            stoch_simul!(context, field)
        elseif field["statementName"] == "verbatim"
            verbatim(field)
        else
            error("Unrecognized statement $(field["statementName"])")
        end
    end
    return context
end

function set_symbol_table!(table::Dict{String, DynareSymbol},
                     modelfile,
                     type::SymbolType)
    count = 0
    for entry in modelfile
        count += 1
        symbol = DynareSymbol(entry["longName"],
                        entry["texName"],
                        type,
                        count)
        table[entry["name"]] = symbol

    end
    return count
end

get_model_info(field) =
    ModelInfo(hcat(field["lead_lag_incidence"]...),
              field["nstatic"],
              field["nfwrd"],
              field["npred"],
              field["nboth"],
              field["nsfwrd"],
              field["nspred"],
              field["ndynamic"],
              field["maximum_endo_lag"],
              field["maximum_endo_lead"],
              field["maximum_exo_lag"],
              field["maximum_exo_lead"],
              field["maximum_exo_det_lag"],
              field["maximum_exo_det_lead"],
              field["maximum_lag"],
              field["maximum_lead"],
              field["orig_maximum_endo_lag"],
              field["orig_maximum_endo_lead"],
              field["orig_maximum_exo_lag"],
              field["orig_maximum_exo_lead"],
              field["orig_maximum_exo_det_lag"],
              field["orig_maximum_exo_det_lead"],
              max(field["orig_maximum_lag"],
                  field["orig_maximum_lag_with_diffs_expanded"]),
              field["orig_maximum_lead"]
              )
                                  
            
function initialize_parameter!(params, field, symboltable)
    s = symboltable[field["name"]]
    k = s.orderintype
    params[k] = eval(Meta.parse(field["value"]))
end

function native_statement(field)
    eval
#    println("NATIVE: $field")
end

function histval!(context, field, symboltable)
    m = context.models[1]
    histval = zeros(m.orig_maximum_lag, m.endogenous_nbr)
    for v in field["vals"]
        k = symboltable[v["name"]].orderintype
        l = m.orig_maximum_lag - v["lag"]
        histval[l, k] = Meta.parse(v["value"])
    end
    @show histval
    context.work.histval = histval
    @show context.work.histval
end

function initval!(field)
end

function shocks!(Sigma, field, symboltable)
    set_variance!(Sigma, field["variance"], symboltable)
    set_stderr!(Sigma, field["stderr"], symboltable)
    set_covariance!(Sigma, field["covariance"], symboltable)
    set_correlation!(Sigma, field["correlation"], symboltable)
end

function set_variance!(Sigma, variance, symboltable)
    for v in variance
        k =  symboltable[v["name"]].orderintype
        Sigma[k, k] = eval(Meta.parse(v["variance"]))
    end
end

function set_stderr!(Sigma, stderr, symboltable)
    for s in stderr
        k =  symboltable[s["name"]].orderintype
        x = eval(Meta.parse(s["stderr"]))
        Sigma[k, k] = x*x
    end
end

function set_covariance!(Sigma, covariance, symboltable)
    for c in covariance
        k1 =  symboltable[c["name"]].orderintype
        k2 =  symboltable[c["name2"]].orderintype
        Sigma[k1, k2] = eval(Meta.parse(c["covariance"]))
        Sigma[k2, k1] = Sigma[k1, k2]
    end
end

function set_correlation!(Sigma, correlation, symboltable)
    for c in correlation
        k1 =  symboltable[c["name"]].orderintype
        k2 =  symboltable[c["name2"]].orderintype
        corr = eval(Meta.parse(c["correlation"]))
        Sigma[k2, k1] = sqrt(Sigma[k1, k1]*Sigma[k2, k2])*corr
    end
end

function display_stoch_simul(x, title, context)
    endogenous_names = get_endogenous_longname(context.symboltable)
    emptyrow = ["" for _= 1:size(x,1)]
    column_header = []
    map(x -> push!(column_header, string(x, "\U0209C")), endogenous_names)
    row_header = [""]
    map(x -> push!(row_header, "ϕ($x)"), endogenous_names[context.models[1].i_bkwrd_b])
    map(x -> push!(row_header, "$x\U0209C"), get_exogenous_longname(context.symboltable))
    data = hcat(row_header,
                vcat(reshape(column_header, 1, length(column_header)),
                     x))
    # Note: ϕ(x) = x_{t-1} - \bar x
    note = string("Note: ϕ(x) = x\U0209C\U0208B\U02081 - ", "\U00305", "x")
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
    compute_stoch_simul!(context)
    compute_variance!(context)
    x = results.linearrationalexpectations.g1
    vx = view(x, :, 1:size(x, 2) - 1)
    display_stoch_simul(vx', "Coefficients of approximate solution function", context)
    if (periods = get(options["stoch_simul"], "periods", 0)) > 0
        simulresults = Matrix{Float64}(undef, periods + 1, model.endogenous_nbr)
        histval = work.histval
        if size(histval, 1) == 0
            y0 = results.endogenous_steady_state
        else
            y0 = view(histval, 1, :)
        end
        C = cholesky(model.Sigma_e)
        x = vcat(zeros(1, model.exogenous_nbr),
                 randn(periods, model.exogenous_nbr)*C.U)
        c = results.endogenous_steady_state
        A = zeros(model.endogenous_nbr, model.endogenous_nbr)
        B = zeros(model.endogenous_nbr, model.exogenous_nbr)
        make_A_B!(A, B, model, results)
        simul_first_order!(simulresults, y0, x, c, A, B, periods)
        if work.model_has_trend
            simulresults .+= transpose(results.endogenous_trend[:,1]) .+ collect(0:size(simulresults, 1)-1) * transpose(results.endogenous_trend[:, 2]) 
        end
        first_period = get(options["stoch_simul"], "first_periods", 1)
        endogenous_names = [Symbol(n) for n in get_endogenous_longname(context.symboltable)]
        data = TimeDataFrame(simulresults, Periods.Undated, first_period, endogenous_names)
        push!(results.simulations, Simulation("", "stoch_simul", options["stoch_simul"], data))
    end
end

function perfect_foresight_setup(context, field)
    context.options["perfect_foresight_setup"] = Dict()
    copy!(context.options["perfect_foresight_setup"], field["options"])
    compute_perfect_foresight_setup(context)
end

function perfect_foresight_solver(context, field)
    context.options["perfect_foresight_solver"] = Dict()
    copy!(context.options["perfect_foresight_solver"], field["options"])
    compute_perfect_foresight_solver(context)
end


function verbatim(field)
#    println("VERBATIM: $field")
end

function check(field)
end

function compute_stoch_simul!(context)
    m = context.models[1]
    results = context.results.model_results[1]
    options = context.options["stoch_simul"]
    options["cyclic_reduction"] = Dict()
    options["generalized_schur"] = Dict()
    work = context.work
    Base.invokelatest(steady_state!, context)
    fill!(results.exogenous_steady_state, 0.0)
    get_jacobian_at_steadystate!(work,
                                 results.endogenous_steady_state,
                                 results.exogenous_steady_state, 
                                 m,
                                 2)
    if isnothing(get(options,"dr_cycle_reduction", nothing))
        algo = "GS"
    else
        algo = "CR"
    end
    ws = LinearRationalExpectationsWs(algo,
                                      m.endogenous_nbr,
                                      m.exogenous_nbr,
                                      m.exogenous_deterministic_nbr,
                                      m.i_fwrd_b,
                                      m.i_current,
                                      m.i_bkwrd_b,
                                      m.i_both,
                                      m.i_static)
    LinearRationalExpectations.remove_static!(work.jacobian, ws)
    if algo == "GS"
        LinearRationalExpectations.get_de!(ws, work.jacobian)
        options["generalized_schur"]["criterium"] = 1 + 1e-6
    else
        LinearRationalExpectations.get_abc!(ws, work.jacobian)
        options["cyclic_reduction"]["tol"] = 1e-8
    end
    LinearRationalExpectations.first_order_solver!(results.linearrationalexpectations,
                                                   algo,
                                                   work.jacobian,
                                                   options,
                                                   ws)
end

function compute_prefect_foresight_setup(context); end
function compute_perfect_foresight_solver(context); end

function get_dynamic_variables!(y::Vector{Float64}, steadystate::Vector{Float64}, lli::Matrix{Int64})
    for i = 1:size(lli,2)
        value = steadystate[i]
        for j = 1:size(lli,1)
            k = lli[j, i]
            if k > 0
                y[k] = value
            end
        end
    end
end

function get_jacobian_at_steadystate!(work::Work, steadystate, exogenous, m::Model, period::Int64)
    lli = m.lead_lag_incidence
    get_dynamic_variables!(work.dynamic_variables, steadystate, lli)
    Base.invokelatest(m.dynamic!.dynamic!,
                      work.temporary_values,
                      work.residuals,
                      work.jacobian,
                      work.dynamic_variables,
                      repeat(exogenous', 3, 1),
                      work.params,
                      steadystate,
                      period)  
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
    vΣy = view(Σy, m.i_bkwrd_b, m.i_bkwrd_b)
    vΣy .= Σ
    n_non_states = m.endogenous_nbr - m.n_states
    B2 = zeros(n_non_states, m.exogenous_nbr)
    vr2 = view(results.linearrationalexpectations.g1_2, m.i_non_states, :)
    B2 .= vr2
    Σ = zeros(m.endogenous_nbr, n_non_states)
    A2 = zeros(n_non_states, m.n_states)
    vr3 = view(results.linearrationalexpectations.g1_1, m.i_non_states, :)
    A2 .= vr3
    tmp1 = zeros(m.endogenous_nbr, m.n_states)
    mul!(tmp1, g1_1, vΣy)
    mul!(Σ, tmp1, A2')
    tmp2 = zeros(m.endogenous_nbr, m.exogenous_nbr)
    mul!(tmp2, g1_2, m.Sigma_e)
    mul!(Σ, tmp2, B2', 1.0, 1.0)
    display(Σ)
    vΣy = view(Σy, :, m.i_non_states)
    vΣy .= Σ
    vΣy = view(Σy, m.i_non_states, m.i_bkwrd_b)
    vΣ = view(Σ, m.i_bkwrd_b, :)
    vΣy .= vΣ'
end
    
function calib_smoother!(context, field, varobs, varobs_ids)
    model = context.models[1]
    options = context.options
    results = context.results.model_results[1]
    options["calib_smoother"] = Dict()
    copy!(options["calib_smoother"], field["options"])
    file = get(options["calib_smoother"], "datafile", "")
    if (file = get(options["calib_smoother"], "datafile", "")) != ""
        df = DataFrame(CSV.File(file))
        if uppercase(names(df)[1]) in ["COLUMN1", "DATE", "DATES",
                                      "PERIOD", "PERIODS"]
            frequency = identify_period_frequency(uppercase(df[1, 1]))
        else
            frequency = Periods.Undated
            firstperiod = 1
        end

        timedataframe = TimeDataFrame(df, frequency, firstperiod)
    end

    ny = length(varobs)
    nobs = size(df, 1) - 1
    start = 1
    last = nobs
    Y = Matrix{Float64}(undef, ny, nobs)
    for (i, v) in enumerate(varobs)
        Y[i, :] .= df[start + 1:last + 1, v] .- results.endogenous_steady_state[varobs_ids[i]]
    end

    statevar_ids = model.i_bkwrd_b
    kalman_statevar_ids = collect(1:model.endogenous_nbr)
    ns = length(kalman_statevar_ids)
    np = model.exogenous_nbr
    kws = KalmanSmootherWs{Float64, Int64}(ny, ns, model.exogenous_nbr, nobs)
    c = zeros(ny)
    k1 = findall(in(varobs_ids), kalman_statevar_ids)
    k2 = findall(in(statevar_ids), kalman_statevar_ids)
    Z = zeros(ny, ns)
    for i in varobs_ids
        Z[i, varobs_ids[1]] = 1.0
    end
    H = zeros(ny, ny)
    d = zeros(ns)
    T = zeros(ns, ns)
    vg1 = view(context.results.model_results[1].linearrationalexpectations.g1_1, kalman_statevar_ids, :)
    T[:, k2] .= vg1
    R = zeros(ns, np)
    vg2 = view(context.results.model_results[1].linearrationalexpectations.g1_2, kalman_statevar_ids, :)
    R .= vg2
    Q = model.Sigma_e
    a0 = zeros(ns)
    alphah = zeros(ns, nobs)
    epsilonh = zeros(ny, nobs)
    etah = zeros(np, nobs)
    P = zeros(ns, ns)
    vv = view(context.results.model_results[1].endogenous_variance, kalman_statevar_ids, kalman_statevar_ids)
    P .= vv 
    Valpha = zeros(ns, ns, nobs)
    Vepsilon = zeros(ny, ny, nobs)
    Veta = zeros(np, np, nobs)
    presample = 0
    data_pattern = Vector{Vector{Integer}}(undef, 0)
    for i = 1:nobs
        push!(data_pattern, findall(Y[:, i] .!= NaN))
    end

    kalman_smoother!(Y, c, Z, H, d, T, R, Q, a0,  P, alphah, epsilonh, etah,
                     Valpha, Vepsilon, Veta, start, last, presample,
                     kws, data_pattern)

    alphah .+= results.endogenous_steady_state
    results.smoother["alphah"] = alphah

    Plots.plot(alphah, layout = 7)
end

function identify_period_frequency(period)::Periods.Frequency
    if !isuppercase(period)
        period = uppercase(period)
    end
    if 'Y' in period
        frequency = Periods.Year
    elseif 'A' in period
        frequency = Periods.Year
    elseif 'S' in period
        frequency = Period.Semester
    elseif 'H' in period
        frequency = Periods.Semester
    elseif 'Q' in period
        frequency = Periods.Quarter
    elseif 'M' in period
        frequency = Periods.Month
    elseif 'W' in period
        frequency = Periods.Week
    elseif 'D' in period
        frequency = Periods.Day
    elseif (cdash = count("-", period)) == 2
        frequency = Periods.Day
    elseif cdash == 1
        frequency = Periods.Month
    else
        throw(ErrorException)
    end
end

function deterministic_trends!(context, field)
    model = context.models[1]
    work = context.work
    trend = context.results.model_results[1].endogenous_trend
    symboltable = context.symboltable
    work.model_has_trend = true
    fill!(trend, 0.0)
    for (key, value) in field["trends"]
        trend[symboltable[key].orderintype,2] = work.params[symboltable[value].orderintype]
    end
end
