function histval!(context::Context, field::Dict{String, Any})
    symboltable = context.symboltable
    m = context.models[1]
    histval = zeros(m.orig_maximum_lag, m.endogenous_nbr)
    for v in field["vals"]
        k = symboltable[v["name"]::String].orderintype::Int64
        l = m.orig_maximum_lag - v["lag"]::Int64
        histval[l, k] = dynare_parse_eval(v["value"]::String, context)
    end
    context.work.histval .= histval
end

function get_date(optionname::String, options)
    if optionname in options
        return ExtendedDates.simpleperiod(options[optionname])
    else
        return nothing
    end
end

"""
   check_predetermined_variables(symboltable, m, datanames)

checks whether the predetermined variables are presents is datanames.
The function raise an exception if it is false. It returns whether
auxiliary variables are in datanames or not. It returns true if there
are no auxiliary variables in the model.
"""
function check_predetermined_variables(symboltable, m, datanames)
    auxiliary_variables_present = true
    for (i, v) in enumerate(get_endogenous(symboltable))
        if i in m.i_bkwrd_b
            if v in datanames
                continue
            elseif i <= m.original_endogenous_nbr
                error("Variable {v} is missing")
            else
                auxiliary_variables_present = false
            end
        end
    end
    return auxiliary_variables_present
end
    
function histval_file!(context::Context, field::Dict{String, Any})
    m = context.models[1]
    options = field["options"]
    symboltable = context.symboltable
    if "datafile" in options
        data = TimeDataFrame(CSV.File(options["datafile"]))
    else
        data = nothing
    end

    first_obs = get_date("first_obs", options)
    last_obs = get_date("last_obs", options)
    first_simulation_period = get_date("first_simulation_period", options)
    nobs = ("nobs" in options) ? options["nobs"] : nothing
    series = ("series" in options) ? options["series"] : nothing

    if first_obs != nothing && first_simulation_period != nothing
        if typeof(first_obs) != typeof(first_simulation_period)
            error("first_obs && first_simulation_period, if both present, must have the same frequency. Note that only one of these options is necessary.")
        elseif first_simulation - first_obs != required_periods
            error("first_obs && first_simulation_period are inconsistent: first_simulation_period should equal first_obs + $required_periods. Note that only one of these options is necessary.")
        end
    end

    if first_obs != nothing && last_obs != nothing
        if typeof(first_obs) != typeof(last_obs)
            error("first_obs && last_obs must have the same frequency")
        elseif last_obs  != first_obs + required_periods - 1
            error("first_obs && last_obs are inconsistent: last_obs must be equal to first_obs + $(required_periods -1)")
        end
    end
    
    if first_simulation_period != nothing && last_obs != nothing
        if typeof(first_simulation_period) != typeof(last_obs)
            error("first_simulation_period && last_obs must have the same frequency")
        elseif first_simulation_period  != last_obs + 1
            error("first_simulation_period && last_obs are inconsistent: last_obs must be equal first_simulation_period - 1.")
        end
    end

    auxiliary_variables_present =
        check_predetermined_variables(context.symboltable, m, names(data))
    required_perdiods = (auxiliary_variables_present) ? 1 : m.orig_maximum_lag
    
    start = stop = nothing
    if first_obs == nothing
        if first_simulation_period != nothing
            start = first_simulation_period - required_periods
            stop = first_simulation_period - 1
        elseif last_obs != nothing
            start = last_obs - required_periods + 1
            stop = last_obs
        else
            start = UndatedDate(1)
            stop = start + required_periods - 1
        end
    else
        start = first_obs
        stop = start + required_periods - 1
    end
    
    data_periods = data.periods
    if typeof(start) != UndatedDate && typeof(start) != typeof(data_periods)
        error("Data frequency is different from frequency used in options")
    end

    istart = value(start) - value(data_periods[1]) + 1
    istop = value(stop) - value(data_periods[1]) + 1
    offset = m.orig_maximum_lag - stop
    data_ = getfield(getfield(data, :data), :data)
    histval = context.work.histval
    let colindex, endogenous_names
        if auxiliary_variable_present
            colindex = m.i_bkwrd_b
            endogenous_names = get_endogenous(symboltable)
        else
            colindex = view(m.i_bkwrd_b, 1:m.orig_endogenous_nbr)
            endogenous_names = view(get_endogenous(symboltable), 1:m.orig_endogenous_nbr)
        end
        columns = (auxiliary_variables_present) ? m.i_bkwrd_b : view(m.i_bkwrd_b, 1:m.orig_endogenous_nbr)
        for (i, vname) in endogenous_names
            colindex = find
            for j in istart:istop
                histval[j + offset, i] = data_[j, colindex] 
            end
        end
    end
end

function param_init!(context::Context, field::Dict{String, Any})
    params = context.work.params
    symboltable::Dict{String, DynareSymbol} = context.symboltable
    s = symboltable[field["name"]::String]
    k = s.orderintype::Int64
    params[k] = Float64(dynare_parse_eval(field["value"]::String, context))
end

function initval!(context::Context, field::Dict{String, Any})
    symboltable = context.symboltable
    m = context.models[1]
    endogenous_steady_state = zeros(m.endogenous_nbr)
    exogenous_steady_state = zeros(m.exogenous_nbr)
    exogenous_det_steady_state = zeros(m.exogenous_deterministic_nbr)
    for v in field["vals"]
        s = symboltable[v["name"]::String]
        typ = s.symboltype
        k = s.orderintype::Int64
        value = dynare_parse_eval(v["value"]::String, context)
        if typ == Endogenous
            endogenous_steady_state[k] = value
        elseif typ == Exogenous
            exogenous_steady_state[k] = value
        elseif typ == ExogenousDeterministic
            exogenous_det_steady_state[k] = value
        else
            throw(error("$(v["name"]) can't be set in INITVAL"))
        end
    end
    params = context.work.params
    if !isempty(m.set_auxiliary_variables!)
        @show m.set_auxiliary_variables!.set_auxiliary_variables!
        Base.invokelatest(
            m.set_auxiliary_variables!.set_auxiliary_variables!(
                endogenous_steady_state,
                exogenous_steady_state,
                params)
        )
    end
    trends = context.results.results_model[1].trends 
    trends.endogenous_steady_state = endogenous_steady_state
    trends.exogenous_steady_state = exogenous_steady_state
    trends.exogenous_det_steady_state = exogenous_det_steady_state
end

function shocks!(context::Context, field::Dict{String, Any})
    Sigma = context.models[1].Sigma_e
    symboltable = context.symboltable
    set_variance!(Sigma, field["variance"], symboltable)
    set_stderr!(Sigma, field["stderr"], symboltable)
    set_covariance!(Sigma, field["covariance"], symboltable)
    set_correlation!(Sigma, field["correlation"], symboltable)
end

function set_variance!(Sigma::Matrix{Float64}, variance::Vector{Any}, symboltable::SymbolTable)
    for v in variance
        k =  symboltable[v["name"]::String].orderintype::Int64
        Sigma[k, k] = dynare_parse_eval(v["variance"]::String, context)::Float64
    end
end

function set_stderr!(Sigma::Matrix{Float64}, stderr::Vector{Any}, symboltable::SymbolTable)
    for s in stderr
        k =  symboltable[s["name"]::String].orderintype::Int64
        x = dynare_parse_eval(s["stderr"]::String, context)::Float64
        Sigma[k, k] = x*x
    end
end

function set_covariance!(Sigma::Matrix{Float64}, covariance::Vector{Any}, symboltable::SymbolTable)
    for c in covariance
        k1 =  symboltable[c["name"]::String].orderintype::Int64
        k2 =  symboltable[c["name2"]::String].orderintype::Int64
        Sigma[k1, k2] = dynare_parse_eval(c["covariance"]::String, context)::Float64
        Sigma[k2, k1] = Sigma[k1, k2]
    end
end

function set_correlation!(Sigma::Matrix{Float64}, correlation::Vector{Any}, symboltable::SymbolTable)
    for c in correlation
        k1 =  symboltable[c["name"]::String].orderintype::Int64
        k2 =  symboltable[c["name2"]::String].orderintype::Int64
        corr = dynare_parse_eval(c["correlation"]::String, context)::Float64
        Sigma[k2, k1] = sqrt(Sigma[k1, k1]*Sigma[k2, k2])*corr
    end
end

function load_params!(context::Context, filename::String)
    param_nbr = context.models[1].parameter_nbr;
    parameters = context.work.params
    symboltable = context.symboltable
    open(filename) do io
        variables = [];
        param_names = get
        for line in readlines(io)
            elem = split(line);
            if elem[1] in keys(symboltable)
                s = symboltable[elem[1]]
                if s.symboltype == Parameter
                    k = s.orderintype::Int64
                    parameters[k] = parse(Float64, elem[2])
                end
            end
        end
    end
end

function load_steadystate!(context::Context, filename::String)
    endogenous = context.results.model_results[1].trends.endogenous_steady_state
    exogenous = context.results.model_results[1].trends.exogenous_steady_state
    exogenous_det = context.results.model_results[1].trends.exogenous_det_steady_state
    m = context.models[1]
    parameters = context.work.params
    symboltable = context.symboltable
    open(filename) do io
        variables = [];
        param_names = get
        for line in readlines(io)
            elem = split(line);
            if elem[1] in keys(symboltable)
                s = symboltable[elem[1]]
                k = s.orderintype::Int64
                if s.symboltype == Endogenous
                    endogenous[k] = parse(Float64, elem[2])
                elseif s.symboltype == Exogenous
                    exogenous[k] = parse(Float64, elem[2])
                elseif s.symboltype == ExogenousDeterministic
                    exogenous_det[k] = parse(Float64, elem[2])
                end
            end
        end
    end
    if !isempty(m.set_auxiliary_variables!)
        Base.invokelatest(m.set_auxiliary_variables!.set_auxiliary_variables!,
                          endogenous,
                          exogenous,
                          parameters)
    end
end

function isempty(m::Module)
    return names(m)[1] == :anonymous
end
