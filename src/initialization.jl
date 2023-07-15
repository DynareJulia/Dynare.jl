function histval!(context::Context, field::Dict{String,Any})
    symboltable = context.symboltable
    m = context.models[1]
    histval = Matrix{Union{Float64, Missing}}(missing, m.orig_maximum_lag, m.endogenous_nbr + m.exogenous_nbr)
    for v in field["vals"]
        k = symboltable[v["name"]::String].orderintype::Int64
        l = m.orig_maximum_lag + v["lag"]::Int64
        histval[l, k] = dynare_parse_eval(v["value"]::String, context)
    end
    aat = AxisArrayTable(
        Matrix{Union{Float64,Missing}}(histval),
        Undated(1):Undated(size(histval,1)),
        [Symbol(s) for s in vcat(get_endogenous(symboltable),
            get_exogenous(symboltable))],
    )
    DFunctions.dynamic_auxiliary_variables!(aat, context.work.params)
    context.work.histval = Matrix(aat)
end

function get_date(optionname::String, options)
    if haskey(options, optionname)
        return periodparse(options[optionname])
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
                error("Variable {$v} is missing")
            else
                auxiliary_variables_present = false
                break
            end
        end
    end
    return auxiliary_variables_present
end

Delta(::Type{YearSE}) = Year
Delta(::Type{SemesterSE}) = Semester
Delta(::Type{QuarterSE}) = Quarter
Delta(::Type{MonthSE}) = Month
Delta(::Type{WeekSE}) = Week
Delta(::Type{DaySE}) = Day
Delta(::Type{Undated}) = Undated
"""
    check_periods_options(options::Options; required_lags::Int64=0, required_leads::Int64=0)
checks the following assertions if the fields exist:
    - all periods must have the same frequency
    - first_simulation_period >= first_obs + required_periods
    - last_obs >= last_simulation_periods + required
    - last_obs == first_obs + nobs
"""
function check_periods_options(options; required_lags::Int64 = 0, required_leads::Int64 = 0)
    first_obs = get_date("first_obs", options)
    tf = typeof(first_obs)
    last_obs = get_date("last_obs", options)
    tl = typeof(last_obs)
    first_simulation_period = get_date("first_simulation_period", options)
    tfs = typeof(first_simulation_period)
    last_simulation_period = get_date("last_simulation_period", options)
    tls = typeof(last_simulation_period)
    nobs_str = get(options, "nobs", nothing)
    nobs = !isnothing(nobs_str) ? parse(Int64, nobs_str) : nothing

    common_type = UndatedDate
    for f in [first_obs, last_obs, first_simulation_period, last_simulation_period]
        if !isnothing(f)
            tf = typeof(f)
            if common_type == UndatedDate
                common_type = tf
            else
                @assert tf == common_type "period fields must have the same frequency"
            end
        end
    end

    if !isnothing(first_obs) && !isnothing(first_simulation_period)
        @assert first_simulation_period >= first_obs + Delta(tf)(required_lags) "first_simulation_period should be greater or equal to first_obs + $required_lags."
    end

    if !isnothing(last_obs) && !isnothing(last_simulation_period)
        @assert last_obs >= last_simulation_period + Delta(tf)(required_leads) "last_obs should be greater or equal to first_simulation_period + $required_leads."
    end

    if !isnothing(first_obs) && !isnothing(last_obs) && !isnothing(nobs)
        @assert first_obs + Delta(tf)(nobs - 1) == last_obs "last_obs must equal first_obs + nobs"
    end

    return (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        common_type,
    )

end

function add_auxiliary_variables(
    aat::AxisArrayTable,
    symboltable::SymbolTable,
    endogenous_nbr::Int64,
    original_endogenous_nbr::Int64,
)
    aa = getfield(aat, :data)
    nrow = size(aa, 1)
    vnames = names(aat)
    for v in get_endogenous(symboltable)[original_endogenous_nbr+1:endogenous_nbr]
        if !(v in vnames)
            insertcols!(df, Symbol(v) => Vector{Union{Float64,Missing}}(missing, nrow))
        end
    end
    for v in get_exogenous(symboltable)
        if !(v in vnames)
            insertcols!(df, Symbol(v) => Vector{Union{Float64,Missing}}(missing, nrow))
        end
    end
end

function histval_file!(context::Context, field::Dict{String,Any})
    m = context.models[1]
    options = field["options"]
    symboltable = context.symboltable

    if haskey(options, "datafile")
        tdf = MyTimeDataFrame(options["datafile"])
        add_auxiliary_variables(
            tdf,
            symboltable,
            m.endogenous_nbr,
            m.original_endogenous_nbr,
        )
    elseif haskey(options, "series")
        tdf = options["series"]
    else
        error("option datafile or series must be provided")
    end

    # if auxiliary variables are present, we need only one period of observations
    auxiliary_variables_present = check_predetermined_variables(symboltable, m, names(tdf))
    required_lags = (auxiliary_variables_present) ? 1 : m.orig_maximum_lag
    # check options consistency
    (first_obs, last_obs, first_simulation_period, last_simulation_period, nobs, P) =
        check_periods_options(options; required_lags)

    tdf_periods = getfield(tdf, :periods)
    if P != UndatedDate && P != typeof(tdf_periods)
        error("Data frequency is different from frequency used in options")
    end

    # DP is type of period interval
    DP = typeof(P(1) - P(0))

    start = stop = nothing
    if isnothing(first_obs)
        if !isnothing(first_simulation_period)
            start = first_simulation_period - DP(required_lags)
            stop = first_simulation_period - DP(1)
        elseif !isnothing(last_obs)
            start = last_obs - DP(required_lags + 1)
            stop = last_obs
        else
            start = P(1)
            stop = start + DP(required_lags - 1)
        end
    else
        start = first_obs
        stop = start + required_lags - 1
    end

    istart = ExtendedDates.value(start) - ExtendedDates.value(tdf_periods[1]) + 1
    istop = ExtendedDates.value(stop) - ExtendedDates.value(tdf_periods[1]) + 1
    data_ = Matrix(getfield(tdf, :data))
    histval = context.work.histval
    let colindex, endogenous_names
        if auxiliary_variables_present
            colindex = m.i_bkwrd_b
            endogenous_names = get_endogenous(symboltable)
        else
            colindex = view(m.i_bkwrd_b, 1:m.orig_endogenous_nbr)
            endogenous_names = view(get_endogenous(symboltable), 1:m.orig_endogenous_nbr)
        end
        columns =
            (auxiliary_variables_present) ? m.i_bkwrd_b :
            view(m.i_bkwrd_b, 1:m.orig_endogenous_nbr)
        dnames = names(tdf)
        for (i, vname) in enumerate(endogenous_names)
            colindex = findfirst(x -> x == vname, dnames)
            for j = istart:istop
                histval[j, i] = tdf[j, colindex]
            end
        end
    end
end

function initval_file!(context::Context, field::Dict{String,Any})
    m = context.models[1]
    work = context.work
    options = field["options"]
    symboltable = context.symboltable

    if haskey(options, "datafile")
        tdf = MyTimeDataFrame("$(options["datafile"])")
    else
        error("option datafile or series must be provided")
    end
    tdf_nrow = TimeDataFrames.nrow(tdf)

    # if auxiliary variables are present, we need only one period of observations
    auxiliary_variables_present =
        check_predetermined_variables(context.symboltable, m, names(tdf))
    if auxiliary_variables_present
        required_lags = 1
    else
        required_lags = m.orig_maximum_lag
        add_auxiliary_variables(
            tdf,
            symboltable,
            m.endogenous_nbr,
            m.original_endogenous_nbr,
        )
        DFunctions.dynamic_auxiliary_variables!(tdf, work.params)
    end

    # check options consistency
    (first_obs, last_obs, first_simulation_period, last_simulation_period, nobs, P) =
        check_periods_options(options; required_lags)

    tdf_periods = getfield(tdf, :periods)
    if P != UndatedDate && P != typeof(tdf_periods)
        error("Data frequency is different from frequency used in options")
    end

    # DP is type of period interval
    DP = typeof(P(1) - P(0))

    start = stop = nothing
    if isnothing(first_obs)
        if isnothing(first_simulation_period)
            start = UndatedDate(1)
        else
            start = first_simulation_period - DP(required_lags)
        end
    end
    if isnothing(last_obs)
        if isnothing(last_simulation_period)
            stop = UndatedDate(tdf_nrow)
        else
            stop = last_simulation_period + DP(required_leads)
        end
    end

    istart = ExtendedDates.value(start) - ExtendedDates.value(tdf_periods[1]) + 1
    istop = ExtendedDates.value(stop) - ExtendedDates.value(tdf_periods[1]) + 1
    nobs = istop - istart + 1
    data_ = Matrix(getfield(tdf, :data))
    let colindex, endogenous_names, exogenous_names, exogenousdeterministic_names
        endogenous_names = get_endogenous(symboltable)
        exogenous_names = get_exogenous(symboltable)
        exogenousdeterministic_names = get_exogenousdeterministic(symboltable)
        dnames = names(tdf)
        work.initval_endogenous = Matrix{Float64}(undef, m.endogenous_nbr, nobs)
        work.initval_exogenous = Matrix{Float64}(undef, m.exogenous_nbr, nobs)
        work.initval_exogenous_deterministic =
            Matrix{Float64}(undef, m.exogenous_deterministic_nbr, nobs)
        lli = m.lead_lag_incidence
        for (i, vname) in enumerate(endogenous_names)
            colindex = findfirst(x -> x == vname, dnames)
            if !isnothing(colindex)
                work.initval_endogenous[i, :] .= tdf[!, colindex][istart:istop]
            end
        end
        for (i, vname) in enumerate(exogenous_names)
            colindex = findfirst(x -> x == vname, dnames)
            if isnothing(colindex)
                @warn(
                    "INITVAL_FILE: Variable $vname doesn't exist in $(options["datafile"]). Initial value set to zero"
                )
            else
                work.initval_exogenous[i, :] .= tdf[!, colindex][istart:istop]
            end
        end
        for (i, vname) in enumerate(exogenousdeterministic_names)
            colindex = findfirst(x -> x == vname, dnames)
            if isnothing(colindex)
                @warn(
                    "INITVAL_FILE: Variable $vname doesn't exist in $(options["datafile"]). Initial value set to zero"
                )
            else
                work.initval_exogenous_deterministic[i, :] .= tdf[!, colindex][istart:istop]
            end
        end
    end
end

function param_init!(context::Context, field::Dict{String,Any})
    params = context.work.params
    symboltable::Dict{String,DynareSymbol} = context.symboltable
    s = symboltable[field["name"]::String]
    k = s.orderintype::Int64
    params[k] = Float64(dynare_parse_eval(field["value"]::String, context))
end

function initval!(context::Context, field::Dict{String,Any})
    symboltable = context.symboltable
    m = context.models[1]
    work = context.work
    work.initval_endogenous = zeros(m.endogenous_nbr, 1)
    work.initval_exogenous = zeros(m.exogenous_nbr, 1)
    work.initval_exogenous_deterministic = zeros(m.exogenous_deterministic_nbr, 1)
    xs = [Endogenous, Exogenous, ExogenousDeterministic]
    xw = [work.initval_endogenous, work.initval_exogenous, work.initval_exogenous_deterministic]
    for v in field["vals"]
        s = symboltable[v["name"]::String]
        typ = s.symboltype
        k = s.orderintype::Int64
        value = dynare_parse_eval(v["value"]::String, context, xs = xs, xw = xw)
        if typ == Endogenous
            work.initval_endogenous[k, 1] = value
        elseif typ == Exogenous
            work.initval_exogenous[k, 1] = value
        elseif typ == ExogenousDeterministic
            work.initval_exogenous_determinisitic[k, 1] = value
        else
            throw(error("$(v["name"]) can't be set in INITVAL"))
        end
    end
    params = context.work.params
    DFunctions.static_auxiliary_variables!(work.initval_endogenous, work.initval_exogenous, params)
end

function endval!(context::Context, field::Dict{String,Any})
    symboltable = context.symboltable
    m = context.models[1]
    work = context.work
    endval_endogenous = zeros(m.endogenous_nbr)
    endval_exogenous = zeros(m.exogenous_nbr)
    endval_exogenous_det = zeros(m.exogenous_deterministic_nbr)
    xs = [Endogenous, Exogenous, ExogenousDeterministic]
    xw = [endval_endogenous, endval_exogenous, endval_exogenous_det]
    work.endval_endogenous = zeros(m.endogenous_nbr,1)
    work.endval_exogenous = zeros(m.exogenous_nbr,1)
    work.endval_exogenous_deterministic = zeros(m.exogenous_deterministic_nbr,1)
    for v in field["vals"]
        s = symboltable[v["name"]::String]
        typ = s.symboltype
        k = s.orderintype::Int64
        value = dynare_parse_eval(v["value"]::String, context, xs = xs, xw = xw)
        if typ == Endogenous
            work.endval_endogenous[k] = value
        elseif typ == Exogenous
            work.endval_exogenous[k] = value
        elseif typ == ExogenousDeterministic
            work.endval_exogenous_determinisitic[k] = value
        else
            throw(error("$(v["name"]) can't be set in ENDVAL"))
        end
    end
    params = context.work.params
    DFunctions.static_auxiliary_variables!(work.endval_endogenous, work.endval_exogenous, params)
end

function shocks!(context::Context, field::Dict{String,Any})
    observed_variables = context.work.observed_variables
    Sigma = context.models[1].Sigma_e
    Sigma_m = context.work.Sigma_m
    shocks = context.work.shocks
    symboltable = context.symboltable
    set_variance!(Sigma, Sigma_m, field["variance"], symboltable, observed_variables)
    set_stderr!(Sigma, Sigma_m, field["stderr"], symboltable, observed_variables)
    set_covariance!(Sigma, Sigma_m, field["covariance"], symboltable, observed_variables)
    set_correlation!(Sigma, Sigma_m, field["correlation"], symboltable, observed_variables)
    if haskey(field, "deterministic_shocks")
        context.work.shocks = 
            set_deterministic_shocks!(
                    field["deterministic_shocks"],
                    symboltable,
                    context.models[1].exogenous_nbr,
                    context.results.model_results[1].trends.exogenous_steady_state,
            )
    end
end

function scenario!(; 
                   name::Symbol,
                   period::PeriodsSinceEpoch,
                   value::Float64,
                   context::Context = context,
                   infoperiod::PeriodsSinceEpoch = Undated(1),
                   exogenous::Symbol = Symbol()
                   )
    scenario = context.work.scenario
    if haskey(scenario, infoperiod)
        if haskey(scenario[infoperiod], period)
            scenario[infoperiod][period][name] = (value => exogenous)
        else
            scenario[infoperiod][period] = Dict(name => (value => exogenous))
        end
    else
        scenario[infoperiod] = Dict(period => Dict(name => (value => exogenous)))
    end
end
               
function set_variance!(
    Sigma::Matrix{Float64},
    Sigma_m::Matrix{Float64},
    variance::Vector{Any},
    symboltable::SymbolTable,
    observed_variables
)
    for v in variance
        value = dynare_parse_eval(v["variance"]::String, context)
        sv = symboltable[v["name"]]
        if sv.symboltype == Exogenous
            k = sv.orderintype
            Sigma[k, k] = value
        elseif sv.symboltype === Endogenous
            k = findfirst(ov -> ov == v["name"], observed_variables)
            Sigma_m[k, k] = value
        end 
    end
end

function set_stderr!(
    Sigma::Matrix{Float64}, 
    Sigma_m::Matrix{Float64},
    stderr::Vector{Any}, 
    symboltable::SymbolTable,
    observed_variables
)
    for s in stderr
        value = dynare_parse_eval(s["stderr"]::String, context)
        sv = symboltable[s["name"]]
        if sv.symboltype == Exogenous
            k = sv.orderintype
            Sigma[k, k] = value * value
        elseif sv.symboltype === Endogenous
            k = findfirst(ov -> ov == s["name"], observed_variables)
            Sigma_m[k, k] = value * value
        end 
    end
end

function set_covariance!(
    Sigma::Matrix{Float64},
    Sigma_m::Matrix{Float64},
    covariance::Vector{Any},
    symboltable::SymbolTable,
    observed_variables
)
    for c in covariance
        sv1 = symboltable[c["name"]]
        sv2 = symboltable[c["name2"]]
        value = dynare_parse_eval(c["covariance"]::String, context)
        if sv1.symboltype == Exogenous
            k1 = sv1.orderintype
            k2 = sv2.orderintype
            Sigma[k1, k2] = value
            Sigma[k2, k1] = value
        elseif sv1.symboltype === Endogenous
            k1 = findfirst(c["name", observed_variables])
            k2 = findfirst(c["name", observed_variables])
            Sigma_m[k1, k2] = value
            Sigma_m[k2, k1] = value
        end 
    end
end

function set_correlation!(
    Sigma::Matrix{Float64},
    Sigma_m::Matrix{Float64},
    correlation::Vector{Any},
    symboltable::SymbolTable,
    observed_variables
)
    for c in correlation
        value = dynare_parse_eval(c["correlation"]::String, context)
        sv1 = symboltable[s["name"]]
        sv2 = symboltable[s["name2"]]
        if sv.symboltype == Exogenous
            k1 = sv1.orderintype
            k2 = sv2.orderintype
            x = sqrt(Sigma[k1, k1] * Sigma[k2, k2]) * value
            Sigma[k1, k2] = x
            Sigma[k2, k1] = x
        elseif sv.symboltype === Endogenous
            k1 = findfirst(c["name", observed_variables])
            k2 = findfirst(c["name2", observed_variables])
            x = sqrt(Sigma_m[k1, k1] * Sigma_m[k2, k2]) * value
            Sigma_m[k1, k2] = x
            Sigma_m[k2, k1] = x
        end 
    end
end

function set_deterministic_shocks!(
    shocks::Vector{Any},
    symboltable::SymbolTable,
    exogenous_nbr::Int64,
    exogenous_steady_state::Vector{Float64},
)
    pmax = maximum((s) -> maximum(p -> p["period2"], s["values"]), shocks)::Int64
    x= Vector{Float64}(undef, pmax * exogenous_nbr)
    x .= repeat(exogenous_steady_state, pmax)
    for s in shocks
        for v in s["values"]
            i = symboltable[s["var"]].orderintype
            d1 = v["period1"]
            d2 = v["period2"]
            val = dynare_parse_eval(v["value"], context)
            set_exogenous(x, i, d1, d2, val, exogenous_nbr)
        end
    end
    return x
end

function set_exogenous(
    x::Vector{Float64},
    i::Int64,
    p1::Int64,
    p2::Int64,
    v::Real,
    exogenous_nbr::Int64,
)
    offset = (p1 - 1) * exogenous_nbr
    for p = p1:p2
        x[offset + i] = v
        offset += exogenous_nbr
    end
end

function load_params!(context::Context, filename::String)
    param_nbr = context.models[1].parameter_nbr
    parameters = context.work.params
    symboltable = context.symboltable
    open(filename) do io
        variables = []
        param_names = get
        for line in readlines(io)
            elem = split(line)
            if haskey(symboltable, elem[1])
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
    parameters = context.work.params
    symboltable = context.symboltable
    open(filename) do io
        variables = []
        param_names = get
        for line in readlines(io)
            elem = split(line)
            if haskey(symboltable, elem[1])
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
    DFunctions.static_auxiliary_variables!(endogenous, exogenous, parameters)
end
#=
import Base: +, -, *, /

+(t::TimeDataFrame, a::Float64) = a .+ t     
-(t::TimeDataFrame, a::Float64) = t .- a     
*(t::TimeDataFrame, a::Float64) = a .* t   
/(t::TimeDataFrame, a::Float64) = t ./ a     
+(a::Float64, t::TimeDataFrame) = a .+ t     
-(a::Float64, t::TimeDataFrame) = a .- t     
*(a::Float64, t::TimeDataFrame) = a .* t     
/(a::Float64, t::TimeDataFrame) = a ./ t     
=#
