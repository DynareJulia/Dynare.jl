## IRFs
"""
    irf(; context = context, model = 1)

returns a dictionary. The keys are the name of the exogenous variables. For each key, the value
is an AxisArrayTable with the impulse response of all endogenous variables to a shock the exogenous variable
in the key.

    irf(shock; context = context, model = 1)

returns an AxisArrayTable with the impulse response of all endogenous variables

    irf(shock, varnames; context = context, model = 1)

returns an AxisArrayTable with the impulse response of `varnames`.

# Arguments
- shock::Union{String, Symbol}: name of the exogenous variable that is shocked
- varnames::Union{String, Symbol, Vector{String}, Vector{Symbol}}: list of endogenous variables

# Keyword arguments
- context::Context: context of the computation (default: context)
- model: model number (must be always 1 for the time being)
"""
irf(; context = context, model = 1) = context.results.model_results[model].irfs
irf(shock; context = context, model = 1) =
    context.results.model_results[model].irfs[Symbol(shock)]# IRF for a given model, shock and variable
irf(shock, varnames; context = context, model = 1) =
    context.results.model_results[model].irfs[Symbol(shock)][Symbol.(varnames)]
#irf(shock::String, variable; context = context, model = 1) =
#    irf(Symbol(shock), Symbol.(variable); context = context, model = model)

## SIMULATIONS
"""
    simulation(; context = context, model = 1)
Returns a vector of AxisArrayTable with all the simulations done

    simulation(simnbr::Int64; context = context, model = 1, firstperiod=simulation(context=context, model=model)[simnbr].firstperiod, lastperiod=simulation(context=context, model=model)[simnbr].lastperiod) 
Returns an AxisArrayTable with simulation number simnbr

    simulation(varname; context = context, model = 1, simnbr = length(context.results.model_results[model].simulations), firstperiod = context.results.model_results[model].simulations[simnbr].firstperiod, lastperiod = context.results.model_results[model].simulations[simnbr].lastperiod)
Returns an AxisArrayTable with variables in varname in simulation number simnbr (by default: last simulation)

# Arguments
- simnbr: number of the simulation
- varname: list of variable names

# Keyword arguments
- context: context of the computation
- firstperiod: first period
- lastperiod: last period
- model: model number (always 1 for the time being)
- simnbr: number of the simulation

"""
simulation(; context = context, model = 1) =
    context.results.model_results[model].simulations

# all simulated variables for a given model and a given simulation
function simulation(simnbr::Int64; 
        context = context, 
        model = 1, 
        firstperiod=simulation(context=context, model=model)[simnbr].firstperiod, 
        lastperiod=simulation(context=context, model=model)[simnbr].lastperiod)
    return simulation(context=context, 
        model=model)[simnbr].data[firstperiod..lastperiod]
end

# one simulated variable for a given model and a given simulation
function simulation(varname; 
        context = context, 
        model = 1, 
        simnbr = length(context.results.model_results[model].simulations), 
        firstperiod = context.results.model_results[model].simulations[simnbr].firstperiod,
        lastperiod = context.results.model_results[model].simulations[simnbr].lastperiod)
    return context.results.model_results[model].simulations[simnbr].data[firstperiod..lastperiod, Symbol.(varname)]
end

## FILTER

function getfilter(; context = context::Context,
    firstperiod = row_labels(context.results.model_results[1].filter)[1],
    lastperiod = row_labels(context.results.model_results[1].filter)[end],
    )
    return context.results.model_results[1].filter[firstperiod:lastperiod]
end

function getfilter(varnames; context = context,
    firstperiod = row_labels(context.results.model_results[1].filter)[1],
    lastperiod = row_labels(context.results.model_results[1].filter)[end],
    )
    return context.results.model_results[1].filter[firstperiod:lastperiod, Symbol.(varnames)]
end

## SMOOTHER
"""
    smoother(; context = context, firstperiod = row_labels(context.results.model_results[1].smoother)[1], lastperiod = row_labels(context.results.model_results[1].smoother)[end])
Return an AxisArrayTable with the smoothed value of the endogenous and exogenous variables

    smoother(varnames; context = context, firstperiod = row_labels(context.results.model_results[1].smoother)[1], lastperiod = row_labels(context.results.model_results[1].smoother)[end], )
Returns an AxisArrayTable with the smoothed value of variables `varname`.

#Arguments
- varnames: a list of variable names. It can be a stringm a symbol or a vector of strings or symbols.

# Keyword arguments
- context: the context of the computation
- firstperiod: first period
- lastperiod: last period
"""
function smoother(; context = context,
    firstperiod = row_labels(context.results.model_results[1].smoother)[1],
    lastperiod = row_labels(context.results.model_results[1].smoother)[end],
    )
    return context.results.model_results[1].smoother[firstperiod:lastperiod]
end

function smoother(varnames; context = context,
    firstperiod = row_labels(context.results.model_results[1].smoother)[1],
    lastperiod = row_labels(context.results.model_results[1].smoother)[end],
    )
    return context.results.model_results[1].smoother[firstperiod:lastperiod, Symbol.(varnames)]
end

## FORECAST
"""
    forecast(; context = context, firstperiod = row_labels(context.results.model_results[1].forecast[1])[1], lastperiod = row_labels(context.results.model_results[1].forecast[1])[end], informationperiod = Undated(typemin(Int)) ) 
 Returns an AxisArrayTable with the forecasted value of the endogenous variables.

    forecast(varnames; context = context, firstperiod = row_labels(context.results.model_results[1].forecast)[1], lastperiod = row_labels(context.results.model_results[1].forecast)[end], informationperiod = Undated(typemin(Int)) )
Returns an AxisArrayTable with the forecasted values of variables `varnames`

# Arguments
- varnames

# Keyword arguments
- firstperiod:
- information period: period when information about future shocks or future values of some endogenous variables is known and taken into account when computing the forecast (default: the last period before the beginning of the forecast)
- lastperiod:
"""

function forecast(; context = context,
    firstperiod = row_labels(context.results.model_results[1].forecast[1])[1],
    lastperiod = row_labels(context.results.model_results[1].forecast[1])[end],
    informationperiod = Undated(typemin(Int))
    )
    forecast_ = context.results.model_results[1].forecast
    if length(forecast_) == 1
        return forecast_[1][firstperiod:lastperiod]
    else
        D = Dict((row_labels(d)[1] => d) for d in forecast_)
        if informationperiod != Undated(typemin(Int))
            return D[informationperiod]
        else
            return D
        end 
    end
end

function forecast(varnames; 
    context = context,
    firstperiod = row_labels(context.results.model_results[1].forecast)[1],
    lastperiod = row_labels(context.results.model_results[1].forecast)[end],
    informationperiod = Undated(typemin(Int))
    )
    forecast_ = context.results.model_results[1].forecast
    if length(forecast_) == 1
        return forecast_[1][firstperiod:lastperiod, Symbol.(varnames)]
    else
        D = Dict((row_labels(d)[1] => d[Symbol(varnames)]) for d in forecast_)
        if informationperiod != Undated(typemin(Int))
            return D[informationperiod][firstperiod:lastperiod]
        else
            return D
        end 
    end
end

