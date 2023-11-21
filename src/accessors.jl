## IRFs
# all the IRFs for a given model
irf(; context = context, model = 1) = context.results.model_results[model].irfs
# all the IRFs for a given model and a given shock
irf(shock::Symbol; context = context, model = 1) =
    context.results.model_results[model].irfs[shock]
irf(shock::String; context = context, model = 1) =
    irf(Symbol(shock); context = context, model = model)
# IRF for a given model, shock and variable
irf(shock::Symbol, variable::Symbol; context = context, model = 1) =
    getproperty(context.results.model_results[model].irfs[shock], variable)
irf(shock::String, variable::String; context = context, model = 1) =
    irf(Symbol(shock), Symbol(variable); context = context, model = model)

## SIMULATIONS
# all simulations for a model
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
#=
function simulation(varname::String; 
        context = context, 
        model = 1, 
        simnbr = length(context.results.model_results[model].simulations), 
        firstperiod=context.results.model_results[model].simulations[simnbr].firstperiod,
        lastperiod=context.results.model_results[model].simulations[simnbr].lastperiod)
    return context.results.model_results[model].simulations[simnbr].data[firstperiod..lastperiod, Symbol(varname)]
end

# subset of simulated variables for a  given model and a given simulation
function simulation(varnames::Vector{Symbol}; 
                    context = context, 
                    model = 1, 
                    simnbr = length(context.results.model_results[model].simulations), 
                    firstperiod=context.results.model_results[model].simulations[simnbr].firstperiod,
                    lastperiod=context.results.model_results[model].simulations[simnbr].lastperiod)
    context.results.model_results[model].simulations[simnbr].data[:, varnames]
end

simulation(varnames::Vector{String};
           context = context,
           model = 1,
           simnbr = length(context.results.model_results[model].simulations), 
           firstperiod=simulation(context=context, model=model).firstperiod, lastperiod=simulation(context=context, model=model).lastperiod) = 
               simulation([Symbol(v) for v in varnames]; context = context, model = model, simnbr = simnbr, firstperiod=firstperiod, lastperiod=lastperiod)

simulation(varnames::Tuple;
           context = context,
           model = 1,
           simnbr = length(context.results.model_results[model].simulations), 
           firstperiod=context.results.model_results[model].simulations[simnbr].firstperiod,
           lastperiod=context.results.model_results[model].simulations[simnbr].lastperiod) =
               simulation([Symbol(v) for v in varnames]; context = context, model = model, simnbr = simnbr, firstperiod=firstperiod, lastperiod=lastperiod)
=#

## SMOOTHER
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
    return context.results.model_results[1].smoother[firstperiod:lastperiod, varnames]
end

## FORECAST
function forecast(; context = context,
    firstperiod = row_labels(context.results.model_results[1].smoother)[1],
    lastperiod = row_labels(context.results.model_results[1].smoother)[end],
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
    firstperiod = row_labels(context.results.model_results[1].smoother)[1],
    lastperiod = row_labels(context.results.model_results[1].smoother)[end],
    informationperiod = Undated(typemin(Int))
    )
    forecast_ = context.results.model_results[1].forecast
    if length(forecast_) == 1
        return forecast_[1][firstperiod:lastperiod]
    else
        D = Dict((row_labels(d)[1] => d[varnames]) for d in forecast_)
        if informationperiod != Undated(typemin(Int))
            return D[informationperiod]
        else
            return D
        end 
    end
end

