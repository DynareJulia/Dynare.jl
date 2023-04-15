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
function simulation(varname::Symbol; 
        context = context, 
        model = 1, 
        simnbr = 1, 
        firstperiod = context.results.model_results[model].simulations[simnbr].firstperiod,
        lastperiod = context.results.model_results[model].simulations[simnbr].lastperiod)
    return context.results.model_results[model].simulations[simnbr].data[firstperiod..lastperiod, varname]
end

function simulation(varname::String; 
        context = context, 
        model = 1, 
        simnbr = 1, 
        firstperiod=context.results.model_results[model].simulations[simnbr].firstperiod,
        lastperiod=context.results.model_results[model].simulations[simnbr].lastperiod)
    return context.results.model_results[model].simulations[simnbr].data[firstperiod..lastperiod, Symbol(varname)]
end

# subset of simulated variables for a  given model and a given simulation
function simulation(varnames::Vector{Symbol}; 
                context = context, 
                model = 1, 
                simnbr = 1,
                firstperiod=context.results.model_results[model].simulations[simnbr].firstperiod,
                lastperiod=context.results.model_results[model].simulations[simnbr].lastperiod)
    context.results.model_results[model].simulations[simnbr].data[:, varnames]
end

simulation(varnames::Vector{String}; context = context, model = 1, simnbr = 1,
firstperiod=simulation(context=context, model=model).firstperiod, lastperiod=simulation(context=context, model=model).lastperiod) = 
simulation([Symbol(v) for v in varnames]; context = context, model = model, simnbr = simnbr, firstperiod=firstperiod, lastperiod=lastperiod)

simulation(varnames::Tuple; context = context, model = 1, simnbr = 1,
firstperiod=context.results.model_results[model].simulations[simnbr].firstperiod,
lastperiod=context.results.model_results[model].simulations[simnbr].lastperiod) =
simulation([Symbol(v) for v in varnames]; context = context, model = model, simnbr = simnbr, firstperiod=firstperiod, lastperiod=lastperiod)
