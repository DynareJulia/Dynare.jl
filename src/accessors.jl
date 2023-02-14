## IRFs
# all the IRFs for a given model
irf(; thiscontext = context, model = 1) = thiscontext.results.model_results[model].irfs
# all the IRFs for a given model and a given shock
irf(shock::Symbol; thiscontext = context, model = 1) =
    thiscontext.results.model_results[model].irfs[shock]
irf(shock::String; thiscontext = context, model = 1) =
    irf(Symbol(shock); thiscontext = context, model = model)
# IRF for a given model, shock and variable
irf(shock::Symbol, variable::Symbol; thiscontext = context, model = 1) =
    getproperty(thiscontext.results.model_results[model].irfs[shock], variable)
irf(shock::String, variable::String; thiscontext = context, model = 1) =
    irf(Symbol(shock), Symbol(variable); thiscontext = context, model = model)

## SIMULATIONS
# all simulations for a model
simulation(; thiscontext = context, model = 1) =
    thiscontext.results.model_results[model].simulations
# all simulated variables for a given model and a given simulation
simulation(simnbr::Int64; thiscontext = context, model = 1) =
    thiscontext.results.model_results[model].simulations[simnbr].data
# one simulated variable for a given model and a given simulation
simulation(varname::Symbol; thiscontext = context, model = 1, simnbr = 1) =
    getproperty(thiscontext.results.model_results[model].simulations[simnbr].data, varname)
simulation(varname::String; thiscontext = context, model = 1, simnbr = 1) =
    simulation(Symbol(varname); thiscontext = context, model = model, simnbr = simnbr)
# subset of simulated variables for a  given model and a given simulation
function simulation(varnames::Vector{Symbol}; thiscontext = context, model = 1, simnbr = 1)
    tdf = thiscontext.results.model_results[model].simulations[simnbr].data
    TimeDataFrame(select(dataframe(tdf), varnames), periods(tdf), iscontinuous(tdf))
end
simulation(varnames::Vector{String}; thiscontext = context, model = 1, simnbr = 1) =
    simulation([Symbol(v) for v in varnames]; thiscontext = context, model = model, simnbr = simnbr)
simulation(varnames::Tuple; thiscontext = context, model = 1, simnbr = 1) =
    simulation([Symbol(v) for v in varnames]; thiscontext = context, model = model, simnbr = simnbr)
