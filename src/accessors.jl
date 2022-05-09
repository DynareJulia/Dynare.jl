## IRFs
# all the IRFs for a given model
irf(; context=context, model=1) = context.results.model_results[model].irfs
# all the IRFs for a given model and a given shock
irf(shock::Symbol; context=context, model=1) = context.results.model_results[model].irfs[shock]
irf(shock::String; context=context, model=1) = irf(Symbol(shock); context=context, model=model) 
# IRF for a given model, shock and variable
irf(shock::Symbol, variable::Symbol; context=context, model=1) = getproperty(context.results.model_results[model].irfs[shock], variable)
irf(shock::String, variable::String; context=context, model=1) = irf(Symbol(shock), Symbol(variable); context=context, model=model) 

## SIMULATIONS
# all simulations for a model
simulation(; context=context, model=1) = context.results.model_results[model].simulations
# all simulated variables for a given model and a given simulation
simulation(simnbr::Int64; context=context, model=1) = context.results.model_results[model].simulations[simnbr].data
# one simulated variable for a given model and a given simulation
simulation(varname::Symbol; context=context, model=1, simnbr=1) = getproperty(context.results.model_results[model].simulations[simnbr].data, varname)
simulation(varname::String; context=context, model=1, simnbr=1) = simulation(Symbol(varname); context=context, model=model, simnbr=simnbr)
# subset of simulated variables for a  given model and a given simulation
function simulation(varnames::Vector{Symbol}; context=context, model=1, simnbr=1)
    tdf = context.results.model_results[model].simulations[simnbr].data
    TimeDataFrame(select(dataframe(tdf), varnames),
                  periods(tdf),
                  iscontinuous(tdf))
end
simulation(varnames::Vector{String}; context=context, model=1, simnbr=1) = simulation([Symbol(v) for v in varnames]; context, model=model, simnbr=simnbr)
simulation(varnames::Tuple; context=context, model=1, simnbr=1) = simulation([Symbol(v) for v in varnames]; context, model=model, simnbr=simnbr)

