
function compute_steady_state!(context::Context)
    model = context.models[1]
    work = context.work
    steadystatemodule = model.steady_state! 
    results = context.results.model_results[1]
    fill!(results.trends.exogenous_steady_state, 0.0)
    if isempty(steadystatemodule)
        fill!(results.trends.endogenous_steady_state, 0.0)    
    else
        Base.invokelatest(steadystatemodule.steady_state!,
                          results.trends.endogenous_steady_state,
                          results.trends.exogenous_steady_state,
                          work.params)
    end
end
	
