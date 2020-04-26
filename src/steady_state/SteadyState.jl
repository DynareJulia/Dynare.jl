
function steady_state!(context)
    model = context.models[1]
    work = context.work
    steadystatemodule = model.steady_state! 
    if isnothing(steadystatemodule)
    else
        results = context.results.model_results[1]
        fill!(results.exogenous_steady_state, 0.0)
        steadystatemodule.steady_state!(results.endogenous_steady_state,
                                        results.exogenous_steady_state,
                                        work.params)
    end
end
	
