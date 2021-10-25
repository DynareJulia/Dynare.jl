
function compute_steady_state!(context::Context)
    model = context.models[1]
    results = context.results.model_results[1]
    work = context.work
    steadystatemodule = model.steady_state!
    if Symbol("steady_state!") in names(steadystatemodule)
        # explicit steady state
        evaluate_steady_state!(results, steadystatemodule, work.params)
    end
end
	
function evaluate_steady_state!(results::ModelResults,
                                static_module::Module,
                                params::AbstractVector{Float64})
    fill!(results.trends.exogenous_steady_state, 0.0)
    Base.invokelatest(static_module.steady_state!,
                      results.trends.endogenous_steady_state,
                      results.trends.exogenous_steady_state,
                      params)
end

function solve_steady_state!(results::ModelResults,
                             static_module::Module,
                             params::AbstractVector{Float64})
    static_function(x::AbstractVector{Float64}) =
        Base.invokelatest(static_module.static!,
                          x,
                          results.trends.exogenous_steady_state,
                          params)
    static_jacobian(x::AbstractVector{Float64}) =
        Base.invokelatest(static_module.static!,
                          x,
                          results.trends.exogenous_steady_state,
                          params)
    
    results.trends.endogenous_steady_state = trustregion(static_function,
                                                         static_jacobian,
                                                         x0::Vector{Float64})
end
