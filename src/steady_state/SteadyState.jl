using NLsolve

mutable struct DynareSteadyStateComputationFailed <: Exception end
Base.showerror(io::IO, e::DynareSteadyStateComputationFailed) =
    print(
"""
Dynare couldn't compute the steady state.
Either there is no solution or the guess values
are too far from the solution
"""
    )

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

function solve_steady_state!(context::Context,
                             x0::Vector{Float64})

    ws = StaticWs(context)
    m = context.models[1]
    w = context.work
    results = context.results.model_results[1]
    exogenous = results.trends.exogenous_steady_state
    
    residual_function(x::AbstractVector{Float64}) =
        get_static_residuals!(ws, w.params, x, exogenous, m)

    jacobian_function(x::AbstractVector{Float64}) =
        get_static_jacobian!(ws, w.params, x, exogenous, m)

    residual_function(x0)
    @show "OK1"
    jacobian_function(x0)
    @show "OK2"
    result = nlsolve(residual_function,
                     jacobian_function,
                     x0::Vector{Float64})
    if converged(result)
        results.trends.endogenous_steady_state .= result.zero
    else
        @debug "Steady state computation failed with\n $result"
        throw(DynareSteadyStateComputationFailed)
    end
end
