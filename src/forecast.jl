@enum ForecastModes histval calibsmoother estimation

function forecast!(; periods, forecast_mode, context=context, datafile="", order=1)
    nendo = context.models[1].endogenous_nbr
    nexo = context.models[1].exogenous_nbr
    Y = Matrix{Float64}(nendo, periods)
    c = context.results.model_results[1].trends.endogenous_steady_state
    A = zeros(model.endogenous_nbr, model.endogenous_nbr)
    B = zeros(model.endogenous_nbr, model.exogenous_nbr)
    if forecast_mode == histval
        y0 = histval[end,:]
    elseif forecast_mode = calibsmoother
        
    make_A_B!(A, B, model, results)
    if order == 1
        foreccast_!(Y, y0, c, A, B, periods)
    end
end
)

function recursive_forecast!
end

function conditional_forecast!
end

function recursive_conditional_forecast!
    end
    

function forecast_!(Y::AbstractMatrix, y0::AbstracVector, c::AbstractVector,
    A::AbstractMatrix, B::AbstractMatrix, periods::Integer)
    x = zeros(size(B, 2), periods)
    simul_first_order!(Y, y0, x, c, A, B, periods)
    return Y
end
