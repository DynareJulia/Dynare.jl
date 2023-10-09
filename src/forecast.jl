@enum ForecastModes histval calibsmoother estimation

function forecasting!(; periods, forecast_mode, context=context, datafile="", order=1)
    model = context.models[1]
    nendo = model.endogenous_nbr
    nexo = model.exogenous_nbr
    results = context.results.model_results[1]
    Y = Matrix{Float64}(undef, periods + 1, nendo)
    c = results.trends.endogenous_steady_state
    A = zeros(model.endogenous_nbr, model.endogenous_nbr)
    B = zeros(model.endogenous_nbr, model.exogenous_nbr)
    if forecast_mode == histval
        y0 = results.trends.endogenous_steady_state
        for (i, v) in enumerate(context.work.histval)
            !ismissing(v) && (y0[i] = v)
        end
    elseif forecast_mode == calibsmoother
        calibsmoother!(context=context, datafile=datafile, first_obs=firstobs, last_obs=lastobs)
        @views y0 = results.smoother[:,end]
    end 
        
    make_A_B!(A, B, model, results)
    if order == 1
        forecast_!(Y, y0, c, A, B, periods)
    end
    results.forecast = [AxisArrayTable(Y, 
                                       Undated(0):Undated(periods), 
                                       [Symbol(v) for v in get_endogenous(context.symboltable)])]
end

function recursive_forecast!
end

function conditional_forecast!
end

function recursive_conditional_forecast!
end
    

function forecast_!(Y::AbstractMatrix, y0::AbstractVector, c::AbstractVector,
    A::AbstractMatrix, B::AbstractMatrix, periods::Integer)
    x = zeros(periods + 1, size(B, 2))
    simul_first_order!(Y, y0, x, c, A, B, periods)
    return Y
end