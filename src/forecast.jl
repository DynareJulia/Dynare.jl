@enum ForecastModes histval calibsmoother estimation

"""
    forecasting!(; periods::Integer,
                   forecast_mode::ForecastModes,
                   context::Context=context,
                   datafile::String="",
                   first_obs::PeriodsSinceEpoch=1,
                   last_obs::PeriodsSinceEpoch=0,
                   first_period::PeriodsSinceEpoch=0,
                   order::Integer=1)
computes an unconditional forecast of the variables of the model

# Keyword arguments
- `periods::Integer`: number of forecasted periods [required]
- `forecast_mode::ForecastModes`: one of `histval` or `calibsmoother` [required]
- `datafile::String`: file with the observations for the smoother
- `first_obs::PeriodsSinceEpoch`: first period used by smoother (default: first observation in the file)  
- `last_obs::PeriodsSinceEpoch`: last period used by smoother  (default: last observation in the file)
- `first_period::PeriodsSinceEpoch`: initial_period for the forecast (default when histval: Undated(0),
                                                                     default when calibsmoother: last period of the smoother)
- `order::Integer`: order of local approximation
"""
function forecasting!(; periods::Integer,
                      forecast_mode::ForecastModes,
                      context::Context=context,
                      datafile::String="",
                      first_obs::PeriodsSinceEpoch=Undated(1),
                      last_obs::PeriodsSinceEpoch=Undated(0),
                      first_period::PeriodsSinceEpoch=Undated(0),
                      order::Integer=1)
    results = context.results.model_results[1]
    Y = forecasting_(; periods, forecast_mode, context=context, datafile=datafile, 
                       first_obs=first_obs, last_obs=last_obs, order=order)
    if (forecast_mode == calibsmoother) && first_period == 0
        first_period = size(results.smoother, 1)
    end                     
    results.forecast = [AxisArrayTable(Y, 
                        first_period .+ (0:Undated(periods)), 
                        [Symbol(v) for v in get_endogenous(context.symboltable)])]
    return results.forecast
end

function forecasting_(; periods, forecast_mode::ForecastModes, context=context, datafile="", first_obs=1, last_obs=0, first_period = 0, order=1)
    model = context.models[1]
    nendo = model.endogenous_nbr
    nexo = model.exogenous_nbr
    results = context.results.model_results[1]
    Y = Matrix{Float64}(undef, periods + 1, nendo)
    if !isempty(results.trends.endogenous_linear_trend)
        c = results.trends.endogenous_steady_state .+ 
            results.trends.endogenous_linear_trend.*transpose(0:periods)
    else
        c = results.trends.endogenous_steady_state
    end
    A = zeros(model.endogenous_nbr, model.endogenous_nbr)
    B = zeros(model.endogenous_nbr, model.exogenous_nbr)
    let y0
        if forecast_mode == histval
            y0 = copy(results.trends.endogenous_steady_state)
            for (i, v) in enumerate(context.work.histval[end,:])
                !ismissing(v) && (y0[i] = v)
            end
        elseif forecast_mode == calibsmoother
            calibsmoother!(context=context, datafile=datafile, first_obs=first_obs, last_obs=last_obs)
            @views y0 = Vector(results.smoother[end, :])
        end 
        
        make_A_B!(A, B, model, results)
        if order == 1
            forecast_!(Y, y0, c, A, B, periods)
        end
    end 
end

"""
function recursive_forecasting!(; periods::Integer,
                                first_period::PeriodsSinceEpoch, 
                                last_period::PeriodsSinceEpoch, 
                                context::Context=context,
                                datafile::String="", 
                                first_obs::PeriodsSinceEpoch=Undated(1), 
                                last_obs::PeriodsSinceEpoch=Undated(0), 
                                order::Integer=1)
computes an unconditional recursive forecast for one variable

# Keyword arguments
- `periods::Integer`: number of forecasted periods [required]
- `first_period::PeriodsSinceEpoch`: initial period of first forecast [required]
- `last_period::PeriodsSinceEpoch`: initial period of last forecast [required]
- `datafile::String`: file with the observations for the smoother
- `first_obs::PeriodsSinceEpoch`: first period used by smoother (default: first observation in the file)  
- `last_obs::PeriodsSinceEpoch`: last period used by smoother  (default: last observation in the file)
- `order::Integer`: order of local approximation
"""
function recursive_forecasting!(; periods::Integer,
                                first_period::PeriodsSinceEpoch, 
                                last_period::PeriodsSinceEpoch, 
                                context::Context=context,
                                datafile::String="", 
                                first_obs::PeriodsSinceEpoch=Undated(1), 
                                last_obs::PeriodsSinceEpoch=Undated(0), 
                                order::Integer=1)
    #@assert last_period <= last_obs
    results = context.results.model_results[1]
    empty!(results.forecast)
    for p = first_period:last_period
        Y = forecasting_(context=context, periods=periods, forecast_mode=calibsmoother, first_obs=first_obs, last_obs=p, datafile=datafile, order=order)
        if p == first_period
            results.initial_smoother = copy(results.smoother)
        end 
        push!(results.forecast, AxisArrayTable(Y, 
                                               Undated(p):Undated(p + periods), 
                                               [Symbol(v) for v in get_endogenous(context.symboltable)]))
    end
end

function conditional_forecast!
end

function recursive_conditional_forecast!
end
    

function forecast_!(Y::AbstractMatrix, y0::AbstractVector, c::AbstractVecOrMat,
    A::AbstractMatrix, B::AbstractMatrix, periods::Integer)
    x = zeros(periods + 1, size(B, 2))
    simul_first_order!(Y, y0, x, c, A, B, periods)
    return Y
end
