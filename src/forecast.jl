@enum ForecastModes histval calibsmoother estimation

"""
    forecasting!(; periods::Integer,
                   forecast_mode::ForecastModes,
                   context::Context=context,
                   datafile::String="",
                   first_obs::PeriodsSinceEpoch=Undated(typemin(Int)),
                   first_period::PeriodsSinceEpoch=Undated(0),
                   last_obs::PeriodsSinceEpoch=Undated(typemin(Int)),
                   order::Integer=1)
computes an unconditional forecast of the variables of the model

# Keyword arguments
- `periods::Integer`: number of forecasted periods [required]
- `forecast_mode::ForecastModes`: one of `histval` or `calibsmoother` [required]
- `datafile::String`: file with the observations for the smoother
- `first_obs::PeriodsSinceEpoch`: first period used by smoother (default: first observation in the file)  
- `first_period::PeriodsSinceEpoch`: initial_period for the forecast (default when histval: Undated(0),
                                                                     default when calibsmoother: last period of the smoother)
- `last_obs::PeriodsSinceEpoch`: last period used by smoother  (default: last observation in the file)
- `order::Integer`: order of local approximation
"""
function forecasting!(; periods::Integer,
                      forecast_mode::ForecastModes,
                      context::Context=context,
                      data::AxisArrayTable=AxisArrayTable(Matrix{Float64}(undef, 0, 0), PeriodsSinceEpoch[], Symbol[]),
                      datafile::String="",
                      first_obs::PeriodsSinceEpoch=Undated(typemin(Int)),
                      last_obs::PeriodsSinceEpoch=Undated(typemin(Int)),
                      nobs::Int=0,
                      first_period::PeriodsSinceEpoch=Undated(0),
                      order::Integer=1)
    results = context.results.model_results[1]
    Y = forecasting_(periods = periods, 
                    forecast_mode = forecast_mode, 
                    context = context, 
                    data = data, 
                    datafile = datafile, 
                    first_obs=first_obs,
                    last_obs=last_obs,
                    nobs = nobs, 
                    first_period = first_period, 
                    order=order)
    if (forecast_mode == calibsmoother) && first_period == 0
        first_period = size(results.smoother, 1)
    end                     
    results.forecast = [AxisArrayTable(Y, 
                        first_period .+ (0:Undated(periods)), 
                        [Symbol(v) for v in get_endogenous(context.symboltable)])]
    return results.forecast
end

function forecasting_(; periods, 
    forecast_mode::ForecastModes, 
    context=context,
    data = AxisArrayTable(Matrix{Float64}(undef, 0, 0), PeriodsSinceEpoch[], Symbol[]), 
    datafile = "", 
    first_obs::PeriodsSinceEpoch = Undated(typemin(Int)), 
    last_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
    nobs::Int = 0, 
    first_period = 0, 
    order = 1)
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
            calibsmoother!(context=context, data = data, datafile=datafile, first_obs=first_obs, last_obs=last_obs, nobs = 0)
            # take only endogenous variables
            @views y0 = Vector(results.smoother[end, 1:nendo])
        end 
        
        make_A_B!(A, B, model, results)
        if order == 1
            forecast_!(Y, y0, c, A, B, periods)
        end
    end 
end

"""
    recursive_forecasting!(; Np::Integer,
                           first_period::PeriodsSinceEpoch=Undated(typemin(Int)), 
                           last_period::PeriodsSinceEpoch=Undated(typemax(Int)), 
                           context::Context=context,
                           datafile::String="", 
                           first_obs::PeriodsSinceEpoch=Undated(typemin(Int)), 
                           last_obs::PeriodsSinceEpoch=Undated(typemax(Int)), 
                           order::Integer=1)
computes an unconditional recursive forecast by adding one period to the
sample used for the smoother before forecasting over `Np` periods.

# Keyword arguments
- `Np::Integer`: number of forecasted periods [required]
- `first_period::PeriodsSinceEpoch`: initial period of first forecast [required]
- `last_period::PeriodsSinceEpoch`: initial period of last forecast [required]
- `datafile::String`: file with the observations for the smoother
- `first_obs::PeriodsSinceEpoch`: first period used by smoother (default: first observation in the file)  
- `last_obs::PeriodsSinceEpoch`: last period used by smoother  (default: last observation in the file)
- `order::Integer`: order of local approximation
"""
function recursive_forecasting!(; Np::Integer,
                                first_period::PeriodsSinceEpoch=Undated(typemin(Int)), 
                                last_period::PeriodsSinceEpoch=Undated(typemin(Int)), 
                                context::Context=context,
                                data::AxisArrayTable = AxisArrayTable(Matrix{Float64}(undef, 0, 0), PeriodsSinceEpoch[], Symbol[]),
                                datafile::String="",
                                variables::Vector{<:Union{Symbol, String}} = Union{Symbol, String}[], 
                                first_obs::PeriodsSinceEpoch=Undated(typemin(Int)), 
                                last_obs::PeriodsSinceEpoch=Undated(typemin(Int)),
                                nobs::Int = 0, 
                                order::Integer=1)
    #@assert last_period <= last_obs
    results = context.results.model_results[1]
    data_ = get_data!(context, datafile, data, variables, first_obs, last_obs, nobs)
    first_period == Undated(typemin(Int)) && (first_period = row_labels(data_)[1]) 
    last_period == Undated(typemin(Int)) && (last_period = row_labels(data_)[end]) 
    T = typeof(first_period).parameters[1]
    empty!(results.forecast)
    for p = first_period:last_period
        Y = forecasting_(context=context, periods=Np, forecast_mode=calibsmoother, first_obs=first_obs, last_obs=p, data = data_, order=order)
        if p == first_period
            results.initial_smoother = copy(results.smoother)
        end
        p1 = p + T(Np)
        push!(results.forecast, AxisArrayTable(Y, 
                                               p:p1, 
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
