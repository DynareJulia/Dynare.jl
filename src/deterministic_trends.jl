function deterministic_trends!(context::Context, field::Dict{String,Any})
    model = context.models[1]
    work = context.work
    trend = context.results.model_results[1].trends.endogenous_linear_trend
    symboltable = context.symboltable
    work.model_has_trend[1] = true
    isempty(trend) && resize!(trend, model.endogenous_nbr)
    fill!(trend, 0.0)
    for (key, value) in field["trends"]
        trend[symboltable[key].orderintype] = work.params[symboltable[value].orderintype]
    end
end

function add_steady_state!(dataout::Any,
                           datain::Any,
                           steady_state::AbstractVector{Float64},
                           dim;
                           start::Int64 = 0,
                           )
    @inbounds for i in axes(datain, 2)
        if dim == 1
            for j in axes(datain, 1)
                dataout[j, i] = datain[j] + steady_state[j]
            end
        else
            for j in axes(datain, 1)
                dataout[j, i] = datain[i] + steady_state[i]
            end
        end
    end
    return dataout
end

function remove_steady_state!(dataout::Any,
                              datain::Any,
                              steady_state::AbstractVector{Float64},
                              dim;
                              start::Int64 = 0,
                              )
    @inbounds for i in axes(datain, 2)
        trend = start
        if dim == 1
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] - steady_state[j]
            end
        else
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] - steady_state[i]
            end
        end
    end
    return dataout
end

function add_linear_trend!(
    dataout::Any,
    datain::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    dim;
    start::Int64 = 0,
)
    linear_trend = start
    @inbounds for i in axes(datain, 2)
        if dim == 1
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] + steady_state[j] + linear_trend_coeffs[j] * linear_trend 
            end
            linear_trend += 1
        else
            linear_trend = start
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] + steady_state[i] + linear_trend_coeffs[i] * linear_trend 
                linear_trend += 1
            end
        end
    end
    return dataout
end

function remove_linear_trend!(
    dataout::Any,
    datain::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    dim;
    start::Int64 = 0, 
)
    linear_trend = start
    @inbounds for i in axes(datain, 2)
        if dim == 1
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] - (steady_state[j] + linear_trend_coeffs[j] * linear_trend)
            end
            linear_trend += 1
        else
            linear_trend = start
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] - (steady_state[i] + linear_trend_coeffs[i] * linear_trend)
                linear_trend += 1
            end
        end
    end
    return dataout
end

function add_quadratic_trend!(
    dataout::Any,
    datain::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64},
    dim;
    start::Int64 = 0,
)
    trend = start
    @inbounds for i in axes(datain, 2)
        if dim == 1
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] + steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2
            end
            trend += 1
        else
            trend = start
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] + steady_state[i] + linear_trend_coeffs[i] * trend + quadratic_trend_coeffs[i] * trend^2
                trend += 1
            end
        end
    end
end

function remove_quadratic_trend!(
    dataout::Any,
    datain::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64},
    dim;
    start = 0,
)
    trend = start
    @inbounds for i in axes(datain, 2)
        if dim == 1
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] - (steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2)
            end
            trend += 1
        else
            trend = start
            for j in axes(datain, 1)
                dataout[j, i] = datain[j, i] - (steady_state[i] + linear_trend_coeffs[i] * trend + quadratic_trend_coeffs[i] * trend^2)
                trend += 1
            end
        end
    end
    return dataout
end

function add_steady_state!(data::Any,
                           steady_state::AbstractVector{Float64},
                           dim;
                           start::Int64 = 0,
                           )
    @inbounds for i in axes(data, 2)
        if dim == 1
            for j in axes(data, 1)
                data[j, i] += steady_state[j]
            end
        else
            for j in axes(data, 1)
                data[j, i] += steady_state[i]
            end
        end
    end
    return data
end

function remove_steady_state!(data::Any,
                              steady_state::AbstractVector{Float64},
                              dim;
                              start::Int64 = 0,
                              )
    @inbounds for i in axes(data, 2)
        trend = start
        if dim == 1
            for j in axes(data, 1)
                data[j, i] -= steady_state[j]
                trend += 1
            end
        else
            for j in axes(data, 1)
                data[j, i] -= steady_state[i]
                trend += 1
            end
        end
    end
    return data
end

function add_linear_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    dim;
    start::Int64 = 0,
)
    linear_trend = start
    @inbounds for i in axes(data, 2)
        if dim == 1
            for j in axes(data, 1)
                data[j, i] += steady_state[j] + linear_trend_coeffs[j] * linear_trend 
            end
            linear_trend += 1
        else
            linear_trend = start
            for j in axes(data, 1)
                data[j, i] += steady_state[i] + linear_trend_coeffs[i] * linear_trend 
                linear_trend += 1
            end
        end
    end
end

function remove_linear_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    dim;
    start::Int64 = 0,
)
    linear_trend = start
    @inbounds for i in axes(data, 2)
        if dim == 1
            for j in axes(data, 1)
                data[j, i] -= steady_state[j] + linear_trend_coeffs[j] * linear_trend
            end
            linear_trend += 1
        else
            linear_trend = start
            for j in axes(data, 1)
                data[j, i] -= steady_state[i] + linear_trend_coeffs[i] * linear_trend
                linear_trend += 1
            end
        end
    end
end

function add_quadratic_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64},
    dim;
    start::Int64 = 0,
)
    trend = start
    @inbounds for i in axes(data, 2)
        if dim == 1
            for j in axes(data, 1)
                data[j, i] += steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2
            end
            trend += 1
        else
            trend = start
            for j in axes(data, 1)
                data[j, i] += steady_state[i] + linear_trend_coeffs[i] * trend + quadratic_trend_coeffs[i] * trend^2
                trend += 1
            end
        end
    end
end

function remove_quadratic_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64},
    dim;
    start = 0,
)
    trend = start
    @inbounds for i in axes(data, 2)
        if dim == 1
            for j in axes(data, 1)
                data[j, i] -= steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2
            end
            trend += 1
        else
            trend = start
            for j in axes(data, 1)
                data[j, i] -= steady_state[i] + linear_trend_coeffs[i] * trend + quadratic_trend_coeffs[i] * trend^2
                trend += 1
            end
        end
    end
end
