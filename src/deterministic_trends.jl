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

function add_linear_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0
)
    linear_trend = start
    @inbounds for i in axes(data_out, 2)
        for j in axes(data_out, 1)
            data_out[j, i] = data_in[j, i] + steady_state[j] + linear_trend_coeffs[j] * linear_trend
        end
        linear_trend += 1
    end
end

function remove_linear_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0
)
    linear_trend = start
    @inbounds for i in axes(data_out, 2)
        for j in axes(data_out, 1)
            data_out[j, i] = data_in[j, i] - steady_state[j] - linear_trend_coeffs[j] * linear_trend
        end
        linear_trend += 1
    end
end

function add_quadratic_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0
)
    trend = start
    @inbounds for i in axes(data_out, 2)
        for j in axes(data_out, 1)
            data_out[j, i] = data_in[j, i] + steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2
        end
        linear_trend += 1
    end
end

function remove_quadratic_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0
)
    trend = start
    @inbounds for i in axes(data_out, 2)
        for j in axes(data_out, 1)
            data_out[j, i] = data_in[j, i] - steady_state[j] - linear_trend_coeffs[j] * trend - quadratic_trend_coeffs[j] * trend^2
        end
        linear_trend += 1
    end
end

function add_linear_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0
)
    linear_trend = start
    @inbounds for i in axes(data, 2)
        for j in axes(data, 1)
            data[j, i] += steady_state[j] + linear_trend_coeffs[j] * linear_trend 
        end
        linear_trend += 1
    end
end

function remove_linear_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0
)
    linear_trend = start
    @inbounds for i in axes(data, 2)
        for j in axes(data, 1)
            data[j, i] -= steady_state[j] + linear_trend_coeffs[j] * linear_trend
        end
        linear_trend += 1
    end
end

function add_quadratic_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64};
    start::Int64 = 0,
)
    trend = start
    @inbounds for i in axes(data, 2)
        for j in axes(data, 1)
            data[j, i] += steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2
    end
    trend += 1
end

function remove_quadratic_trend!(
    data::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64};
    start = 0,
)
    trend = start
    @inbounds for i in axes(data, 2)
        for j in axes(data, 1)
            data[j, i] -= steady_state[j] + linear_trend_coeffs[j] * trend + quadratic_trend_coeffs[j] * trend^2
        end
        trend += 1
    end
end

end
