function deterministic_trends!(context::Context, field::Dict{String,Any})
    model = context.models[1]
    work = context.work
    trend = context.results.model_results[1].trends.endogenous_linear_trend
    symboltable = context.symboltable
    work.model_has_trend[1] = true
    fill!(trend, 0.0)
    for (key, value) in field["trends"]
        trend[symboltable[key].orderintype] = work.params[symboltable[value].orderintype]
    end
end

function remove_linear_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    row::Int64 = 1,
)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    data_out .= data_in .- steady_state .- linear_trend_coeffs .* transpose(linear_trend)
end

function add_linear_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64};
    row::Int64 = 1,
)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    data_out .= data_in .+ steady_state .+ linear_trend_coeffs .* transpose(linear_trend)
end

function remove_quadratic_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64};
    row = 1,
)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    quadratic_trend = collect(row - 1 .+ (1:n)) .^ 2

    data_out = (
        data_in .- steady_state .- linear_trend_coeffs .* transpose(linear_trend) .-
        quadratic_trend_coeffs .* transpose(quadratic_trend)
    )
end

function add_quadratic_trend!(
    data_out::Any,
    data_in::Any,
    steady_state::AbstractVector{Float64},
    linear_trend_coeffs::AbstractVector{Float64},
    quadratic_trend_coeffs::AbstractVector{Float64};
    row::Int64 = 1,
)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    quadratic_trend = collect(row - 1 .+ (1:n)) .^ 2

    data_out = (
        data_in .+ steady_state .+ linear_trend_coeffs .* transpose(linear_trend) .+
        quadratic_trend_coeffs .* transpose(quadratic_trend)
    )
end
