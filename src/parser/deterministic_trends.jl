struct Trends
    endogenous_steady_state::Vector{Float64}
    endogenous_linear_trend::Vector{Float64}
    endogenous_quadratic_trend::Vector{Float64}
    exogenous_steady_state::Vector{Float64}
    exogenous_linear_trend::Vector{Float64}
    exogenous_quadratic_trend::Vector{Float64}
    exogenous_det_steady_state::Vector{Float64}
    exogenous_det_linear_trend::Vector{Float64}
    exogenous_det_quadratic_trend::Vector{Float64}
    function Trends(ny, nx, nxd)
        endogenous_steady_state = Vector{Float64}(undef, ny)
        endogenous_linear_trend = Vector{Float64}(undef, ny)
        endogenous_quadratic_trend = Vector{Float64}(undef, ny)
        exogenous_steady_state = Vector{Float64}(undef, nx)
        exogenous_linear_trend = Vector{Float64}(undef, nx)
        exogenous_quadratic_trend = Vector{Float64}(undef, nx)
        exogenous_det_steady_state = Vector{Float64}(undef, nxd)
        exogenous_det_linear_trend = Vector{Float64}(undef, nxd)
        exogenous_det_quadratic_trend = Vector{Float64}(undef, nxd)
        new(endogenous_steady_state, endogenous_linear_trend,
            endogenous_quadratic_trend, exogenous_steady_state,
            exogenous_linear_trend, exogenous_quadratic_trend,
            exogenous_det_steady_state, exogenous_det_linear_trend,
            exogenous_det_quadratic_trend)
    end
end

function remove_linear_trend!(data_out, data_in, steady_state, linear_trend_coeffs; row = 1)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    data_out .= data_in .- steady_state .- linear_trend_coeffs.*transpose(linear_trend)
end
    
function add_linear_trend!(data_out, data_in, steady_state, linear_trend_coeffs; row = 1)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    data_out = data_in .+ steady_state .+ linear_trend_coeffs.*transpose(linear_trend)
end
    
function remove_quadratic_trend!(data_out, data_in, steady_state, linear_trend_coeffs,
                                 quadratic_trend_coeffs; row = 1)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    quadratic_trend = collect(row - 1 .+ (1:n)).^2
    
    data_out = (data_in .- steady_state .- linear_trend_coeffs.*transpose(linear_trend)
                .- quadratic_trend_coeffs.*transpose(quadratic_trend))
end
    
function add_quadratic_trend!(data_out, data_in, steady_state, linear_trend_coeffs,
                              quadratic_trend_coeffs; row = 1)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    quadratic_trend = collect(row - 1 .+ (1:n)).^2
    
    data_out = (data_in .+ steady_state .+ linear_trend_coeffs.*transpose(linear_trend)
                .+ quadratic_trend_coeffs.*transpose(quadratic_trend))
end
    

    
    
