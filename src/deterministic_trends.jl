function remove_linear_trend!(data_out, data_in, steady_state, linear_trend_coeffs; row = 1)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    data_out .= data_in .- steady_state .- linear_trend_coeffs.*transpose(linear_trend)
end
    
function add_linear_trend!(data_out, data_in, steady_state, linear_trend_coeffs; row = 1)
    n = size(data_in, 2)
    linear_trend = collect(row - 1 .+ (1:n))
    data_out .= data_in .+ steady_state .+ linear_trend_coeffs.*transpose(linear_trend)
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
    

    
    
