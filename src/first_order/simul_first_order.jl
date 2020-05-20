using LinearAlgebra

function simul_first_order!(results, initial_values, x, c, A, B, periods)
    y_1 = initial_values
    for t = 1:periods
        r = view(results, t, :)
        e = view(x, t, :)
        mul!(r, B, e)
        mul!(r, A, y_1, 1.0, 1.0)
        r .=  c
    end
end
