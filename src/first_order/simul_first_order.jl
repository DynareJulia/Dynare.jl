using LinearAlgebra

function simul_first_order!(results, initial_values, x, c, A, B, periods)
    r_1 = view(results, 1, :)
    r_1 .= initial_values .- c
    for t = 2:periods + 1
        r = view(results, t, :)
        e = view(x, t, :)
        mul!(r, B, e)
        mul!(r, A, r_1, 1.0, 1.0)
        r_1 .+=  c
        r_1 = r
    end
    r_1 .+= c
end
