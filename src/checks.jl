function is_jacobian_full_rank(jacobian::AbstractMatrix{Float64})
    return rank(jacobian) == size(jacobian, 1)
end
