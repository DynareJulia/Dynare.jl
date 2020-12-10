function is_jacobian_full_rank(jacobian)
    return rank(jacobian) == size(jacobian, 1)
end
