function beta_specification(μ, σ2; lb = 0, ub = 1, name = "")
    
    if μ < lb
        error("Prior for $(name): the expectation $(μ) cannot be smaller than $(lb), the lower bound of the Beta distribution!")
    end
    
    if μ > ub
        error("Prior for $(name): the expectation $(μ) cannot be greater than $(ub), the upper bound of the Beta distribution!")
    end
    
    len = ub - lb
    
    μ = (μ - lb)/len
    σ2 = σ2/(len*len)
    
    if σ2 > (1 - μ)*μ
        error("Beta prior for $(name): given the declared prior expectation, prior lower and upper bounds, the prior std. has to be smaller than $(sqrt((1-μ)*μ))!")
    end
    
    α = (1-μ)*μ*μ/σ2 - μ
    β = α*(1/μ - 1) 
    return α, β
end


function gamma_specification(μ, σ2)
    θ = μ/σ2
    α = μ*θ
    return α, θ
end
