using Roots
using SpecialFunctions

function beta_specification(μ, σ2; lb = 0, ub = 1, name = "")

    if μ < lb
        error(
            "Prior for $(name): the expectation $(μ) cannot be smaller than $(lb), the lower bound of the Beta distribution!",
        )
    end

    if μ > ub
        error(
            "Prior for $(name): the expectation $(μ) cannot be greater than $(ub), the upper bound of the Beta distribution!",
        )
    end

    len = ub - lb

    μ = (μ - lb) / len
    σ2 = σ2 / (len * len)

    if σ2 > (1 - μ) * μ
        error(
            "Beta prior for $(name): given the declared prior expectation, prior lower and upper bounds, the prior std. has to be smaller than $(sqrt((1-μ)*μ))!",
        )
    end

    α = (1 - μ) * μ * μ / σ2 - μ
    β = α * (1 / μ - 1)
    return α, β
end

function gamma_specification(μ, σ2)
    θ = μ / σ2
    α = μ * θ
    return α, θ
end

function inverse_gamma_1_specification(μ, σ2)
    μ2 = μ * μ
    f(x) = (x - 1) * (σ2 + μ2) - μ2 * gamma(x) / gamma(x - 0.5)
    α = find_zero(f, (1, 100))
    θ = (α - 1) * (σ2 + μ2)
    return α, θ
end

function inverse_gamma_2_specification(μ, σ2)
    μ2 = μ * μ
    α = 2 + μ2 / σ2
    θ = (α - 1) * μ
    return α, θ
end


function uniform_specification(μ, σ2)
    b = sqrt(12 * σ2) / 2 + μ
    a = 2 * μ - b
    return a, b
end

function weibull_specification(μ, σ2)
    f(x) = σ2 + μ * μ * (1 + gamma(1 + 2 / x) / gamma(1 + 1 / x)^2)
    α = find_zero(f, (0, 100))
    θ = μ / gamma(1 + 1 / α)
    return α, θ
end
