using Test

include("../src/estimation/smc.jl")

function data_generation!(y, ϕ1, ϕ2, ϕ3)
    s1 = 0
    s2 = 0
    for i in eachindex(y)
        s2 = ϕ2*s1 + ϕ3*s2
        s1 = ϕ1*s1 + rand()
        y[i] = s1 + s2
    end
    return y
end

initial_P(ϕ) = 1/(1 - ϕ[1]^2)*[1 ϕ[1]*ϕ[2]/(1 - ϕ[1]*ϕ[3]); ϕ[1]*ϕ[2]/(1 - ϕ[1]*ϕ[3])  ϕ[2]^2*(1 + 2*ϕ[1]*ϕ[3]/(1 - ϕ[1]*ϕ[3]))/(1 - ϕ[3]^2)]

function loglikelihood(y, ϕ1, ϕ2, ϕ3)
    s = zeros(2)
    for i in eachindex(y)
        ν = y[i] - s[1] - s[2]
    end
end

function reduced_form!(ϕ, θ)
    θ1sq = θ[1]*θ[1]
    ϕ[1] = θ1sq
    ϕ[2] = 1 - θ1sq - θ[1]*θ[2]
    ϕ[3] = 1 - θ1sq
    return ϕ
end

θ = [0.45, 0.45]
ϕ = zeros(3)

reduced_form!(ϕ, θ)

T(ϕ) = [ϕ[1] 0; ϕ[2] ϕ[3]]

P0 = initial_P(ϕ)
@test P0 ≈ T(ϕ)*P0*T(ϕ)' + [1 0; 0 0]

