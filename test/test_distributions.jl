using Distributions
using Test

include("../src/distribution_parameters.jl")
include("../src/distributions/inversegamma1.jl")

ig1 = InverseGamma1(1.3, 2.5)
m = mean(ig1)
s2 = var(ig1)
α, θ =  inverse_gamma_1_specification(m, s2)
@test α ≈ 1.3
@test θ ≈ 2.5

ig2 = InverseGamma(2.3, 2.5)
m = mean(ig2)
s2 = var(ig2)
α, θ =  inverse_gamma_2_specification(m, s2)
@test α ≈ 2.3
@test θ ≈ 2.5

