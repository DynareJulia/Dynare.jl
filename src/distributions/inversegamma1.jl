"""
    InverseGamma1(α, θ)

The *inverse Gamma distribution* with shape parameter `α` and scale `θ` has probability
density function

```math
f(x; \\alpha, \\theta) = \\frac{\\theta^\\alpha x^{-(\\alpha + 1)}}{\\Gamma(\\alpha)}
e^{-\\frac{\\theta}{x}}, \\quad x > 0
```

It is related to the [`Gamma`](@ref) distribution: if ``X \\sim \\operatorname{Gamma}(\\alpha, \\beta)``, then ``1 / \\sqrt{X} \\sim \\operatorname{InverseGamma1}(\\alpha, \\beta^{-1})``.

```julia
InverseGamma1()        # Inverse Gamma distribution with unit shape and unit scale, i.e. InverseGamma1(1, 1)
InverseGamma1(α)       # Inverse Gamma distribution with shape α and unit scale, i.e. InverseGamma1(α, 1)
InverseGamma1(α, θ)    # Inverse Gamma distribution with shape α and scale θ

params(d)        # Get the parameters, i.e. (α, θ)
shape(d)         # Get the shape parameter, i.e. α
scale(d)         # Get the scale parameter, i.e. θ
```

External links

* [Inverse gamma distribution on Wikipedia](http://en.wikipedia.org/wiki/Inverse-gamma_distribution)
"""
struct InverseGamma1{T<:Real} <: ContinuousUnivariateDistribution
    invd::Gamma{T}
    θ::T
    InverseGamma1{T}(α::T, θ::T) where {T<:Real} = new{T}(Gamma(α, inv(θ), check_args=false), θ)
end

function InverseGamma1(α::T, θ::T; check_args::Bool=true) where {T <: Real}
    @check_args InverseGamma1 (α, α > zero(α)) (θ, θ > zero(θ))
    return InverseGamma1{T}(α, θ)
end

InverseGamma1(α::Real, θ::Real; check_args::Bool=true) = InverseGamma1(promote(α, θ)...; check_args=check_args)
InverseGamma1(α::Integer, θ::Integer; check_args::Bool=true) = InverseGamma1(float(α), float(θ); check_args=check_args)
InverseGamma1(α::Real; check_args::Bool=true) = InverseGamma1(α, one(α); check_args=check_args)
InverseGamma1() = InverseGamma1{Float64}(1.0, 1.0)

@distr_support InverseGamma1 0.0 Inf

#### Conversions
convert(::Type{InverseGamma1{T}}, α::S, θ::S) where {T <: Real, S <: Real} = InverseGamma1(T(α), T(θ))
function Base.convert(::Type{InverseGamma1{T}}, d::InverseGamma1) where {T<:Real}
    return InverseGamma1{T}(T(shape(d)), T(d.θ))
end
Base.convert(::Type{InverseGamma1{T}}, d::InverseGamma1{T}) where {T<:Real} = d

#### Parameters

shape(d::InverseGamma1) = shape(d.invd)
scale(d::InverseGamma1) = d.θ
rate(d::InverseGamma1) = scale(d.invd)

params(d::InverseGamma1) = (shape(d), scale(d))
partype(::InverseGamma1{T}) where {T} = T


#### Parameters

mean(d::InverseGamma1{T}) where {T} = ((α, θ) = params(d); α  > 1 ? θ / (α - 1) : T(Inf))

mode(d::InverseGamma1) = scale(d) / (shape(d) + 1)

function var(d::InverseGamma1{T}) where T<:Real
    (α, θ) = params(d)
    α > 2 ? θ^2 / ((α - 1)^2 * (α - 2)) : T(Inf)
end

function skewness(d::InverseGamma1{T}) where T<:Real
    α = shape(d)
    α > 3 ? 4sqrt(α - 2) / (α - 3) : T(NaN)
end

function kurtosis(d::InverseGamma1{T}) where T<:Real
    α = shape(d)
    α > 4 ? (30α - 66) / ((α - 3) * (α - 4)) : T(NaN)
end

function entropy(d::InverseGamma1)
    (α, θ) = params(d)
    α + loggamma(α) - (1 + α) * digamma(α) + log(θ)
end

function kldivergence(p::InverseGamma1, q::InverseGamma1)
    # We can reuse the implementation of Gamma
    return kldivergence(p.invd, q.invd)
end


#### Evaluation

function logpdf(d::InverseGamma1, x::Real)
    (α, θ) = params(d)
    log(2) + α * log(θ) - loggamma(α) - (0.5*α + 1) * log(x) - θ / (x .* x)
end

function zval(::InverseGamma1, x::Real) = inv(sqrt(max(x, 0)))

cdf(d::InverseGamma1, x::Real) = ccdf(d.invd, zval(d, x))
ccdf(d::InverseGamma1, x::Real) = cdf(d.invd, zval(d, x))
logcdf(d::InverseGamma1, x::Real) = logccdf(d.invd, zval(d, x))
logccdf(d::InverseGamma1, x::Real) = logcdf(d.invd, zval(d, x))

quantile(d::InverseGamma1, p::Real) = inv(cquantile(d.invd, p))
cquantile(d::InverseGamma1, p::Real) = inv(quantile(d.invd, p))
invlogcdf(d::InverseGamma1, p::Real) = inv(invlogccdf(d.invd, p))
invlogccdf(d::InverseGamma1, p::Real) = inv(invlogcdf(d.invd, p))

function mgf(d::InverseGamma1{T}, t::Real) where T<:Real
    (a, b) = params(d)
    t == zero(t) ? one(T) : 2(-b*t)^(0.5a) / gamma(a) * besselk(a, sqrt(-4*b*t))
end

function cf(d::InverseGamma1{T}, t::Real) where T<:Real
    (a, b) = params(d)
    t == zero(t) ? one(T)+zero(T)*im : 2(-im*b*t)^(0.5a) / gamma(a) * besselk(a, sqrt(-4*im*b*t))
end


#### Evaluation

rand(rng::AbstractRNG, d::InverseGamma1) = 1 / sqrt(rand(rng, d.invd))
