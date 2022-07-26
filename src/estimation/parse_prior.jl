include("../distributions/inversegamma1.jl")

mutable struct Prior
    boundaries::String
    domain::Any
    init::Any
    interval::Any
    jscale::Any
    mean::Any
    median::Any
    mode::Any
    name::Any
    name1::Any
    name2::Any
    regimes::Any
    shape::Any
    shift::Any
    stdev::Any
    subsample::Any
    truncate::Any
    variance::Any
    Prior() = new()
end

function parse_prior!(context, field)
    p = Prior()
    _parse!(p, field)
    ep = context.work.estimated_parameters
    push!(ep.name, p.name)
    push!(ep.prior, make_prior_distribution(p, Val(Symbol(p.shape))))
end

function _parse!(p, field)
    for (k, v) in field
        if k == "statementName"
            nothing
        elseif k == "options"
            for (k1, v1) in v
                setfield!(p, Symbol(k1), v1)
            end
        else
            setfield!(p, Symbol(k), v)
        end
    end
    return p
end

function make_prior_distribution(p::Prior, ::Val{Symbol("beta")})
    isdefined(p, :stdev) && !isdefined(p, :variance) && (p.variance = p.stdev * p.stdev)
    α, β = beta_specification(p.mean, p.variance)
    return Beta(α, β)
end

function make_prior_distribution(p::Prior, ::Val{Symbol("gamma")})
    isdefined(p, :stdev) && !isdefined(p, :variance) && (p.variance = p.stdev * p.stdev)
    α, θ = gamma_specification(p.mean, p.variance)
    return Gamma(α, θ)
end

function make_prior_distribution(p::Prior, ::Val{Symbol("inv_gamma")})
    isdefined(p, :stdev) && !isdefined(p, :variance) && (p.variance = p.stdev * p.stdev)
    α, θ = inverse_gamma_1_specification(p.mean, p.variance)
    return InverseGamma1(α, θ)
end

function make_prior_distribution(p::Prior, ::Val{Symbol("inv_gamma2")})
    isdefined(p, :stdev) && !isdefined(p, :variance) && (p.variance = p.stdev * p.stdev)
    α, θ = inverse_gamma_2_specification(p.mean, p.variance)
    return InverseGamma(α, θ)
end

function make_prior_distribution(p::Prior, ::Val{Symbol("normal")})
    !isdefined(p, :mean) && isdefined(p, :median) && (p.mean = p["median"])
    !isdefined(p, :mean) && isdefined(p, :mode) && (p.mean = p["mode"])
    !isdefined(p, :stdev) && isdefined(p, :variance) && (σ = sqrt(p.variance))
    return Normal(p.mean, p.stdev)
end

function make_prior_distribution(p::Prior, ::Val{Symbol("uniform")})
    if !isdefined(p, :domain)
        !isdefined(p, :mean) && isdefined(p, :median) && (p.mean = p["median"])
        !isdefined(p, :variance) && isdefined(p, :stdev) && (p.variance = p.stdev * p.stdev)
        a, b = uniform_specification(p.mean, p.variance)
    else
        a, b = p.domain
    end
    return Uniform(a, b)
end

function make_prior_distribution(p::Prior, ::Val{Symbol("weibull")})
    !isdefined(p, :variance) && isdefined(p, :stdev) && (p.variance = p.stdev * p.stdev)
    α, θ = weibull_specification(p.mean, p.variance)
    return Weibull(α, θ)
end

#=
testcase = [
Dict("statementName" => "prior", "name" => "alpha", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.356, "stdev" => 0.02)), 
Dict("statementName" => "prior", "name" => "beta", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.993, "stdev" => 0.002)), 
Dict("statementName" => "prior", "name" => "rho", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.129, "stdev" => 0.223)), 
Dict("statementName" => "prior", "name" => "delta", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.01, "stdev" => 0.005)), 
Dict("statementName" => "prior", "name" => "theta", "subsample" => "", "shape" => "normal", "options" => Dict("mean" => 3.0, "stdev" => 1.0)), 
Dict("statementName" => "prior", "name" => "tau", "subsample" => "", "shape" => "beta", "options" => Dict("mean" => 0.03, "stdev" => 0.01))
]

for c in testcase
    @show parse_prior([], c)
end  
=#
