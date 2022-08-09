include("../distributions/inversegamma1.jl")

mutable struct Prior
    boundaries::Any
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
    Prior() = new(repeat([missing], 18)...)
end

function parse_prior!(context, field)
    p = Prior()
    _parse!(p, field)
    @assert !ismissing(p.name) || (!ismissing(p.name1) && !ismissing(p.name2))
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    index, name = get_index_name(p, symboltable)
    push!(ep.index, index)
    push!(ep.name, name)
    push!(ep.initialvalue, p.init)
    if !ismissing(p.name)
        parametertype = symboltable[p.name].symboltype
    else
        parametertype = symboltable[p.name1].symboltype
    end
    push!(ep.parametertype, parametertype)
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

function get_index_name(p, symboltable)
    if !ismissing(p.name1)
        return (
            Pair(symboltable[p.name1].orderintype, symboltable[p.name2].orderintype),
            Pair(p.name1, p.name2),
        )
    else
        return (symboltable[p.name].orderintype, p.name)
    end
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("beta")}
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, β = beta_specification(p.mean, p.variance)
    return Beta(α, β)
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("gamma")}
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, θ = gamma_specification(p.mean, p.variance)
    return Gamma(α, θ)
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("inv_gamma")}
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, θ = inverse_gamma_1_specification(p.mean, p.variance)
    return InverseGamma1(α, θ)
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("inv_gamma2")}
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, θ = inverse_gamma_2_specification(p.mean, p.variance)
    return InverseGamma(α, θ)
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("normal")}
)
    ismissing(p.mean) && !ismissing(p.median) && (p.mean = p["median"])
    ismissing(p.mean) && !ismissing(p.mode) && (p.mean = p["mode"])
    ismissing(p.stdev) && !ismissing(p.variance) && (σ = sqrt(p.variance))

    return Normal(p.mean, p.stdev)
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("uniform")}
)
    if ismissing(p.domain)
        ismissing(p.mean) && !ismissing(p.median) && (p.mean = p["median"])
        ismissing(p.variance) && !ismissing(p.stdev) && (p.variance = p.stdev * p.stdev)
        a, b = uniform_specification(p.mean, p.variance)
    else
        a, b = p.domain
    end
    return Uniform(a, b)
end

function make_prior_distribution(
    p::Prior,
    ::Val{Symbol("weibull")}
)
    ismissing(p.variance) && !ismissing(p.stdev) && (p.variance = p.stdev * p.stdev)
    α, θ = weibull_specification(p.mean, p.variance)
    return Weibull(α, θ)
end
