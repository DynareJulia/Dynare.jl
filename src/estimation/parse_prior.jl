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
    if !ismissing(p.name)
        parametertype = symboltable[p.name].symboltype
    else
        parametertype = symboltable[p.name1].symboltype
    end
    index, name = get_index_name(p, symboltable)
    make_prior_distribution!(ep, p, Val(Symbol(p.shape)), index, name, parametertype)
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

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("beta")},
    index,
    name,
    parametertype,
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, β = beta_specification(p.mean, p.variance)
    push!(
        ep.prior_01,
        (
            index = index,
            initialvalue = p.init,
            name = name,
            prior = Beta(α, β),
            parametertype = parametertype,
        ),
    )
    return nothing
end

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("gamma")},
    index,
    name,
    parametertype,
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, θ = gamma_specification(p.mean, p.variance)
    push!(
        ep.prior_Rplus,
        (
            index = index,
            initialvalue = p.init,
            name = name,
            prior = Gamma(α, θ),
            parametertype = parametertype,
        ),
    )
    return nothing
end

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("inv_gamma")},
    index,
    name,
    parametertype,
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, θ = inverse_gamma_1_specification(p.mean, p.variance)
    push!(
        ep.prior_Rplus,
        (
            index = index,
            initialvalue = p.init,
            name = name,
            prior = InverseGamma1(α, θ),
            parametertype = parametertype,
        ),
    )
    return nothing
end

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("inv_gamma2")},
    index,
    name,
    parametertype,
)
    !ismissing(p.stdev) && ismissing(p.variance) && (p.variance = p.stdev * p.stdev)
    α, θ = inverse_gamma_2_specification(p.mean, p.variance)
    push!(
        ep.prior_Rplus,
        (
            index = index,
            initialvalue = p.init,
            name = name,
            prior = InverseGamma(α, θ),
            parametertype = parametertype,
        ),
    )
    return nothing
end

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("normal")},
    index,
    name,
    parametertype,
)
    ismissing(p.mean) && !ismissing(p.median) && (p.mean = p["median"])
    ismissing(p.mean) && !ismissing(p.mode) && (p.mean = p["mode"])
    ismissing(p.stdev) && !ismissing(p.variance) && (σ = sqrt(p.variance))

    push!(
        ep.prior_R,
        (
            index = index,
            initialvalue = p.init,
            name = name,
            prior = Normal(p.mean, p.stdev),
            parametertype = parametertype,
        ),
    )
    return nothing
end

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("uniform")},
    index,
    name,
    parametertype,
)
    if ismissing(p.domain)
        ismissing(p.mean) && !ismissing(p.median) && (p.mean = p["median"])
        ismissing(p.variance) && !ismissing(p.stdev) && (p.variance = p.stdev * p.stdev)
        a, b = uniform_specification(p.mean, p.variance)
    else
        a, b = p.domain
    end
    push!(
        ep.prior_AB,
        (
            index = index,
            initialvalue = p.init,
            name = name,
            prior = Uniform(a, b),
            parametertype = parametertype,
        ),
    )
    return nothing
end

function make_prior_distribution!(
    ep::EstimatedParameters,
    p::Prior,
    ::Val{Symbol("weibull")},
    index,
    name,
    parametertype,
)
    ismissing(p.variance) && !ismissing(p.stdev) && (p.variance = p.stdev * p.stdev)
    α, θ = weibull_specification(p.mean, p.variance)
    push!(
        ep.prior_Rplus,
        (
            index = index,
            initialvalue = init,
            name = p.name,
            prior = Weibull(α, θ),
            parametertype = parametertype,
        ),
    )
    return nothing
end
