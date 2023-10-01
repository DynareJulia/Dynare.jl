import Base: rand, findfirst
using Distributions
import Distributions: pdf, logpdf, cdf

struct Parameters
    symbols::Vector{Symbol}
    calibrations::Vector{Float64}
    distributions::Vector{Distribution}
    index::Vector{Int}
    initialvalues::Vector{Float64}
    mh_scale::Vector{Float64}
    optim_lb::Vector{Float64}
    optim_ub::Vector{Float64}
    posteriormodes::Vector{Float64}
    posteriormeans::Vector{Float64}
    posteriormedian::Vector{Float64}
    posteriorHPI::Matrix{Float64}
    transformfrombasic::Vector{Function}
    transformtobasic::Vector{Function}
    function Parameters()
        symbols = Vector{Symbol}(undef, 0)
        calibrations = Vector{Float64}(undef, 0)
        distributions = Vector{Distribution}(undef, 0)
        index = Vector{Int}(undef, 0)
        initialvalues = Vector{Float64}(undef, 0)
        mh_scale = Vector{Float64}(undef, 0)
        optim_lb = Vector{Float64}(undef, 0)
        optim_ub = Vector{Float64}(undef, 0)
        posteriormodes = Vector{Float64}(undef, 0)
        posteriormeans = Vector{Float64}(undef, 0)
        posteriormedian = Vector{Float64}(undef, 0)
        posteriorHPI = Matrix{Float64}(undef, 0, 2)
        transformfrombasic = Vector{Function}(undef, 0)
        transformtobasic = Vector{Function}(undef, 0)
        new(
            symbols,
            calibrations,
            distributions,
            index,
            initialvalues,
            mh_scale,
            optim_lb,
            optim_ub,
            posteriormodes,
            posteriormeans,
            posteriormedian,
            posteriorHPI,
            transformfrombasic,
            transformtobasic,
        )
    end
end

function parse_estimated_parameters!(context::Context, fields::Dict{String,Any})
    parameters = context.work.estimated_parameters
    symbol_table = context.symboltable
    for p in fields["params"]
        if "param" in keys(p)
            push!(parameters.name, p["param"])
            push!(parameters.index, symbol_table[p["param"]].orderintype)
        elseif "var" in keys(p)
            push!(parameters.name, p["var"])
            push!(parameters.index, symbol_table[p["var"]].orderintype)
        else
            error("Estimation needs a datafile or an AxisArrayTable")
        end
        initval = dynare_parse_eval(p["init_val"], context)
        if isnan(initval)
            push!(parameters.initialvalue, missing)
        else
            push!(parameters.initialvalue, initval)
        end
            #=
        push!(parameters.optim_lb, dynare_parse_eval(p["lower_bound"], context))
        push!(parameters.optim_ub, dynare_parse_eval(p["upper_bound"], context))
        push!(
            parameters.transformfrombasic,
            get_transformfrombasic(Val(p["prior_distribution"]), p["3"], p["4"]),
        )
        push!(
            parameters.transformtobasic,
            get_transformtobasic(Val(p["prior_distribution"]), p["3"], p["4"]),
        )
        push!(parameters.mh_scale, dynare_parse_eval(p["scale"], context))
        =#
        push!(
            parameters.prior,
            parse_prior_distribution(Val(p["prior_distribution"]), p),
        )
    end
end

function Base.findfirst(pattern::String, collection::Vector{Union{String, Pair{String, String}}})
    for (i, element) in enumerate(collection)
        if pattern == element
            return i
        end 
    end 
    return nothing
end 

function parse_estimated_parameter_init!(context::Context, fields::Dict{String,Any})
    parameters = context.work.estimated_parameters
    np = length(parameters)
    for p in fields["params"]
        if "param" in keys(p)
            k = Base.findfirst(p["param"], parameters.name)
        elseif "var" in keys(p)
            k = Base.findfirst(p["var"], parameters.name)
        end
        if isnothing(k)
            error("estimated_params_init: parameter $(p["param"]) hasn't been declared as an estimated parameter")
        else
            parameters.initialvalue[k] = dynare_parse_eval(p["init_val"], context)
        end
    end
end 

function get_basic_parameters(p::Dict{String,Any})
    if "param" in keys(p)
        name = p["param"]
    elseif "var" in keys(p)
        name = p["var"]
    else
        error("Estimation needs a datafile or an AxisArrayTable")
    end
    return (
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
        p3 = dynare_parse_eval(p["p3"], context),
        p4 = dynare_parse_eval(p["p4"], context),
        lb = dynare_parse_eval(p["lower_bound"], context),
        ub = dynare_parse_eval(p["upper_bound"], context),
        name = name,
    )
end

function parse_prior_distribution(::Val{1}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, β = beta_specification(μ, σ * σ, name = name)
    d = Beta(α, β)
end

function parse_prior_distribution(::Val{2}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, θ = gamma_specification(μ, σ*σ)
    d = Gamma(α, θ)
end

function parse_prior_distribution(::Val{3}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    d = Normal(μ, σ)
end

# Inverse gamma 1
function parse_prior_distribution(::Val{4}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, θ = inverse_gamma_1_specification(μ, σ*σ)
    d = InverseGamma1(α, θ)
end

function parse_prior_distribution(::Val{5}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    @show μ, σ, p3, p4
    if isnan(μ) && isnan(σ)
        a = p3
        b = p4
    else
        a, b = uniform_specification(μ, σ*σ)
    end
    d = Uniform(a, b)
 end

# Inverse gamma 2
function parse_prior_distribution(::Val{6}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, θ = inverse_gamma_2_specification(μ, σ*σ)
    d = InverseGamma(α, θ)
end

# Weibull
function parse_prior_distribution(::Val{8}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, θ = weibull_specification(μ, σ*σ)
    d = Weibull(α, θ)
end

function get_transformfrom(::Val{1}, p3, p4)
    if p3 == 0.0 && p4 == 1.0
        return (identity, (x) -> 1)
    else
        return ((x) -> (x - p3) / (p4 - p3), (x) -> p4 - p3)
    end
end

function get_transformto(::Val{1}, p3, p4)
    if p3 == 0.0 && p4 == 1.0
        return identity
    else
        return (x) -> x * (p4 - p3) + p3
    end
end

function get_transformfrom(::Val{2}, p3, p4)
    if p3 == 0.0 && p4 == 1.0
        return identity
    else
        return (x) -> x - p3
    end
end

function get_transformto(::Val{2}, p3, p4)
    if p3 == 0.0 && p4 == 1.0
        return identity
    else
        return (x) -> x + p3
    end
end
