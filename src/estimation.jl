import Base: rand
using Distributions
import Distributions: pdf, logpdf, cdf

struct Parameters
    symbols::Vector{Symbol}
    calibrations::Vector{Float64}
    distributions::Vector{Distribution}
    initialvalues::Vector{Float64}
    mh_scale::Vector{Float64}
    optim_lb::Vector{Float64}
    optim_ub::Vector{Float64}
    posteriormodes::Vector{Float64}
    posteriormeans::Vector{Float64}
    posteriormedian::Vector{Float64}
    posteriorHPI::Matrix{Float64}
    transformfrombasic::Function
    transformtobasic::Function
    function Parameters()
        symbols = Vector{Symbol}(undef, 0)
        calibrations = Vector{Float64}(undef, 0)
        distributions = Vector{Distribution}(undef, 0)
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
            initialvalues,
            mh_scale,
            optim_lb,
            optim_ub,
            posteriormodes,
            posteriormeans,
            posteriormedian,
            posteriorHPI,
            transfrombasic,
            transformtobasic,
        )
    end
end

dynare_parse_eval(s, c) = Dynare.dynare_parse_eval(s, c)

function parse_estimated_parameters!(parameters, context, fields::Dict{String,Any})
    #    parameters = context.results.model_results[1].parameters
    for p in fields["params"]
        push!(parameters.symbols, Symbol(p["param"]))
        push!(parameters.initialvalues, dynare_parse_eval(p["init_val"], context))
        @show dynare_parse_eval(p["lower_bound"], context)
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
        #        push!(parameters.mh_scale, dynare_parse_eval(p["scale"], context))
        push!(
            parameters.distributions,
            parse_prior_distribution(Val(p["prior_distribution"]), p),
        )
    end
end

function get_basic_parameters(p::Dict{String,Any})
    return (
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
        p3 = dynare_parse_eval(p["p3"], context),
        p4 = dynare_parse_eval(p["p4"], context),
        lb = dynare_parse_eval(p["lower_bound"], context),
        ub = dynare_parse_eval(p["upper_bound"], context),
        name = p["param"],
    )
end

function parse_prior_distribution(::Val{1}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, β = beta_specification(μ, σ * σ, name = name)
    d = Beta(α, β)
end

function parse_prior_distribution(::Val{2}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    d = Gamma(
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
    )
end

function parse_prior_distribution(::Val{3}, p)
    d = Normal(
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
    )
end

# Inverse gamma 1
function parse_prior_distribution(::Val{4}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, θ = inverse_gamma_1_specification(μ, σ)
    d = InverseGamma1(
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
    )
end

function parse_prior_distribution(::Val{5}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    a, b = uniform_specification(μ, σ, p3, p4)
    d = Uniform(
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
    )
end

# Inverse gamma 2
function parse_prior_distribution(::Val{6}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    α, θ = inverse_gamma_2_specification(μ, σ)
    d = InverseGamma(
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
    )
end

# Weibull
function parse_prior_distribution(::Val{8}, p)
    d = Weibull(
        μ = dynare_parse_eval(p["mean"], context),
        σ = dynare_parse_eval(p["std"], context),
    )
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



#pdf(p::Parameters, i::Int64, x::Float64) = pdf(p.distributions[i], p.

parameters = Parameters()
fields = Dict(
    "statementName" => "estimated_params",
    "params" => [
        Dict(
            "param" => "alp",
            "init_val" => "NaN",
            "lower_bound" => "(-Inf)",
            "upper_bound" => "Inf",
            "prior_distribution" => 1,
            "mean" => "0.356",
            "std" => "0.02",
            "p3" => "NaN",
            "p4" => "NaN",
            "jscale" => "NaN",
        ),
        Dict(
            "param" => "bet",
            "init_val" => "NaN",
            "lower_bound" => "(-Inf)",
            "upper_bound" => "Inf",
            "prior_distribution" => 1,
            "mean" => "0.993",
            "std" => "0.002",
            "p3" => "NaN",
            "p4" => "NaN",
            "jscale" => "NaN",
        ),
    ],
)
context = @dynare "test/models/estimation/fs2000.mod"

parameters = Parameters()
parse_estimated_parameters!(parameters, context, fields)
