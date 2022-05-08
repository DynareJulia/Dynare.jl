using Distributions
using Dynare

struct Parameters
    symbols::Vector{Symbol}
    calibrations::Vector{Float64}
    distributions::Vector{Distribution}
    initialvalues::Vector{Float64}
    mh_scale::Vector{Float64}
    optim_lb::Vector{Float64}
    optin_ub::Vector{Float64}
    posteriormodes::Vector{Float64}
    posteriormeans::Vector{Float64}
    posteriormedian::Vector{Float64}
    posteriorHPI::Matrix{Float64}
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
        new(symbols, calibrations, distributions, initialvalues, mh_scale,
            optim_lb, optim_ub, 
            posteriormodes, posteriormeans, posteriormedian, posteriorHPI)
    end
end

function parse_estimated_parameters!(context::Context, fields::Dict{String, Any})
    parameters = context.results.model_results[1].parameters
    for p in fields["params"]
        push!(parameters.symbols, Symbol(p["param"]))
        push!(parameters.initialvalues, dynare_parse_eval(p["init_val"]))
        push!(parameters.optim_lb, dynare_parse_eval(p["lower_bound"])),
        push!(parameters.optim_ub, dynare_parse_eval(p["upper_bound"])),
        push!(parameters.mh_scale, dynare_parse_eval(p["scale"])),
        push!(parameters.distributions, parse_prior_distribution(Val(p["prior_distribution"]), p))
    end
end

function get_basic_parameters(p::Dict{String, Any})
    return (μ = dynare_parse_eval(p["mean"]),
            σ = dynare_parse_eval(p["std"]),
            p3 = dynare_parse_eval(p["p3"]),
            p4 = dynare_parse_eval(p["p4"]),
            name = p["param"])
end

function parse_prior_distribution(::Val{1}, p)
    μ, σ, p3, p4, lb, ub, name = get_basic_parameters(p)
    
    α, β = beta_specification(μ, σ*σ, lb, ub, name)
    d = Beta()
end

function parse_prior_distribution(::Val{2}, p)
    d = Gamma(μ = dynare_parse_eval(p["mean"]), σ = dynare_parse_eval(p["std"]))
end

function parse_prior_distribution(::Val{3}, p)
    d = Normal(μ = dynare_parse_eval(p["mean"]), σ = dynare_parse_eval(p["std"]))
end

# Inverse gamma 1
function parse_prior_distribution(::Val{4}, p)
    d = Normal(μ = dynare_parse_eval(p["mean"]), σ = dynare_parse_eval(p["std"]))
end

function parse_prior_distribution(::Val{5}, p)
    d = Uniform(μ = dynare_parse_eval(p["mean"]), σ = dynare_parse_eval(p["std"]))
end

# Inverse gamma 2
function parse_prior_distribution(::Val{6}, p)
    d = Uniform(μ = dynare_parse_eval(p["mean"]), σ = dynare_parse_eval(p["std"]))
end

function parse_prior_distribution(::Val{8}, p)
    d = Weibull(μ = dynare_parse_eval(p["mean"]), σ = dynare_parse_eval(p["std"]))
end

parameters = Parameters()
fields = Dict("statementName" => "estimated_params", "params" => [Dict("param" => "alp", "init_val" => "NaN", "lower_bound" => "(-Inf)", "upper_bound" => "Inf", "prior_distribution" => 1, "mean" => "0.356", "std" => "0.02", "p3" => "NaN", "p4" => "NaN", "jscale" => "NaN"),
                                                         Dict("param" => "bet", "init_val" => "NaN", "lower_bound" => "(-Inf)", "upper_bound" => "Inf", "prior_distribution" => 1, "mean" => "0.993", "std" => "0.002", "p3" => "NaN", "p4" => "NaN", "jscale" => "NaN")])
context = @dynare "dynare.jl/test/models/estimation/fs2000.mod"

parse_estimated_parameters!(context, fields)
