import Base: findfirst
using Distributions
import Distributions: pdf, logpdf, cdf

include("../distributions/inversegamma1.jl")

# prior!() function

# alias for InverseGamma distribution
InverseGamma2 = InverseGamma

struct corr
    s1::Symbol
    s2::Symbol
    function corr(s1_arg, s2_arg)
        s1 = s1_arg
        s2 = s2_arg
        new(s1, s2)
    end 
end 

struct stdev
    s::Symbol
    function stdev(s_arg::Symbol)
        s = s_arg
        new(s)
    end
end

struct variance
    s::Symbol
    function stdev(s_arg::Symbol)
        s = s_arg
        new(s)
    end
end

"""
     prior!(s::Symbol; shape::{<:Distributions}, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
     prior!(s::stdev; shape::{<:Distributions}, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
     prior!(s::variance; shape::{<:Distributions}, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
     prior!(s::corr; shape::{<:Distributions}, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)

generates a prior for a symbol of a parameter, the standard deviation (`stdev`) or the variance (`variance`) of an exogenous variable or an endogenous variable (measurement error) or the correlation (`corr`) between 2 endogenous or exogenous variables

# Keyword arguments
- `shape <: Distributions`: the shape of the prior distribution (`Beta`, `InvertedGamma`, `InvertedGamma1`, `Gamma`, `Normal`, `Uniform`, `Weibull`) [required]
- `context::Context=context`: context in which the prior is declared
- `domain::Vector{<:Real}=Float64[]`: domain for a uniform distribution
- `initialvalue::Union{Real,Missing}=missing`: initialvalue for mode finding or MCMC iterations
- `mean::Union{Real,Missing}=missing`: mean of the prior distribution
- `stdev::Union{Real,Missing}=missing`: stdev of the prior distribution
- `variance::Union{Real,Missing}=missing`: variance of the prior distribution
"""
    prior!(s; 
           shape, 
           initialvalue::Union{Real,Missing}=missing, 
           mean::Union{Real,Missing}=missing, 
           stdev::Union{Real,Missing}=missing, 
           domain=[], variance ::Union{Real,Missing}=missing, 
           context::Context=context
           )

function prior!(s::Symbol; 
    shape, initialvalue::Union{Real,Missing}=missing, 
    mean::Union{Real,Missing}=missing, 
    stdev::Union{Real,Missing}=missing, 
    domain=[], variance ::Union{Real,Missing}=missing, 
    context::Context=context
    )
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    name = string(s)
    push!(ep.index, symboltable[name].orderintype)
    push!(ep.name, name)
    push!(ep.initialvalue, initialvalue)
    push!(ep.parametertype, EstParameter)
    prior_(s, shape, mean, stdev, domain, variance, ep)
end
    
function prior!(s::stdev; shape, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    name = string(s.s)
    push!(ep.name, name)
    push!(ep.initialvalue, initialvalue)
    if symboltable[name].symboltype == Exogenous
        index = symboltable[name].orderintype
        push!(ep.index, (index => index))
        push!(ep.parametertype, EstSDShock)
    elseif symboltable[name].symboltype == Endogenous
        index = findfirst(name .== varobs)
        push!(ep.index, (index => index))
        push!(ep.parametertype, EstSDMeasurement)
    else
        error("Wrong symboltype")
    end 
    prior_(s, shape, mean, stdev, domain, variance, ep)
end

function prior!(s::variance; shape, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    name = string(s.s)
    push!(ep.name, name)
    push!(ep.initialvalue, initialvalue)
    if symboltable[name].symboltype == Exogenous
        index = symboltable[name].orderintype
        push!(ep.index, (index => index))
        push!(ep.parametertype, EstVarShock)
    elseif symboltable[name].symboltype == Endogenous
        index = findfirst(name .== varobs)
        push!(ep.index, (index => index))
        push!(ep.parametertype, EstVarMeasurement)
    else
        error("Wrong symboltype")
    end 
    prior_(s, shape, mean, stdev, domain, variance, ep)
end

function prior!(s::corr; shape, initialvalue::Union{Real,Missing}=missing, mean::Union{Real,Missing}=missing, stdev::Union{Real,Missing}=missing, domain=[], variance ::Union{Real,Missing}=missing, context::Context=context)
    ep = context.work.estimated_parameters
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    name1 = string(s.s1)
    name2 = string(s.s2)
    push!(ep.name, (name1 => name2))
    push!(ep.initialvalue, initialvalue)
    if symboltable[name1].symboltype == Exogenous
        index1 = symboltable[name1].orderintype
        index2 = symboltable[name2].orderintype
        push!(ep.index, (index1 =>index2))
        push!(ep.parametertype, EstCorrShock)
    elseif symboltable[name1].symboltype == Endogenous
        index1 = findfirst(name1 .== varobs)
        index2 = findfirst(name2 .== varobs)
        push!(ep.index, (index1 =>index2))
        push!(ep.parametertype, EstCorrMeasurement)
    else
        error("Wrong symboltype")
    end 
    prior_(s, shape, mean, stdev, domain, variance, ep)
end

function prior_(s, shape, mean, stdev, domain, variance, ep)
    if ismissing(variance) && !ismissing(stdev)
        variance = stdev*stdev
    end 
    if shape == Beta
        α, β = beta_specification(mean, variance)
        push!(ep.prior, Beta(α, β))
    elseif shape == Gamma
        α, β = gamma_specification(mean, variance) 
        push!(ep.prior, Gamma(α, β))
    elseif shape == InverseGamma1
        α, β = inverse_gamma_1_specification(mean, variance) 
        push!(ep.prior, InverseGamma1(α, β))
    elseif shape == InverseGamma2
        α, β = inverse_gamma_2_specification(mean, variance)
        push!(ep.prior, InverseGamma(α, β))
    elseif shape == Normal
        push!(ep.prior, Normal(mean, stdev))
    elseif shape == Uniform
        if isempty(domain)
            a, b = uniform_specification(mean, variance)
        else
            a, b = domain
        end 
        push!(ep.prior, Uniform(a, b))
    elseif shape == Weibull
        α, β = weibull_specification(mean, stdev)
        push!(ep.prior, Weibull(α, β))
    else
        error("Unknown prior shape")
    end 
end         

#=
function get_index_name(s::Symbol, symboltable::SymbolTable)
    name = String(s)
    index = symboltable[name].orderintype
    return (index, name)
end 
=#

"""
    plot_priors(; context::Context = context, n_points::Int = 100)

plots prior density

# Keyword arguments
 - `context::Context = context`: context in which to take the date to be ploted
 - `n_points::Int = 100`: number of points used for a curve
"""
function plot_priors(; context::Context = context, n_points = 100)
    ep = context.work.estimated_parameters
    @assert length(ep.prior) > 0 "There is no defined prior"

    path = "$(context.modfileinfo.modfilepath)/graphs/"
    mkpath(path)
    filename = "$(path)/Priors"
    nprior = length(ep.prior)
    (nbplt, nr, nc, lr, lc, nstar) = pltorg(nprior)
    ivars = collect(1:nr*nc)
    names = get_parameter_names(ep)
    for p = 1:nbplt
        filename1 = "$(filename)_$(p).png"
        p == nbplt && (ivars = ivars[1:nstar])
        plot_panel_priors(
            ep.prior,
            names,
            ivars,
            nr,
            nc,
            filename1
        )
        ivars .+= nr * nc
    end
end

function plot_panel_priors(
    prior,
    ylabels,
    ivars,
    nr,
    nc,
    filename;
    kwargs...
    )
    sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]
    for (i, j) in enumerate(ivars)
        sp[i] = Plots.plot(prior[j], title = ylabels[j], labels = "Prior", kwargs...)
    end
    
    pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title = "Prior distributions")
    graph_display(pl)
    savefig(filename)
end

"""
    plot_prior_posterior(chains; context::Context=context)

plots priors posterios and mode if computed on the same plots

# Keyword arguments
- `context::Context=context`: context used to get the estimation results

# Output
- the plots are saved in `./<modfilename>/Graphs/PriorPosterior_<x>.png`
"""
function plot_prior_posterior(chains; context::Context=context)
    ep = context.work.estimated_parameters
    mode = context.results.model_results[1].estimation.posterior_mode
    names = get_parameter_names(ep)
    @assert length(ep.prior) > 0 "There is no defined prior"
    @assert length(chains) > 0 "There is no MCMC chain"

    path = "$(context.modfileinfo.modfilepath)/graphs/"
    mkpath(path)
    filename = "$(path)/PriorPosterior"
    nprior = length(ep.prior)
    (nbplt, nr, nc, lr, lc, nstar) = pltorg(nprior)
    ivars = collect(1:nr*nc)
    for p = 1:nbplt-1
        filename1 = "$(filename)_$(p).png"
        plot_panel_prior_posterior(
            ep.prior,
            chains,
            mode,
            names,
            ivars,
            nr,
            nc,
            filename1
        )
        ivars .+= nr * nc
    end
    ivars = ivars[1:nstar]
    filename = "$(filename)_$(nbplt).png"
    plot_panel_prior_posterior(
        ep.prior,
        chains,
        mode,
        names,
        ivars,
        lr,
        lc,
        filename
    )
end

function plot_panel_prior_posterior(
    prior,
    chains,
    mode,
    ylabels,
    ivars,
    nr,
    nc,
    filename;
    kwargs...
)
sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]
for (i, j) in enumerate(ivars)
    posterior_density = kde(vec(get(chains, Symbol(ylabels[j]))[1].data))
    
    sp[i] = Plots.plot(prior[j], title = ylabels[j], labels = "Prior", kwargs...)
    plot!(posterior_density, labels = "Posterior")
end

pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title = "Prior and posterior distributions")
graph_display(pl)
savefig(filename)
end



# Dynare estimated_parameters statement
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
    symboltable = context.symboltable
    for p in fields["params"]
        if "param" in keys(p)
            push!(parameters.name, p["param"])
            push!(parameters.index, symboltable[p["param"]].orderintype)
            push!(parameters.parametertype, EstParameter)
        elseif "var" in keys(p)
            push!(parameters.name, p["var"])
            if is_endogenous(p["var"], symboltable)
                k = findfirst(p["var"] .== context.work.observed_variables)
                push!(parameters.index, Pair(k, k))
                push!(parameters.parametertype, EstSDMeasurement)
            elseif is_exogenous(p["var"], symboltable)
                k = symboltable[p["var"]].orderintype
                push!(parameters.index, Pair(k, k))
                push!(parameters.parametertype, EstSDShock)
            else
                error("unknown variable type: $(p["var"])")
            end
        elseif "var1" in keys(p)
            push!(parameters.name, (p["var1"] => p["var2"]))
            if is_endogenous(p["var1"], symboltable)
                k1 = findfirst(p["var1"] .== context.work.observed_variables)
                k2 = findfirst(p["var2"] .== context.work.observed_variables)
                push!(parameters.index, Pair(k1, k2))
                push!(parameters.parametertype, EstCorrMeasurement)
            elseif is_exogenous(p["var1"], symboltable)
                k1 = symboltable[p["var1"]].orderintype
                k2 = symboltable[p["var2"]].orderintype
                push!(parameters.index, Pair(k1, k2))
                push!(parameters.parametertype, EstCorrShock)
            else
                error("unknown variable type: $(p["var1"])")
            end
        else
            error("Unrecognized parameter name")
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
    elseif "var1" in keys(p)
        name = (p["var1"] => p["var2"])
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
