"""
NonLinearSolveAlgos - enumerate

Nonlinear system of equations available algorithms:
- mcp: mixed complementarity problem algorithm fro NLsolve package
- trustregion: trust-region algorithm from NLsolve package
"""
@enum PerfectForesightAlgos mcp, trustregion 

struct PerfectForesightOptions
    algo::PerfectForesightAlgos
    datafile::String
    display::Bool
    maxit::Int64
    periods::Int64
    tolf::Float64
    tolx::Float64
    function PerfectForesightOptions(options::Dict{String,Any})
        algo = trustregion
        datafile = ""
        maxit = 50
        tolf = 1e-5
        tolx = 1e-5
        print = true
        for (k, v) in pairs(options)
            if k == "algo"
                algo = v::PerfectForesightAlgos
            elseif k == "datafile"
                datafile = v::String
            elseif k == "maxit"
                maxit = v::Int64
            elseif k == "periods"
                periods = v::Int64
            elseif k == "tolf"
                tolf = v::Float64
            elseif k == "tolx"
                tolx = v::Float64
            else
                error("Unknown option $k")
            end
        end
        new(algo, datafile, display, maxit, periods,
            tolf, tolx)
    end
end

struct PerfectForesightWs
    y::Vector{Float64}
    x::Vector{Float64}
    function PerfectForesightWs(endogenous_nbr,
                                exogenous_nbr,
                                periods)
        y = Vector{Float64}(undef, endogenous_nbr*periods)
        x = Vector{Float64}(undef, exogenous_nbr*periods)
        new(y, x)
    end
end
    
struct PerfectForesightMcpWs
    y::Vector{Float64}
    x::Vector{Float64}
    lb::Vector{Float64}
    ub::Vector{Float64}
    function PerfectForesightWs(endogenous_nbr,
                                exogenous_nbr,
                                periods)
        y = Vector{Float64}(undef, endogenous_nbr*periods)
        x = Vector{Float64}(undef, exogenous_nbr*periods)
        lb = Vector{Float64}(undef, endogenous_nbr*periods)
        ub = Vector{Float64}(undef, endogenous_nbr*periods)
        new(y, x, lb, ub)
    end
end
    
function perfect_foresight_setup!(context::Context, field::Dict{String,Any})
    if !("options" in keys(field)) || !("periods" in keys(field["options"]))
        error("perfect_foresight_setup must contain option 'periods'")
    end
    context.options["perfect_foresight_setup"] = Dict()
    copy!(context.options["perfect_foresight_setup"], field["options"])
    
end

function perfect_foresight_solver!(context::Context, field::Dict{String,Any})
    context.options["perfect_foresight_solver"] = Dict()
    copy!(context.options["perfect_foresight_solver"], field["options"])
    compute_perfect_foresight_solver(context)
end

function compute_prefect_foresight_setup(context::Context)
end

function compute_perfect_foresight_solver(context::Context)
end
