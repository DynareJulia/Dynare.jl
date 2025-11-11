# traits used for PathSolver extension
abstract type NonLinearSolver end
struct PathNLS <: NonLinearSolver end
struct DefaultNLS <: NonLinearSolver end

# Dummy definition for PathSolver extension

# perfectforesight
function mcp_perfectforesight_core!(::DefaultNLS,
                                    perfect_foresight_ws::PerfectForesightWs,
                                    context::Dynare.Context,
                                    periods::Int64,
                                    guess_values::Vector{Float64},
                                    initialvalues::Vector{Float64},
                                    terminalvalues::Vector{Float64},
                                    dynamic_ws::Dynare.DynamicWs;
                                    maxit = maxit,
                                    tolf = tolf,
                                    tolx = tolx
                                    )
end

# perfectforesight with information about future endogenous variables
function mcp_perfectforesight_core!(::DefaultNLS,
                                    perfect_foresight_ws::PerfectForesightWs,
                                    context::Dynare.Context,
                                    periods::Int64,
                                    guess_values::Vector{Float64},
                                    initialvalues::Vector{Float64},
                                    terminalvalues::Vector{Float64},
                                    dynamic_ws::Dynare.DynamicWs,
                                    flipinfo::FlipInformation,
                                    infoperiod;
                                    maxit= maxit,
                                    tolf = tolf,
                                    tolx = tolx
                                    )
end

# sparsegrids
function mcp_sg_core!(::DefaultNLS)
end

function mcp_solve!(::DefaultNLS)
end
