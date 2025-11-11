include("blocks.jl")
include("block_firstorder.jl")
include("monomial.jl")
include("smolyak_gh.jl")
include("SG.jl")
include("DDSG.jl")
include("helper_functions.jl")

using AxisArrays: AxisArray
using DualNumbers
using Distributions
using FiniteDiff
using LinearAlgebra
import NLsolve
using NonlinearSolve
using Roots
using Tasmanian

export sparsegridapproximation, simulate!

@enum SGSolver NLsolver NonlinearSolver PATHSolver

struct SGModel
    dyn_endogenous::Vector{Float64} # dynamic endogenous variables in period t-1, t, t+1 
    exogenous::Vector{Float64}      # exogenous variables in period t
    endogenous_nbr::Int             # number of endogenous variables in period t
    exogenous_nbr::Int              # number of exogenous variables
    parameters::Vector{Float64}     # model parameters    
    steadystate::Vector{Float64}    # steady state of the model
end

struct SGIndices
    forward_in_system::Vector{Int} # indices of forward looking variables in system variables
    state_in_system::Vector{Int}  # indices of variables that belong both to state and system variables
end 

function SGIndices(endogenous_nbr, forward_variables, state_variables, system_variables)
    forward_in_system = findall(in(forward_variables), system_variables .- endogenous_nbr)
    state_in_system = findall(in(system_variables .- endogenous_nbr), state_variables)
    return SGIndices(
        forward_in_system,
        state_in_system
    )
end 

"""
    simulation_approximation_error(; context, grid::TasmanianSG, drawsnbr=10000, quantile_probability=0.999, sgws)

Computes the approximation error of the simulation by evaluating system equations at simulated state values.

## Arguments:
- `context::Context`: The simulation context containing model parameters, results, and settings.
- `grid::TasmanianSG`: The sparse grid used for approximating the policy function.
- `drawsnbr::Int` (default: `10000`): The number of simulation draws (excluding the first initialization period).
- `quantile_probability::Float64` (default: `0.999`): The quantile probability level for reporting error statistics.
- `sgws::SparseGridWorkspace`: A workspace containing system equations, dynamic state variables, and steady-state values.

## Process:
1. **Simulation of the System**:
   - Calls `simulate!` to generate `drawsnbr + 1` periods of simulated data.
   - Allocates memory for error calculations.

2. **Error Computation**:
   - Iterates over each simulation period (excluding the first) to compute errors.
   - Splits state variables into those evaluated at `t-1` and `t`.
   - Evaluates system equations using `sysOfEqs!` and stores errors.

3. **Error Analysis**:
   - Computes and logs the quantile and mean errors per equation.
   - Computes overall mean and quantile errors.

4. **Results Storage**:
   - Saves error statistics in `context.results.model_results[1].sparsegrids`.

## Returns:
- `Matrix{Float64}`: A matrix containing approximation errors for each equation across all periods.

## Notes:
- Uses `quantile(abs.(errors[i, :]), quantile_probability)` for error reporting.
- Modifies `context.results.model_results[1].sparsegrids` in place.
"""
function sparsegridapproximation(; context::Context=context,
                                 dimRef = -1,
                                 ftol = 1e-5,
                                 iterRefStart = 25,
                                 gridDepth = 2,
                                 gridOrder = 1,
                                 gridRule = "localp",
                                 maxiter = 300,
                                 maxRef = 1,
                                 maxRefLevel = gridDepth + maxRef,
                                 mcp = false,
                                 method = NewtonRaphson(),
                                 numstart = 0,
                                 savefreq = 10,
                                 scaleCorrInclude = [],
                                 scaleCorrExclude = [],
                                 show_trace = false,
                                 solver = nothing,
                                 surplThreshold= 1e-3,
                                 tol_ti = 1e-4,
                                 TT = 10000,
                                 typeRefinement = "classic",
                                 )
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    model = context.models[1]
    endogenous_nbr = model.endogenous_nbr
    exogenous_nbr = model.exogenous_nbr
    work = context.work
    (state_variables, predetermined_variables, system_variables,
     forward_equations_nbr, other_equations_nbr,
     forward_expressions_eqs, other_expressions_eqs,
     backward_block, forward_block, preamble_block) = make_block_functions(context)
    lb = Vector{Float64}(undef, 0)
    ub = Vector{Float64}(undef, 0)
    bmcps = Vector{Vector{Int}}(undef, 0)
    # must be consistent with SysOfEqs()
    block_eqs = union(forward_expressions_eqs, other_expressions_eqs)
    # computes lb, ub, bmcps
    block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables)
    equation_xref_list, variable_xref_list = xref_lists(context)
    params = context.work.params
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    nstates = length(state_variables)

    sgmodel = SGModel(repeat(steadystate, 3),
                      zeros(exogenous_nbr),
                      endogenous_nbr,
                      exogenous_nbr,
                      params,
                      steadystate
                      )

    ids = SGIndices(
        endogenous_nbr,
        model.i_fwrd_b,
        state_variables,
        system_variables
    )
    forward_system_variables = model.i_fwrd_b
    setdiff!(forward_system_variables, predetermined_variables .- endogenous_nbr)

    # Initialize error storage
    num_equations = length(system_variables)
    errorMat = zeros(num_equations, drawsnbr, replications)

    # Identify variable indices
    state_vars_tminus1 = filter(x -> x < endogenous_nbr, dynamic_state_variables)
    state_vars_t = filter(x -> x > endogenous_nbr, dynamic_state_variables)
    state_vars_t_adjusted = state_vars_t .- endogenous_nbr  # precompute once

    num_state_vars = length(dynamic_state_variables)
    num_state_vars_tminus1 = length(state_vars_tminus1)

    # Allocate reusable buffers
    x = Vector{Float64}(undef, num_state_vars)
    policy = Vector{Float64}(undef, num_equations)
    error = Vector{Float64}(undef, num_equations)

    # Compute equation errors across all simulation periods
    for r in 1:replications
        for t in 2:size(Y, 1)
            @inbounds for i in 1:num_state_vars_tminus1
                x[i] = Y[t-1, state_vars_tminus1[i], r]
            end
            @inbounds for i in 1:(num_state_vars - num_state_vars_tminus1)
                x[num_state_vars_tminus1 + i] = Y[t, state_vars_t_adjusted[i], r]
            end
            @inbounds for i in 1:num_equations
                policy[i] = Y[t, system_variables[i], r]
            end
            sysOfEqs!(error, policy, x, grid, sgws)
            errorMat[:,t-1,r] .= error
        end
    end

    return errorMat
end

function DDSGapproximation(opts::DDSGOptions; context::Context = context)
    # Extract all option fields
    @unpack k_max, ftol, gridDepth, gridOrder, gridRule, maxiter, maxRef, mcp, method, savefreq, scaleCorrInclude, scaleCorrExclude, show_trace, solver, surplThreshold, tol_ti, polUpdateWeight, maxIterEarlyStopping, drawsnbr, typeRefinement, initialPolGuess,
    quadrature, smolyak_level = opts

    # Extract model information
    model, endogenous_nbr, exogenous_nbr, params, steadystate,
    dynamic_state_variables, endogenous_variables, predetermined_variables, system_variables,
    backward_block, forward_block, preamble_block, system_block = initialize_sparsegrid_context(context)

    # Compute MCP constraints (bounds & complementarity conditions)
    lb, ub, bmcps = initialize_mcp_constraints(context, forward_block, backward_block, system_variables)

    #=
    endogenous = zeros(3*endogenous_nbr)
    exogenous_nbr = model.exogenous_nbr
    exogenous = zeros(exogenous_nbr)
    parameters = work.params

    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    =#
    
    polGuess = guess_policy(context, aNum, nPols, aPoints, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, lb, ub, sgmodel)
    
    loadNeededPoints!(grid0, polGuess)

    # params = context.work.params

    (nodes, weights) = monomial_power(exogenous_nbr)

    # Scale correction in the refinement process:
    scaleCorr = make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables])

    # Set the solver configuration options for sparse grid approximation.
    sgsolveroptions = SGSolverOptions(mcp, method, solver, ftol, show_trace)

    # Initialize sparse grid workspace
    sgws_ = initialize_sparsegrid_workspace(monomial, dynamic_state_variables, system_variables,
                                            bmcps, ids, lb, ub, backward_block, forward_block,
                                            preamble_block, system_block, n_nodes, gridDim,
                                            gridOut, endogenous_nbr, sgmodel, sgsolveroptions)

        # Index of current grid level to control the number of refinements
        ilev = gridDepth
        while ((getNumNeeded(grid1) > 0) && (ilev <=  maxRefLevel))
            grid1 = ti_step(grid1, polGuess1, pol, grid0, nPols, nodes, weights,
                            forward_equations_nbr, state_variables, system_variables, bmcps, dynamicws, model, 
                            predetermined_variables, forward_system_variables,
                            forward_expressions_eqs, other_expressions_eqs, ids, lb, ub,
                            backward_block, forward_block, preamble_block,
                            mcp, method, solver, ftol, show_trace, sgmodel)
            # We start the refinement process after a given number of iterations
            if (iter0 - 1 > iterRefStart)
                grid1, polGuess1 = refine(grid1, gridOut, scaleCorr, surplThreshold, dimRef, typeRefinement)
            end
            # Track the grid level
            ilev += 1
        end
        
        # Calculate (approximate) errors on tomorrow's policy grid
        metric, polGuess, grid0 = policy_update(grid0, grid1, length(system_variables))
 
        println("Iteration: $iter0, Grid pts: $(getNumPoints(grid0)), Level: $ilev, Metric: $metric")

        if (mod(iter0 + 1, savefreq) == 0)
            #            save_grid(grid0,iter0)
        end

        if (metric < tol_ti)
            break
        end
    else
        sgws[1] = sgws_
    end
#    Y = evaluateBatch(grid0, test_points)
#    display(Y')
    return (grid0, state_variables, system_variables)
end

function simulate!(context, grid, periods, policy_variables, state_variables, replications=1)
    model = context.models[1]
    params = context.work.params
    results = context.results.model_results[1]
    trends = results.trends
    perfect_foresight_ws = PerfectForesightWs(context, periods)
    shocks = perfect_foresight_ws.shocks
    initial_values = get_dynamic_initialvalues(context)
    endogenous_nbr = model.endogenous_nbr
    exogenous_nbr = model.exogenous_nbr
    steadystate = trends.endogenous_steady_state
    Sigma_e = model.Sigma_e
    chol_sigma_e_L = cholesky(Sigma_e).L
    x = reshape(perfect_foresight_ws.x, exogenous_nbr, periods)
        display(x)
    
    Y = Array{Float64}(undef, periods, endogenous_nbr + exogenous_nbr, replications)
    y = Vector{Float64}(undef, 3*endogenous_nbr)
    z = Vector{Float64}(undef, endogenous_nbr)
    random_shocks = Matrix{Float64}(undef, exogenous_nbr, periods)
    sv_buffer = Vector{Float64}(undef, length(state_variables))
    pv_buffer = Vector{Float64}(undef, length(policy_variables))
    @inbounds for r in 1:replications
        mul!(random_shocks, chol_sigma_e_L, randn(exogenous_nbr, periods))
        if !isempty(shocks)
            random_shocks .+= x
        end 
        y[1:endogenous_nbr] .= initial_values
        for p in 1:periods
            @views begin
                preamble_block(y, random_shocks[:, p], params, steadystate)
                sv_buffer .= y[state_variables]
                evaluateBatch!(pv_buffer, grid, sv_buffer)
                y[policy_variables] .= pv_buffer
                Y[p, 1:endogenous_nbr, r] .= y[endogenous_nbr .+ (1:endogenous_nbr)]
                Y[p, endogenous_nbr .+ (1:exogenous_nbr), r] .= random_shocks[:, p]

            end  
            circshift!(y, -endogenous_nbr)
        end
    end
    symboltable = context.symboltable
    varnames = vcat(Symbol.(get_endogenous(symboltable)),
                    Symbol.(get_exogenous(symboltable)))
    return AxisArray(Y, Undated(1):Undated(periods), varnames, 1:replications) 
end

"""
    Build initial grid with user defined parameter
"""
function set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain)
    # Generate the grid structure
    grid = Tasmanian.TasmanianSG(gridDim, gridOut, gridDepth)
    makeLocalPolynomialGrid!(grid, iOrder = gridOrder, sRule = gridRule)
    # Transform the domain
    setDomainTransform!(grid, gridDomain)
    # Get the points that require function values
    aPoints = getPoints(grid)
    # Get the number of points that require function values
    aNum = getNumPoints(grid)
    return grid, aPoints, aNum, gridDim, gridDepth, gridDomain
end

#=
# waiting for the linear approximation of a single block
function policy_guess(context, aNum, nPols, aPoints, i_polv, i_state)
    results = context.results.model_results[1]
    model = context.models[1]
    
    steadystate = results.trends.endogenous_steady_state
    y0 = zeros(model.endogenous_nbr)
    simulresults = Matrix{Float64}(undef, 2, model.endogenous_nbr)
    x = zeros(2, model.exogenous_nbr)
    A = zeros(model.endogenous_nbr, model.endogenous_nbr)
    B = zeros(model.endogenous_nbr, model.exogenous_nbr)
    Dynare.make_A_B!(A, B, model, results)
    polGuess = zeros(aNum, nPols)
    for (i,r) in enumerate(eachrow(aPoints))
        view(y0, i_state) .= r
        Dynare.simul_first_order!(simulresults, y0, x, steadystate, A, B, 1)
        @views polGuess[i, :] .= simulresults[2, i_polv]
    end
    return polGuess
end
# TODO insure that guess is in variable domain
=#

"""
    given the state variables, find a solution that is constant in t and t+1 in absence of future shocks
"""
function make_guess_system(residuals, T, state, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    fill!(exogenous, 0.0)
    dyn_endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(x)
        dyn_endogenous[system_variables] .= x
        dyn_endogenous[system_variables .+ endogenous_nbr] .= x
        other_block_([], res1, dyn_endogenous, exogenous, parameters, steadystate)
        forward_block_([], res2, dyn_endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function make_guess_system1(residuals, T, state, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    dyn_endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(residuals, x, params)
        dyn_endogenous[system_variables] .= x
        dyn_endogenous[system_variables .+ endogenous_nbr] .= x
        DFunctions.SparseDynamicResidTT!(T, dyn_endogenous, exogenous, parameters, steadystate)
        other_block_(T, res1, dyn_endogenous, exogenous, parameters, steadystate)
        forward_block_(T, res2, dyn_endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function guess_policy(context, aNum, nPols, aPoints, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, lb, ub, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    n =  length(system_variables)
    guess_values = zeros(n, aNum)
    @views x0 = copy(sgmodel.steadystate[system_variables .- endogenous_nbr])
    residuals = similar(x0)
    ws = DynamicWs(context)
    T = []
    for i in axes(aPoints, 2)
        @views state = aPoints[:, i]
        f1 = make_guess_system1(residuals, T, state, state_variables, system_variables,
                              forward_equations_nbr, other_equations_nbr, sgmodel)
        problem = NonlinearProblem(f1, x0, parameters)
        result = NonlinearSolve.solve(problem, NewtonRaphson(autodiff=AutoFiniteDiff()))
        if result.retcode == ReturnCode.Success
            @views guess_values[:, i] .= result.u
        else
            error("guess_policy failed")
        end
    end
    return guess_values
end

"""
    Evaluates residual of forward looking equations for each integration node
"""
function forward_looking_equation_residuals(x, T, state, grid, nodes, state_variables, system_variables, backward_block, forward_block, preamble_block, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    numNodes = size(nodes, 1)
    residuals = zeros(length(forward_block.variables), numNodes)
    y = zeros(length(dyn_endogenous), numNodes)
    evalPt = zeros(getNumDimensions(grid), numNodes)
    for i in 1:numNodes
        @views y[state_variables, i] .= state
        # 1) Determine t+1  variables using preamble
        @views preamble_block.set_endogenous_variables!([], y[endogenous_nbr + 1: end, i], nodes[i,:], parameters, steadystate)
        @views y[system_variables, i] .= x
        # 2) Extract next period's state variables
        @views evalPt[:, i] .= y[endogenous_nbr .+ state_variables, i]
    end
    # 3) Determine relevant variables within the expectations operator
    X = evaluateBatch(grid, evalPt)

    fill!(exogenous, 0.0)
    for i in 1:numNodes
        # setting future variables
        @views begin
            y[endogenous_nbr .+ system_variables, i] .= X[:, i]
            forward_block.get_residuals!([], residuals[:, i], y[:, i], exogenous, parameters, steadystate)
        end
    end
    return residuals
end

function expectation_equations!(expected_residuals, x, T, state, grid, nodes, weights, state_variables, system_variables, backward_block, forward_block, preamble_block, sgmodel)
    # get numNodes x numForwardEquations
    Integrand = forward_looking_equation_residuals(x, T, state, grid, nodes, state_variables, system_variables, backward_block, forward_block, preamble_block, sgmodel)
    for i in axes(expected_residuals, 1)
        @views expected_residuals[i] = dot(Integrand[i, :], weights)
    end
    return expected_residuals
end

function sysOfEqs(policy, T, state, grid, nPols, nodes, weights, state_variables, system_variables, bmcps, backward_block, forward_block, preamble_block, sgmodel)
    @unpack dyn_endogenous, exogenous, parameters, steadystate = sgmodel
    @views begin
        dyn_endogenous[state_variables] .= state
        dyn_endogenous[system_variables] .= policy 
    end
    nfwrd = length(forward_block.variables)
    res = zeros(length(policy))
    @views begin
        expectation_equations!(res[1:nfwrd], policy, T, state, grid, nodes, weights, state_variables, system_variables, backward_block, forward_block, preamble_block, sgmodel)
        backward_block.get_residuals!([], res[nfwrd + 1:end], dyn_endogenous, exogenous, parameters, steadystate)
    end
    reorder!(res, bmcps)
    return res
end

function sysOfEqs_derivatives_update!(J, Jexpected, Jother, policy, y, state, grid, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, bmcps)
    @unpack dyn_endogenous, exogenous, parameters, steadystate = sgmodel
    @views begin
        dyn_endogenous[state_variables] .= state
        dyn_endogenous[system_variables] .= policy 
    end 
    res = zeros(length(policy))

    fill!(exogenous, 0.0)
    @views begin
        expectation_equations_derivatives!(Jexpected, policy, state, grid, nodes, weights, forward_equations_nbr, state_variables, system_variables)
        other_block_derivatives(Jother, dyn_endogenous, exogenous, parameters, steadystate)
    end
    reorder!(res, bmcps)
    return res
end

function make_scaleCorr(scaleCorrInclude, scaleCorrExclude, policy_variables)
    @assert isempty(scaleCorrInclude) || isempty(scaleCorrExclude) "scaleCorrInclude and scaleCorrExclude are mutually exclusive"

    n = length(policy_variables)
    if !isempty(scaleCorrInclude)
        scaleCorr = [any(occursin.(v, scaleCorrInclude)) for v in policy_variables]
        return scaleCorr
    elseif !isempty(scaleCorrExclude)
        scaleCorr = [any(.!occursin.(v, scaleCorrExclude)) for v in policy_variables]
        return scaleCorr
    else
        return ones(n)
    end
end 

"""
   Time iteration step
"""
function ti_step(grid, X0, pol, gridZero, nPols, nodes, weights, forward_equations_nbr, 
                 state_variables, system_variables, bmcps, dynamicws, model,
                 predetermined_variables, forward_system_variables,
                 forward_expressions_eqs, other_expressions_eqs, ids, lb, ub,
                 backward_block, forward_block, preamble_block,
                 mcp, method, solver, ftol, show_trace, sgmodel)

    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    # Get the points that require function values
    aPoints1 = Tasmanian.getNeededPoints(grid)
    # Get the number of points that require function values
    aNumAdd = Tasmanian.getNumNeeded(grid)
    # Array for intermediate update step
    n = length(system_variables)
    Y = zeros(n, nstates)
    # TODO: check the case with non-autocorrelated shocks
    fill!(exogenous, 0.0)
    T = dynamicws.temporary_values
    J = zeros(n, n)
    fx = zeros(n)

    state = Vector{Float64}(undef, size(states, 1))
    x0 =  Vector{Float64}(undef, n)
    
    f(x) = sysOfEqs(x, T, state, gridZero, nPols, nodes, weights,
                    state_variables, system_variables, bmcps,
                    backward_block, forward_block, preamble_block, sgmodel)
    JA1!(J, x) = sysOfEqs_derivatives!(J, x, state, gridZero, nPols, nodes, weights, forward_equations_nbr,
                                       state_variables, system_variables, bmcps, dynamicws, model, 
                                       predetermined_variables, forward_system_variables,
                                       forward_expressions_eqs, other_expressions_eqs, ids, preamble_block,
                                       sgmodel)                   
    
    function f1!(fx, x)
        fx .= f(x)
    end

    function f2!(fx, x, p)
        fx .= f(x)
    end

    function JA2!(J, x, p)
        JA1!(J, x)
    end
    
    if isnothing(solver)
        solver = mcp ? NLsolver : NonlinearSolver 
    end
    
    # Time Iteration step
    if solver == NonlinearSolver
        NonLinearSolver_solve!(Y, f2!, JA2!, J, X0, states, state, x0, parameters, method, ftol, show_trace )
    elseif solver == NLsolver
        df = OnceDifferentiable(f1!, JA1!, X0, fx, J)
        NLsolve_solve!(Y, df, X0, lb, ub, states, X0, ftol, show_trace)
    elseif solver == PATHSolver
        PATHsolver_solve!(Y, f1!, JA1!, X0, lb, ub, states, X0, ftol, show_trace)
    else
        error("Sparsegrids: unknown solver")
    end

    # Add the new function values to grid1
    Tasmanian.loadNeededPoints!(grid, polInt)

    return grid
end

"""
    Grid refinement
"""
function refine(grid, nPols, scaleCorr, surplThreshold, dimRef, typeRefinement)

    # Get the points that require function values
    aNumLoad = Tasmanian.getNumLoaded(grid)
    # Scaling to only allow for those policies that are supposed to
    # determine the refinement process (in this case only the capital policies)
    scaleCorrMat = zeros(nPols, aNumLoad )
    scaleCorrMat .= repeat(scaleCorr, 1, aNumLoad)
    
    # Refine the grid based on the surplus coefficients
    Tasmanian.setSurplusRefinement!(grid, surplThreshold, iOutput=dimRef, sCriteria=typeRefinement, llfScaleCorrection=scaleCorrMat)
    if Tasmanian.getNumNeeded(grid) > 0

	# Get the new points and the number of points
	nwpts = Matrix(Tasmanian.getNeededPoints(grid))
	aNumNew = Tasmanian.getNumNeeded(grid)
        
	# We assign (for now) function values through interpolation#
	pol_guess = zeros(aNumNew, nPols)
	pol_guess = Tasmanian.evaluateBatch(grid, nwpts)
    else
        
	pol_guess = []

    end
    return grid, pol_guess
end

"""
   New grid construction
"""
function fresh_grid(gridDim, gridOut, gridDepth, gridDomain, gridOrder, gridRule)

    # Generate the grid structure
    grid = Tasmanian.TasmanianSG(gridDim, gridOut, gridDepth)
    Tasmanian.makeLocalPolynomialGrid!(grid, iOrder = gridOrder, sRule = gridRule)
    # Transform the domain
    Tasmanian.setDomainTransform!(grid, gridDomain)

    return grid
end

"""
   Checking convergence and updating policies
"""
function policy_update(gridOld, gridNew, nPols)

    # Get the points and the number of points from grid1
    aPoints2 = Tasmanian.getPoints(gridNew)
    aNumTot = Tasmanian.getNumPoints(gridNew)

    # Evaluate the grid points on both grid structures
    polGuessTr1 = Tasmanian.evaluateBatch(gridNew, Matrix(aPoints2))
    polGuessTr0 = Tasmanian.evaluateBatch(gridOld, Matrix(aPoints2))

    # 1) Compute the Sup-Norm

    metricAux = zeros(nPols)

    for imet in 1:nPols
        @views metricAux[imet] = maximum(abs.(polGuessTr0[imet, :] - polGuessTr1[imet, :]))
    end
    metricSup = maximum(metricAux)

    # 2) Compute the L2-Norm

    metricL2 = 0.0

    for imetL2 in 1:nPols
        @views metricL2 += sum(abs.(polGuessTr0[imetL2, :] - polGuessTr1[imetL2, :]).^2)
    end
    
    metricL2 = (metricL2/(aNumTot*nPols))^0.5

    metric = min(metricL2, metricSup)

    # Now update pol_guess and grid

    polGuess = zeros(nPols, aNumTot)
    for iupd in 1:nPols
        @views polGuess[iupd, :] = 0.5*polGuessTr0[iupd, :] + 0.5*polGuessTr1[iupd, :]
    end

    gridOld = Tasmanian.copyGrid(gridNew)

    return metric, polGuess, gridOld
end

"""
   Grid storage
"""
function save_grid(grid, iter)

    grid.write(data_location_smooth + "grid_iter_$i.txt")

    return
end

"""
    monomial nodes and weights
"""
function monomial_power(nShocks)
    # Number of integration nodes
    numNodes = 2*nShocks^2 + 1

    z0 = zeros(numNodes, nShocks)

    # Deviations in one dimension (note that the origin is row one)
    for i1 in 1:nShocks
        z0[(i1 - 1)*2 + 2,i1] = 1.0
        z0[(i1 - 1)*2 + 3,i1] = -1.0
    end
    
    i0 = 0
    # Deviations in two dimensions
    for i1 in 1:nShocks
        for i2 in i1+1:nShocks
            z0[2*nShocks+2+i0*4, i1] = 1.0
            z0[2*nShocks+3+i0*4, i1] = 1.0
            z0[2*nShocks+4+i0*4, i1] = -1.0
            z0[2*nShocks+5+i0*4, i1] = -1.0
            z0[2*nShocks+2+i0*4, i2] = 1.0
            z0[2*nShocks+3+i0*4, i2] = -1.0
            z0[2*nShocks+4+i0*4, i2] = 1.0
            z0[2*nShocks+5+i0*4, i2] = -1.0
            i0 += 1
        end
    end
    # Nodes
    IntNodes = zeros(numNodes, nShocks)
    IntNodes[2:nShocks*2+1, :] = z0[2:nShocks*2+1, :]*sqrt(2.0 + nShocks)
    IntNodes[nShocks*2+2:end, :] = z0[nShocks*2+2:end, :]*sqrt((2.0+nShocks)/2.0)

    # Weights
    IntWeights = zeros(numNodes)

    IntWeights[1] = 2.0/(2.0+nShocks)
    IntWeights[2:nShocks*2+1] .= (4-nShocks)/(2*(2+nShocks)^2)
    IntWeights[nShocks*2+2:end] .= 1.0/(nShocks+2)^2

    return IntNodes, IntWeights
end

function block_mcp!(lb, ub, bmcps, context, block_eqs, block_vars)
    block_eqs_ = copy(block_eqs)
    n = length(block_eqs_)
    resize!(lb, n)
    resize!(ub, n)
    fill!(lb, -Inf)
    fill!(ub, Inf)
    for m in context.models[1].mcps
        (eqn, var, op, expr) = m
        bvar = findfirst(var .== block_vars .- context.models[1].endogenous_nbr)
        beqn = findfirst(eqn .== block_eqs_)
        block_eqs_[[bvar, beqn]] .= block_eqs_[[beqn, bvar]]
        push!(bmcps, [beqn, bvar])
        boundary = Dynare.dynare_parse_eval(String(expr), context)
        if op[1] == '<'
            ub[bvar] = boundary
        elseif op[1] == '>'
            lb[bvar] = boundary
        else
            error("MCP operator must be <, <=, > or >=")
        end
    end
end

function forward_looking_equation_derivatives!(J, x, state, dynamicws, grid, nodes, weight,
                                               model, forward_equations_nbr, state_variables,
                                               system_variables, forward_system_variables, forward_expressions_eqs,
                                               predetermined_variables, ids, preamble_block, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    iDimensions = getNumDimensions(grid)
    iOutputs = getNumOutputs(grid)
    evalPt = zeros(iDimensions)
    T = dynamicws.temporary_values
    @views begin
        dyn_endogenous[state_variables] .= state
        preamble_block.set_endogenous_variables!(T, dyn_endogenous[endogenous_nbr + 1: end], nodes, parameters, steadystate)
        dyn_endogenous[system_variables] .= x
        evalPt .= dyn_endogenous[endogenous_nbr .+ state_variables]
    end

    # policy function Jacobian
    policy_j = zeros(iDimensions, Tasmanian.getNumOutputs(grid))
    Tasmanian.differentiate!(policy_j, grid, evalPt)
    # f(x) = evaluateBatch(grid, x)
    # policy_j = transpose(FiniteDiff.finite_difference_jacobian(f, evalPt[:,1]))

    # future policy Jacobian = policy Jacobian x derivatives of predetermined variables w.r. state
    future_policy_j = zeros(length(ids.state_in_system), length(forward_system_variables))
    copy_from_submatrix!(future_policy_j, policy_j, ids.state_in_system, ids.forward_in_system)
    X = zeros(iOutputs)
    # replace with Tasmanian.evaluate()
    X = evaluateBatch(grid, evalPt)
    jacobian = Matrix(get_dynamic_jacobian!(
        dynamicws,
        parameters,
        dyn_endogenous,
        exogenous,
        steadystate,
        model,
        2,
    ))
    n_forward_system_variables = length(forward_system_variables)
    forward_j = Matrix{Float64}(undef, forward_equations_nbr, n_forward_system_variables)
    copy_from_submatrix!(forward_j, jacobian, forward_expressions_eqs, forward_system_variables .+ 2*endogenous_nbr)
    M1 = Matrix{Float64}(undef, forward_equations_nbr, length(ids.state_in_system))
    mul!(M1, forward_j, transpose(future_policy_j))
    jacobian[forward_expressions_eqs, state_variables[ids.state_in_system] .+ endogenous_nbr] .+= M1
    @views lmul!(weight, jacobian[forward_expressions_eqs, system_variables])
    J .+= jacobian[forward_expressions_eqs, system_variables]
    return J
end

function forward_equations_derivatives(x, nodes, weights, dynamicws, model, 
                                       forward_equations_nbr, state, grid, state_variables,
                                       system_variables, forward_system_variables, forward_expressions_eqs, 
                                       predetermined_variables, ids, preamble_block, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    jacobian = get_dynamic_jacobian!(
        dynamicws,
        parameters,
        dyn_endogenous,
        nodes[1, :],
        steadystate,
        model,
        2,
    )
    
    J = zeros(forward_equations_nbr, length(system_variables))
    for i in axes(nodes, 1)
        forward_looking_equation_derivatives!(J, x, state, dynamicws, grid, nodes[i, :], weights[i],
                                              model, forward_equations_nbr, state_variables,
                                              system_variables, forward_system_variables, forward_expressions_eqs,
                                              predetermined_variables, ids, preamble_block, sgmodel)
    end
    return J
end 
    
function sysOfEqs_derivatives!(J, x, state, grid,
                               nPols, nodes, weights, forward_equations_nbr,
                               state_variables, system_variables, bmcps, dynamicws, model, 
                               predetermined_variables, forward_system_variables,
                               forward_expressions_eqs, other_expressions_eqs, ids, preamble_block, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    @views begin
        dyn_endogenous[state_variables] .= state
        dyn_endogenous[system_variables] .= x 
    end 
    
    Jfwrd = forward_equations_derivatives(x, nodes, weights, dynamicws, model, 
                                  forward_equations_nbr, state, grid, state_variables,
                                  system_variables, forward_system_variables, forward_expressions_eqs, 
                                  predetermined_variables, ids, preamble_block, sgmodel)
    copy_to_submatrix!(J, 1:forward_equations_nbr, 1:length(system_variables), Jfwrd)

    jacobian = get_dynamic_jacobian!(
        dynamicws,
        parameters,     
        dyn_endogenous,
        exogenous,
        steadystate,
        model,
        2,
    )
    copy_submatrix(J, forward_equations_nbr .+ (1:length(other_expressions_eqs)), 1:length(system_variables), 
                   jacobian, other_expressions_eqs, system_variables)
    reorder_rows!(J, bmcps)
end

function copy_submatrix(dest, rowsd, colsd, src, rowss, colss)
    @assert length(rowsd) == length(rowss)
    @assert length(colsd) == length(colss)
    for (j1, j2) in enumerate(colsd)
        k = colss[j1]
        for (i1, i2) in enumerate(rowsd)
            dest[i2, j2] = src[rowss[i1], k]
        end
    end
    return dest
end 

function copy_from_submatrix!(dest, src, rows, cols)
    for j in axes(dest, 2)
        k = cols[j]
        for i in axes(dest, 1)
            dest[i, j] = src[rows[i], k]
        end
    end
    return dest
end 

function copy_to_submatrix!(dest, rows, cols, src)
    for j in axes(src, 2)
        k = cols[j]
        for i in axes(src, 1)
            dest[rows[i], k] = src[i, j]
        end
    end
    return dest
end 

function reorder_rows!(x::AbstractMatrix, permutations)
    for c in axes(x, 2)
        @views reorder!(x[:, c], permutations)
    end
end

function add_to_submatrix!(dest, rows, cols, src)
    for j in axes(src, 2)
        k = cols[j]
        for i in axes(src, 1)
            dest[rows[i], k] += src[i, j]
        end
    end
    return dest
end 

function PATHsolver_solve!(Y, f!, JA!, X0, lb, ub, states, x0, ftol, show_trace)
    for i in axes(states, 2)
        @views begin
            state = states[:, i]
            x0 .= X0[:, i]
        end
        (status, results, info) = mcp_solve!(PathNLS(), f!, JA!, JJ, lb, ub, pol, silent=true, convergence_tolerance=1e-4)
        if status != 1
            @show status
            @show res.info
            error("sparsegrids: solution update failed")
        end
        @views Y[:, i] .= results
    end
end

function NLsolve_solve!(Y, df, X0, lb, ub, states, x0, ftol, show_trace)
    for i in axes(states, 2)
        @views begin
            state = states[:, i]
            x0 .= X0[:, i]
        end
        res = NLsolve.mcpsolve(df, lb, ub, x0, ftol = ftol, show_trace = show_trace)
        if !res.f_converged
            @show res.f_converged
            error("sparsegrids: solution update failed")
        end
        @views Y[:, i] .= res.zero
    end
end

function NonLinearSolver_solve!(Y, f, JA, J, X0, states, state, x0, parameters, method, ftol, show_trace )
    nlf = NonlinearFunction(f, jac = JA, jac_prototype=J)

    for i in axes(states, 2)
        @views begin 
            state .= states[:, i]
            x0 .= X0[:, i]
        end
        nlp = NonlinearProblem(nlf, x0, parameters)
        res = NonlinearSolve.solve(nlp, method, show_trace = Val(show_trace))
        @views Y[:, i] .= res.u
    end
end
