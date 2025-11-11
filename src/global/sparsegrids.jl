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
function simulation_approximation_error!(;
    context::Context,
    grid::Union{TasmanianSG,DDSG},
    drawsnbr::Int = 11000,
    replications::Int = 100,
    sgws::SparsegridsWs
)
    # Extract workspace parameters
    @unpack dynamic_state_variables, system_block, system_variables = sgws
    @unpack endogenous_nbr = sgws.sgmodel

    # Simulate system dynamics
    Y = simulate!(context = context, grid = grid, periods = drawsnbr + 1, replications = replications, sgws = sgws)

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

    # Initialize SGModel & SGIndices
    sgmodel = initialize_sgmodel(steadystate, endogenous_nbr, exogenous_nbr, params)
    ids, forward_system_variables = initialize_sg_indices(model, dynamic_state_variables, system_variables, predetermined_variables)

    # Compute grid dimensions
    gridDim, gridOut = get_grid_dimensions(dynamic_state_variables, system_variables)

    # Compute domain limits for the grid
    gridDomain = get_grid_domain(context.work.limits, ids, context, gridDim)

    # Constructs a monomial quadrature rule
    monomial = (quadrature == :smolyakgh) ?
               SmolyakGHIntegration(exogenous_nbr, smolyak_level) :
               MonomialPowerIntegration(exogenous_nbr)
    n_nodes = length(monomial.nodes)

    # Scale correction in the refinement process:
    scaleCorr = make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables])

    # Set the solver configuration options for sparse grid approximation.
    sgsolveroptions = SGSolverOptions(mcp, method, solver, ftol, show_trace)

    # Initialize sparse grid workspace
    sgws_ = initialize_sparsegrid_workspace(monomial, dynamic_state_variables, system_variables,
                                            bmcps, ids, lb, ub, backward_block, forward_block,
                                            preamble_block, system_block, n_nodes, gridDim,
                                            gridOut, endogenous_nbr, sgmodel, sgsolveroptions)

    # Handle multithreading
    sgws=Vector{SparsegridsWs}(undef, Threads.nthreads())
    if Threads.nthreads() > 1
        BLAS.set_num_threads(1)
        for i in 1:Threads.nthreads()
            sgws[i] = deepcopy(sgws_)
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
function make_guess_system(x, residuals, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
    endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(x)
        endogenous[system_variables] .= x
        endogenous[system_variables .+ endogenous_nbr] .= x
        other_block(res1, endogenous, exogenous, parameters, steadystate)
        forward_block(res2, endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function guess_policy(context, aNum, nPols, aPoints, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
    n =  length(system_variables)
    guess_values = zeros(n, aNum)
    @views x = copy(steadystate[system_variables .- endogenous_nbr])
    residuals = similar(x)
    for i in axes(aPoints, 2)
        @views state = aPoints[:, i]
        f = make_guess_system(x, residuals, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
        result = Dynare.nlsolve(f, x; method = :robust_trust_region, show_trace = false, ftol = cbrt(eps()), iterations = 50)
        @views guess_values[:, i] .= result.zero
    end
    return guess_values
end

"""
    Evaluates residual of forward looking equations for each integration node
"""
function forward_looking_equation_residuals(x, state, params, grid, nodes, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
    numNodes = size(nodes, 1)
    residuals = zeros(numNodes, forward_equations_nbr)
    y = zeros(numNodes, 3*endogenous_nbr)
    evalPt = zeros(getNumDimensions(grid), numNodes)
    for i in 1:numNodes
        @views y[i, state_variables] .= state
        # 1) Determine t+1  variables using preamble
        @views preamble_block(y[i, endogenous_nbr + 1: end], nodes[i,:], params, steadystate)
        @views y[i, system_variables] .= x
        # 2) Determine next period's state variables
        @views evalPt[:, i] .= y[i, endogenous_nbr .+ state_variables]
    end
    # 3) Determine relevant variables within the expectations operator
    X = evaluateBatch(grid, evalPt)
    
    resid = zeros(forward_equations_nbr)
    for i in 1:numNodes
        @views y[i, endogenous_nbr .+ system_variables] .= X[:, i]
        forward_block(resid, y[i, :], zeros(exogenous_nbr), params, steadystate)
        residuals[i,:] .= resid
    end
    return residuals
end

function expectation_equations!(expected_residuals, x, state, params, grid0, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
    # get numNodes x numForwardEquations
    Integrand = forward_looking_equation_residuals(x, state, params, grid0, nodes, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
    for i in axes(expected_residuals, 1)
        @views expected_residuals[i] = dot(Integrand[:, i], weights)
    end
    return expected_residuals
end

function sysOfEqs(policy, y, state, grid, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables)
    @views begin
        y[state_variables] .= state
        y[system_variables] .= policy 
    end 
    res = zeros(length(policy))

    @views begin
        expectation_equations!(res[1:forward_equations_nbr], policy, state, params, grid, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
        other_block(res[forward_equations_nbr + 1:end], y, zeros(exogenous_nbr), params, steadystate)
    end 
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
function ti_step(grid, pol_guess, gridZero, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, y)

    # Get the points that require function values
    aPoints1 = Tasmanian.getNeededPoints(grid)
    # Get the number of points that require function values
    aNumAdd = Tasmanian.getNumNeeded(grid)
    # Array for intermediate update step
    polInt = zeros(length(system_variables), aNumAdd)
    let state
        fnl(x) = sysOfEqs(x, y, state, gridZero, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables)
        
        # Time Iteration step
        for ii1 in 1:aNumAdd
            @views state = aPoints1[:, ii1]
            @views pol = pol_guess[:, ii1]
            res = Dynare.nlsolve(fnl, pol, method = :robust_trust_region, show_trace = false)
            polInt[:, ii1] .= res.zero
        end
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

    # Save the results
    results = context.results.model_results[1].sparsegrids
    results.iteration_number = iter
    results.average_iteration_time = average_time
    results.drawsnbr = drawsnbr
    results.ftol = ftol
    results.gridDepth = gridDepth
    results.gridOrder = gridOrder
    results.gridRule = gridRule
    results.maxRef = maxRef
    results.mcp = mcp
    results.method = method
    results.solver = solver
    results.surplThreshold = surplThreshold
    return (ddsg, sgws[1])
end
