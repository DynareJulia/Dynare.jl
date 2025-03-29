include("blocks.jl")
include("block_firstorder.jl")
include("monomial.jl")
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

export SGapproximation, simulate!, simulation_approximation_error!, DDSGapproximation, DDSGOptions, SGOptions, UserPolicyGuess

"""
    SGapproximation(; context::Context, ...)

Performs adaptive sparse grid approximation for solving high-dimensional 
dynamic stochastic models.

# Arguments
- `context::Context`: Dynare model context with system equations.

"""
function SGapproximation(opts::SGOptions; context::Context=context)
    # Extract all option fields
    @unpack dimRef, ftol, iterRefStart, gridDepth, gridOrder, gridRule, maxiter, maxRef, mcp, method, savefreq, scaleCorrInclude, scaleCorrExclude, show_trace, solver, surplThreshold, tol_ti, polUpdateWeight, maxIterEarlyStopping, drawsnbr, typeRefinement, initialPolGuess = opts

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

    # Create the initial sparse grid (merged function)
    grid0, aPoints, aNum = initialize_sparse_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain)

    # Constructs a monomial quadrature rule
    monomial = MonomialPowerIntegration(exogenous_nbr)
    n_nodes = length(monomial.nodes)
    
    # Make a deep copy of the initial grid to perform the time iteration on
    grid = copyGrid(grid0)
 
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

    # Get a guess for the policy function using first-order perturbation
    polGuess = initialize_policy_approximation(context, sgws[1], aNum, aPoints, initialPolGuess)

    # Load the values of the policy function guess on the grid
    loadNeededPoints!(grid, polGuess)

    # Sparse time-iteration
    grid, polGuess, average_time, iter = sparse_grid_time_iteration!(
        grid, grid0, polGuess, sgws, scaleCorr, surplThreshold, dimRef,
        typeRefinement, maxiter, maxRef, iterRefStart, maxIterEarlyStopping,
        tol_ti, polUpdateWeight,savefreq
    )

    # Save the results
    results = context.results.model_results[1].sparsegrids
    results.iteration_number = iter
    results.average_iteration_time = average_time
    results.drawsnbr = drawsnbr
    results.ftol = ftol
    results.grid = grid
    results.gridDepth = gridDepth
    results.gridOrder = gridOrder
    results.gridRule = gridRule
    results.iterRefStart = iterRefStart
    results.maxRef = maxRef
    results.mcp = mcp
    results.method = method
    results.solver = solver
    results.surplThreshold = surplThreshold

    return (grid, sgws[1])
end

"""
    simulate!(; context, grid::TasmanianSG, periods=1000, replications=1, sgws)

Simulates the dynamic system over a given number of periods and replications using a sparse grid approximation.

## Arguments:
- `context::Context`: The simulation context containing model parameters, results, and settings.
- `grid::TasmanianSG`: The sparse grid used for approximating the policy function.
- `periods::Int` (default: `1000`): The number of time periods for which the simulation is run.
- `replications::Int` (default: `1`): The number of Monte Carlo replications.
- `sgws::SparseGridWorkspace`: A workspace structure containing model variables, steady state values, and system settings.

## Process:
1. **Model and Workspace Extraction**:
   - Extracts the first model from `context.models` and retrieves relevant parameters.
   - Extracts variables related to the sparse grid workspace (`sgws`), including steady-state values and system variables.
   - Retrieves trends from `context.results.model_results[1]`.

2. **Initialization**:
   - Constructs a `PerfectForesightWs` object to handle exogenous shock simulations.
   - Computes the initial state values using `get_dynamic_initialvalues(context)`.
   - Precomputes the Cholesky decomposition of the exogenous shock covariance matrix `Sigma_e`.

3. **Memory Allocation**:
   - Preallocates arrays for storing simulation results:
     - `Y`: Stores the simulated values of endogenous and exogenous variables.
     - `y`: Stores current and lagged endogenous states.
     - `z`: Temporary storage for endogenous variables.
     - `random_shocks`: Holds the generated exogenous shocks.
     - `sv_buffer` and `pv_buffer`: Buffers for state variables and policy function evaluations.

4. **Monte Carlo Simulation**:
   - For each replication:
     - Generates exogenous shocks via a Cholesky transformation.
     - Adds deterministic shocks if present.
     - Initializes the endogenous state vector.
     - Iterates over `periods`, updating state variables and computing policy function evaluations using `interpolate!`.
     - Stores results in `Y`.

5. **Post-processing**:
   - Constructs variable names from the model's symbol table.
   - Returns an `AxisArray` containing the simulated time series.

## Returns:
- `AxisArray{Float64,3}`
"""
function simulate!(;
    context::Context = context, 
    grid::Union{TasmanianSG,DDSG} = grid, 
    periods::Int = 1000,
    replications::Int = 1,
    sgws::SparsegridsWs = sgws
)
    model = context.models[1]
    params = context.work.params
    
    # Extract relevant model parameters
    @unpack endogenous_nbr, exogenous_nbr, Sigma_e = model
    @unpack dynamic_state_variables, preamble_block, system_variables = sgws
    @unpack steadystate = sgws.sgmodel

    # Precompute Cholesky of innovation covariance
    chol_sigma_e_L = cholesky(Sigma_e).L

    # Retrieve the user-specified innovation paths
    exogenous_shocks = reshape(PerfectForesightWs(context, periods).x, exogenous_nbr, periods)
    shock_flag = !isempty(exogenous_shocks)

    # Useful buffer arrays
    y = Vector{Float64}(undef, 3 * endogenous_nbr)
    sv_buffer = Vector{Float64}(undef, length(dynamic_state_variables))
    pv_buffer = Vector{Float64}(undef, length(system_variables))
    random_shocks = Matrix{Float64}(undef, exogenous_nbr, periods)

    # Allocate the output array
    Y = Array{Float64}(undef, periods, endogenous_nbr + exogenous_nbr, replications)

    # Retrieve useful variables for future function calls
    ws = DynamicWs(context)
    T = ws.temporary_values
    initial_values = get_dynamic_initialvalues(context)

    for r in 1:replications
        # Generate the full shock path
        randn!(random_shocks)
        mul!(random_shocks, chol_sigma_e_L, random_shocks)
        if shock_flag
            @. random_shocks += exogenous_shocks
        end

        # Initialize endogenous state
        copyto!(y, 1, initial_values, 1, endogenous_nbr)

        for p in 1:periods
            # Compute endogenous variables
            preamble_block.set_endogenous_variables!(T, y, view(random_shocks,:,p), params, steadystate)
            @views copyto!(sv_buffer, y[dynamic_state_variables])

            # Evaluate policy function
            interpolate!(pv_buffer, grid, sv_buffer)

            # Shift state
            circshift!(y, -endogenous_nbr)
            @views begin
                y[system_variables] .= pv_buffer
                # Store results
                Y[p, 1:endogenous_nbr, r] .= y[1:endogenous_nbr]
                Y[p, endogenous_nbr .+ (1:exogenous_nbr), r] .= random_shocks[:, p]
            end
        end
    end

    # Retrieve variable names
    varnames = vcat(Symbol.(get_endogenous(context.symboltable)), Symbol.(get_exogenous(context.symboltable)))

    return AxisArray(Y, Undated(1):Undated(periods), varnames, 1:replications)
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
    errors = zeros(num_equations, drawsnbr, replications)

    # Identify variable indices
    state_vars_tminus1 = filter(x -> x < endogenous_nbr, dynamic_state_variables)
    state_vars_t = filter(x -> x > endogenous_nbr, dynamic_state_variables)
    state_vars_t_adjusted = state_vars_t .- endogenous_nbr  # precompute once

    num_state_vars = length(dynamic_state_variables)
    num_state_vars_tminus1 = length(state_vars_tminus1)

    # Allocate reusable buffers
    x = Vector{Float64}(undef, num_state_vars)
    policy = Vector{Float64}(undef, num_equations)

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
            sysOfEqs!(errors[:, t-1, r], policy, x, grid, sgws)
        end
    end

    return errors
end

function DDSGapproximation(opts::DDSGOptions; context::Context = context)
    # Extract all option fields
    @unpack k_max, ftol, gridDepth, gridOrder, gridRule, maxiter, maxRef, mcp, method, savefreq, scaleCorrInclude, scaleCorrExclude, show_trace, solver, surplThreshold, tol_ti, polUpdateWeight, maxIterEarlyStopping, drawsnbr, typeRefinement, initialPolGuess = opts

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
    monomial = MonomialPowerIntegration(exogenous_nbr)
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

    # Guess for the policy function
    M, N, _ = block_first_order_approximation(context, sgws[1])
    polGuess = X -> guess_policy(context, size(X,2), X, sgws[1], M, N, initialPolGuess)

    # Initialization
    ddsg = DDSG(gridDim, gridOut, gridDepth, gridDepth+maxRef, k_max; order=gridOrder, rule=gridRule, domain=gridDomain)
    DDSG_init!(ddsg)
    DDSG_build!(ddsg, polGuess, ddsg.centroid; refinement_type=typeRefinement, refinement_tol=surplThreshold, scale_corr_vec=scaleCorr)

    # Time-iteration
    ddsg, average_time, iter = ddsg_time_iteration!(
        ddsg,
        scaleCorr,
        surplThreshold,
        typeRefinement,
        maxiter,
        maxIterEarlyStopping,
        tol_ti,
        polUpdateWeight,
        savefreq,
        context,
        sgws
    )

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