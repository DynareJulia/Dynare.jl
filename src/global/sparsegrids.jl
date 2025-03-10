include("blocks.jl")
include("block_firstorder.jl")
include("monomial.jl")
include("helper_functions.jl")
include("SG.jl")
include("DDSG.jl")

using AxisArrays: AxisArray
using DualNumbers
using Distributions
using FiniteDiff
using LinearAlgebra
import NLsolve
using NonlinearSolve
using Roots
using Tasmanian

export sparsegridapproximation, simulate!, simulation_approximation_error!, DDSGapproximation

"""
    sparsegridapproximation(; context::Context, ...)

Performs adaptive sparse grid approximation for solving high-dimensional 
dynamic stochastic models.

# Arguments
- `context::Context`: Dynare model context with system equations.
- `dimRef = -1`: Outputs that are considered in the refinement process (-1 implies that all outputs are considered)
- `ftol = 1e-5`: Convergence tolerance for the nonlinear solver,
- `iterRefStart = 25`: Iteration at which the refinement starts,
- `gridDepth = 2`: Initial sparse grid depth,
- `gridOrder = 1`: Polynomial order for interpolation,
- `gridRule = "localp": Type of base functions`,
- `maxiter = 300`: Maximum iterations for time iteration,
- `maxRef = 1`: Number of maximum refinements
- `maxRefLevel = gridDepth + maxRef`: Maximum Level of ASG,
- `mcp = false`: Enables handling of occasionally binding constraints using a Mixed Complementarity Problem (MCP) solver,
- `method = NewtonRaphson()`: Solver for nonlinear system,
- `savefreq = 10`: Frequency for grid saving,
- `scaleCorrInclude = []`: Specifies which policy variables should have scale correction applied,
- `scaleCorrExclude = []`: Specifies which policy variables should be excluded from scale correction,
- `show_trace = false`: Enables detailed solver iteration logs for debugging nonlinear convergence issues,
- `solver = mcp ? NLsolver : NonlinearSolver`: Specifies the nonlinear solver for equilibrium equations, supporting standard (NonlinearSolve) or mixed complementarity (NLsolve) problems,
- `surplThreshold= 1e-3`: Threshold for grid refinement,
- `tol_ti = 1e-4`: Convergence criterion for time iteration,
- `polUpdateWeight = 0.5` : Weight of the current-step policy function when computing the updated policy function. The weight of the previous-step policy function is `1-polUpdateWeight`
- `maxIterEarlyStopping = 1` : Number of iterations after which TI stops after TI convergence measure starts increasing
- `drawsnbr = 10000`: Number of random draws for the error computation,
- `typeRefinement = "classic"`,
- `initialPolGuess::UserPolicyGuess = `
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
                                 savefreq = 10,
                                 scaleCorrInclude = [],
                                 scaleCorrExclude = [],
                                 show_trace = false,
                                 solver = mcp ? NLsolver : NonlinearSolver,
                                 surplThreshold= 1e-3,
                                 tol_ti = 1e-4,
                                 polUpdateWeight = 0.5,
                                 maxIterEarlyStopping = 0,
                                 drawsnbr = 10000,
                                 typeRefinement = "classic",
                                 initialPolGuess::UserPolicyGuess = UserPolicyGuess(),
                                 )

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
    sgoptions = SGOptions(mcp, method, solver, ftol, show_trace)

    # Initialize sparse grid workspace
    sgws_ = initialize_sparsegrid_workspace(monomial, dynamic_state_variables, system_variables,
                                            bmcps, ids, lb, ub, backward_block, forward_block,
                                            preamble_block, system_block, n_nodes, gridDim,
                                            gridOut, endogenous_nbr, sgmodel, sgoptions)

    # Handle multithreading
    if Threads.nthreads() > 1
        BLAS.set_num_threads(1)
        sgws = [deepcopy(sgws_) for _ in 1:Threads.nthreads()]
    else
        sgws = [sgws_]
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
    grid::TasmanianSG = grid, 
    periods = 1000, 
    replications = 1, 
    sgws = sgws
)
    model = context.models[1]
    params = context.work.params
    
    # Extract relevant model parameters
    @unpack endogenous_nbr, exogenous_nbr, Sigma_e = model
    @unpack dynamic_state_variables, preamble_block, system_variables = sgws
    @unpack steadystate = sgws.sgmodel

    # Initialize necessary structures
    perfect_foresight_ws = PerfectForesightWs(context, periods)
    shocks = perfect_foresight_ws.shocks
    initial_values = get_dynamic_initialvalues(context)
    chol_sigma_e_L = cholesky(Sigma_e).L

    # Allocate arrays
    exogenous_shocks = reshape(perfect_foresight_ws.x, exogenous_nbr, periods)
    Y = Array{Float64}(undef, periods, endogenous_nbr + exogenous_nbr, replications)
    y = Vector{Float64}(undef, 3 * endogenous_nbr)
    random_shocks = Matrix{Float64}(undef, exogenous_nbr, periods)
    sv_buffer = Vector{Float64}(undef, length(dynamic_state_variables))
    pv_buffer = Vector{Float64}(undef, length(system_variables))

    ws = DynamicWs(context)
    T = ws.temporary_values

    @inbounds for r in 1:replications
        # Generate shocks
        mul!(random_shocks, chol_sigma_e_L, randn(exogenous_nbr, periods))
        if !isempty(shocks)
            random_shocks .+= exogenous_shocks
        end

        # Initialize state
        y[1:endogenous_nbr] .= initial_values

        for p in 1:periods
            @views begin
                # Compute endogenous variables
                preamble_block.set_endogenous_variables!(T, y, random_shocks[:, p], params, steadystate)

                # Store dynamic state variables
                sv_buffer .= y[dynamic_state_variables]

                # Evaluate policy function
                interpolate!(pv_buffer, grid, sv_buffer)

                # Update state variables
                circshift!(y, -endogenous_nbr)
                y[system_variables] .= pv_buffer

                # Store results
                Y[p, 1:endogenous_nbr, r] .= y[1:endogenous_nbr]
                Y[p, endogenous_nbr .+ (1:exogenous_nbr), r] .= random_shocks[:, p]
            end
        end
    end

    # Retrieve variable names
    symboltable = context.symboltable
    varnames = vcat(Symbol.(get_endogenous(symboltable)), Symbol.(get_exogenous(symboltable)))

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
    grid::TasmanianSG, 
    drawsnbr=11000, 
    burnin=1000,
    quantile_probability=0.999, 
    sgws
)
    # Extract workspace parameters
    @unpack dynamic_state_variables, system_block, system_variables = sgws
    @unpack endogenous_nbr = sgws.sgmodel
    equations = system_block.equations

    # Simulate system dynamics
    Y = simulate!(context = context, grid = grid, periods = drawsnbr + 1, replications = 1, sgws = sgws)

    # Initialize error storage
    num_equations = length(system_variables)
    errors = zeros(num_equations, drawsnbr)

    # Identify variable indices
    state_vars_tminus1 = filter(x -> x < endogenous_nbr, dynamic_state_variables)
    state_vars_t = filter(x -> x > endogenous_nbr, dynamic_state_variables)
    
    num_state_vars = length(dynamic_state_variables)
    num_state_vars_tminus1 = length(state_vars_tminus1)

    # Allocate buffer for state variable inputs
    x = zeros(num_state_vars)

    # Compute equation errors across all simulation periods
    @inbounds for i in 2:size(Y, 1)  # Skip first period (used for initialization)
        @views begin
            x[1:num_state_vars_tminus1] .= Y[i-1, state_vars_tminus1]  # Use t-1 values
            x[num_state_vars_tminus1 + 1:num_state_vars] .= Y[i, state_vars_t .- endogenous_nbr]  # Use t values
            
            sysOfEqs!(errors[:, i-1], Y[i, system_variables], x, grid, sgws)
        end
    end

    # Compute per-equation quantile and mean errors
    equation_quantile_errors = zeros(num_equations)
    equation_average_errors = zeros(num_equations)

    println("\nMaximum approximation error by equation:")
    @inbounds for i in 1:num_equations
        equation_quantile_errors[i] = quantile(abs.(errors[i, :]), quantile_probability)
        println("Equation $(equations[i]): absolute error $(100 * quantile_probability)% quantile: $(equation_quantile_errors[i])")
    end

    println("\nMean approximation error by equation:")
    @inbounds for i in 1:num_equations
        equation_average_errors[i] = mean(abs.(errors[i, :]))
        println("Equation $(equations[i]): mean absolute error: $(equation_average_errors[i])")
    end

    # Compute overall statistics
    average_error = mean(abs.(errors[:,burnin+1:end]))
    quantile_error = quantile(abs.(vec(errors[:,burnin+1:end])), quantile_probability)
    max_error = maximum(abs.(errors[:,burnin+1:end]))

    println("Overall absolute error $(100 * quantile_probability)% quantile: $(quantile_error)")
    println("Overall mean absolute error: $(average_error)")
    println("Overall max error: $max_error")

    # Store results in the context
    results = context.results.model_results[1].sparsegrids
    results.average_error = average_error
    results.max_error = max_error
    results.drawsnbr = drawsnbr
    results.burnin = burnin
    results.equation_average_errors = equation_average_errors
    results.equation_quantile_errors = equation_quantile_errors
    results.quantile_probability = quantile_probability
    results.quantile_error = quantile_error

    return errors
end

"""
    DDSGapproximation(; 
        context::Context=context,
        k_max = 1,
        ftol = 1e-5,
        gridDepth = 2,
        gridOrder = 1,
        gridRule = "localp",
        maxiter = 300,
        maxRef = 0,
        maxRefLevel = gridDepth + maxRef,
        mcp = false,
        method = NewtonRaphson(),
        savefreq = 10,
        scaleCorrInclude = [],
        scaleCorrExclude = [],
        show_trace = false,
        solver = mcp ? NLsolver : NonlinearSolver,
        surplThreshold= 1e-3,
        tol_ti = 1e-4,
        polUpdateWeight = 0.5,
        maxIterEarlyStopping = 0,
        drawsnbr = 10000,
        typeRefinement = "classic",
        initialPolGuess::UserPolicyGuess = UserPolicyGuess()
    )

Performs adaptive sparse grid approximation using **Dimensionally Decomposed Sparse Grids (DDSG)** 
for solving high-dimensional dynamic stochastic models.

# Arguments
- `context::Context`: Dynare model context with system equations.
- `k_max::Int`: Maximum order of interaction terms in DDSG decomposition.
- `ftol::Float64`: Function tolerance for solver convergence.
- `gridDepth::Int`: Initial depth of the sparse grid.
- `gridOrder::Int`: Polynomial interpolation order.
- `gridRule::String`: Type of basis functions used in sparse grids (e.g., `"localp"`).
- `maxiter::Int`: Maximum number of iterations for the time iteration process.
- `maxRef::Int`: Maximum number of refinement levels for the sparse grid.
- `maxRefLevel::Int`: Effective maximum grid refinement level (`gridDepth + maxRef`).
- `mcp::Bool`: Whether to use a **Mixed Complementarity Problem (MCP)** solver for constraints.
- `method`: Nonlinear solver method (e.g., `NewtonRaphson()`).
- `savefreq::Int`: Frequency for saving grid states during iterations.
- `scaleCorrInclude::Vector{Int}`: Indices of policy variables for which scale correction is applied.
- `scaleCorrExclude::Vector{Int}`: Indices of policy variables to be excluded from scale correction.
- `show_trace::Bool`: Enables verbose solver output for debugging.
- `solver`: Solver type (`NLsolver` for MCP problems, `NonlinearSolver` otherwise).
- `surplThreshold::Float64`: Surplus threshold for grid refinement.
- `tol_ti::Float64`: Convergence tolerance for time iteration.
- `polUpdateWeight::Float64`: Weight factor for updating the policy function during iteration.
- `maxIterEarlyStopping::Int`: Number of iterations before early stopping if convergence slows.
- `drawsnbr::Int`: Number of Monte Carlo draws for error estimation.
- `typeRefinement::String`: Type of refinement strategy applied to the sparse grid.
- `initialPolGuess::UserPolicyGuess`: Initial policy function guess.

# Returns
- `ddsg::DDSG`: The DDSG sparse grid approximation of the policy function.
- `sgws::SparsegridsWs`: The workspace storing model variables, parameters, and solver configurations.

# Process
1. **Model Extraction**  
   - Extracts **state-space representation**, endogenous/exogenous variables, and constraints.
   - Computes bounds (`lb`, `ub`) and **Mixed Complementarity Problem (MCP) conditions**.

2. **Sparse Grid Initialization**  
   - Constructs an **initial sparse grid** using model-defined domain boundaries.
   - Uses **monomial quadrature** to approximate integrals in stochastic models.

3. **First-Order Approximation & Policy Function Initialization**  
   - Computes an initial policy function guess using **first-order perturbation**.
   - Constructs `DDSG` grid representation and initializes the approximation.

4. **Time Iteration Process**  
   - Calls `ddsg_time_iteration!`, iteratively refining the policy function until convergence.
   - Performs **adaptive refinement** based on the surplus error threshold.

5. **Results Storage**  
   - Stores **iteration count, grid depth, solver tolerance, and refinement parameters** 
     in `context.results.model_results[1].sparsegrids`.

# Usage Context
This function is used for **solving high-dimensional dynamic models** where standard **sparse grids** 
are computationally expensive. **DDSG decomposes interactions between dimensions**, significantly reducing 
the number of required interpolation nodes while maintaining accuracy.
"""
function DDSGapproximation(;context::Context=context,
    k_max = 1,
    ftol = 1e-5,
    gridDepth = 2,
    gridOrder = 1,
    gridRule = "localp",
    maxiter = 300,
    maxRef = 0,
    maxRefLevel = gridDepth + maxRef,
    mcp = false,
    method = NewtonRaphson(),
    savefreq = 10,
    scaleCorrInclude = [],
    scaleCorrExclude = [],
    show_trace = false,
    solver = mcp ? NLsolver : NonlinearSolver,
    surplThreshold= 1e-3,
    tol_ti = 1e-4,
    polUpdateWeight = 0.5,
    maxIterEarlyStopping = 0,
    drawsnbr = 10000,
    typeRefinement = "classic",
    initialPolGuess::UserPolicyGuess = UserPolicyGuess(),
)
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
    sgoptions = SGOptions(mcp, method, solver, ftol, show_trace)

    # Initialize sparse grid workspace
    sgws_ = initialize_sparsegrid_workspace(monomial, dynamic_state_variables, system_variables,
                                            bmcps, ids, lb, ub, backward_block, forward_block,
                                            preamble_block, system_block, n_nodes, gridDim,
                                            gridOut, endogenous_nbr, sgmodel, sgoptions)

    # Handle multithreading
    if Threads.nthreads() > 1
        BLAS.set_num_threads(1)
        sgws = [deepcopy(sgws_) for _ in 1:Threads.nthreads()]
    else
        sgws = [sgws_]
    end

    # Guess for the policy function 
    M, N, _ = block_first_order_approximation(context, sgws[1])
    polGuess = X -> guess_policy(context, size(X,2), X, sgws[1], M, N, initialPolGuess)

    # Initialization
    ddsg = DDSG(gridDim, gridOut, gridDepth, maxRefLevel, k_max; order=gridOrder, rule=gridRule, domain=gridDomain)
    DDSG_init!(ddsg)
    DDSG_build!(ddsg, polGuess, ddsg.centroid; refinement_type=typeRefinement, refinement_tol=surplThreshold, scale_corr_vec=scaleCorr)

    # Time-iteration
    ddsg, average_time, iter = ddsg_time_iteration!(
        ddsg, polGuess, sgws, scaleCorr, surplThreshold, typeRefinement,
        maxiter, maxIterEarlyStopping, tol_ti, polUpdateWeight, savefreq
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

# function plot_policy_function(plot_variables, vstate, state, grid, state_variables, system_variables; N=50, context=context)
#     sv_buffer = repeat(state, 1, N)
#     k = findfirst(state_variables .== context.symboltable[vstate].orderintype)
#     kp = [findfirst(context.symboltable[v].orderintype .== system_variables) for v in plot_variables]
#     gridDomain = getDomainTransform(grid)
#     x = range(gridDomain[k,1], gridDomain[k, 2], N)
#     sv_buffer[k, :] .= x
#     Y = evaluateBatch(grid, sv_buffer)
#     return x, Y
# end
