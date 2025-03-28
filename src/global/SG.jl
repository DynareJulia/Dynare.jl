"""
    SGOptions

Struct to hold the options for the sparse-grid approximation

# Arguments
- `dimRef = -1`: Outputs that are considered in the refinement process (-1 implies that all outputs are considered)
- `ftol = 1e-5`: Convergence tolerance for the nonlinear solver,
- `iterRefStart = 25`: Iteration at which the refinement starts,
- `gridDepth = 2`: Initial sparse grid depth,
- `gridOrder = 1`: Polynomial order for interpolation,
- `gridRule = "localp": Type of base functions`,
- `maxiter = 300`: Maximum iterations for time iteration,
- `maxRef = 0`: Maximum number of refinements
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
- `initialPolGuess::UserPolicyGuess = UserPolicyGuess()`
"""
Base.@kwdef struct SGOptions
    dimRef::Int = -1
    ftol::Float64 = 1e-5
    iterRefStart::Int = 25
    gridDepth::Int = 2
    gridOrder::Int = 1
    gridRule::String = "localp"
    maxiter::Int = 300
    maxRef::Int = 0
    mcp::Bool = false
    method::NonlinearSolve.GeneralizedFirstOrderAlgorithm = NewtonRaphson()
    savefreq::Int = 10
    scaleCorrInclude::Vector{String} = Vector{String}()
    scaleCorrExclude::Vector{String} = Vector{String}()
    show_trace::Bool = false
    solver::Dynare.SGSolver = NonlinearSolver
    surplThreshold::Float64 = 0.
    tol_ti::Float64 = 1e-4
    polUpdateWeight::Float64 = 0.5
    maxIterEarlyStopping::Int = 0
    drawsnbr::Int = 10000
    typeRefinement::String = "classic"
    initialPolGuess::UserPolicyGuess = UserPolicyGuess()
end

"""
    initialize_sparse_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain)

Creates and initializes a sparse grid with the specified parameters.

# Arguments
- `gridDim` : Number of state variables (grid dimensions).
- `gridOut` : Number of policy function outputs.
- `gridDepth` : Initial sparse grid depth.
- `gridOrder` : Order of polynomial interpolation.
- `gridRule` : Rule for sparse grid basis functions.
- `gridDomain` : Boundaries of the state space.

# Returns
- `grid` : Initialized sparse grid.
- `aPoints` : Initial grid points that require function evaluation.
- `aNum` : Number of grid points needing function values.
"""
function initialize_sparse_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain)
    # Generate sparse grid
    grid = Tasmanian.TasmanianSG(gridDim, gridOut, gridDepth)
    makeLocalPolynomialGrid!(grid, order=gridOrder, rule=gridRule)

    # Apply domain transformation
    setDomainTransform!(grid, gridDomain)

    # Extract grid points requiring evaluation
    aPoints = getPoints(grid)
    aNum = getNumPoints(grid)

    return grid, aPoints, aNum
end

"""
    initialize_policy_approximation(context, sgws, aNum, gridDim, aPoints)

Computes the first-order perturbation approximation and initializes the policy function.

# Arguments
- `context::Context` : Dynare context storing the model structure.
- `sgws::SparsegridsWs` : Workspace struct storing grid and solver data.
- `aNum::Int` : Number of grid points requiring function values.
- `aPoints` : Sparse grid evaluation points.
- `initialPolGuess::UserPolicyGuess` : User-provided policy guess

# Returns
- `polGuess::Matrix{Float64}` : Initial policy function approximation.
"""
function initialize_policy_approximation(context, sgws, aNum, aPoints, initialPolGuess)
    # First-order approximation
    M, N, _ = block_first_order_approximation(context, sgws)

    # Compute initial policy guess
    polGuess = guess_policy(context, aNum, aPoints, sgws, M, N, initialPolGuess)

    return polGuess
end

"""
    sparse_grid_time_iteration(grid, grid0, polGuess, sgws, scaleCorr, surplThreshold, dimRef, 
                               typeRefinement, maxiter, maxRef, iterRefStart, tol_ti, savefreq, timings, X_sample)

Performs the time iteration loop to refine the sparse grid and compute the policy function until convergence.

# Arguments
- `grid` : Sparse grid structure.
- `grid0` : Initial sparse grid.
- `polGuess::Matrix{Float64}` : Initial policy function guess.
- `sgws::Vector{SparsegridsWs}` : Workspace structures for each thread.
- `scaleCorr::Vector{Float64}` : Scale correction vector.
- `surplThreshold::Float64` : Surplus threshold for adaptive refinement.
- `dimRef::Int` : Refinement dimension indicator.
- `typeRefinement::String` : Type of refinement (e.g., "classic").
- `maxiter::Int` : Maximum number of iterations.
- `maxRef::Int` : Maximum refinement level.
- `iterRefStart::Int` : Iteration at which refinement starts.
- `maxIterEarlyStopping::Int` : Number of iterations after which TI stops after TI convergence measure starts increasing
- `tol_ti::Float64` : Convergence criterion.
- `polUpdateWeight::Float64` : Weight of the current-step policy function when computing the updated policy function. The weight of the previous-step policy function is `1-polUpdateWeight`
- `savefreq::Int` : Frequency for saving the grid.
- `timings::Vector{Millisecond}` : Storage for iteration timings.

# Returns
- `grid` : Updated sparse grid after refinement.
- `polGuess` : Updated policy function.
- `average_time::Float64` : Average iteration time.
"""
function sparse_grid_time_iteration!(
    grid, grid0, polGuess, sgws, scaleCorr, surplThreshold, dimRef,
    typeRefinement, maxiter, maxRef, iterRefStart, maxIterEarlyStopping,
    tol_ti, polUpdateWeight, savefreq
)
    # Initialize timing variables
    total_time = 0
    context.timings["sparsegrids"] = Vector{Millisecond}()
    timings = context.timings["sparsegrids"]
    iter = 1
    iterEarlyStopping = 0
    previousMetric = Inf
    scaleCorrMat = repeat(scaleCorr, 1, Tasmanian.getNumLoaded(grid))

    while iter <= maxiter && iterEarlyStopping <= maxIterEarlyStopping
        tin = now()
        polGuess1 = copy(polGuess)
        newgrid = copyGrid(grid0)

        # Refinement loop
        ilev = 0
        while (getNumNeeded(newgrid) > 0)
            if (ilev <= maxRef)
                # Reset Jacobian for each workspace
                map(s -> fill!(s.J, 0.0), sgws)

                # Compute new policy function values
                ti_step!(newgrid, grid, polGuess1, sgws)

                # Perform adaptive refinement
                if iter >= iterRefStart && ilev < maxRef
                    polGuess1 = refine!(newgrid, polGuess1, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)
                    # Track refinement level
                    ilev += 1
                end
            else
                # If refinement is done, load function values
                loadNeededPoints!(newgrid, polGuess1)
            end
        end

        # Compute error and update policy function
        metric, polGuess, grid = policy_update!(grid, newgrid, polGuess, polGuess1, length(sgws[1].system_variables), polUpdateWeight)
        iteration_walltime = now() - tin

        # Verify that the metric decreased w.r.t the previous iteration and
        # adjust early-stopping indicators accordingly
        if previousMetric > metric
            iterEarlyStopping = 0
        else
            iterEarlyStopping += 1
        end
        previousMetric = metric

        # Store iteration timing
        iter > 1 && (total_time += iteration_walltime.value)
        push!(timings, iteration_walltime)

        # Logging iteration details
        println("Iteration: $iter, Grid pts: $(getNumPoints(grid)), Refinement level: $(ilev), Metric: $metric, Computing time: $(iteration_walltime)")

        # Save grid state periodically
        if (mod(iter + 1, savefreq) == 0)
            # save_grid(grid0, iter)
        end

        # Convergence check
        if metric < tol_ti
            break
        end

        iter += 1
    end

    # Compute average iteration time
    average_time = total_time / max(iter - 1, 1)  # Avoid division by zero
    println("Last grid points: $(getNumPoints(grid))")
    println("Average iteration time (except first one): $average_time")

    return grid, polGuess, average_time, iter
end

"""
    ti_step!(newgrid, oldgrid, polGuess, sgws)

Updates grid points by solving the nonlinear problem at newly needed grid locations.

# Arguments
- `newgrid` : The sparse grid that needs new function values.
- `oldgrid` : The previous iteration's grid.
- `polGuess::Matrix{Float64}` : Matrix to store the computed function values.
- `sgws::Vector{SparsegridsWs}` : The workspace containing all model and solver settings.

# Effect
- Modifies `newgrid` by loading computed function values into it.
"""
function ti_step!(newgrid, oldgrid, polGuess, sgws)
    # Extract necessary data from sgws
    ws = sgws[1]
    @unpack system_variables, lb, ub, sgsolveroptions, J, fx, sgmodel = ws
    @unpack mcp, method, solver, ftol, show_trace = sgsolveroptions
    @unpack exogenous, parameters = sgmodel

    # Ensure solver is set
    solver = isnothing(solver) ? (mcp ? NLsolver : NonlinearSolver) : solver

    # Retrieve points requiring function values
    states = getNeededPoints(newgrid)
    nstates = size(states, 2)

    # Initialize exogenous variables (assuming no autocorrelated shocks)
    fill!(exogenous, 0.0)

    # Solve the nonlinear problem
    SG_NLsolve!(polGuess, lb, ub, fx, J, states, oldgrid, sgws, solver, method, ftol, show_trace)

    # Load new function values into the sparse grid
    @views loadNeededPoints!(newgrid, polGuess[:, 1:nstates])
end

"""
    refine!(grid, polguess, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)

Refines the sparse grid based on surplus coefficients and updates function values.

# Arguments
- `grid` : The sparse grid to refine.
- `polguess::Matrix{Float64}` : Current policy function values.
- `scaleCorr::Vector{Float64}` : Scale correction vector.
- `scaleCorrMat::Matrix{Float64}` : Scale correction matrix.
- `surplThreshold::Float64` : Surplus threshold for refinement.
- `dimRef::Int` : Output dimensions considered in refinement.
- `typeRefinement::String` : Type of refinement to apply.

# Effect
- Updates `grid` and `polguess` in place.

# Returns
- `polguess::Matrix{Float64}` : Updated policy function values.
"""
function refine!(grid, polguess, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)
    # Ensure `scaleCorrMat` is correctly allocated
    num_loaded = getNumLoaded(grid)
    if size(scaleCorrMat, 2) < num_loaded
        scaleCorrMat = repeat(scaleCorr, 1, num_loaded)
    end

    # Perform grid refinement based on surplus coefficients
    setSurplusRefinement!(grid, surplThreshold, output=dimRef, refinement_type=typeRefinement, scale_correction=scaleCorrMat')

    # Check if new points were added
    num_needed = getNumNeeded(grid)
    if num_needed > 0
        # Get new required points
        nwpts = getNeededPoints(grid)
        if size(polguess,1) == num_needed
            # Update policy function values in place
            interpolate!(polguess, grid, nwpts)
        else
            polguess = interpolate(grid, nwpts)
        end
    end

    return polguess
end

"""
    policy_update(gridOld, gridNew, polGuessOld, polGuessNew, gridOut)

Checks convergence and updates policy function values.

# Arguments
- `gridOld` : Previous sparse grid.
- `gridNew` : Updated sparse grid.
- `polGuessOld::Matrix{Float64}` : Previous policy function values.
- `polGuessNew::Matrix{Float64}` : Updated policy function values.
- `gridOut::Int` : Number of policy variables.
- `polUpdateWeight::Float64` : Weight of the current-step policy function when computing the updated policy function. The weight of the previous-step policy function is `1-polUpdateWeight`

# Effect
- Updates `polGuessNew` in place.

# Returns
- `metric::Float64` : Convergence metric (min of L2 and Sup-norm).
- `polGuessNew::Matrix{Float64}` : Updated policy function values.
- `gridNew` : Copied grid.
"""
function policy_update!(gridOld, gridNew, polGuessOld, polGuessNew, gridOut, polUpdateWeight)

    # Get the points and the number of points from grid1
    aPoints = interpolationPoints(gridNew)
    aNumTot = size(aPoints,2)

    # Evaluate the grid points on both grid structures
    if size(polGuessNew, 2) == aNumTot
        interpolate!(polGuessNew, gridNew, Matrix(aPoints))
        
    else    
        polGuessNew = interpolate(gridNew, Matrix(aPoints))
    end 
    # Load the polGuessOld points on gridOld
    loadNeededPoints!(gridOld, polGuessOld)
    if size(polGuessOld, 2) == aNumTot
        interpolate!(polGuessOld, gridOld, Matrix(aPoints))
    else    
        polGuessOld = interpolate(gridOld, Matrix(aPoints))
    end 

    # 1) Compute the Sup-Norm
    metricAux = zeros(gridOut)

    for imet in 1:gridOut
        @views metricAux[imet] = maximum(abs.(polGuessOld[imet, :] - polGuessNew[imet, :]))
    end
    metricSup = maximum(metricAux)

    # 2) Compute the L2-Norm
    metricL2 = 0.0

    for imetL2 in 1:gridOut
        @views metricL2 += sum(abs.(polGuessOld[imetL2, :] - polGuessNew[imetL2, :]).^2)
    end
    metricL2 = (metricL2/(aNumTot*gridOut))^0.5
    metric = min(metricL2, metricSup)

    # Now updating pol_guess and grid
    @views @. polGuessNew[1:gridOut, :] = polUpdateWeight * polGuessNew[1:gridOut, :] + (1 - polUpdateWeight) * polGuessOld[1:gridOut, :]

    # Load the corresponding points on the grid
    loadNeededPoints!(gridNew, polGuessNew)

    #return metric, polGuessNew, copyGrid(gridNew)
    return metric, polGuessNew, makeCopy(gridNew)
end

"""
    SG_NLsolve!(X, lb, ub, fx, J, states, parameters, oldgrid, sgws, solver, method, ftol, show_trace)

Solves the nonlinear system using the specified solver.

# Arguments
- `polGuess::Matrix{Float64}` : Solution matrix to update.
- `lb::Vector{Float64}` : Lower bounds.
- `ub::Vector{Float64}` : Upper bounds.
- `fx::Vector{Float64}` : Function evaluation vector.
- `J::Matrix{Float64}` : Jacobian matrix.
- `states::Matrix{Float64}` : State variables for solving.
- `parameters` : Model parameters.
- `oldgrid` : Previous iteration’s grid.
- `sgws::Vector{SparsegridsWs}` : Workspace structure for each thread.
- `solver` : The selected nonlinear solver.
- `method` : Solver method (e.g., Newton-Raphson).
- `ftol::Float64` : Function tolerance for solver convergence.
- `show_trace::Bool` : Enables solver debugging output.

# Effect
- Updates `polGuess` with the computed solution values.

# Returns
- `polGuess::Matrix{Float64}` : Updated solution matrix.
"""
function SG_NLsolve!(polGuess, lb, ub, fx, J, states, oldgrid, sgws, solver, method, ftol, show_trace)
    ns = size(states, 2)  # Number of state variable instances

    # Get BLAS to work on a single thread to avoid conflicts with the non-linear
    # solvers
    previous_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    # Partition workload among available threads
    chunks = Iterators.partition(1:ns, (ns ÷ Threads.nthreads()) + 1) |> collect
    tasks = []
    for (i, chunk) in enumerate(chunks)
        push!(tasks, Threads.@spawn begin
            if solver == PATHSolver
                # PATHsolver is not thread-safe
                Threads.nthreads() == 1 || error("PATHSolver is not thread-safe! Run with `JULIA_NUM_THREADS=1`.")
                PATHsolver_solve!(polGuess, lb, ub, states, fx, J, oldgrid, sgws[1], show_trace, chunk)
            elseif solver == NonlinearSolver
                NonlinearSolver_solve!(polGuess, states, J, oldgrid, method, sgws[i], ftol, show_trace, chunk)
            elseif solver == NLsolver
                NLsolve_solve!(polGuess, lb, ub, fx, J, states, oldgrid, sgws[i], ftol, show_trace, chunk)
            else
                error("Unknown non-linear solver! The available options are PATHSolver, NonLinearSolver and NLsolve")
            end
        end)
    end

    fetch.(tasks)

    # Get BLAS to use the number of threads before SG_NLsolve! call
    BLAS.set_num_threads(previous_blas_threads)
    return polGuess
end
