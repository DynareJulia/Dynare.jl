using Combinatorics
using AxisArrays: AxisArray
using DualNumbers
using Distributions
using FiniteDiff
using LinearAlgebra
import NLsolve
using NonlinearSolve
using Roots
using Tasmanian

"""
    interpolate(grid, X)

Interpolates the function values at the given points using the specified sparse grid.

# Arguments
- `grid`: The sparse grid structure used for interpolation.
- `X`: The input points for which interpolation is performed. Can be a single point or a batch of points.

# Returns
- Interpolated function values at `X`.

# Behavior
- If the grid has the property `is_ddsg`, the function delegates to `DDSG_evaluate(grid, X)`.
- If `X` has the same number of dimensions as the grid, it calls `evaluate(grid, X)`.
- Otherwise, it processes `X` as a batch using `evaluateBatch(grid, X)`.
"""
function interpolate(grid,X)
    if hasproperty(grid, :is_ddsg) 
        return DDSG_evaluate(grid,X)
    else
        if length(X)==getNumDimensions(grid)
            return evaluate(grid,X)
        else
            return evaluateBatch(grid,X) 
        end
    end
end

"""
    interpolate!(Y, grid, X)

Computes the interpolated values at the given points `X` and stores them in `Y`.

# Arguments
- `Y`: Preallocated storage for the interpolated results.
- `grid`: The sparse grid structure used for interpolation.
- `X`: The input points for which interpolation is performed.

# Effect
- The interpolated values are stored in `Y`.

# Behavior
- If the grid has the property `is_ddsg`, it calls `DDSG_evaluate!(Y, grid, X)`.
- Otherwise, it uses `evaluateBatch!(Y, grid, X)`.
"""
function interpolate!(Y,grid,X)
    if hasproperty(grid, :is_ddsg) 
        DDSG_evaluate!(Y,grid,X)
    else
        evaluateBatch!(Y,grid,X) 
    end
end

"""
    derivate!(J, grid, x)

Computes the derivatives of the interpolated function at the given point `x` and stores them in `J`.

# Arguments
- `J`: Preallocated storage for the computed derivatives.
- `grid`: The sparse grid structure used for differentiation.
- `x`: The input point at which differentiation is performed.

# Effect
- The computed derivatives are stored in `J`.

# Behavior
- If the grid has the property `is_ddsg`, it delegates to `DDSG_differentiate!(grid, J=J, x=x)`.
- Otherwise, it calls `differentiate!(J, grid, x)`.
"""
function derivate!(J,grid,x)
    if hasproperty(grid, :is_ddsg) 
        DDSG_differentiate!(J,grid,x)
    else
        differentiate!(J,grid,x) 
    end
end

"""
    makeCopy(grid)

Creates a copy of the given sparse grid.

# Arguments
- `grid`: The sparse grid structure to be copied.

# Returns
- A new grid structure that is a copy of `grid`.

# Behavior
- If the grid has the property `is_ddsg`, it returns `DDSG(grid)`.
- Otherwise, it calls `copyGrid(grid)`.
"""
function makeCopy(grid)
    if hasproperty(grid, :is_ddsg) 
        return DDSG(grid)
    else
        return copyGrid(grid) 
    end
end

"""
    interpolationPoints(grid)

Retrieves the interpolation nodes of the given sparse grid.

# Arguments
- `grid`: The sparse grid structure.

# Returns
- A set of interpolation nodes.

# Behavior
- If the grid has the property `is_ddsg`, it returns `DDSG_nodes(grid)`.
- Otherwise, it calls `getPoints(grid)`.
"""
function interpolationPoints(grid)
    if hasproperty(grid, :is_ddsg) 
        return DDSG_nodes(grid)
    else
        return getPoints(grid) 
    end
end

"""
    countPoints(grid)

Computes the total number of interpolation points in the given sparse grid.

# Arguments
- `grid`: The sparse grid structure.

# Returns
- The number of interpolation points.

# Behavior
- If the grid has the property `is_ddsg`, it sums up `grid.grid_points`.
- Otherwise, it returns `getPoints(grid)`.
"""
function countPoints(grid)
    if hasproperty(grid, :is_ddsg) 
        return sum(grid.grid_points)
    else
        return getPoints(grid) 
    end
end

"""
    initialize_sparsegrid_context(context)

Extracts model information, system variables, and parameter values.
"""
function initialize_sparsegrid_context(context)
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    model = context.models[1]

    endogenous_nbr = model.endogenous_nbr
    exogenous_nbr = model.exogenous_nbr
    params = context.work.params
    steadystate = context.results.model_results[1].trends.endogenous_steady_state

    # Identify system variables and equation blocks
    dynamic_state_variables, predetermined_variables, system_variables,
    backward_block, forward_block, preamble_block, system_block = make_block_functions(context)

    return model, endogenous_nbr, exogenous_nbr, params, steadystate,
           dynamic_state_variables, endogenous_variables, predetermined_variables, system_variables,
           backward_block, forward_block, preamble_block, system_block
end

"""
    block_mcp!(lb, ub, bmcps, context, block_eqs, block_vars, limits)

Computes Mixed Complementarity Problem (MCP) constraints, setting up lower and upper bounds for policy variables.

# Arguments
- `lb::Vector{Float64}` : Lower bounds (modified in place).
- `ub::Vector{Float64}` : Upper bounds (modified in place).
- `bmcps::Vector{Vector{Int}}` : Complementarity conditions (modified in place).
- `context` : Dynare context storing model definitions.
- `block_eqs::Vector{Symbol}` : System equations.
- `block_vars::Vector{Symbol}` : Variables corresponding to equations.
- `limits` : Constraint limits.

# Effect
- Updates `lb`, `ub`, and `bmcps` in place.

# Returns
- `nothing`
"""
function block_mcp!(lb, ub, bmcps, context, block_eqs, block_vars, limits)
    # Create a local copy of equations to avoid modifying original
    block_eqs_ = copy(block_eqs)

    # Ensure lb and ub are properly sized and initialized
    n = length(block_eqs_)
    resize!(lb, n)
    resize!(ub, n)
    fill!(lb, -Inf)
    fill!(ub, Inf)

    # Extract MCP constraints from the model
    for (eqn, var, op, expr) in context.models[1].mcps
        bvar = findfirst(==(var), block_vars)
        beqn = findfirst(==(eqn), block_eqs_)

        # Swap equation and variable positions
        block_eqs_[[bvar, beqn]] .= block_eqs_[[beqn, bvar]]
        push!(bmcps, [beqn, bvar])

        # Evaluate boundary condition
        boundary = Dynare.dynare_parse_eval(String(expr), context)

        # Apply bound based on operator
        if startswith(op, "<")
            ub[bvar] = boundary
        elseif startswith(op, ">")
            lb[bvar] = boundary
        else
            error("Invalid MCP operator: $op. Must be one of '<', '<=', '>', '>='")
        end
    end

    return nothing
end

"""
    initialize_mcp_constraints(context, forward_block, backward_block, system_variables)

Computes MCP constraints, extracting lower/upper bounds and complementarity constraints.
"""
function initialize_mcp_constraints(context, forward_block, backward_block, system_variables)
    lb = Float64[]
    ub = Float64[]
    bmcps = Vector{Vector{Int}}()

    block_eqs = union(forward_block.equations, backward_block.equations)
    block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables, context.work.limits)

    return lb, ub, bmcps
end

"""
    initialize_sgmodel(steadystate, endogenous_nbr, exogenous_nbr, params)

Creates an SGModel with repeated steady state values and initialized zero vectors.
"""
function initialize_sgmodel(steadystate, endogenous_nbr, exogenous_nbr, params)
    return SGModel(
        repeat(steadystate, 3),  # State variable initialization
        zeros(2 * endogenous_nbr),  # Placeholder for policy variables
        zeros(exogenous_nbr),  # Placeholder for exogenous shocks
        endogenous_nbr,
        exogenous_nbr,
        params,
        steadystate,
        []
    )
end

"""
    initialize_sg_indices(model, dynamic_state_variables, system_variables, predetermined_variables)

Creates an SGIndices structure and determines forward system variables.
"""
function initialize_sg_indices(model, dynamic_state_variables, system_variables, predetermined_variables)
    ids = SGIndices(
        model.endogenous_nbr,
        model.i_fwrd_b,
        dynamic_state_variables,
        system_variables
    )

    # Determine forward system variables (excluding predetermined)
    forward_system_variables = setdiff(model.i_fwrd_b, predetermined_variables .- model.endogenous_nbr)

    return ids, forward_system_variables
end

"""
    get_grid_dimensions(dynamic_state_variables, system_variables)

Computes the grid dimensions based on the number of state and policy variables.
"""
function get_grid_dimensions(dynamic_state_variables, system_variables)
    gridDim = length(dynamic_state_variables)
    gridOut = length(system_variables)
    return gridDim, gridOut
end

"""
    get_grid_domain(limits, ids, context, gridDim)

Computes the domain limits for the sparse grid based on the model's constraints.
"""
function get_grid_domain(limits, ids, context, gridDim)
    gridDomain = zeros(gridDim, 2)
    for l in limits
        i = findfirst(x -> x == context.symboltable[string(l[1])].orderintype, ids.state_variables)
        if !isnothing(i)
            gridDomain[i, 1] = l[2].min
            gridDomain[i, 2] = l[2].max
        end
    end
    return gridDomain
end

"""
    initialize_sparsegrid_workspace(monomial, dynamic_state_variables, system_variables, 
                                    bmcps, ids, lb, ub, backward_block, forward_block, 
                                    preamble_block, system_block, grid, nodesnbr, 
                                    gridDim, gridOut, dynamic_state_variables, 
                                    endogenous_nbr, sgmodel, sgoptions)

Initializes all matrices, vectors, and workspace structures needed for sparse grid approximation.

# Arguments
- `monomial::MonomialPowerIntegration` : Monomial quadrature rule.
- `dynamic_state_variables::Vector{Int}` : Indices of dynamic state variables.
- `system_variables::Vector{Int}` : Indices of all system variables.
- `bmcps::Vector{Vector{Int}}` : Mixed Complementarity Problem (MCP) constraints.
- `ids::SGIndices` : Stores variable indexing.
- `lb::Vector{Float64}` : Lower bounds for policy variables.
- `ub::Vector{Float64}` : Upper bounds for policy variables.
- `backward_block::BackwardBlock` : Backward-looking equations.
- `forward_block::ForwardBlock` : Forward-looking equations.
- `preamble_block::AssignmentBlock` : Initialization equations.
- `system_block::SimultaneousBlock` : Nonlinear equation system.
- `grid` : Sparse grid.
- `nodesnbr::Int` : Number of monomial integration nodes.
- `gridDim::Int` : Number of state variables.
- `gridOut::Int` : Number of policy function outputs.
- `endogenous_nbr::Int` : Number of endogenous variables.
- `sgmodel::SGModel` : Sparse grid model.
- `sgoptions::SGOptions` : Solver configuration.

# Returns
- `SparsegridsWs` instance containing all required structures.

"""
function initialize_sparsegrid_workspace(monomial, dynamic_state_variables, system_variables,
                                         bmcps, ids, lb, ub, backward_block,
                                         forward_block, preamble_block,
                                         system_block,  n_nodes, gridDim,
                                         gridOut, endogenous_nbr,
                                         sgmodel, sgoptions)

    # Precompute lengths
    n_system = length(system_variables)
    n_forward = length(forward_block.equations)
    n_forward_in_system = length(ids.forward_in_system)
    n_state_in_system = length(ids.state_in_system)
    n_backward = length(backward_block.equations)

    # Jacobian Matrices
    backward_variable_jacobian = Matrix{Float64}(undef, n_backward, n_system)
    forward_jacobian = Matrix{Float64}(undef, n_forward, n_system)
    forward_variable_jacobian = Matrix{Float64}(undef, n_forward, n_forward_in_system)
    J = zeros(n_system, n_system)
    future_policy_jacobian = zeros(n_state_in_system, n_forward_in_system)
    policy_jacobian = zeros(gridDim, n_system)

    # Residuals and Function Evaluations
    residuals = Vector{Float64}(undef, n_system)
    forward_residuals = Matrix{Float64}(undef, n_forward, n_nodes)
    fx = zeros(n_system)

    # Function Evaluation Points
    evalPt = Vector{Float64}(undef, gridDim)
    forward_points = Matrix{Float64}(undef, gridDim, n_nodes)

    # Storage for Endogenous and Policy Functions
    dyn_endogenous_vector = [zeros(3 * endogenous_nbr) for _ in 1:n_nodes]
    policyguess = zeros(gridOut, n_nodes)

    # Temporary Storage Variables
    M1 = Matrix{Float64}(undef, n_forward, length(ids.state_in_system))
    tmp_state_variables = Vector{Float64}(undef, gridDim)

    # Create Sev struct
    sev = Sev(zeros(3 * endogenous_nbr), zeros(gridDim))

    # Create SparsegridsWs workspace struct
    return SparsegridsWs(monomial, dynamic_state_variables, system_variables, 
                         bmcps, ids, lb, ub, backward_block, backward_variable_jacobian, 
                         forward_block, preamble_block, system_block,
                         residuals, forward_jacobian, forward_points, 
                         forward_residuals, forward_variable_jacobian, J, fx,
                         policy_jacobian, evalPt, 
                         future_policy_jacobian, 
                         dyn_endogenous_vector, M1, policyguess, tmp_state_variables,
                         sev, sgmodel, sgoptions)
end

"""
    make_scaleCorr(scaleCorrInclude, scaleCorrExclude, policy_variables) -> Vector{Float64}

Constructs a scale correction vector ensuring proper weighting in sparse grid refinement.

# Arguments
- `scaleCorrInclude::Vector{String}`: Variables to apply scale correction.
- `scaleCorrExclude::Vector{String}`: Variables to exclude from scale correction.
- `policy_variables::Vector{String}`: List of all policy variables.

# Returns
- `scaleCorr::Vector{Float64}`: A vector where:
  - `1.0` → Variables **included** in scale correction.
  - `0.0` → Variables **excluded** from scale correction.
  - **Defaults to `1.0` for all variables** if no exclusions or inclusions are specified.

# Behavior
- If `scaleCorrInclude` is **not empty**, only variables that **partially match** (via `occursin`) are assigned `1.0`, and others get `0.0`.
- If `scaleCorrExclude` is **not empty**, variables that **partially match** are assigned `0.0`, and others get `1.0`.
- If both are empty, **default is `1.0` for all policy variables**.
"""
function make_scaleCorr(scaleCorrInclude, scaleCorrExclude, policy_variables)
    @assert isempty(scaleCorrInclude) || isempty(scaleCorrExclude) "scaleCorrInclude and scaleCorrExclude are mutually exclusive"

    n = length(policy_variables)
    scaleCorr = ones(Float64, n)  # Default: all variables get scale correction (1.0)

    if !isempty(scaleCorrInclude)
        scaleCorr .= [any(occursin.(v, scaleCorrInclude)) ? 1.0 : 0.0 for v in policy_variables]
    elseif !isempty(scaleCorrExclude)
        scaleCorr .= [any(occursin.(v, scaleCorrExclude)) ? 0.0 : 1.0 for v in policy_variables]
    end

    return scaleCorr
end

"""
    set_boundaries!(y, lb, ub)

Ensures that the values in `y` remain within specified lower and upper bounds.

# Arguments
- `y::Matrix{Float64}` : Policy matrix to be constrained.
- `lb::Vector{Float64}` : Lower bounds for each policy variable.
- `ub::Vector{Float64}` : Upper bounds for each policy variable.

# Effect
- Modifies `y` **in place** to ensure all values satisfy `lb <= y <= ub`.

"""
function set_boundaries!(y::Matrix{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
    @inbounds for j in axes(y, 2)  # Iterate over columns (policy evaluations)
        for i in axes(y, 1)  # Iterate over rows (policy variables)
            y[i, j] = clamp(y[i, j], lb[i] + eps(), ub[i] - eps())
        end
    end
end

"""
    guess_policy(context, aNum, aPoints, sgws, M, N)

Computes an initial guess for the policy function using first-order perturbation.

# Arguments
- `context::Context` : Dynare model context.
- `aNum::Int` : Number of grid points requiring function values.
- `aPoints::Matrix{Float64}` : Sparse grid evaluation points.
- `sgws::SparsegridsWs` : Solver workspace.
- `M::Matrix{Float64}` : First-order coefficient matrix for predetermined variables.
- `N::Matrix{Float64}` : First-order coefficient matrix for state variables.
- `initialPolGuess::UserPolicyGuess` : Dictionnary containing proposals for the policy function, e.g. {"c" => }

# Returns
- `polGuess::Matrix{Float64}` : Initial policy function guess.
"""
function guess_policy(context, aNum, aPoints, sgws, M, N, initialPolGuess)
    # Unpack relevant components from the solver workspace
    @unpack ids, dynamic_state_variables, system_variables, preamble_block, 
            lb, ub, sgmodel = sgws
    state_variables = preamble_block.variables

    # Retrieve model details and steady-state values
    model = context.models[1]
    steadystate = context.results.model_results[1].trends.endogenous_steady_state

    # Extract relevant steady-state components
    ssy = steadystate[system_variables]  # Steady-state of system variables
    ssi = steadystate[filter(x -> x <= model.endogenous_nbr, dynamic_state_variables)]  # Lagged endogenous
    sst = steadystate[state_variables]  # State variables' steady-state

    # Initialize variables
    nsv = length(system_variables)
    polGuess = zeros(nsv, aNum)
    y0 = zeros(nsv)

    # Identify indices for state and lagged endogenous variables
    istate = findall(dynamic_state_variables .> sgmodel.endogenous_nbr)
    ilagendo = findall(dynamic_state_variables .≤ sgmodel.endogenous_nbr)
    iy = findall(in(dynamic_state_variables[ilagendo]), system_variables) 

    # Get the names of the policy variables
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    system_variable_names = endogenous_variables[system_variables]

    # Get the set of policy variables that have a corresponding user-provided
    # policy guess function. The others will be assigned a guess using
    # first-order perturbation
    user_guessed = findall(in(initialPolGuess.outputs), system_variable_names)
    dynare_guessed = findall(!in(initialPolGuess.outputs), system_variable_names)

    # Compute the initial policy function guess for the variables specified
    # in initialPolGuess
    if !isempty(user_guessed)
        # Get the names of the endogenous state variables and exogenous shocks
        # (In SGU's terminology). These are the inputs of the policy functions
        # in the order of `aPoints` rows.  `aPoints` is a matrix with `nsv` rows
        # and `aNum` points, where `nsv` is the number of policy variables, and
        # `aNum` is the number of points on the grid. The first rows contain the
        # coordinates of the endogenous lagged variables (state variables in
        # SGU's terminology), whose indices are in `iy`. The last rows contain
        # the coordinates of the exogenous processes (in SGU's terminology),
        # whose indices in endogenous_variables are `istate`.
        input_names = endogenous_variables[[iy;istate]]
        input_order = findfirst.(isequal.(input_names), Ref(initialPolGuess.inputs))
        user_guessed_order = findfirst.(isequal.(system_variable_names[user_guessed]), Ref(initialPolGuess.outputs))
        polGuess[user_guessed, :] = reduce(hcat, [initialPolGuess.polFun(col[input_order])[user_guessed_order] for col in eachcol(aPoints)])
    end

    # Compute the initial policy function guess for the other variables
    dynare_guess = zeros(nsv)
    @views for (i, c) in enumerate(eachcol(aPoints))
        x = c[istate] .- sst  # Deviation from steady-state for state variables
        y0[iy] .= c[ilagendo] .- ssi  # Deviation for lagged endogenous variables
        dynare_guess .= ssy .+ M * y0 .+ N * x
        polGuess[dynare_guessed, i] = dynare_guess[dynare_guessed]
    end

    # Ensure the policy guess stays within bounds
    set_boundaries!(polGuess, lb, ub)

    return polGuess
end


"""
    set_dyn_endogenous_vector!(dyn_endogenous_vector, policy, state, grid, sgws)

Updates `dyn_endogenous_vector` by computing state transitions for each integration node.

# Arguments
- `dyn_endogenous_vector::Vector{Vector{Float64}}` : Output vector to store updated endogenous variables.
- `policy::Vector{Float64}` : Policy function values.
- `state::Vector{Float64}` : Current state variables.
- `grid` : Sparse grid for evaluation.
- `sgws::SparsegridsWs` : Solver workspace.

# Effect
- Modifies `dyn_endogenous_vector` in place with computed values.

# Returns
- `nothing`
"""
function set_dyn_endogenous_vector!(dyn_endogenous_vector, policy, state, grid, sgws)
    # Unpack necessary variables
    @unpack monomial, dynamic_state_variables, system_variables, forward_block, 
            preamble_block, sgmodel, forward_points, forward_residuals, policyguess = sgws
    @unpack shift_dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate, tempterms = sgmodel

    nodes = monomial.nodes
    nodesnbr = length(nodes)

    for i in 1:nodesnbr
        # Local reference to prevent race conditions
        dyn_endogenous_vector_i = dyn_endogenous_vector[i]
        fill!(dyn_endogenous_vector_i, 0.0)

        # Copy state variables
        for (j, k) in enumerate(dynamic_state_variables)
            dyn_endogenous_vector_i[k] = state[j]
        end

        # Copy to local buffer for next-period computations
        copyto!(shift_dyn_endogenous, 1, dyn_endogenous_vector_i, endogenous_nbr + 1, 2 * endogenous_nbr)

        # Compute t+1 endogenous variables using preamble_block
        preamble_block.set_endogenous_variables!(tempterms, shift_dyn_endogenous, nodes[i], parameters, steadystate)

        # Copy computed next-period values back
        copyto!(dyn_endogenous_vector_i, endogenous_nbr + 1, shift_dyn_endogenous, 1, 2 * endogenous_nbr)

        # Update system variables with policy values
        for (j, sv) in enumerate(system_variables)
            dyn_endogenous_vector_i[sv + endogenous_nbr] = policy[j]
        end

        # Extract next-period state variables
        for (j, dsv) in enumerate(dynamic_state_variables)
            forward_points[j, i] = dyn_endogenous_vector_i[dsv + endogenous_nbr]
        end
    end

    # Determine relevant variables within the expectations operator
    interpolate!(policyguess, grid, forward_points)

    # Update `dyn_endogenous_vector` with computed policies
    for (i, k) in enumerate(system_variables), j in 1:nodesnbr
        dyn_endogenous_vector[j][2 * endogenous_nbr + k] = policyguess[i, j]
    end

    return nothing
end

"""
    expectation_equations!(expected_residuals, policy, state, grid, sgws)

Computes the expected residuals by integrating forward-looking equation errors.

# Arguments
- `expected_residuals::Vector{Float64}` : Output vector to store computed expected residuals.
- `policy::Vector{Float64}` : Policy function values.
- `state::Vector{Float64}` : Current state variables.
- `grid` : Sparse grid for evaluation.
- `sgws::SparsegridsWs` : Solver workspace.

# Effect
- Modifies `expected_residuals` in place.

# Returns
- `nothing`
"""
function expectation_equations!(expected_residuals, policy, state, grid, sgws)
    # Unpack necessary variables
    @unpack monomial, forward_residuals = sgws

    # Compute residuals for forward-looking equations
    forward_looking_equation_residuals!(forward_residuals, policy, state, grid, sgws)

    # Compute expected residuals using monomial weights
    mul!(expected_residuals, forward_residuals, monomial.weights)

    return nothing
end

"""
    copy_from_submatrix!(dest, src, rows, cols)

Copies selected elements from `src` into `dest` using **element-wise mapping**.

# Arguments
- `dest::Matrix{Float64}` : Destination matrix (modified in place).
- `src::Matrix{Float64}` : Source matrix to copy from.
- `rows::Vector{Int}` : Row indices for element-wise selection.
- `cols::Vector{Int}` : Column indices for element-wise selection.

# Effect
- Modifies `dest` in place with element-wise mapped values.

# Returns
- `dest::Matrix{Float64}` : Updated destination matrix.
"""
function copy_from_submatrix!(dest, src, rows, cols)
    @inbounds @views for j in eachindex(cols)
        k = cols[j]  # Get column index
        for i in eachindex(rows)
            dest[i, j] = src[rows[i], k]  # Map row `rows[i]` and column `cols[j]`
        end
    end
    return dest
end

"""
    forward_looking_equation_derivatives!(J, grid, dyn_endogenous, weight, sgws)

Computes the Jacobian matrix for forward-looking equations.

# Arguments
- `J::Matrix{Float64}` : Jacobian matrix to update.
- `grid` : Sparse grid for differentiation.
- `dyn_endogenous::Vector{Float64}` : Current endogenous state variables.
- `weight::Float64` : Quadrature weight for integration.
- `sgws::SparsegridsWs` : Solver workspace.

# Effect
- Modifies `J` in place with computed derivatives.

# Returns
- `J::Matrix{Float64}` : Updated Jacobian matrix.
"""
function forward_looking_equation_derivatives!(J, grid, dyn_endogenous, weight, sgws)
    # Unpack required variables
    @unpack dynamic_state_variables, tmp_state_variables, system_variables, ids,
            forward_block, forward_variable_jacobian, sgmodel,
            policy_jacobian, future_policy_jacobian, evalPt, M1 = sgws
    @unpack endogenous_nbr, parameters, steadystate, tempterms = sgmodel

    # Compute policy function Jacobian
    tmp_state_variables .= dyn_endogenous[dynamic_state_variables .+ endogenous_nbr]
    derivate!(policy_jacobian, grid, tmp_state_variables)

    # Compute future policy Jacobian
    copy_from_submatrix!(future_policy_jacobian, policy_jacobian, ids.state_in_system, ids.forward_in_system)

    # Update forward-block Jacobian
    forward_block.update_jacobian!(tempterms, forward_block.jacobian.nzval, dyn_endogenous, sgmodel.exogenous, parameters, steadystate)

    # Extract forward-looking variable Jacobian entries
    @views begin
        for (j, k) in enumerate(ids.forward_in_system)
            k1 = 2 * endogenous_nbr + system_variables[k]
            forward_variable_jacobian[:, j] .= forward_block.jacobian[:, k1]
        end
    end

    # Compute expected future policy effect
    mul!(M1, forward_variable_jacobian, future_policy_jacobian')

    # Update final Jacobian matrix
    for i in axes(J, 1)
        for (j, k) in enumerate(system_variables)
            J[i, j] += weight * forward_block.jacobian[i, k + endogenous_nbr]
        end
        for (j, k) in enumerate(ids.system_in_state)
            J[i, k] += weight * M1[i, j]
        end
    end

    return J
end

"""
    forward_looking_equation_residuals!(forward_residuals, policy, state, grid, sgws)

Computes residuals of forward-looking equations for each integration node.

# Arguments
- `forward_residuals::Matrix{Float64}` : Output matrix to store residuals.
- `policy::Vector{Float64}` : Policy function values.
- `state::Vector{Float64}` : Current state variables.
- `grid` : Sparse grid for evaluation.
- `sgws::SparsegridsWs` : Solver workspace.

# Effect
- Modifies `forward_residuals` in place.

# Returns
- `nothing`
"""
function forward_looking_equation_residuals!(forward_residuals, policy, state, grid, sgws)
    # Update `dyn_endogenous_vector` with the next-period state
    set_dyn_endogenous_vector!(sgws.dyn_endogenous_vector, policy, state, grid, sgws)

    # Unpack necessary variables
    @unpack monomial, forward_block, dyn_endogenous_vector, sgmodel = sgws
    @unpack dyn_endogenous, exogenous, parameters, steadystate, tempterms = sgmodel

    # Reset variables
    fill!(exogenous, 0.0)
    fill!(forward_residuals, 0.0)

    # Number of integration nodes
    n_nodes = length(monomial.nodes)

    for i in 1:n_nodes
        try

            # Compute residuals for forward-looking equations
            @views forward_block.get_residuals!(
                tempterms, forward_residuals[:, i], dyn_endogenous_vector[i], exogenous, parameters, steadystate
            )

        catch err
            @error "Forward-looking residual computation failed at node $i." 
                    exception=err
                    tempterms=tempterms 
                    forward_residuals=forward_residuals[:, i] 
                    dyn_endogenous_vector=dyn_endogenous 
                    exogenous=exogenous 
                    parameters=parameters 
                    steadystate=steadystate

            # Optionally, rethrow the error or recover
            rethrow()
        end
    end

    return nothing
end


"""
    expectation_equations_derivatives(forward_jacobian, policy, state, grid, sgws)

Computes derivatives of expectation equations by integrating forward-looking equation derivatives.

# Arguments
- `forward_jacobian::Matrix{Float64}` : Output matrix to store computed Jacobian derivatives.
- `grid` : Sparse grid for evaluation.
- `sgws::SparsegridsWs` : Solver workspace.

# Effect
- Modifies `forward_jacobian` in place.

# Returns
- `nothing`
"""
function expectation_equations_derivatives!(forward_jacobian, grid, sgws)
    # Unpack required variables
    @unpack monomial, dyn_endogenous_vector = sgws

    # Number of integration nodes
    n_nodes = length(monomial.nodes)

    # Reset the Jacobian matrix
    fill!(forward_jacobian, 0.0)

    for i in 1:n_nodes
        # Compute forward-looking equation derivatives
        forward_looking_equation_derivatives!(
            forward_jacobian, grid, dyn_endogenous_vector[i], monomial.weights[i], sgws
        )
    end

    return nothing
end

"""
    sysOfEqs!(residuals, policy, state, grid, sgws)

Computes the residuals of the system of equations.

# Arguments
- `residuals::Vector{Float64}` : Output vector to store equation residuals.
- `policy::Vector{Float64}` : Policy function values.
- `state::Vector{Float64}` : Current state variables.
- `grid` : The sparse grid for evaluation.
- `sgws::SparsegridsWs` : The solver workspace.

# Effect
- Modifies `residuals` in place to store equation residuals.

# Returns
- `nothing`
"""
function sysOfEqs!(residuals, policy, state, grid, sgws)
    # Unpack required variables
    @unpack dynamic_state_variables, system_variables, bmcps, backward_block, forward_block, 
            dyn_endogenous_vector, sgmodel = sgws
    @unpack dyn_endogenous, endogenous_nbr, exogenous, parameters, steadystate, tempterms = sgmodel

    # Set state and policy variables in dyn_endogenous
    dyn_endogenous[dynamic_state_variables] .= state
    dyn_endogenous[system_variables .+ endogenous_nbr] .= policy

    # Number of forward-looking equations
    nfwrd = length(forward_block.variables)

    # Compute residuals using the local copy of dyn_endogenous
    @views begin
        expectation_equations!(residuals[1:nfwrd], policy, state, grid, sgws)
        backward_block.get_residuals!(tempterms, residuals[nfwrd + 1:end], dyn_endogenous, exogenous, parameters, steadystate)
    end

    # Reorder residuals according to complementarity constraints
    reorder!(residuals, bmcps)

    return nothing
end

"""
    reorder_rows!(x::AbstractMatrix, permutations)

Reorders the rows of matrix `x` according to `permutations`.

# Arguments
- `x::AbstractMatrix`: The matrix whose rows will be reordered.
- `permutations::Vector{Int}`: A permutation vector specifying the new row order.

# Effect
- Modifies `x` in place by reordering its rows.
"""
function reorder_rows!(x::AbstractMatrix, permutations)
    for c in axes(x, 2)
        @views reorder!(x[:, c], permutations)
    end
end

"""
    sysOfEqs_derivatives_update!(J, policy, state, grid, sgws)

Computes and updates the Jacobian matrix (`J`) for the system of equations.

# Arguments
- `J::Matrix{Float64}` : Jacobian matrix to update.
- `policy::Vector{Float64}` : Policy function values.
- `state::Vector{Float64}` : Current state variables.
- `grid` : The sparse grid for evaluation.
- `sgws::SparsegridsWs` : The solver workspace.

# Effect
- Updates `J` in place with computed derivatives.

# Returns
- `nothing`
"""
function sysOfEqs_derivatives_update!(J, policy, state, grid, sgws)
    # Unpack required variables
    @unpack dynamic_state_variables, system_variables, bmcps,
            backward_block, forward_block, sgmodel = sgws
    @unpack dyn_endogenous, endogenous_nbr, exogenous, parameters, steadystate, tempterms = sgmodel

    # Set state and policy variables in local dyn_endogenous copy
    dyn_endogenous[dynamic_state_variables] .= state
    dyn_endogenous[system_variables .+ endogenous_nbr] .= policy

    # Reset exogenous shocks to zero
    fill!(exogenous, 0.0)

    # Number of forward and backward equations
    nfwrd = length(forward_block.equations)
    nback = length(backward_block.equations)

    # Compute derivatives of expectation equations
    set_dyn_endogenous_vector!(sgws.dyn_endogenous_vector, policy, state, grid, sgws)
    @views expectation_equations_derivatives!(J[1:nfwrd, :], grid, sgws)

    # Update backward-block Jacobian
    @views begin
        backward_block.update_jacobian!(tempterms, backward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
        J[nfwrd .+ (1:nback), :] .= backward_block.jacobian[:, system_variables .+ endogenous_nbr]
    end

    # Reorder Jacobian rows for complementarity constraints
    reorder_rows!(J, bmcps)

    return nothing
end

"""
    PATHsolver_solve!(X, lb, ub, states, fx, J, oldgrid, sgws, show_trace, chunk)

Solves a nonlinear complementarity problem (NCP) using PATHSolver.

# Arguments
- `polGuess::Matrix{Float64}` : Solution matrix to update.
- `lb::Vector{Float64}` : Lower bounds.
- `ub::Vector{Float64}` : Upper bounds.
- `states::Matrix{Float64}` : State variables for solving.
- `fx::Vector{Float64}` : Function evaluation vector.
- `J::Matrix{Float64}` : Jacobian matrix.
- `oldgrid` : Previous iteration’s grid.
- `sgws::Vector{SparsegridsWs}` : Workspace for the current thread.
- `show_trace::Bool` : Enables debugging output.
- `chunk::UnitRange{Int}` : The subset of grid points to process.

# Effect
- Updates `X` in place with computed solutions.
"""
function PATHsolver_solve!(polGuess, lb, ub, states, fx, J, oldgrid, sgws, show_trace, chunk)
    # Define system function and Jacobian
    f!(r, x, state) = sysOfEqs!(r, x, state, oldgrid, sgws)
    JA!(Jx, x, state) = sysOfEqs_derivatives_update!(Jx, x, state, oldgrid, sgws)

    for i in chunk
        # View to avoid unnecessary copies
        @views state = states[:, i]
        @views fx .= polGuess[:, i]

        # Solve NCP using PATHSolver
        status, results, info = mcp_solve!(
            PathNLSolver(), 
            (r, x) -> f!(r, x, state), 
            (Jx, x) -> JA!(Jx, x, state), 
            J, lb, ub, fx, 
            silent=true, convergence_tolerance=1e-4
        )

        # Error handling
        if status != 1
            @warn "PATHSolver failed at iteration $i. Status: $status. Info: $info"
            error("SparseGrids: solution update failed at index $i")
        end

        # Store solution in X
        @views polGuess[:, i] .= results
    end

    return nothing
end

"""
    NonlinearSolver_solve!(polGuess, states, J, oldgrid, method, sgws, show_trace, chunk)

Solves the nonlinear system using `NonlinearSolve.jl` for each state in `chunk`.

# Arguments
- `polGuess::Matrix{Float64}` : Solution matrix to update.
- `states::Matrix{Float64}` : State variables for solving.
- `J::Matrix{Float64}` : Jacobian matrix.
- `oldgrid` : Previous iteration’s grid.
- `method` : Solver method (e.g., Newton-Raphson).
- `sgws::Vector{SparsegridsWs}` : Workspace structure for each thread.
- `show_trace::Bool` : Enables debugging output.
- `chunk::UnitRange{Int}` : The subset of grid points to process.

# Effect
- Updates `polGuess` in place with computed solutions.

# Returns
- nothing
"""
function NonlinearSolver_solve!(polGuess, states, J, oldgrid, method, sgws, ftol, show_trace, chunk)
    # Define system function and Jacobian once
    f!(r, x, state) = sysOfEqs!(r, x, state, oldgrid, sgws)
    JA!(Jx, x, state) = sysOfEqs_derivatives_update!(Jx, x, state, oldgrid, sgws)

    # Create reusable NonlinearFunction object
    nlf = NonlinearFunction(f!, jac=JA!, jac_prototype=J)

    # Set show_trace as a `Val` to optimize performance
    show_trace = Val(show_trace)

    for i in chunk
        @views state = states[:, i]
        @views x = polGuess[:, i]

        resid = zeros(size(polGuess,1))
        f!(resid,x,state)
        JA!(J,x,state)

        # Define and solve the nonlinear problem
        nlp = NonlinearProblem(nlf, x, state)
        res = NonlinearSolve.solve(nlp, method, show_trace=show_trace, abstol=ftol)

        # Store the computed solution
        @views polGuess[:, i] .= res.u
    end

    return nothing
end

"""
    NLsolve_solve!(polGuess, lb, ub, fx, J, states, oldgrid, sgws, ftol, show_trace, chunk)

Solves a nonlinear system with bounds using `NLsolve.jl` for each state in `chunk`.

# Arguments
- `X::Matrix{Float64}` : Solution matrix to update.
- `lb::Vector{Float64}` : Lower bounds.
- `ub::Vector{Float64}` : Upper bounds.
- `fx::Vector{Float64}` : Function evaluation vector.
- `J::Matrix{Float64}` : Jacobian matrix.
- `states::Matrix{Float64}` : State variables for solving.
- `oldgrid` : Previous iteration’s grid.
- `sgws::Vector{SparsegridsWs}` : Workspace for each thread.
- `ftol::Float64` : Function tolerance for solver convergence.
- `show_trace::Bool` : Enables debugging output.
- `chunk::UnitRange{Int}` : The subset of grid points to process.

# Effect
- Updates `X` in place with computed solutions.

# Returns
- nothing
"""
function NLsolve_solve!(X, lb, ub, fx, J, states, oldgrid, sgws, ftol, show_trace, chunk)
    # Define system function and Jacobian outside the loop
    f!(r, x, state) = sysOfEqs!(r, x, state, oldgrid, sgws)
    JA(Jx, x, state) = sysOfEqs_derivatives_update!(Jx, x, state, oldgrid, sgws)

    for i in chunk
        @views state = states[:, i]
        @views x = polGuess[:, i]

        # Wrap functions in `OnceDifferentiable`
        df = OnceDifferentiable(f!, JA!, x, fx, J)

        try
            # Solve using NLsolve
            res = NLsolve.mcpsolve(df, lb, ub, x, ftol=ftol, show_trace=show_trace)

            if !res.f_converged
                error("SparseGrids: NLsolve did not converge at index $i")
            end

            x .= res.zero

        catch e
            @error "NLsolve failed at index $i. Recovering with interpolate." exception=e state=state lb=lb ub=ub x=x
            x .= interpolate(oldgrid, state)  # Fallback evaluation
        end
    end

    return nothing
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
    chunk = 1:ns
    i = 1
    # for (i, chunk) in enumerate(chunks)
        # push!(tasks, Threads.@spawn begin
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
        # end)
    # end

    fetch.(tasks)

    # Get BLAS to use the number of threads before SG_NLsolve! call
    BLAS.set_num_threads(previous_blas_threads)
    
    return polGuess
end
