using Combinatorics

"""
    mutable struct DDSG

A data structure representing a high-dimensional sparse grid, incorporating dimensional decomposition and sparse grid methods.

# Fields
- `dim::Int` : The number of dimensions in the grid.
- `dof::Int` : Degrees of freedom for the grid's output.
- `l_min::Int` : The minimum level of the sparse grid refinement.
- `l_max::Int` : The maximum level of the sparse grid refinement.
- `k_max::Int` : The maximum length of combination terms used in the decomposition.
- `order::Int` : The polynomial order for local basis functions in the sparse grid.
- `rule::String` : The rule used for constructing the sparse grid.
- `domain::Matrix{Float64}` : A two-column matrix defining the lower and upper bounds of the domain in each dimension.
- `centroid::Vector{Float64}` : The central point of the domain, computed as the midpoint of each dimension.
- `lookup::Vector{DDSGInxMap}` : A vector of dictionaries storing grid structures associated with different decomposition terms.
- `rlookup::Vector{DDSGInx}` : A vector of tuples representing reverse lookup mappings.
- `grid::Vector{TasmanianSG}` : A vector storing the sparse grid structures.
- `coeff::Vector{Float64}` : A vector storing decomposition coefficients associated with the grid.
- `grid_points::Vector{Int}` : A vector tracking the number of grid points in each component function.
- `X0::Vector{Float64}` : The anchor point in the domain for function evaluation.
- `Y0::Vector{Float64}` : The function values at the anchor point.
- `coeff0::Float64` : The coefficient corresponding to the zeroth-order decomposition component.
- `is_ddsg::Bool` : A flag indicating whether the grid follows a dimensional decomposition structure (`true`) or is a standard sparse grid (`false`).

# Constructors
- `DDSG(; kwargs...)` : Creates a `DDSG` object with user-defined parameters. If `domain` is not provided, it defaults to `[0,1]` for each dimension.
- `DDSG(other::DDSG)` : Creates a copy of another `DDSG` object, ensuring proper handling of mutable fields through shallow or deep copies where necessary.
"""
mutable struct DDSG
    dim::Int                   # Dimension of the grid
    dof::Int                   # Degrees of freedom for grid output
    l_min::Int                 # Minimum sparse grid level
    l_max::Int                 # Maximum sparse grid level
    k_max::Int                 # Maximum combination length
    order::Int                 # Polynomial order for local basis
    rule::String               # Rule for constructing the sparse grid
    domain::Matrix{Float64}    # Domain as a 2-column matrix
    centroid::Vector{Float64}  # Center of the domain as a vector
    lookup::Dict{Int, Dict{Vector{Int}, Int}} # Vector of dictionaries to hold grid structures
    grid::Vector{TasmanianSG}  # Vector to hold grid structures
    coeff::Vector{Float64}     # Vector to hold grid DD coefficients
    grid_points::Vector{Int}   # Vector to hold grid points
    X0::Vector{Float64}        # Anchor point X0
    Y0::Vector{Float64}        # Values at anchor point
    coeff0::Float64            # Coefficient for the zeroth DD order component function
    is_ddsg::Bool              # is ddsg grid

    # Define constructor for the struct
    function DDSG( 
        dim,
        dof,
        l_min,
        l_max,
        k_max;
        order = 1,
        rule = "localp",
        domain = nothing,
        grid = Vector{TasmanianSG}(),
        coeff = Vector{Float64}(),
        grid_points = Vector{Int}(),
        X0 = Vector{Float64}(),
        Y0 = Vector{Float64}(),
        coeff0 = 0.0,
        is_ddsg = true
    )
        # Set the domain limits for each dimension
        if isnothing(domain)
            domain = zeros(dim,2)
            domain[:,2].=1.0
        end

        # Set center as the midpoint of each dimension in the domain
        centroid = vec(@. 0.5*(domain[:, 1] + domain[:, 2]))

        # Initialize empty dictionaries for the resized lookup
        lookup = Dict()
        for k in 1:k_max
            lookup[k] = Dict{Vector{Int},Int}()
        end

        # If the maximum combination length coincides with the number of
        # dimensions, we are back to a standard sparse grid framework
        if dim==k_max
            is_ddsg = false
        end

        return new(
            dim, 
            dof, 
            l_min, 
            l_max, 
            k_max, 
            order, rule, domain, centroid, 
            lookup, grid, coeff, grid_points,
            X0, Y0, coeff0, is_ddsg
        )
    end

    # Copy constructor
    function DDSG(other::DDSG)
        # Immutable-type instances are automatically copied by the compiler. No
        # need to deepcopy them. This includes fields dim, dof, l_min, l_max,
        # k_max, order, rule, coeff0 and is_ddsg. For some mutable-type
        # instances with immutable-type fields, a shallow copy is enough. This
        # includes domain, centroid, coeff, grid_points, X0 and Y0. lookup and
        # grid fields need a deepcopy.
        return new(
            other.dim,
            other.dof,
            other.l_min,
            other.l_max,
            other.k_max,
            other.order,
            other.rule,
            copy(other.domain),
            copy(other.centroid),
            deepcopy(other.lookup),
            [copyGrid(grid) for grid in other.grid],
            copy(other.coeff),
            copy(other.grid_points),
            copy(other.X0),
            copy(other.Y0),
            other.coeff0,
            other.is_ddsg
        )
    end
end

"""
    DDSG_init_u!(self, u, coeff)

Initializes a Tasmanian sparse grid for a given combination of dimensions in the DDSG structure.

# Arguments
- `self::DDSG` : The DDSG instance being initialized.
- `u::Vector{Int}` : The subset of dimensions for which the sparse grid is created.
- `coeff::Float64` : The initial coefficient value associated with the combination.

# Description
This function initializes a sparse grid using Tasmanian's `TasmanianSG` for the given subset of dimensions `u`. The grid is configured with a local polynomial basis according to the DDSG parameters.

- The domain is transformed to match the selected dimensions.
- A unique `DDSGInx` key is generated for `u` and stored in the lookup structures.
- The grid, coefficients, and metadata are updated within the `DDSG` instance.

"""
function DDSG_init_u!(self, u, coeff)
    # Initialize a Tasmanian sparse grid for the current combination
    grid = Tasmanian.TasmanianSG(length(u), self.dof, self.l_min)
    makeLocalPolynomialGrid!(grid, order=self.order, rule=self.rule)

    # Match grid to domain dimensions
    setDomainTransform!(grid, self.domain[u,:])

    # Assign default values to vectors and map vector indices using lookup
    push!(self.grid, grid)
    push!(self.coeff, coeff)
    push!(self.grid_points, 0)
    push!(self.lookup[length(u)], u => length(self.grid))
    nothing
end

"""
    DDSG_init!(self::DDSG)

Constructs and initializes the DDSG structure, setting up the sparse grid components.

# Arguments
- `self::DDSG` : The DDSG instance to be initialized.

# Description
This function populates the DDSG structure by iterating over possible dimension combinations up to `k_max`.
- If `is_ddsg` is `true`, it iterates through all possible subset sizes `k`, generating grids for each subset of dimensions.
- If `is_ddsg` is `false`, the DDSG reduces to a standard sparse grid covering all dimensions.

Each subset of dimensions is passed to `DDSG_init_u!`, which initializes the corresponding sparse grid and stores the necessary references.

"""
function DDSG_init!(self::DDSG)
    # Otherwise, loop over each possible combination length up to k_max and
    # fill the DDSG accordingly
    if self.is_ddsg
        map(k->map(u->DDSG_init_u!(self, u, 0.0), combinations(1:self.dim,k)), 1:self.k_max)
    # If the DDSG instance is equivalent to a classic sparse grid case, load the
    # sole combination u = 1:self.dim
    else
        DDSG_init_u!(
            self,
            1:self.dim,
            1.0
        )
    end
    nothing
end

# # Method to print DDSG instance details
# function DDSG_print(self::DDSG)
#     println("\n=== DDSG Object Details START ===")
#     println("Dimension (dim): ", self.dim)
#     println("Degrees of freedom (dof): ", self.dof)
#     println("Level min (l_min): ", self.l_min)
#     println("Level max (l_max): ", self.l_max)
#     println("Polynomial order (order): ", self.order)
#     println("Rule (rule): ", self.rule)
#     println("Domain (domain): ", self.domain)
#     println("Domain centroid (centroid): ", self.centroid)
#     println("Maximum expansion order (k_max): ", self.k_max)
#     println("Anchor point (X0): ", self.X0)
#     println("Function values at anchor point (Y0): ", self.Y0)
#     println("DD coefficient for 0-order (coeff0): ", self.coeff0)
#     println("Is DDSG (is_ddsg): ", self.is_ddsg)
#     println("Grid details:")
#     for k in 1:length(self.lookup)
#         println("- Order (k): ", k)
#         for (u, inx) in self.lookup[k]
#             println("-- Index (inx): ", inx)
#             println("-- Index (u): ", u)
#             println("-- DD coefficient (coeff[lookup[k][u]]): ", self.coeff[inx])
#             println("-- Grid points (grid_points[lookup[k][u]]): ", self.grid_points[inx])
#             println(self.grid[inx])
#         end
#     end
#     println("=== DDSG Object Details END ===")
# end


# # Function to print matrix with fixed number of decimal places
# function ppnice(matrix, digits=4)
#     for row in eachrow(matrix)
#         for element in row
#             @printf("%.*f  ", digits, element)
#         end
#         println()  # Newline after each row
#     end
# end

"""
    DDSG_refine!(grid, refinement_tol, refinement_type, scale_corr_vec)

Performs adaptive sparse grid refinement based on surplus coefficients.

# Arguments
- `grid` : The Tasmanian sparse grid to be refined.
- `refinement_tol::Float64` : The tolerance threshold for surplus-based refinement.
- `refinement_type::String` : The refinement strategy to use (e.g., `"classic"`).
- `scale_corr_vec::Union{Vector{Float64}, Nothing}` : An optional scale correction vector for adjusting refinement.

# Description
This function refines the given sparse grid using a **surplus-based refinement strategy**:
- If `scale_corr_vec` is `nothing`, the grid is refined using standard surplus-based selection.
- If `scale_corr_vec` is provided, it applies a scale correction to modify the refinement behavior.

The function modifies the `grid` **in place**, ensuring that refinement follows the specified tolerance and refinement type.

# Notes
- The `output=-1` parameter ensures that all outputs are considered for refinement.
- `scale_corr_mat` is constructed by repeating `scale_corr_vec` across all loaded points to apply a consistent scaling across dimensions.
"""
function DDSG_refine!(grid, refinement_tol, refinement_type, scale_corr_vec)
    if isnothing(scale_corr_vec)
        setSurplusRefinement!(grid, refinement_tol, output=-1, refinement_type=refinement_type)
    else
        scale_corr_mat = repeat(scale_corr_vec, 1, getNumLoaded(grid))
        setSurplusRefinement!(grid, refinement_tol, output=-1, refinement_type=refinement_type, scale_correction=scale_corr_mat')
    end
end

"""
    evaluate_at_̄x_xᵥ(F, X0, X, v_inx)

Evaluates the function `F` at modified grid points by replacing selected dimensions in the anchor point `X0`.

# Arguments
- `F::Function` : The function to evaluate, mapping input points to output values.
- `X0::Vector{Float64}` : The anchor point in the domain, representing the reference point for decomposition.
- `X::Matrix{Float64}` : The matrix of needed points to be substituted into `X0` for evaluation.
- `v_inx::Vector{Int}` : The indices of the dimensions in `X0` that should be replaced with `X`.

# Description
This function constructs the modified input matrix `X_full` by:
1. Initializing `X_full` as a repetition of `X0`, ensuring the correct shape.
2. Replacing the dimensions indexed by `v_inx` with corresponding columns from `X`.
3. Evaluating the function `F` at the modified grid points and returning the result.

# Returns
- A vector or matrix containing function evaluations `F(̄x - xᵥ)`, where `xᵥ` represents the replaced values.

# Notes
- The function is typically used in the DDSG framework to compute higher-order decomposition terms.
- The input `X` should have the same number of columns as the intended number of evaluation points.
"""
function evaluate_at_̄x_xᵥ(F, X0, X, v_inx)
    # Build x = ̄x \ xᵥ
    # Initialize x = ̄x
    N = size(X,2)
    X_full = repeat(X0, 1, N)
    # Replace the v-indexed values in x
    @inbounds X_full[v_inx, :] .= X
    # Return F(̄x \ xᵥ)
    return F(X_full)
end

@doc raw"""
    update_ddsg_coefficients!(self::DDSG)

Updates the DDSG decomposition coefficients using the inclusion-exclusion principle.

# Arguments
- `self::DDSG` : The DDSG instance whose decomposition coefficients are being updated.

# Description
This function updates the coefficients stored in `self.coeff` and `self.coeff0` by applying the **inclusion-exclusion principle**. It iterates over all stored combinations of dimensions in `self.lookup` and modifies the corresponding coefficients in the hierarchical decomposition.

### **Steps**
1. **Iterate over all stored combinations**:
   - Loops through `self.lookup`, which contains all `k`-dimensional subsets `u`.
2. **For each subset `u`**, update coefficients for its non-empty subsets `v ⊆ u`:
   - Iterates over all non-empty subsets `v` of `u`.
   - Uses **alternating signs** based on the length difference `|u| - |v|`.
   - Updates `self.coeff[inx]` accordingly.
3. **Update the coefficient for the empty subset `∅`**:
   - The base coefficient `self.coeff0` is adjusted using the same alternating-sign rule.
4. **Final adjustment for `u = ∅`**:
   - Ensures `self.coeff0` is updated properly.

The coefficient update follows the formula:
\[
\text{coeff}[v] += (-1)^{|u| - |v|}
\]
which ensures a **hierarchical decomposition** where higher-order terms account for lower-order interactions.

# Notes
- This function **modifies `self.coeff` and `self.coeff0` in place**.
- It ensures that the **DDSG decomposition remains consistent** with the dimensional decomposition framework.
- The final `nothing` return value explicitly signifies that the function modifies the DDSG instance in place without returning a value.

"""
function update_ddsg_coefficients!(self)
    # For all dimensions up to the cut
    for len_u in 1:self.k_max
        # For all k-dimension tuples u, k ≥ 1
        for (u, _) in self.lookup[len_u]
            for len_v in len_u:-1:1
                # For all v ⊆ u, v ≠ ∅, update the relevant coefficients
                # stored in ddsg.coeff
                for v in combinations(u,len_v)
                    inx = self.lookup[len_v][v]
                    self.coeff[inx] += (-1)^(len_u-len_v)
                end
            end
            # For v = ∅, update the relevant coefficient coeff0
            self.coeff0 += (-1)^len_u
        end
    end
    # For u = ∅
    self.coeff0 += (-1)^0
    nothing
end

"""
    DDSG_build!(
        self::DDSG; 
        F::Function, 
        X0::Vector{Float64},
        refinement_tol::Float64,
        refinement_type::String="classic",
        scale_corr_vec::Union{Vector{Float64},Nothing}=nothing
    )

Constructs the DDSG interpolant by evaluating the function at required points and refining the sparse grids.

# Arguments
- `self::DDSG` : The DDSG instance being built.
- `F::Function` : The target function to approximate, mapping input points to output values.
- `X0::Vector{Float64}` : The anchor point in the domain, used for dimensional decomposition.
- `refinement_tol::Float64` : The refinement tolerance for adaptive sparse grid refinement.
- `refinement_type::String="classic"` : The refinement strategy used to refine the sparse grid.
- `scale_corr_vec::Union{Vector{Float64}, Nothing}=nothing` : An optional scale correction vector for refining the sparse grid.

# Description
This function constructs the DDSG approximation by:
1. **Computing the anchor evaluation** `F(X0)`, which serves as the base component of the decomposition.
2. **Evaluating all component functions** `F(X̄ \\ xᵥ)`, corresponding to different sparse grid levels.
3. **Refining the sparse grids** iteratively based on surplus coefficients.
4. **Updating DDSG coefficients** using the combinatorial inclusion-exclusion principle to ensure proper decomposition.
5. **Handling the case where DDSG reduces to a standard sparse grid**, setting appropriate coefficients.

The process involves:
- Iterating over all dimensional combinations stored in `lookup`, fetching grid structures, and computing function values at needed points.
- Replacing the corresponding dimensions in `X0` to evaluate function values for subspaces.
- Refining the sparse grid adaptively based on the surplus refinement strategy.
- Updating decomposition coefficients in `coeff` using combinatorial rules.

# Notes
- If `self.is_ddsg` is `true`, the function applies the hierarchical decomposition approach.
- If `self.is_ddsg` is `false`, the DDSG reduces to a standard sparse grid and assigns the last coefficient to `1.0`.
"""
function DDSG_build!(
    self::DDSG,
    F::Function,
    X0::Vector{Float64};
    refinement_tol=0.,
    refinement_type="classic",
    scale_corr_vec=Vector{Float64}()
)
    # Set ̄x
    self.X0 = copy(X0)
    # Compute F(̄x), the first term of the summation
    self.Y0 = vec(F(X0))
    # Compute all the higher-dimension terms F(̄x \ xᵥ) for all possible v ⊆
    # {1,…, k_max}.
    for k in 1:self.k_max
        # For all v, inx contains the position of the associated grid in
        # self.grid, self.grid_points and self.coeff.
        for (v, inx) in self.lookup[k]
            # Load the grid
            grid = self.grid[inx]
            self.grid_points[inx] = 0
            # Refinement based on surplus coefficients
            for l in (self.l_min):(self.l_max)
                ((getNumNeeded(grid) < 1) && break)
                # Get and prepare needed points for the v-grid
                X = getNeededPoints(grid)
                # Compute F(̄x \ xᵥ)
                Y_val = evaluate_at_̄x_xᵥ(F, self.X0, X, v)
                # Load the associated values on the v-grid.
                loadNeededPoints!(grid, Y_val)
                # Increment the number of the v-grid points
                self.grid_points[inx]+=size(X,2)
                # Adaptive sparse refinement
                DDSG_refine!(grid, refinement_tol, refinement_type, scale_corr_vec)
            end
        end
    end
    if self.is_ddsg
        update_ddsg_coefficients!(self)
    else
        self.coeff0=0.0
        self.coeff[end]=1.0
    end
    nothing
end

"""
    DDSG_evaluate!(ddsg::DDSG; Y::AbstractVecOrMat{Float64}, X::AbstractVecOrMat{Float64})

Evaluates the DDSG interpolant at given input points `X` and stores the result in `Y`.

# Arguments
- `ddsg::DDSG` : The DDSG instance to be evaluated.
- `Y::AbstractVecOrMat{Float64}` : The output storage for the evaluation results.
- `X::AbstractVecOrMat{Float64}` : The input points where the DDSG interpolant is evaluated.

# Description
This function evaluates the DDSG approximation at specified input points `X` and stores the results in `Y`. It supports **both single-point evaluation and batch evaluation**:

1. **Detects whether `X` is a batch (matrix) or a single-point vector**.
2. **Handles the case where DDSG reduces to a standard sparse grid**:
   - If `ddsg.is_ddsg == false`, it evaluates the function directly using the last sparse grid stored in `ddsg.grid`.
3. **Evaluates the full DDSG decomposition**:
   - Initializes `Y` with the function value at the anchor point `ddsg.Y0`, scaled by `ddsg.coeff0`.
   - Iterates over all component functions in `ddsg.rlookup`, evaluating them at their respective dimensions.
   - Aggregates contributions using the DDSG decomposition coefficients.

# Behavior
- If `X` is a **matrix**, `Y` is updated column-wise using batch evaluations.
- If `X` is a **vector**, `Y` is updated using standard evaluation.
- Uses `evaluateBatch!` for batched evaluations to optimize performance.
- Uses `evaluate(grid, vec(X_partial))` for single-point evaluation.

# Notes
- The function modifies `Y` **in place**, avoiding unnecessary allocations.
- Ensures compatibility with both **standard sparse grids** and **DDSG hierarchical decompositions**.
- The coefficient `ddsg.coeff0` is applied to the base function evaluation before summing contributions from hierarchical terms.
"""
function DDSG_evaluate!(Y::AbstractVecOrMat{Float64}, ddsg::DDSG, X::AbstractVecOrMat{Float64})
    is_batch = X isa AbstractMatrix
    if !ddsg.is_ddsg
        if is_batch
            evaluateBatch!(Y, ddsg.grid[end], X) 
        else
            Y .= evaluate(ddsg.grid[end], X)
        end
    else
        if is_batch
            for j in 1:size(Y, 2)
                Y[:,j] .= ddsg.Y0
            end
        else
            Y .= ddsg.Y0
        end
        Y .= Y .* ddsg.coeff0
        Y_buffer = similar(Y)
        for k in 1:ddsg.k_max
            for (u, inx) in ddsg.lookup[k]
                if is_batch
                    evaluateBatch!(Y_buffer, ddsg.grid[inx], X[u,:])
                else
                    Y_buffer .= evaluate(ddsg.grid[inx], X[u])
                end
                @. Y += Y_buffer * ddsg.coeff[inx]
            end
        end
    end
end

"""
    DDSG_evaluate(ddsg::DDSG; X::AbstractVecOrMat{Float64}) -> Matrix{Float64}

Evaluates the DDSG interpolant at given input points `X` and returns the computed values.

# Arguments
- `ddsg::DDSG` : The DDSG instance to be evaluated.
- `X::AbstractVecOrMat{Float64}` : The input points where the DDSG interpolant is evaluated.

# Returns
- `Y::Matrix{Float64}` : The evaluation results, where each column corresponds to the function evaluation at a given input point.

# Description
This function computes the DDSG approximation at specified input points `X` and returns the results in a newly allocated matrix `Y`. It serves as a **non-mutating wrapper** around `DDSG_evaluate!`, which performs in-place evaluations.

### **Steps:**
1. Initializes an output matrix `Y` of shape `(ddsg.dof, size(X,2))`, where:
   - `ddsg.dof` is the number of degrees of freedom (output variables).
   - `size(X,2)` is the number of input points to evaluate.
2. Calls `DDSG_evaluate!` to populate `Y` with the evaluation results.
3. Returns `Y`.

# Notes
- This function is useful when an explicit return value is needed instead of modifying a preallocated output.
- Internally, it delegates to `DDSG_evaluate!`, ensuring consistency in how evaluations are performed.
"""
function DDSG_evaluate(ddsg::DDSG, X::AbstractVecOrMat{Float64})
    Y = zeros(ddsg.dof,size(X,2))
    DDSG_evaluate!(Y,ddsg,X)
    return Y
end

function DDSG_differentiate!(J::Matrix{Float64}, ddsg::DDSG, x::Vector{Float64})
    if !ddsg.is_ddsg
        differentiate!(J,ddsg.grid[end],x)
    else
        for k in 1:ddsg.k_max
            for (u,inx) in ddsg.lookup[k]
                grid = ddsg.grid[inx]
                J_buffer = similar(J, k, ddsg.dof)
                x_partial = x[u]
                differentiate!(J_buffer,grid,x_partial)
                J[u,:] .+= J_buffer * ddsg.coeff[inx]
            end
        end
    end
end

function DDSG_nodes(ddsg::DDSG)
    N = sum(ddsg.grid_points)
    # Initialize an empty matrix to store the cumulative horizontal concatenation
    X_loc_total = Float64[]  # Start with an empty matrix
    for k in 1:ddsg.k_max
        for (u, inx) in ddsg.lookup[k]
            # Load the grid
            grid = ddsg.grid[inx]
            X_partial = getPoints(grid)
            # Repeat ddsg.X0 to match the size of X_partial
            X_loc = repeat(ddsg.X0, 1, size(X_partial, 2))
            X_loc[u, :] .= X_partial
            # Concatenate X_loc horizontally with X_loc_total
            if length(X_loc_total) == 0
                X_loc_total = X_loc
            else
                X_loc_total = hcat(X_loc_total, X_loc)
            end
        end
    end
    return X_loc_total  # Return the concatenated result
end

"""
    ddsg_ti_step!(newgrid, oldgrid, sgws)

Updates grid points by solving the nonlinear problem at newly needed grid locations.

# Arguments
- `newgrid` : The ddsg grid that needs new function values.
- `oldgrid` : The previous iteration's grid.
- `sgws::Vector{SparsegridsWs}` : The workspace containing all model and solver settings.

# Effect
- Returns the policy function values.
"""
function ddsg_ti_step!(states, oldPol, oldgrid, sgws)
    # Extract necessary data from sgws
    ws = sgws[1]
    @unpack system_variables, lb, ub, sgoptions, J, fx, sgmodel = ws
    @unpack mcp, method, solver, ftol, show_trace = sgoptions
    @unpack exogenous, parameters = sgmodel

    # Ensure solver is set
    solver = isnothing(solver) ? (mcp ? NLsolver : NonlinearSolver) : solver

    # Retrieve points requiring function values
    nstates = size(states, 2)

    # Initialize exogenous variables (assuming no autocorrelated shocks)
    fill!(exogenous, 0.0)

    # Solve the nonlinear problem
    polGuess = oldPol(states)
    SG_NLsolve!(polGuess, lb, ub, fx, J, states, oldgrid, sgws, solver, method, ftol, show_trace)
    return polGuess
end

make_updatedPol(f,g,w) = x->(1-w)*f(x)+w*g(x,f)

"""
    ddsg_time_iteration(oldDDSG, oldPol, sgws, scaleCorr, surplThreshold, dimRef, 
                               typeRefinement, maxiter, tol_ti, savefreq, timings, X_sample)

Performs the time iteration loop to refine the sparse grid and compute the policy function until convergence.

# Arguments
- `oldDDSG` : Sparse grid structure.
- `oldPol` : Policy function.
- `sgws::Vector{SparsegridsWs}` : Workspace structures for each thread.
- `scaleCorr::Vector{Float64}` : Scale correction vector.
- `surplThreshold::Float64` : Surplus threshold for adaptive refinement.
- `typeRefinement::String` : Type of refinement (e.g., "classic").
- `maxiter::Int` : Maximum number of iterations.
- `maxIterEarlyStopping::Int` : Number of iterations after which TI stops after TI convergence measure starts increasing
- `tol_ti::Float64` : Convergence criterion.
- `polUpdateWeight::Float64` : Weight of the current-step policy function when computing the updated policy function. The weight of the previous-step policy function is `1-polUpdateWeight`
- `savefreq::Int` : Frequency for saving the grid.
- `timings::Vector{Millisecond}` : Storage for iteration timings.

# Returns
- `grid::DDSG` : Updated DDSG.
- `average_time::Float64` : Average iteration time.
- `iter::Int` : Number of iterations
"""
function ddsg_time_iteration!(
    oldDDSG, oldPol, sgws, scaleCorr, surplThreshold,
    typeRefinement, maxiter, maxIterEarlyStopping,
    tol_ti, polUpdateWeight, savefreq
)
    # Initialize timing variables
    total_time = 0
    context.timings["sparsegrids"] = Vector{Millisecond}()
    timings = context.timings["sparsegrids"]
    iter = 1
    iterEarlyStopping = 0
    previousMetric = Inf

    while iter <= maxiter && iterEarlyStopping <= maxIterEarlyStopping
        tin = now()

        # Reset Jacobian for each workspace
        map(s -> fill!(s.J, 0.0), sgws)

        # Set new policy function
        newPol = (X,f) -> ddsg_ti_step!(X, f, oldDDSG, sgws)

        # Build the associated grid
        newDDSG = DDSG(oldDDSG.dim, oldDDSG.dof, oldDDSG.l_min, oldDDSG.l_max, oldDDSG.k_max; order=oldDDSG.order, rule=oldDDSG.rule, domain=copy(oldDDSG.domain))
        DDSG_init!(newDDSG)
        DDSG_build!(newDDSG, X->newPol(X,oldPol), newDDSG.centroid;refinement_type=typeRefinement, refinement_tol=surplThreshold, scale_corr_vec=scaleCorr)

        # Compute error
        metric = get_ddsg_metric(oldDDSG, newDDSG, length(sgws[1].system_variables))
        # Update the policy function
        # updatedPol = make_updatedPol(oldPol, newPol, polUpdateWeight)
        # updatedDDSG = DDSG(newDDSG.dim, newDDSG.dof, newDDSG.l_min, newDDSG.l_max, newDDSG.k_max; order=newDDSG.order, rule=newDDSG.rule, domain=copy(newDDSG.domain))
        # DDSG_init!(updatedDDSG)
        # DDSG_build!(updatedDDSG, updatedPol, newDDSG.centroid;refinement_type=typeRefinement, refinement_tol=surplThreshold, scale_corr_vec=scaleCorr)
        oldPol = make_updatedPol(oldPol, newPol, polUpdateWeight)
        updatedDDSG = DDSG(newDDSG.dim, newDDSG.dof, newDDSG.l_min, newDDSG.l_max, newDDSG.k_max; order=newDDSG.order, rule=newDDSG.rule, domain=copy(newDDSG.domain))
        DDSG_init!(updatedDDSG)
        DDSG_build!(updatedDDSG, oldPol, newDDSG.centroid;refinement_type=typeRefinement, refinement_tol=surplThreshold, scale_corr_vec=scaleCorr)

        oldDDSG = makeCopy(updatedDDSG)
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
        println("Iteration: $iter, Grid pts: $(countPoints(oldDDSG)), Metric: $metric, Computing time: $(iteration_walltime)")

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
    println("Last grid points: $(countPoints(oldDDSG))")
    println("Average iteration time (except first one): $average_time")

    return oldDDSG, average_time, iter
end

function get_ddsg_metric(oldDDSG, newDDSG, gridOut)
    # Get the points and the number of points from grid1
    aPoints = interpolationPoints(newDDSG)
    aNumTot = size(aPoints,2)

    # Evaluate the grid points on both grid structures
    polGuessNew = interpolate(newDDSG, Matrix(aPoints))
    polGuessOld = interpolate(oldDDSG, Matrix(aPoints))

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

    return metric 
end