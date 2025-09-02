@enum SGSolver NLsolver NonlinearSolver PATHSolver

@doc raw"""
    struct MonomialPowerIntegration

Represents a monomial quadrature rule for numerical integration in sparse grid approximation.

# Fields
- `nodes::Vector{Vector{Float64}}`  
  A list of integration nodes, where each node is a `d`-dimensional vector representing a point in the state space.
- `weights::Vector{Float64}`  
  A vector of corresponding integration weights, used to approximate expected values of functions over the state space.

# Purpose
The monomial quadrature rule is used to approximate integrals of the form:

\[
I[f] = \int_{\mathbb{R}^d} f(x) d\mu(x) \approx \sum_{i} w_i f(x_i)
\]

where:
- \( f(x) \) is the function being integrated (e.g., future policy function values).
- \( x_i \) are the quadrature nodes (`nodes`).
- \( w_i \) are the integration weights (`weights`).
- \( d\mu(x) \) is the probability measure (e.g., Gaussian distribution for stochastic shocks).

# Use Case
- In **dynamic stochastic models**, `MonomialPowerIntegration` is used for **expectation operators** when computing forward-looking equations.
- It replaces traditional Monte Carlo sampling by using **deterministic quadrature points** for efficient approximation.

# Example Usage
```julia
# Define quadrature nodes (e.g., 3 integration points in a 2D state space)
nodes = [[-0.5, 0.0], [0.0, 0.5], [0.5, 0.0]]

# Define corresponding integration weights
weights = [0.3, 0.4, 0.3]

# Create a MonomialPowerIntegration object
integration_rule = MonomialPowerIntegration(nodes, weights)

# Compute an expected value approximation
function compute_expectation(f, integration_rule)
    return sum(w * f(x) for (x, w) in zip(integration_rule.nodes, integration_rule.weights))
end

# Define a simple function to integrate
f(x) = sum(x)^2

# Compute the integral approximation
expected_value = compute_expectation(f, integration_rule)
println("Approximate integral: ", expected_value)
"""
struct MonomialPowerIntegration
    nodes::Vector{Vector{Float64}}
    weights::Vector{Float64}
end

struct SmolyakGHIntegration
    nodes   :: Vector{Vector{Float64}} # d × M (columns are nodes)
    weights :: Vector{Float64} # length M
end

struct BlockIndices
    equation_pointers::Vector{Int}
    variable_pointers::Vector{Int}
    expressions::Vector{Expr}
end

struct Block{F1 <: Function, F2 <: Function, F3 <: Function}
    assignment::Bool
    forward::Bool
    jacobian::SparseMatrixCSC{Float64, Int}
    assigment_fcn::F1
    jacobian_fcn::F2
    residual_fcn::F3
    indices::BlockIndices
end

abstract type AbstractBlock end
abstract type AbstractPreambleBlock <: AbstractBlock end

struct AssignmentBlock <: AbstractPreambleBlock
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
    set_endogenous_variables!::Function
    is_linear::Bool
end

struct PreambleBlock <: AbstractPreambleBlock
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
    get_residuals!::Function
    update_jacobian!::Function
end

struct ForwardBlock <: AbstractBlock
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
    get_residuals!::Function
    update_jacobian!::Function
end

struct BackwardBlock <: AbstractBlock    
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
    get_residuals!::Function
    update_jacobian!::Function
end

struct SimultaneousBlock <: AbstractBlock    
    equations::Vector{Int}
    variables::Vector{Int}
    expressions::Vector{Expr}
    jacobian::SparseMatrixCSC{Float64, Int}
    get_residuals!::Function
    update_jacobian!::Function
end

struct SGModel
    dyn_endogenous::Vector{Float64}       # dynamic endogenous variables in period t-1, t, t+1 
    shift_dyn_endogenous::Vector{Float64} # dynamic endogenous variables in period t, t+1 
    exogenous::Vector{Float64}            # exogenous variables in period t
    endogenous_nbr::Int                   # number of endogenous variables in period t
    exogenous_nbr::Int                    # number of exogenous variables
    parameters::Vector{Float64}           # model parameters    
    steadystate::Vector{Float64}          # steady state of the model
    tempterms::Vector{Float64}            # temporary terms
end

struct SGIndices
    forward_in_system::Vector{Int} # indices of forward looking variables in system variables
    state_in_system::Vector{Int}   # indices of system variables that belong to state
    state_variables::Vector{Int} 
    system_in_state::Vector{Int}   # indices of state variables that belong to system
 end 

function SGIndices(endogenous_nbr::Int, forward_variables, dynamic_state_variables, system_variables)
    forward_in_system = findall(in(forward_variables), system_variables)
    state_in_system = findall(in(system_variables), dynamic_state_variables)
    state_variables = [s > endogenous_nbr ? s - endogenous_nbr : s for s in dynamic_state_variables] 
    system_in_state = findall(in(state_variables), system_variables)
    return SGIndices(
        forward_in_system,
        state_in_system,
        state_variables,
        system_in_state
    )
end 

"""
    SGSolverOptions(mcp, method, solver, ftol, show_trace)

Stores solver configuration options for sparse grid approximation.

# Fields
- `mcp::Bool`  
  Whether to use a **Mixed Complementarity Problem (MCP)** solver (for occasionally binding constraints).
- `method`  
  Nonlinear solver method (e.g., `NewtonRaphson()`).
- `solver`  
  Solver type (e.g., `NLsolve` or `NonlinearSolve`).
- `ftol::Float64`  
  Function tolerance for solver convergence (stopping criterion).
- `show_trace::Bool`  
  Whether to print detailed solver logs for debugging.

# Purpose
Encapsulates solver settings in a structured way, ensuring cleaner function calls and easier debugging.
"""
struct SGSolverOptions
    mcp::Bool          # Use Mixed Complementarity Problem (MCP) solver?
    method             # Nonlinear solver method (e.g., Newton-Raphson)
    solver             # Solver type (e.g., NLsolve or NonlinearSolve)
    ftol::Float64      # Function tolerance (convergence criterion)
    show_trace::Bool   # Show detailed solver output?
end

"""
    struct Sev

Stores state-dependent variables used in sparse grid evaluations.

# Fields
- `dyn_endogenous_variable::Vector{Float64}`  
  A vector storing values of endogenous variables at the current node.
- `node::Vector{Float64}`  
  A vector representing the current state node in the sparse grid.

# Purpose
The `Sev` structure is used to store **state-dependent** and **policy-relevant** values for each grid node during sparse grid interpolation. It provides an efficient way to track **dynamic endogenous variables** as they evolve through the iteration process.
"""
struct Sev
    dyn_endogenous_variable::Vector{Float64}
    node::Vector{Float64}
end

"""
    struct UserPolicyGuess

Structure for user-provided policy function guesses

# Fields
- polFun::Function
  A vector-valued policy function.
  ```julia
  function (x)
    Z = x[1]
    K = x[2]
    return [
      α*β*exp(Z)*K^α*ss_l^(1-α),
      (1-α*β)*exp(Z)*K^α*ss_l^(1-α),
      c_pol(K,Z)+k_pol(K,Z),
      (1-α)*y_pol(K,Z)/K,
    ]
  ```
- `inputs::Vector{String}`
  Ordered list of the state variables used as inputs in the policy function field `function`, e.g.
  ```julia
  ["z","k"]
  ```
- `outputs`:: Vector{String}
  Ordered list of the policy function field `function` output
  ```julia
  ["k","c","y","rk"]
  ```
"""
struct UserPolicyGuess
    polFun::Function
    inputs::Vector{String}
    outputs::Vector{String}
    UserPolicyGuess() = new(x->nothing,Vector{String}(),Vector{String}())
    UserPolicyGuess(f,i,o) = new(f,i,o)
end

"""
    struct SparsegridsWs

Workspace structure for sparse grid approximation, storing essential model variables, Jacobians, and solver configurations.

# Fields
- `monomial::MonomialPowerIntegration`  
  Monomial quadrature rule used for numerical integration.
- `dynamic_state_variables::Vector{Int}`  
  Indices of dynamic state variables in the model.
- `system_variables::Vector{Int}`  
  Indices of all system variables used in approximation.
- `bmcps::Vector{Vector{Int}}`  
  Mixed Complementarity Problem (MCP) constraints, defining which variables have inequality constraints.
- `ids::SGIndices`  
  Struct storing variable indexing and forward-looking system details.
- `lb::Vector{Float64}`  
  Lower bounds for policy variables.
- `ub::Vector{Float64}`  
  Upper bounds for policy variables.
- `backward_block::BackwardBlock`  
  Block of equations defining backward-looking dynamics.
- `backward_variable_jacobian::Matrix{Float64}`  
  Jacobian matrix for backward equations.
- `forward_block::ForwardBlock`  
  Block of equations defining forward-looking dynamics.
- `preamble_block::AssignmentBlock`  
  Block storing initialization and parameter assignment equations.
- `system_block::SimultaneousBlock`  
  Block containing the system of nonlinear equations.
- `residuals::Vector{Float64}`  
  Residual vector for the full system of equations.
- `forward_jacobian::Matrix{Float64}`  
  Jacobian matrix for forward equations.
- `forward_points::Matrix{Float64}`  
  Integration nodes used in forward-looking expectations.
- `forward_residuals::Matrix{Float64}`  
  Residual matrix for forward integration.
- `forward_variable_jacobian::Matrix{Float64}`  
  Jacobian matrix for forward-looking variables.
- `J::Matrix{Float64}`  
  Global Jacobian matrix used in solver iterations.
- `fx::Vector{Float64}`  
  Function evaluation vector.
- `policy_jacobian::Matrix{Float64}`  
  Jacobian of policy function values with respect to state variables.
- `evalPt::Vector{Float64}`  
  Evaluation points for function approximation.
- `future_policy_jacobian::Matrix{Float64}`  
  Jacobian matrix for future policy functions.
- `dyn_endogenous_vector::Vector{Vector{Float64}}`  
  Stores dynamic endogenous variables at multiple integration nodes.
- `M1::Matrix{Float64}`  
  Temporary storage matrix for forward system computations.
- `policyguess::Matrix{Float64}`  
  Initial guess for policy function values.
- `tmp_state_variables::Vector{Float64}`  
  Temporary storage vector for state variables.
- `sev::Sev`  
  Struct holding the state-dependent values used in evaluations.
- `sgmodel::SGModel`  
  Sparse grid model containing parameters and steady-state values.
- `sgsolveroptions::SGSolverOptions`
  Configuration settings for nonlinear solver and sparse grid approximation.

# Purpose
The `SparsegridsWs` structure acts as a **workspace** for sparse grid-based time iteration, holding:
1. **Variable indexing** and **constraint handling**.
2. **Jacobian matrices** and **residual vectors** for nonlinear solving.
3. **Storage for forward-looking expectations and policy function approximations**.
4. **Configuration settings for solvers and sparse grids**.
"""
struct SparsegridsWs
    monomial::Union{MonomialPowerIntegration,SmolyakGHIntegration}
    dynamic_state_variables::Vector{Int}
    system_variables::Vector{Int} 
    bmcps::Vector{Vector{Int}}
    ids::SGIndices
    lb::Vector{Float64}
    ub::Vector{Float64}
    backward_block::BackwardBlock
    backward_variable_jacobian::Matrix{Float64}
    forward_block::ForwardBlock 
    preamble_block::AssignmentBlock 
    system_block::SimultaneousBlock
    residuals::Vector{Float64}
    forward_jacobian::Matrix{Float64}
    forward_points::Matrix{Float64} 
    forward_residuals::Matrix{Float64}    
    forward_variable_jacobian::Matrix{Float64}
    J::Matrix{Float64}
    fx::Vector{Float64}
    policy_jacobian::Matrix{Float64} 
    evalPt::Vector{Float64} 
    future_policy_jacobian::Matrix{Float64} 
    dyn_endogenous_vector::Vector{Vector{Float64}} 
    M1::Matrix{Float64} 
    policyguess::Matrix{Float64}
    tmp_state_variables::Vector{Float64}
    sev::Sev
    sgmodel::SGModel
    sgsolveroptions::SGSolverOptions
end

mutable struct SparsegridsResults
    average_error::Float64
    max_error::Float64
    iteration_number::Int
    average_iteration_time::Float64
    drawsnbr::Int
    burnin::Int
    equation_average_errors::Vector{Float64}
    equation_quantile_errors::Vector{Float64}
    ftol::Float64
    grid::TasmanianSG
    gridDepth::Int
    gridOrder::Int
    gridRule::String
    iterRefStart::Int
    maxRef::Int
    mcp::Bool
    method::NonlinearSolve.GeneralizedFirstOrderAlgorithm
    quantile_error::Float64
    quantile_probability::Float64
    solver::SGSolver
    surplThreshold::Float64
    function SparsegridsResults()
        average_error = 0
        max_error = 0
        iteration_number = 0
        average_iteration_time = 0
        drawsnbr = 0
        burnin = 0
        equation_average_errors = Vector{Float64}(undef, 0)
        equation_quantile_errors = Vector{Float64}(undef, 0)
        ftol = 0
        grid = TasmanianSG(0, 0, 0)
        gridDepth = 0
        gridOrder = 0
        gridRule = ""
        iterRefStat = 0
        maxRef = 0
        mcp = false
        method = NewtonRaphson()
        quantile_error = 0
        quantile_probability = 0
        solver = NonlinearSolver
        surplThreshold = 0
        new(average_error, iteration_number, max_error, average_iteration_time,
            drawsnbr, burnin, equation_average_errors,
            equation_quantile_errors, ftol, grid, gridDepth, gridOrder,
            gridRule, iterRefStat, maxRef, mcp, method, quantile_error,
            quantile_probability, solver, surplThreshold)
    end
end    