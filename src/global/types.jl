@enum SGSolver NLsolver NonlinearSolver PATHSolver

struct MonomialPowerIntegration
    nodes::Vector{Vector{Float64}}
    weights::Vector{Float64}
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

struct SGOptions
    mcp::Bool
    method
    solver
    ftol::Float64
    show_trace::Bool
end

struct Sev
    dyn_endogenous_variable::Vector{Float64}
    node::Vector{Float64}
end

struct SparsegridsWs
    monomial::MonomialPowerIntegration
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
    sgoptions::SGOptions
end

mutable struct SparsegridsResults
    average_error::Float64
    average_iteration_time::Float64
    drawsnbr::Int
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
        average_iteration_time = 0
        drawsnbr = 0
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
        new(average_error, average_iteration_time, drawsnbr,
            equation_average_errors, equation_quantile_errors, ftol,
            grid, gridDepth, gridOrder, gridRule, iterRefStat, maxRef,
            mcp, method, quantile_error, quantile_probability, solver,
            surplThreshold)
    end
end    

