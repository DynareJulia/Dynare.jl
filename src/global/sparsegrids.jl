include("blocks.jl")

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
    # Number of dimensions (capital stock and tfp for each country)
    gridDim = nstates
    # Number of outputs (capital policy & multiplier for each country + ARC multiplier)
    gridOut = nPols
    # Grid level (we start with a sparse grid of level 3)
    gridDepth = 2
    # 1=linear, 2=quadratic, 3=cubic
    gridOrder = 1
    # Type of base functions
    gridRule = "localp"
    # Surplus threshold
    surplThreshold = 1e-3
    # Number of maximum refinements
    maxRef = 1
    # Maximum Level of ASG
    maxRefLevel = gridDepth + maxRef
    # Outputs that are considered in the refinement process (-1 implies that all outputs are considered)
    dimRef = -1
    # Type of grid refinements
    typeRefinement = "classic"
    # Iteration to start at (start from scratch -- numstart =0; restart: numstart >0)
    numstart = 0
    # Maximum number of iterations
    maxiter = 300
    # Iteration at which the refinement starts
    iterRefStart = 25
    # Convergence criterion in time iteration
    tol_ti = 1e-4
    # Number of random draws for the error computation
    TT = 10000
    # Number of burn-in periods for SIMULATION
    burnin = 1000
    # Frequency of saving grid
    savefreq = 10
    ftol = 1e-5
    show_trace = false
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
    work = context.work
    (state_variables, predetermined_variables, system_variables,
     forward_equations_nbr, other_equations_nbr,
     forward_expressions_eqs, other_expressions_eqs, preamble_block) = make_block_functions(context)
    lb = Vector{Float64}(undef, 0)
    ub = Vector{Float64}(undef, 0)
    bmcps = Vector{Vector{Int}}(undef, 0)
    # must be consistent with SysOfEqs()
    block_eqs = union(forward_expressions_eqs, other_expressions_eqs)
    # computes lb, ub, bmcps
    block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables)
    equation_xref_list, variable_xref_list = xref_lists(context)
    params = context.work.params
    nstates = length(state_variables)

    ids = SGIndices(
        endogenous_nbr,
        model.i_fwrd_b,
        state_variables,
        system_variables
    )
    forward_system_variables = model.i_fwrd_b
    setdiff!(forward_system_variables, predetermined_variables .- endogenous_nbr)

    gridDim = nstates
    gridOut = length(system_variables)
    limits = context.work.limits
    nl = length(limits)
    gridDomain = zeros(gridDim, 2)
    sv = [ s > endogenous_nbr ? s - endogenous_nbr : s for s in state_variables ]
    for l in limits
        i = findfirst(x -> x == context.symboltable[string(l[1])].orderintype, sv) 
        if !isnothing(i)
            gridDomain[i,1] = l[2].min
            gridDomain[i,2] = l[2].max
        end
    end

    dynamicws = DynamicWs(context)

    # rows = union(forward_expressions_eqs, other_expressions_eqs)
    #J = make_sparse_submatrix(colptr, rowval, colval, rows, system_variables) 

    #Jexpected = make_sparse_submatrix(colptr, rowval, colval, forward_expressions_eqs, system_variables) 
    #Jother = make_sparse_submatrix(colptr, rowval, colval, other_equations_eqs, system_variables) 
    grid0, aPoints, aNum, nPols, =
        set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain);

    ################################################################################
    #                        Adaptivity parameters                                 #
    ################################################################################

    endogenous = zeros(3*endogenous_nbr)
    exogenous_nbr = model.exogenous_nbr
    exogenous = zeros(exogenous_nbr)
    parameters = work.params

    steadystate = context.results.model_results[1].trends.endogenous_steady_state

    polGuess = guess_policy(context, aNum, nPols, aPoints, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr, lb, ub)
    
    loadNeededPoints!(grid0, polGuess)

    params = context.work.params

    (nodes, weights) = monomial_power(exogenous_nbr)

    # Scale correction in the refinement process:
    @views scaleCorr = make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables .- endogenous_nbr])

    polGuess1 = copy(polGuess)
    pol = Vector{Float64}(undef, size(polGuess1, 1))
    for iter0 in 1:maxiter

        polGuess1 = copy(polGuess)
        grid1 = fresh_grid(gridDim, gridOut, gridDepth, gridDomain, gridOrder, gridRule)

        # Index of current grid level to control the number of refinements
        ilev = gridDepth
        while ((getNumNeeded(grid1) > 0) && (ilev <=  maxRefLevel))
            grid1 = ti_step(grid1, polGuess1, pol, grid0, nPols, nodes, weights, exogenous, params, 
                            steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, 
                            state_variables, system_variables, endogenous, bmcps, dynamicws, model, 
                            predetermined_variables, forward_system_variables,
                            forward_expressions_eqs, other_expressions_eqs, ids, lb, ub, preamble_block,
                            mcp, method, solver, ftol, show_trace)
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

    ws = DynamicWs(context)
    T = ws.temporary_values
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
                preamble_block!(T, y, random_shocks[:, p], params, steadystate)
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
    makeLocalPolynomialGrid!(grid, order = gridOrder, sRule = gridRule)
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
function make_guess_system(residuals, T, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
    endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(x)
        endogenous[system_variables] .= x
        endogenous[system_variables .+ endogenous_nbr] .= x
        DFunctions.SparseDynamicResidTT!(T, endogenous, exogenous, parameters, steadystate)
        other_block_(T, res1, endogenous, exogenous, parameters, steadystate)
        forward_block_(T, res2, endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function make_guess_system1(residuals, T, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
    endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(residuals, x, params)
        endogenous[system_variables] .= x
        endogenous[system_variables .+ endogenous_nbr] .= x
        DFunctions.SparseDynamicResidTT!(T, endogenous, exogenous, parameters, steadystate)
        other_block_(T, res1, endogenous, exogenous, parameters, steadystate)
        forward_block_(T, res2, endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function make_guess_system2(residuals, T, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
    endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(residuals, x)
        endogenous[system_variables] .= x
        endogenous[system_variables .+ endogenous_nbr] .= x
        DFunctions.SparseDynamicResidTT!(T, endogenous, exogenous, parameters, steadystate)
        other_block_(T, res1, endogenous, exogenous, parameters, steadystate)
        forward_block_(T, res2, endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    function ff(residuals, x)
        try
            f(residuals, x)
        catch
            residuals = Inf*ones(length(residuals))
        end
    end 
    return f
end

function make_guess_system(residuals, T, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr, permutations)
    endogenous[state_variables] .= state
    res1 = zeros(other_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(x, params)
        endogenous[system_variables] .= x
        endogenous[system_variables .+ endogenous_nbr] .= x
        DFunction.SparseDynamicResidTT!(T, endogenous, exogenous, parameters, steadystate)
        other_block!(T, res1, endogenous, exogenous, parameters, steadystate)
        forward_block!(T, res2, endogenous, exogenous, parameters, steadystate)
        residuals[1:other_equations_nbr] .= res1
        residuals[other_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function guess_policy(context, aNum, nPols, aPoints, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr, lb, ub)
    n =  length(system_variables)
    guess_values = zeros(n, aNum)
    @views x0 = copy(steadystate[system_variables .- endogenous_nbr])
    residuals = similar(x0)
    ws = DynamicWs(context)
    T = Vector{Dual{Float64}}(undef, length(ws.temporary_values))
#    f1 = make_guess_system2(residuals, T, aPoints[:, 1], endogenous, exogenous, parameters,
#                          steadystate, state_variables, system_variables,
#                          forward_equations_nbr, other_equations_nbr, endogenous_nbr)
#    problem = NonlinearProblem(f1, x0, params)
    for i in axes(aPoints, 2)
        @views state = aPoints[:, i]
        f1 = make_guess_system1(residuals, T, state, endogenous, exogenous, parameters,
                              steadystate, state_variables, system_variables,
                              forward_equations_nbr, other_equations_nbr, endogenous_nbr)
        f2 = make_guess_system2(residuals, T, state, endogenous, exogenous, parameters,
                              steadystate, state_variables, system_variables,
                              forward_equations_nbr, other_equations_nbr, endogenous_nbr)
#        res = NLsolve.mcpsolve(f, lb, ub, x0, show_trace=true, autodiff=:central)
#        !res.f_converged && throw(Error("sparsegrids: couldn't find the initial policy"))
        #        @views guess_values[:, i] .= res.zero
        problem = NonlinearProblem(f1, x0, params)
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
function forward_looking_equation_residuals(x, T, state, params, grid, nodes, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables, preamble_block)
    numNodes = size(nodes, 1)
    residuals = zeros(numNodes, forward_equations_nbr)
    y = zeros(numNodes, 3*endogenous_nbr)
    evalPt = zeros(getNumDimensions(grid), numNodes)
    for i in 1:numNodes
        @views y[i, state_variables] .= state
        # 1) Determine t+1  variables using preamble
        @views preamble_block.set_endogenous_variables!([], y[i, endogenous_nbr + 1: end], nodes[i,:], params, steadystate)
        @views y[i, system_variables] .= x
        # 2) Extract next period's state variables
        @views evalPt[:, i] .= y[i, endogenous_nbr .+ state_variables]
    end
    # 3) Determine relevant variables within the expectations operator
    X = evaluateBatch(grid, evalPt)
    
    resid = zeros(forward_equations_nbr)
    YY = []
    for i in 1:numNodes
        # setting future variables
        @views begin
            y[i, endogenous_nbr .+ system_variables] .= X[:, i]
            forward_block!(T, resid, y[i, :], zeros(exogenous_nbr), params, steadystate)
            residuals[i,:] .= resid
        end
    end
    return residuals
end

function expectation_equations!(expected_residuals, x, T, state, params, grid0, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables, preamble_block)
    # get numNodes x numForwardEquations
    Integrand = forward_looking_equation_residuals(x, T, state, params, grid0, nodes, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables, preamble_block)
    for i in axes(expected_residuals, 1)
        @views expected_residuals[i] = dot(Integrand[:, i], weights)
    end
    return expected_residuals
end

function sysOfEqs(policy, T, y, exogenous, state, grid, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, bmcps, preamble_block)
    @views begin
        y[state_variables] .= state
        y[system_variables] .= policy 
    end 
    res = zeros(length(policy))

    @views begin
        expectation_equations!(res[1:forward_equations_nbr], policy, T, state, params, grid, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables, preamble_block)
        other_block!(T, res[forward_equations_nbr + 1:end], y, exogenous, params, steadystate)
    end
    reorder!(res, bmcps)
    return res
end

function sysOfEqs_derivatives_update!(J, Jexpected, Jother, policy, y, state, grid, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, bmcps)
    @views begin
        y[state_variables] .= state
        y[system_variables] .= policy 
    end 
    res = zeros(length(policy))

    @views begin
        expectation_equations_derivatives!(Jexpected, policy, state, params, grid, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
        other_block_derivatives(Jother, y, zeros(exogenous_nbr), params, steadystate)
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
function ti_step(grid, X0, pol, gridZero, nPols, nodes, weights, exogenous, params, steadystate, forward_equations_nbr, 
                 endogenous_nbr, exogenous_nbr, state_variables, system_variables, y, bmcps, dynamicws, model,
                 predetermined_variables, forward_system_variables,
                 forward_expressions_eqs, other_expressions_eqs, ids, lb, ub, preamble_block, mcp, method, solver, ftol, show_trace)

    # Get the points that require function values
    states = Tasmanian.getNeededPoints(grid)
    # Get the number of points that require function values
    nstates = Tasmanian.getNumNeeded(grid)
    # Array for intermediate update step
    n = length(system_variables)
    Y = zeros(n, nstates)
    # TODO: check the case with non-autocorrelated shocks
    exogenous = zeros(exogenous_nbr)
    T = dynamicws.temporary_values
    J = zeros(n, n)
    fx = zeros(n)

    state = Vector{Float64}(undef, size(states, 1))
    x0 =  Vector{Float64}(undef, n)
    
    f(x) = sysOfEqs(x, T, y, exogenous, state, gridZero, nPols, nodes, weights,
                    params, steadystate, forward_equations_nbr, endogenous_nbr,
                    exogenous_nbr, state_variables, system_variables, bmcps, preamble_block)
    JA1!(J, x) = sysOfEqs_derivatives!(J, x, y, state, gridZero,
                                  nPols, exogenous, nodes, weights, params, steadystate,
                                  forward_equations_nbr, endogenous_nbr, exogenous_nbr,
                                  state_variables, system_variables, bmcps, dynamicws, model, 
                                  predetermined_variables, forward_system_variables,
                                  forward_expressions_eqs, other_expressions_eqs, ids, preamble_block)                   
    
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
        NonLinearSolver_solve!(Y, f2!, JA2!, J, X0, states, state, x0, params, method, ftol, show_trace )
    elseif solver == NLsolver
        df = OnceDifferentiable(f1!, JA1!, X0, fx, J)
        NLsolve_solve!(Y, df, X0, lb, ub, states, X0, ftol, show_trace)
    elseif solver == PATHSolver
        PATHsolver_solve!(Y, f1!, JA1!, X0, lb, ub, states, X0, ftol, show_trace)
    else
        error("Sparsegrids: unknown solver")
    end
    #=            
#            mcp_sg_core!(PathNLS(), polInt, grid, pol_guess, pol, gridZero, nPols, nodes, weights, exogenous, params, steadystate, forward_equations_nbr, 
            mcp_sg_core!(polInt, grid, pol_guess, pol, gridZero, nPols, nodes, weights, exogenous, params, steadystate, forward_equations_nbr, 
                      endogenous_nbr, exogenous_nbr, state_variables, system_variables, y, bmcps, dynamicws, model,
                      predetermined_variables, forward_system_variables,
                      forward_expressions_eqs, other_expressions_eqs, ids, aNumAdd, aPoints1, J, lb, ub, preamble_block)
        else
            f2!(x) = sysOfEqs(x, T,y, exogenous, state, gridZero, nPols, nodes,
                             weights, params, steadystate, forward_equations_nbr,
                             endogenous_nbr, exogenous_nbr, state_variables,
                             system_variables, bmcps, preamble_block)
#=
            JA2!(x) = sysOfEqs_derivatives!(J, x, y, exogenous, state, gridZero, nPols,
                                            nodes, weights, params, steadystate,
                                            forward_equations_nbr, endogenous_nbr, exogenous_nbr,
                                            state_variables, system_variables, bmcps, dynamicws,
                                            predetermined_variables, forward_expressions_eqs, other_expressions_eqs)
=#        
            for ii1 in 1:aNumAdd
                @views state = aPoints1[:, ii1]
                @views pol = pol_guess[:, ii1]
                res = Dynare.nlsolve(f2!, pol, method = :trust_region, show_trace = false)
                polInt[:, ii1] .= res.zero
            end
        end
=#

    # Add the new function values to grid1
    Tasmanian.loadNeededPoints!(grid, Y)

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
    Tasmanian.setSurplusRefinement!(grid, surplThreshold, output=dimRef, refinement_type=typeRefinement, scale_correction=scaleCorrMat')
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
    Tasmanian.makeLocalPolynomialGrid!(grid, order = gridOrder, sRule = gridRule)
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

function forward_looking_equation_derivatives!(J, y, x, exogenous, state, params, dynamicws, grid, nodes, weight, steadystate, 
                                               model, forward_equations_nbr, state_variables, endogenous_nbr, 
                                               exogenous_nbr, system_variables, forward_system_variables, forward_expressions_eqs,
                                               predetermined_variables, ids, preamble_block)
    iDimensions = getNumDimensions(grid)
    iOutputs = getNumOutputs(grid)
    evalPt = zeros(iDimensions)
    T = dynamicws.temporary_values
    @views begin
        y[state_variables] .= state
        preamble_block.set_endogenous_variables!(T, y[endogenous_nbr + 1: end], nodes, params, steadystate)
        y[system_variables] .= x
        evalPt .= y[endogenous_nbr .+ state_variables]
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
        params,
        y,
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

function forward_equations_derivatives(x, y, nodes, weights, dynamicws, params, exogenous, 
                                       steadystate, model, 
                                       forward_equations_nbr, state, grid, state_variables, endogenous_nbr, 
                                       exogenous_nbr, system_variables, forward_system_variables, forward_expressions_eqs, 
                                       predetermined_variables, ids, preamble_block)
    jacobian = get_dynamic_jacobian!(
        dynamicws,
        params,
        y,
        nodes[1, :],
        steadystate,
        model,
        2,
    )
    
    J = zeros(forward_equations_nbr, length(system_variables))
    for i in axes(nodes, 1)
        forward_looking_equation_derivatives!(J, y, x, exogenous, state, params, dynamicws, grid, nodes[i, :], weights[i], steadystate, 
                                              model, forward_equations_nbr, state_variables, endogenous_nbr, 
                                              exogenous_nbr, system_variables, forward_system_variables, forward_expressions_eqs,
                                              predetermined_variables, ids, preamble_block)
    end
    return J
end 
    
function sysOfEqs_derivatives!(J, x, y, state, grid,
                               nPols, exogenous, nodes, weights, params, steadystate,
                               forward_equations_nbr, endogenous_nbr, exogenous_nbr,
                               state_variables, system_variables, bmcps, dynamicws, model, 
                               predetermined_variables, forward_system_variables,
                               forward_expressions_eqs, other_expressions_eqs, ids, preamble_block)
    @views begin
        y[state_variables] .= state
        y[system_variables] .= x 
    end 
    
    Jfwrd = forward_equations_derivatives(x, y, nodes, weights, dynamicws, params, exogenous, 
                                  steadystate, model, 
                                  forward_equations_nbr, state, grid, state_variables, endogenous_nbr, 
                                  exogenous_nbr, system_variables, forward_system_variables, forward_expressions_eqs, 
                                  predetermined_variables, ids, preamble_block)
    copy_to_submatrix!(J, 1:forward_equations_nbr, 1:length(system_variables), Jfwrd)

    jacobian = get_dynamic_jacobian!(
        dynamicws,
        params,     
        y,
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

#=
function mcp_sg_core!(polInt, grid, pol_guess, pol, gridZero, nPols, nodes, weights, exogenous, params, steadystate, forward_equations_nbr, 
                      endogenous_nbr, exogenous_nbr, state_variables, system_variables, y, bmcps, dynamicws, model,
                      predetermined_variables, forward_system_variables,
                      forward_expressions_eqs, other_expressions_eqs, ids, aNumAdd, aPoints1, JJ, lb, ub, preamble_block)
    T = dynamicws.temporary_values
    let state
        function f1!(fx, x)
            fx .= Dynare.sysOfEqs(x, T, y, exogenous, state, gridZero, nPols, nodes, weights,
                                  params, steadystate, forward_equations_nbr, endogenous_nbr,
                                  exogenous_nbr, state_variables, system_variables, bmcps, preamble_block)
        end
        
        function JA1!(J, x)
            Dynare.sysOfEqs_derivatives!(J, x, y, state, gridZero,
                                         nPols, exogenous, nodes, weights, params, steadystate,
                                         forward_equations_nbr, endogenous_nbr, exogenous_nbr,
                                         state_variables, system_variables, bmcps, dynamicws, model, 
                                         predetermined_variables, forward_system_variables,
                                         forward_expressions_eqs, other_expressions_eqs, ids, preamble_block)                    
        end

        if solve_algo == PATHsolver
            for ii1 in 1:aNumAdd
                @views begin
                    state = aPoints1[:, ii1]
                    pol .= pol_guess[:, ii1]
                end
                (status, results, info) = mcp_solve!(PathNLS(), f1!, JA1!, JJ, lb, ub, pol, silent=true, convergence_tolerance=1e-4)
                if status != 1
                    @show status
                    @show info
                    error("PATHSolver failed")
                end
                @views polInt[:, ii1] .= results
            end
        elseif mcp_solve_algo == 
            fx = zeros(size(JJ, 1))
            df = OnceDifferentiable(f1!, JA1!, pol, fx, JJ)
       
    end
end
=#

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

function NonLinearSolver_solve!(Y, f, JA, J, X0, states, state, x0, params, method, ftol, show_trace )
    nlf = NonlinearFunction(f, jac = JA, jac_prototype=J)

    for i in axes(states, 2)
        @views begin 
            state .= states[:, i]
            x0 .= X0[:, i]
        end
        nlp = NonlinearProblem(nlf, x0, params)
        res = NonlinearSolve.solve(nlp, method, show_trace = Val(show_trace))
        @views Y[:, i] .= res.u
    end
end
