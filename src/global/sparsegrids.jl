include("blocks.jl")

using AxisArrays: AxisArray
using Distributions
using LinearAlgebra
using Roots
using Tasmanian

export sparsegridapproximation, simulate!

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
"""
function sparsegridapproximation(; context::Context=context,
                                 dimRef = -1,
                                 iterRefStart = 25,
                                 gridDepth = 2,
                                 gridOrder = 1,
                                 gridRule = "localp",
                                 maxiter = 300,
                                 maxRef = 1,
                                 maxRefLevel = gridDepth + maxRef,
                                 numstart = 0,
                                 savefreq = 10,
                                 scaleCorrInclude = [],
                                 scaleCorrExclude = [],
                                 surplThreshold= 1e-3,
                                 tol_ti = 1e-4,
                                 TT = 10000,
                                 typeRefinement = "classic",
                                 )
                                 
    model = context.models[1]
    work = context.work
    (state_variables, predetermined_variables, system_variables,
     forward_equations_nbr, other_equations_nbr,
     forward_expressions_eqs, system_expressions_eqs) = make_block_functions(context)
    lb = Vector{Float64}(undef, 0)
    ub = Vector{Float64}(undef, 0)
    bmcps = Vector{Vector{Int}}(undef, 0)
    block_eqs = union(forward_expressions_eqs, system_expressions_eqs)
    block_mcp!(lb, ub, bmcps, context, block_eqs,
               system_variables)
    @show lb
    @show ub
    @show bmcps
    
    equation_xref_list, variable_xred_list = xref_lists(context)

    params = context.work.params
    nstates = length(state_variables)

    model = context.models[1]

    gridDim = nstates
    gridOut = length(system_variables)
    limits = context.work.limits
    nl = length(limits)
    gridDomain = zeros(nl, 2)
    endogenous_nbr = model.endogenous_nbr
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    for i in 1:nl
        k = state_variables[i]
        k > endogenous_nbr && (k -= endogenous_nbr)
        vname = Symbol(endogenous_variables[k])
        gridDomain[i,1] = limits[vname].max
        gridDomain[i,2] = limits[vname].min
    end
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

    polGuess = guess_policy(context, aNum, nPols, aPoints, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr)
    loadNeededPoints!(grid0, polGuess)

    params = context.work.params

    (nodes, weights) = monomial_power(exogenous_nbr)

    # Scale correction in the refinement process:
    @views scaleCorr = make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables .- endogenous_nbr])

    polGuess1 = copy(polGuess)
    for iter0 in 1:maxiter

        polGuess1 = copy(polGuess)
        grid1 = fresh_grid(gridDim, gridOut, gridDepth, gridDomain, gridOrder, gridRule)

        # Index of current grid level to control the number of refinements
        ilev = gridDepth
        while ((getNumNeeded(grid1) > 0) && (ilev <=  maxRefLevel))
            grid1 = ti_step(grid1, polGuess1, grid0, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, endogenous, bmcps)
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

function make_guess_system(x, residuals, state, endogenous, exogenous, parameters, steadystate, state_variables, system_variables, forward_equations_nbr, other_equations_nbr, endogenous_nbr, permutations)
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

function sysOfEqs(policy, y, state, grid, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, bmcps)
    @views begin
        y[state_variables] .= state
        y[system_variables] .= policy 
    end 
    res = zeros(length(policy))

    @views begin
        expectation_equations!(res[1:forward_equations_nbr], policy, state, params, grid, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
        other_block(res[forward_equations_nbr + 1:end], y, zeros(exogenous_nbr), params, steadystate)
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
function ti_step(grid, pol_guess, gridZero, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, y, bmcps)

    # Get the points that require function values
    aPoints1 = Tasmanian.getNeededPoints(grid)
    # Get the number of points that require function values
    aNumAdd = Tasmanian.getNumNeeded(grid)
    # Array for intermediate update step
    polInt = zeros(length(system_variables), aNumAdd)
    let state
        fnl(x) = sysOfEqs(x, y, state, gridZero, nPols, nodes, weights, params, steadystate, forward_equations_nbr, endogenous_nbr, exogenous_nbr, state_variables, system_variables, bmcps)
        
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


function block_mcp!(lb, ub, bmcps, context, block_eqs, block_vars)
    n = length(block_eqs)
    resize!(lb, n)
    resize!(ub, n)
    fill!(lb, -Inf)
    fill!(ub, Inf)
    for m in context.models[1].mcps
        (eqn, var, op, expr) = m
        @show m
        @show block_vars  .- context.models[1].endogenous_nbr
        bvar = findfirst(var .== block_vars .- context.models[1].endogenous_nbr)
        @show bvar
        @show block_eqs
        beqn = findfirst(eqn .== block_eqs)
        @show beqn
        push!(bmcps, [beqn, bvar])
        @show bmcps
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
    
