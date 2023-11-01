module SparseGrids

include("blocks.jl")

#using Blocks
using Distributions
#using Dynare
using LinearAlgebra
using Roots
using Tasmanian

include("time_iteration.jl")

function set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain)
    #=
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

    # Set the grid domain to [kmin,kmax]^n x [amin,amax]^n
    gridDomain = zeros(gridDim,2)
    for i in 1:nstates
        
    
    sigE = 0.01
    rhoZ = 0.95
    # Lower bound for capital
    kMin = 0.8
    # Upper bound for capital
    kMax = 1.2
    # Lower bound for TFP
    aMin = -0.8*sigE/(1.0-rhoZ)
    # Upper bound for TFP
    aMax = 0.8*sigE/(1.0-rhoZ)

    i_k = 1:2
    i_a = 3:4
    gridDomain[i_k, 1] .= kMin
    gridDomain[i_k, 2] .= kMax
    gridDomain[i_a, 1] .= aMin
    gridDomain[i_a, 2] .= aMax
    =#
    
    # Generate the grid structure
    grid = Tasmanian.TasmanianSG(gridDim, gridOut, gridDepth)
    Tasmanian.makeLocalPolynomialGrid!(grid, iOrder = gridOrder, sRule = gridRule)
    # Transform the domain
    Tasmanian.setDomainTransform!(grid, gridDomain)
    # Get the points that require function values
    aPoints = Tasmanian.getPoints(grid)
    # Get the number of points that require function values
    aNum = Tasmanian.getNumPoints(grid)
    return grid, aPoints, aNum, gridDim, gridDepth, gridDomain
end

# TODO insure that guess is in variable domain
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

#=
function  ARC_zero(lamb, gridPt, delta, t, gamma, kappa, phi, A, F)
    
    res = 0.0
    for i1 in 1:2
        res += (F(gridPt[i1], gridPt[2 + i1]) - (delta*phi/2.0)^2 - (lamb/t[i1])^(-gamma[i1]))
    end
    
    return res
end

function policy_guess(aPoints, delta, t, gamma, kappa, phi, A, aNum, nPols, F)
    polGuess = zeros(aNum, nPols)
    @views polGuess[:, 1:2] .= (1 - delta) .* aPoints[:, 1:2] 
    for i in axes(aPoints, 1)
        f(x) = ARC_zero(x, aPoints[i, :], delta, t, gamma, kappa, phi, A, F)
        polGuess[i, 3] = find_zero(f, [0.1, 1.5], Bisection(), atol=1e-8)
    end
    return polGuess
end
=#

"""
Evaluates residual of forward looking equations for each integration node
"""
function ExpectFOC(x, state, params, grid, nodes, nCountries, steadystate)
    numNodes = size(nodes, 1)
    residuals = zeros(numNodes, nCountries)
    y = zeros(numNodes, 15)
    evalPt = zeros(numNodes, nCountries*2)
    for i in 1:numNodes
        @views y[i, i_state] .= state
        yy = copy(y[i,:])
        # 1) Determine t+1  tfp states using preamble between t-1 and 6
        @views preamble_block!(y, nodes[i,:], params, steadystate)
        y[i, :] .= yy
        @views y[i, [12, 14]] .= y[i, [7, 9]]
        @views y[i, [7, 9]] .= y[i, [2, 4]]
        # 2) Determine next period's state variables
        @views evalPt[i, 1:nCountries] .= x[1:nCountries]
        @views evalPt[i, nCountries + 1:end] .= y[i, [12, 14]]
    end
    # 3) Determine relevant variables within the expectations operator
    X = Tasmanian.evaluateBatch(grid, evalPt)
    
    capPrPr = X[:, 1:2]
    lambPr = X[:, 3]


    resid = zeros(2)
    for i in 1:numNodes
        #=
        if typeIRBC=='non-smooth':
        gzAlphaPr = grid.evaluateBatch(evalPt)[:,nCountries+1:]
        gzAplusPr = np.maximum(0.0,gzAlphaPr)
        
        
        # Compute tomorrow's marginal productivity of capital
        MPKtom = np.zeros((numNodes,nCountries))
        for impk in range(nCountries):
        MPKtom[:,impk] = 1.0 - delta + Fk(ktemp[impk],newstate[:,impk]) - AdjCost_k(ktemp[impk],capPrPr[:,impk])
        
        
        #Compute Density
        if typeInt=='GH-quadrature':
        density = np.pi**(-(nCountries+1) * 0.5)
        else:
        density = 1.0
        =#
        #Specify Integrand
        y[i, [6,8]] .= x[1:2]
        y[i, 10] = x[3]
        #yp1[:, 2:2:4] .= newstate
        y[i, [11,13]] .= capPrPr[i, :]
        y[i, 15] = lambPr[i]
        # residuals numNodes x nbr forward equations
        forwardblock!(resid, y[i, :], zeros(2), params)
        residuals[i,:] .= resid
        #=
        if typeIRBC=='non-smooth':
        
        for iexp in range(nCountries):
        ExpectFOC[:,iexp] = (MPKtom[:,iexp]*lambPr - (1.0-delta)*gzAplusPr[:,iexp]) * density
        
        else:
        
        for iexp in range(nCountries):
        ExpectFOC[:,iexp] = MPKtom[:,iexp]*lambPr * density
        
        =#
    end
    return residuals
end

#=
function forward_looking_equations(y, yp1, ym1, params)
    Y = zeros(21)
    @views Y[1:4] .= ym1
    @views Y[6:10] .= y
    x = zeros(2)
    residuals = zeros(size(yp1,1), 2)
    @views for i in axes(yp1, 1)
        Y[11:15] .= yp1[i, :]
        forwardblock!(residuals[i, :], Y, x, params)
    end
    return residuals
end
=#

function monomial_power(nShocks, sigE)
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

function expectation_equations!(expected_residuals, x, state, params, grid0, nodes, weights, nCountries)
    # get numNodes x numForwardEquations
    Integrand = ExpectFOC(x, state, params, grid0, nodes, nCountries)
    for i in axes(expected_residuals, 1)
        expected_residuals[i] = dot(Integrand[:, i], weights)
    end
    return expected_residuals
end

function sysOfEqs(policy, state, grid, F, AdjCost, nCountries, nPols, nodes, weights, params)
    # State variables
    capStates = state[[1, 3]]
    tfpStates = state[[2, 4]]

    # Policy values
    capPolicies = x[1:nCountries]
    lamb = x[nCountries + 1]

    res = zeros(nPols)

    expectation_equations!(res, policy, state, params, grid, nodes, weights, nCountries)
    # Aggregate resource constraint
    gamma = params[1:3:6]
    t = params[3:3:6]
    kappa, beta, delta, phi, rho, A, sigE = params[7:13]
    for ires2 in 1:nCountries
        res[nCountries + 1] += (F(capStates[ires2],tfpStates[ires2]) + (1.0 - delta)*capStates[ires2] - capPolicies[ires2]
                                - AdjCost(capStates[ires2],capPolicies[ires2]) - (lamb/t[ires2])^(-gamma[ires2]))
#        @show ((lamb/t[ires2])^(-gamma[ires2]))
    end
#    @show (res[nCountries + 1])
    return res
end

"""
    #=
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
                                 iterrRefStart = 25,
                                 gridDepth = 2,
                                 gridOrder = 1,
                                 gridRule = "localp",
                                 maxiter = 300,
                                 maxRef = 1,
                                 maxRefLevel = gridDepth + maxRef,
                                 numstart = 0,
                                 savefreq = 10,
                                 scaleCorr = nothing,
                                 surplThreshold= 1e-3,
                                 tol_ti = 1e-4,
                                 TT = 10000,
                                 typeRefinement = "classic",
                                 )
                                 
    model = context.models[1]
    states, predetermined_variables, system_variables = make_block_functions(context)
    equation_xref_list, variable_xred_list = xref_lists(context)

    params = context.work.params
    nstates = length(states)

    model = context.models[1]

    gridDim = nstates
    gridOut = length(system_variables)
    limits = context.work.limits
    nl = length(limits)
    gridDomain = zeros(nl, 2)
    endogenous_nbr = model.endogenous_nbr
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    for i in 1:nl
        k = states[i]
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

    # Scale correction in the refinement process:
    isnothing(scaleCorr) && (scaleCorr = zeros(nPols))
    #=
    # We only let the capital policies of each country determine the addition of grid points
    scaleCorr = zeros(nPols)
    scaleCorr[1:nCountries] .= 1
    =#
    
    #=
    gamma = params[1:3:6]
    t = params[3:3:6]
    

#    polGuess = policy_guess(context, aNum, nPols, aPoints, [1, 3, 5], model.i_bkwrd_b)    
    polGuess = policy_guess(aPoints, delta, t, gamma, kappa, phi, A, aNum, nPols, F)
    Tasmanian.loadNeededPoints!(grid0, polGuess)

#    steadystate = context.results.model_results[1].trends.endogenous_steady_state

    state = zeros(4)
#    state[1:2:4] = steadystate[1:2:4]
#    state[2:2:4] = steadystate[2:2:4]
#    ktemp = steadystate[2:2:4]
    state = [1, 0, 1, 0]
    ktemp = [1, 1]
    
    params = context.work.params
    nshocks = 3

    (nodes, weights) = monomial_power(nshocks, sigE)

    iter0 = 0
    iterRefStart = 25
    #polGuess1 = copy(polGuess)
    gridOrder = 1
    gridRule = "localp"
    gridDim = length(state)
    gridOut = nPols

    maxiter = 300
    
    for iter0 in 1:maxiter

        polGuess1 = copy(polGuess)
        grid1 = fresh_grid(gridDim, gridOut, gridDepth, gridDomain, gridOrder, gridRule)

        # Index of current grid level to control the number of refinements
        ilev = gridDepth
        
        while ((Tasmanian.getNumNeeded(grid1) > 0) && (ilev <=  maxRefLevel))
        
            grid1 = ti_step(grid1, polGuess1, grid0, F, AdjCost, nPols, nCountries, nodes, weights, params)
        
            # We start the refinement process after a given number of iterations
            if (iter0 - 1 > iterRefStart)
                grid1, polGuess1 = refine(grid1, nPols, scaleCorr, surplThreshold, dimRef, typeRefinement)
                #grid1 = time_iter.refine(grid1)
            end
            
            # Track the grid level
            ilev += 1
        end
        
        ## Calculate (approximate) errors on tomorrow's policy grid
        metric, polGuess, grid0 = policy_update(grid0, grid1, nPols)
 
        println("Iteration: $iter0, Grid pts: $(Tasmanian.getNumPoints(grid0)), Level: $ilev, Metric: $metric")

        if (mod(iter0 + 1, savefreq) == 0)
            #            save_grid(grid0,iter0)
        end
        
        if (metric < tol_ti)
            break
        end
    end
    =#
end

function test_sparsegrids()
    context = @dynare "irbc_small"  "notmpterms" "stoponerror" "savemacro";
    sparsegridapproximation(context=context)
end


end # module SparseGrids

SparseGrids.test_sparsegrids()
