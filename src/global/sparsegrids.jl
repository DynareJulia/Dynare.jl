include("blocks.jl")
include("block_firstorder.jl")
include("monomial.jl")

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
    drawsnbr = 10000
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
                                 solver = mcp ? NLsolver : NonlinearSolver,
                                 surplThreshold= 1e-3,
                                 tol_ti = 1e-4,
                                 drawsnbr = 10000,
                                 typeRefinement = "classic",
                                 )
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    model = context.models[1]
    endogenous_nbr = model.endogenous_nbr
    exogenous_nbr = model.exogenous_nbr
    work = context.work
    (dynamic_state_variables, predetermined_variables, system_variables,
     backward_block, forward_block, preamble_block, system_block) = make_block_functions(context)
    lb = Vector{Float64}(undef, 0)
    ub = Vector{Float64}(undef, 0)
    bmcps = Vector{Vector{Int}}(undef, 0)
    # must be consistent with SysOfEqs()
    block_eqs = union(forward_block.equations, backward_block.equations)
    # computes lb, ub, bmcps
    limits = context.work.limits
    block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables, limits)
    equation_xref_list, variable_xref_list = xref_lists(context)
    params = context.work.params
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    nstates = length(dynamic_state_variables)

    sgmodel = SGModel(repeat(steadystate, 3),
                      zeros(2*endogenous_nbr),
                      zeros(exogenous_nbr),
                      endogenous_nbr,
                      exogenous_nbr,
                      params,
                      steadystate,
                      []
                      )

    ids = SGIndices(
        endogenous_nbr,
        model.i_fwrd_b,
        dynamic_state_variables,
        system_variables
    )
    forward_system_variables = model.i_fwrd_b
    setdiff!(forward_system_variables, predetermined_variables .- endogenous_nbr)

    gridDim = nstates
    gridOut = length(system_variables)
    nl = length(limits)
    gridDomain = zeros(gridDim, 2)
    for l in limits
        i = findfirst(x -> x == context.symboltable[string(l[1])].orderintype, ids.state_variables) 
        if !isnothing(i)
            gridDomain[i,1] = l[2].min
            gridDomain[i,2] = l[2].max
        end
    end

    grid0, aPoints, aNum, nPols =
        set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain);

    ################################################################################
    #                        Adaptivity parameters                                 #
    ################################################################################

    grid = copyGrid(grid0)

    monomial = MonomialPowerIntegration(exogenous_nbr)
    nodesnbr = length(monomial.nodes)
 
    # Scale correction in the refinement process:
    @views scaleCorr = make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables])

    n = length(system_variables)
    backward_variable_jacobian = Matrix{Float64}(undef, length(backward_block.equations), length(system_variables)) 
    forward_equations_nbr = length(forward_block.equations)
    sgoptions = SGOptions(mcp, method, solver, ftol, show_trace)
    residuals = Vector{Float64}(undef, n)
    forward_jacobian = Matrix{Float64}(undef, forward_equations_nbr, n)
    forward_points = Matrix{Float64}(undef, getNumDimensions(grid), nodesnbr)
    forward_residuals = Matrix{Float64}(undef, length(forward_block.equations), nodesnbr)
    forward_variable_jacobian = Matrix{Float64}(undef, length(forward_block.equations), length(ids.forward_in_system))
    J = zeros(n, n)
    fx = zeros(n)
    evalPt = Vector{Float64}(undef, nstates)
    future_policy_jacobian = zeros(length(ids.state_in_system), length(ids.forward_in_system))
    policy_jacobian = Matrix{Float64}(undef, gridDim, n)
    dyn_endogenous_vector = [zeros(3*endogenous_nbr) for i in 1:nodesnbr]
    shift_dyn_endogenous_vector = zeros(2*endogenous_nbr, nodesnbr)
    forward_equations_nbr = length(forward_block.equations)
    M1 = Matrix{Float64}(undef, forward_equations_nbr, length(ids.state_in_system))
    policyguess = zeros(gridOut, nodesnbr)
    tmp_state_variables = Vector{Float64}(undef, length(dynamic_state_variables))
    sev = Sev(zeros(3*endogenous_nbr), zeros(nstates))
    sgws_ = SparsegridsWs(monomial, dynamic_state_variables, system_variables, 
                         bmcps, ids, lb, ub, backward_block, backward_variable_jacobian, 
                         forward_block, preamble_block, system_block,
                         residuals, forward_jacobian, forward_points, 
                         forward_residuals, forward_variable_jacobian, J, fx,
                         policy_jacobian, evalPt, 
                         future_policy_jacobian, 
                         dyn_endogenous_vector, M1, policyguess, tmp_state_variables,
                         sev, sgmodel, sgoptions)
    if Threads.nthreads() > 1
        BLAS.set_num_threads(1)
        sgws = [deepcopy(sgws_) for i in 1:Threads.nthreads()]
    else
        sgws = [sgws_]
    end
    
    M, N, P = block_first_order_approximation(context, sgws[1])
    polGuess = guess_policy(context, aNum, nPols, aPoints, sgws[1], M, N, P)
#    polGuess = guess_policy(context, aNum, nPols, aPoints, dynamic_state_variables, system_variables, lb, ub, backward_block, forward_block, sgmodel)
    loadNeededPoints!(grid, polGuess)

    scaleCorrMat = repeat(Float64.(scaleCorr), 1, Tasmanian.getNumLoaded(grid))

    total_time = 0
    context.timings["sparsegrids"] = Vector{Millisecond}(undef, 0)
    timings = context.timings["sparsegrids"]
    iter = 1
    while iter <= maxiter
        tin = now()
        # new policy guess to be computed
        polGuess1 = copy(polGuess)
        newgrid = copyGrid(grid0)
        # Index of current grid level to control the number of refinements
        ilev = gridDepth
        while (getNumNeeded(newgrid) > 0)
            if (ilev <=  maxRefLevel)
                map(s -> fill!(s.J, 0.0), sgws)
                ti_step!(newgrid, grid, polGuess1, sgws)
                # We start the refinement process after a given number of iterations
                if iter >= iterRefStart && ilev < maxRefLevel
                    polGuess1 = refine!(newgrid, polGuess1, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)
                end
            else
                loadNeededPoints!(newgrid, polGuess1)
            end
            # Track the grid level
            ilev += 1
        end
        # Calculate (approximate) errors on tomorrow's policy grid
        metric, polGuess, grid = policy_update(grid, newgrid, polGuess, polGuess1, length(system_variables))
        iteration_walltime = now() - tin
        iter > 1 && (total_time += iteration_walltime.value)
        push!(timings, iteration_walltime)
        println("Iteration: $iter, Grid pts: $(getNumPoints(grid)), Level: $ilev, Metric: $metric, Computing time: $(now() -tin)")
        if (mod(iter + 1, savefreq) == 0)
            #            save_grid(grid0,iter0)
        end

        if (metric < tol_ti)
            break
        end
        iter += 1
    end
    average_time = total_time/(iter - 1)
    println("Last grid points: $(getNumPoints(grid))")
    println("Average iteration time (except first one): $average_time")

    results = context.results.model_results[1].sparsegrids
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

function simulate!(; context = context, grid::TasmanianSG = grid, periods = 1000, replications=1, sgws = sgws)
    model = context.models[1]
    @unpack endogenous_nbr, exogenous_nbr, Sigma_e = model
    @unpack dynamic_state_variables, preamble_block, system_variables = sgws
    @unpack steadystate = sgws.sgmodel
    params = context.work.params
    results = context.results.model_results[1]
    trends = results.trends
    perfect_foresight_ws = PerfectForesightWs(context, periods)
    shocks = perfect_foresight_ws.shocks
    initial_values = get_dynamic_initialvalues(context)
    chol_sigma_e_L = cholesky(Sigma_e).L
    x = reshape(perfect_foresight_ws.x, exogenous_nbr, periods)

    ws = DynamicWs(context)
    T = ws.temporary_values
    Y = Array{Float64}(undef, periods, endogenous_nbr + exogenous_nbr, replications)
    y = Vector{Float64}(undef, 3*endogenous_nbr)
    z = Vector{Float64}(undef, endogenous_nbr)
    random_shocks = Matrix{Float64}(undef, exogenous_nbr, periods)
    sv_buffer = Vector{Float64}(undef, length(dynamic_state_variables))
    pv_buffer = Vector{Float64}(undef, length(system_variables))
    @inbounds for r in 1:replications
        mul!(random_shocks, chol_sigma_e_L, randn(exogenous_nbr, periods))
        if !isempty(shocks)
            random_shocks .+= x
        end 
        y[1:endogenous_nbr] .= initial_values
        for p in 1:periods
            @views begin
                preamble_block.set_endogenous_variables!(T, y, random_shocks[:, p], params, steadystate)
                sv_buffer .= y[dynamic_state_variables]
                evaluateBatch!(pv_buffer, grid, sv_buffer)
                circshift!(y, -endogenous_nbr)
                y[system_variables] .= pv_buffer
                Y[p, 1:endogenous_nbr, r] .= y[1:endogenous_nbr]
                Y[p, endogenous_nbr .+ (1:exogenous_nbr), r] .= random_shocks[:, p]
            end
        end
    end
    symboltable = context.symboltable
    varnames = vcat(Symbol.(get_endogenous(symboltable)),
                    Symbol.(get_exogenous(symboltable)))
    return AxisArray(Y, Undated(1):Undated(periods), varnames, 1:replications) 
end

function simulation_approximation_error(; context = context, grid::TasmanianSG = grid, drawsnbr = 10000, quantile_probability = 0.999, sgws = sgws)
    @unpack dynamic_state_variables, system_block, system_variables = sgws
    @unpack endogenous_nbr = sgws.sgmodel

    # the first simulation period is used for initialization
    Y = simulate!(context = context, grid = grid, periods = drawsnbr + 1, replications = 1, sgws = sgws)
    errors = zeros(length(system_variables), drawsnbr)
    ir1 = filter(x -> x .< endogenous_nbr, dynamic_state_variables)
    n = length(dynamic_state_variables)
    n1 = length(ir1)
    ir2 = filter(x -> x .> endogenous_nbr, dynamic_state_variables)
    x = zeros(length(dynamic_state_variables))
    @views begin
        for i in 1:size(Y, 1) 
            # initial values
            i == 1 && continue
            x[1:n1] .= Y[i-1, ir1]
            x[n1 .+ 1:n] .= Y[i, ir2 .- endogenous_nbr]
            y = Y[i, system_variables]
            sysOfEqs!(errors[:, i - 1], y, x, grid, sgws)
        end
    end
    equations = system_block.equations
    println("\nMaximum approximation error by equation:\n")
    equation_quantile_errors = zeros(length(equations))
    equation_average_errors = zeros(length(equations))
    @views begin
        for i in axes(errors, 1)
            equation_quantile_errors[i] = quantile(abs.(errors[i, :]), quantile_probability)
            println("Equation $(equations[i]): absolute error $(100*quantile_probability) quantile: $(quantile(abs.(errors[i, :]), quantile_probability))")
        end
    end
    println("\nMean approximation error by equation:\n")
    @views begin
        for i in axes(errors, 1)
            equation_average_errors[i] = mean(abs.(errors[i, :]))
            println("Equation $(equations[i]): mean absolute error: $(equation_average_errors[i])")
        end
    end
    average_error = mean(abs.(errors))
    quantile_error = quantile(abs.(vec(errors)), quantile_probability)
    println("\nOverall absolute error $(100*quantile_probability)% quantile: $(quantile_error)" ) 
    println("Overall mean absolute error: $(average_error)\n" ) 
    results = context.results.model_results[1].sparsegrids
    results.average_error = average_error
    results.drawsnbr = drawsnbr
    results.equation_average_errors = equation_average_errors
    results.equation_quantile_errors = equation_quantile_errors
    results.quantile_probability = quantile_probability
    results.quantile_error = quantile_error
    return errors
end


function plot_policy_function(plot_variables, vstate, state, grid, state_variables, system_variables; N=50, context=context)
    sv_buffer = repeat(state, 1, N)
    k = findfirst(state_variables .== context.symboltable[vstate].orderintype)
    kp = [findfirst(context.symboltable[v].orderintype .== system_variables) for v in plot_variables]
    gridDomain = getDomainTransform(grid)
    x = range(gridDomain[k,1], gridDomain[k, 2], N)
    sv_buffer[k, :] .= x
    Y = evaluateBatch(grid, sv_buffer)
    return x, Y
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

function guess_policy(context, aNum, nPols, aPoints, sgws, M, N, P)
    @unpack ids, dynamic_state_variables, system_variables, preamble_block, lb, ub, backward_block, forward_block, sgmodel = sgws
    state_variables = preamble_block.variables

    results = context.results.model_results[1]
    model = context.models[1]
    
    steadystate = results.trends.endogenous_steady_state
    ssy = steadystate[system_variables]
    ssi = steadystate[filter(x -> x <= model.endogenous_nbr, dynamic_state_variables)]
    sst = steadystate[state_variables]

    nsv = length(system_variables)
    y0 = zeros(nsv)
    polGuess = zeros(nsv, aNum)
    istate = findall(dynamic_state_variables .> sgmodel.endogenous_nbr)
    ilagendo = findall(dynamic_state_variables .<= sgmodel.endogenous_nbr)
    iy = findall(in(dynamic_state_variables[ilagendo]), system_variables) 
    @views for (i, c) in enumerate(eachcol(aPoints))
        x = c[istate] .- sst
        y0[iy] .= c[ilagendo] .- ssi
        polGuess[:, i] .= ssy .+ M*y0 .+ N*x
    end
    set_boundaries!(polGuess, lb, ub)
    return polGuess
end

function set_boundaries!(y, lb, ub)
    for c in eachcol(y)
        for (i, x) in enumerate(c)
            if x <= lb[i]
                x = lb[i] + eps()
            elseif x >= ub[i]
                x = ub[i] - eps()
            end
        end
    end
end
            


"""
    given the state variables, find a solution that is constant in t and t+1 in absence of future shocks
"""

function make_guess_system1(residuals, T, state, dynamic_state_variables, system_variables, backward_block, forward_block, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    dyn_endogenous[dynamic_state_variables] .= state
    backward_equations_nbr = length(backward_block.equations)
    forward_equations_nbr = length(forward_block.equations)
    res1 = zeros(backward_equations_nbr)
    res2 = zeros(forward_equations_nbr)
    function f(residuals, x, params)
        dyn_endogenous[system_variables .+ endogenous_nbr] .= x
        dyn_endogenous[system_variables .+ 2*endogenous_nbr] .= x
        DFunctions.SparseDynamicResidTT!(T, dyn_endogenous, exogenous, parameters, steadystate)
        other_block_(T, res1, dyn_endogenous, exogenous, parameters, steadystate)
        forward_block_(T, res2, dyn_endogenous, exogenous, parameters, steadystate)
        residuals[1:backward_equations_nbr] .= res1
        residuals[backward_equations_nbr .+ (1:forward_equations_nbr)] .= res2
        return residuals
    end
    return f
end

function guess_policy(context, aNum, nPols, aPoints, dynamic_state_variables, system_variables, lb, ub, backward_block, forward_block, sgmodel)
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
    n =  length(system_variables)
    guess_values = zeros(n, aNum)
    @views x0 = copy(steadystate[system_variables])
    residuals = similar(x0)
    ws = DynamicWs(context)
    T = []
    for i in axes(aPoints, 2)
        @views state = aPoints[:, i]
        f1 = make_guess_system1(residuals, T, state, dynamic_state_variables, system_variables,
                                backward_block, forward_block, sgmodel)
        problem = NonlinearProblem(f1, x0, parameters)
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
    ti_step!(newgrid, oldgrid, X, sgws)
updates grid points    
"""
function ti_step!(newgrid, oldgrid, X, sgws)
    @unpack system_variables, lb, ub, sgoptions, J, fx,
            sgmodel = sgws[1]
    @unpack mcp, method, solver, ftol, show_trace = sgoptions
    @unpack exogenous, parameters = sgmodel
    
    # Get the points that require function values
    states = getNeededPoints(newgrid)
    nstates = size(states, 2)
    
    # Get the number of points that require function values
    # Array for intermediate update step
    n = length(system_variables)
    # TODO: check the case with non-autocorrelated shocks
    fill!(exogenous, 0.0)

    state = view(states, :, 1)

    # Default solver
    if isnothing(solver)
        solver = mcp ? NLsolver : NonlinearSolver 
    end

    # Solving nonlinear problem
    SG_NLsolve!(X, lb, ub, fx, J, states, parameters, oldgrid, sgws, solver, method, ftol, show_trace )
    
    # Add the new function values to newgrid
    loadNeededPoints!(newgrid, view(X, :, 1:nstates))
end

""" 
    Grid refinement
"""
function refine!(grid, polguess, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)

    # Refine the grid based on the surplus coefficients

    if size(scaleCorrMat, 2) < getNumLoaded(grid)
        scaleCorrMat = repeat(Float64.(scaleCorr), 1, getNumLoaded(grid))
    end
    setSurplusRefinement!(grid, surplThreshold, output=dimRef, refinement_type=typeRefinement, scale_correction=scaleCorrMat')
    
    if getNumNeeded(grid) > 0
	# Get the new points and the number of points
	nwpts = getNeededPoints(grid)
	aNumNew = getNumNeeded(grid)
        
	# We assign (for now) function values through interpolation#
	#polguess = zeros(aNumNew, gridOut)
        if size(polguess, 1) == aNumNew
            evaluateBatch!(polguess, grid, nwpts)
        else           
            polguess = evaluateBatch(grid, nwpts)
        end
    end
    return polguess
end

"""
   Checking convergence and updating policies
"""
function policy_update(gridOld, gridNew, polGuessOld, polGuessNew, gridOut)

    # Get the points and the number of points from grid1
    aPoints2 = getPoints(gridNew)
    aNumTot =  getNumPoints(gridNew)

    # Evaluate the grid points on both grid structures
    if size(polGuessNew, 2) == aNumTot
        evaluateBatch!(polGuessNew, gridNew, Matrix(aPoints2))
    else    
        polGuessNew = evaluateBatch(gridNew, Matrix(aPoints2))
    end 
    if size(polGuessOld, 2) == aNumTot
        evaluateBatch!(polGuessOld, gridOld, Matrix(aPoints2))
    else    
        polGuessOld = evaluateBatch(gridOld, Matrix(aPoints2))
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

    # Now update pol_guess and grid

    for iupd in 1:gridOut
        @views polGuessNew[iupd, :] .= 0.5*(polGuessOld[iupd, :] .+ polGuessNew[iupd, :])
    end

    return metric, polGuessNew, copyGrid(gridNew)
end

"""
   Grid storage
"""
function save_grid(grid, iter)

    grid.write(data_location_smooth + "grid_iter_$i.txt")

    return
end

function block_mcp!(lb, ub, bmcps, context, block_eqs, block_vars, limits)
    block_eqs_ = copy(block_eqs)
    n = length(block_eqs_)
    resize!(lb, n)
    resize!(ub, n)
    fill!(lb, -Inf)
    fill!(ub, Inf)
    for m in context.models[1].mcps
        (eqn, var, op, expr) = m
        bvar = findfirst(var .== block_vars)
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

function sysOfEqs!(residuals, policy, state, grid, sgws)
    @unpack dynamic_state_variables, system_variables, bmcps, backward_block, forward_block,
            dyn_endogenous_vector, sgmodel = sgws
    @unpack dyn_endogenous, endogenous_nbr, exogenous, parameters, steadystate, tempterms = sgmodel
    @views begin
        dyn_endogenous[dynamic_state_variables] .= state
        dyn_endogenous[system_variables .+ endogenous_nbr] .= policy 
    end
    nfwrd = length(forward_block.variables)
#    set_dyn_endogenous_vector!(dyn_endogenous_vector, policy, state, grid, sgws)
    @views begin
        expectation_equations!(residuals[1:nfwrd], policy, state, grid, sgws)
        backward_block.get_residuals!(tempterms, residuals[nfwrd  + 1:end], dyn_endogenous, exogenous, parameters, steadystate)
    end
    reorder!(residuals, bmcps)
    return nothing
end

function sysOfEqs_derivatives!(J, policy, state, grid, sgws)
    @unpack dynamic_state_variables, system_variables, bmcps,
    backward_block, forward_block, sgmodel = sgws
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate, tempterms = sgmodel
    forward_equations_nbr = length(forward_block.equations)
    @views begin
        dyn_endogenous[dynamic_state_variables] .= state
        dyn_endogenous[system_variables .+ endogenous_nbr] .= policy 
    end 

    Jfwrd = forward_looking_equation_derivatives!(policy, state, grid, sgws)
    copy_to_submatrix!(J, 1:forward_equations_nbr, 1:length(system_variables), Jfwrd)

    backward_block.update_jacobian!(tempterms, forward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
    @views J[length(forward_equations_nbr) .+ (1:length(backward_block.equations)), :] .=
        backward_block.jacobian[:, system_variables .+ endogenous_nbr]
    reorder_rows!(J, bmcps)
end

function sysOfEqs_derivatives_update!(J, policy, state, grid, sgws)
    @unpack dynamic_state_variables, system_variables, bmcps,
        backward_block, forward_block, sgmodel = sgws
    @unpack dyn_endogenous, endogenous_nbr, exogenous, parameters, steadystate, tempterms = sgmodel
    @views begin
        dyn_endogenous[dynamic_state_variables] .= state
        dyn_endogenous[system_variables .+ endogenous_nbr] .= policy 
    end 

    fill!(exogenous, 0.0)
    nfwrd = length(forward_block.equations)
    set_dyn_endogenous_vector!(sgws.dyn_endogenous_vector, policy, state, grid, sgws)
    @views begin
        expectation_equations_derivatives(J[1:nfwrd, :], policy, state, grid, sgws)
        backward_block.update_jacobian!(tempterms, backward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
        J[nfwrd .+ (1:length(backward_block.equations)), :] .= backward_block.jacobian[:, system_variables .+ endogenous_nbr]
    end
    reorder_rows!(J, bmcps)
    return nothing
end

function sysOfEqs_with_jacobian!(residuals, J, policy, state, grid, sgws)
    @unpack dynamic_state_variables, dyn_endogenous_vector, system_variables, bmcps, backward_block, forward_block,
            sgmodel = sgws
    @unpack dyn_endogenous, endogenous_nbr, exogenous, parameters, steadystate, tempterms = sgmodel
    nfwrd = length(forward_block.variables)
    for (i, k) in enumerate(dynamic_state_variables)
        dyn_endogenous[k] = state[i]
    end
    for (i, k) in enumerate(system_variables)
        dyn_endogenous[k + endogenous_nbr] = policy[i]
    end
    
    @views begin
        expectation_equations!(residuals[1:nfwrd], policy, state, grid, sgws)
        backward_block.get_residuals!(tempterms, residuals[nfwrd  + 1:end], dyn_endogenous, exogenous, parameters, steadystate)
    end
    reorder!(residuals, bmcps)

    fill!(exogenous, 0.0)
    set_dyn_endogenous_vector!(dyn_endogenous_vector, policy, state, grid, sgws)
    @views begin
        expectation_equations_derivatives(J[1:nfwrd, :], policy, state, grid, sgws)
        backward_block.update_jacobian!(tempterms, backward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
        J[nfwrd .+ (1:length(backward_block.equations)), :] .= Matrix(backward_block.jacobian[:, system_variables .+ endogenous_nbr])
    end
    reorder_rows!(J, bmcps)
end

function expectation_equations!(expected_residuals, policy, state, grid, sgws)
    @unpack monomial, forward_residuals = sgws
    forward_looking_equation_residuals!(forward_residuals, policy, state, grid, sgws)
    mul!(expected_residuals, forward_residuals, monomial.weights)
    return expected_residuals
end

function expectation_equations_derivatives(forward_jacobian, policy, state, grid, sgws)
    @unpack monomial, dynamic_state_variables, system_variables, ids, dyn_endogenous_vector,
        backward_block, forward_block, preamble_block, sgmodel = sgws
    @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel

    fill!(forward_jacobian, 0.0)
    for i in 1:length(monomial.nodes)
        forward_looking_equation_derivatives!(forward_jacobian,
                                              grid, dyn_endogenous_vector[i],
                                              monomial.weights[i],
                                              sgws)
    end
    nothing
end 
    
"""
    Evaluates residual of forward looking equations for each integration node
"""
function set_dyn_endogenous_vector!(dyn_endogenous_vector, policy, state, grid, sgws)
    @unpack monomial, dynamic_state_variables, system_variables, forward_block, 
            preamble_block, sgmodel, forward_points, forward_residuals, policyguess, sgmodel = sgws
    @unpack shift_dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate, tempterms = sgmodel
    nodes = monomial.nodes
    nodesnbr = length(nodes)
    
    for i in 1:nodesnbr
        dyn_endogenous_vector_i = dyn_endogenous_vector[i]
        fill!(dyn_endogenous_vector_i, 0.0)
        for (j, k)  in enumerate(dynamic_state_variables)
            dyn_endogenous_vector_i[k] = state[j]
        end
        copyto!(shift_dyn_endogenous, 1, dyn_endogenous_vector_i, endogenous_nbr + 1, 2*endogenous_nbr)
        # 1) Determine t+1  variables using preamble
        preamble_block.set_endogenous_variables!(tempterms, shift_dyn_endogenous, nodes[i], parameters, steadystate)
        copyto!(dyn_endogenous_vector_i, endogenous_nbr + 1, shift_dyn_endogenous, 1, 2*endogenous_nbr)
        for (j, sv)  in enumerate(system_variables)
            dyn_endogenous_vector_i[sv + endogenous_nbr] = policy[j]
        end
        for (j, dsv) in enumerate(dynamic_state_variables)
            # 2) Extract next period's state variables
            forward_points[j, i] = dyn_endogenous_vector_i[dsv + endogenous_nbr]
        end
    end
    # 3) Determine relevant variables within the expectations operator
    evaluateBatch!(policyguess, grid, forward_points)
    for (i, k) in enumerate(system_variables), j in 1:nodesnbr
        dyn_endogenous_vector[j][2*endogenous_nbr + k] = policyguess[i, j]
    end
    nothing
end

function forward_looking_equation_residuals!(forward_residuals, policy, state, grid, sgws)
    set_dyn_endogenous_vector!(sgws.dyn_endogenous_vector, policy, state, grid, sgws)
    
    @unpack monomial, system_variables, forward_block, dyn_endogenous_vector, sgmodel = sgws
    @unpack dyn_endogenous, exogenous, parameters, steadystate, tempterms = sgmodel
    
    fill!(exogenous, 0.0)
    fill!(forward_residuals, 0.0)
    n = length(monomial.nodes)
    for i in 1:n
        # setting future variables
        @views forward_block.get_residuals!(tempterms, forward_residuals[:, i], dyn_endogenous_vector[i], exogenous, parameters, steadystate)
    end
    nothing
end

function forward_looking_equation_derivatives!(J, grid, dyn_endogenous, weight, sgws)
    @unpack dynamic_state_variables, tmp_state_variables, system_variables, ids,
      backward_block, forward_block, forward_variable_jacobian, preamble_block, sgmodel,
      policy_jacobian, future_policy_jacobian, evalPt, M1 = sgws
    @unpack exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate, tempterms = sgmodel

    # policy function Jacobian
    for (i, k) in enumerate(dynamic_state_variables)
        tmp_state_variables[i] = dyn_endogenous[k + endogenous_nbr]
    end
    differentiate!(policy_jacobian, grid, tmp_state_variables)
    # future policy Jacobian = policy Jacobian x derivatives of predetermined variables w.r. state
    copy_from_submatrix!(future_policy_jacobian, policy_jacobian, ids.state_in_system, ids.forward_in_system)
    forward_block.update_jacobian!(tempterms, forward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
    for (j, k) in enumerate(ids.forward_in_system)
        k1 = 2*endogenous_nbr .+ system_variables[k]
        for i in 1:length(forward_block.equations)
            forward_variable_jacobian[i, j] = forward_block.jacobian[i, k1]
        end
    end
    mul!(M1, forward_variable_jacobian, transpose(future_policy_jacobian))
    for i in axes(J, 1)
        for (j, k) in enumerate(system_variables)
            J[i, j] += weight*forward_block.jacobian[i, k + endogenous_nbr]
        end
        for (j, k) in enumerate(ids.system_in_state)
            J[i, k] += weight*M1[i, j]
        end
    end
    return J
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

function SG_NLsolve!(X, lb, ub, fx, J, states, parameters, oldgrid, sgws, solver, method, ftol, show_trace)
    ns = size(states, 2)
    nx = size(X, 2)
    if solver == PATHSolver
        # PATHsolver is not thread-safe
        PATHsolver_solve!(X, lb, ub, states, fx, J, oldgrid, sgws[1], show_trace, 1:ns)
    elseif Threads.nthreads() == 1
        if solver == NonlinearSolver
            NonlinearSolver_solve!(X, states, J, oldgrid, method, sgws[1], show_trace, 1:ns)
        elseif solver == NLsolver
            NLsolve_solve!(X, lb, ub, fx, J, states, oldgrid, sgws[1], ftol, show_trace, 1:ns)
        else
            error("Sparsegrids: unknown solver")
        end
    else
        previousblasnumthreads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        chunks = collect(Iterators.partition(1:ns, (ns ÷ Threads.nthreads()) + 1))
        if solver == NonlinearSolver
            tasks = map(enumerate(chunks)) do (i, chunk)
                Threads.@spawn NonlinearSolver_solve!(X, states, J, oldgrid, method, sgws[i], show_trace, chunk)
            end
        elseif solver == NLsolver
            tasks = map(enumerate(chunks)) do (i, chunk)
                Threads.@spawn NLsolve_solve!(X, lb, ub, fx, J, states, oldgrid, sgws[i], ftol, show_trace, chunk)
            end
        else
            error("Sparsegrids: unknown solver")
        end
        fetch.(tasks)
        BLAS.set_num_threads(previousblasnumthreads)
    end
    return X
end

function PATHsolver_solve!(X, lb, ub, states, fx, J, oldgrid, sgws, show_trace, chunk)
    for i in chunk
        @views begin 
            state = states[:, i]
            fx .= X[:, i]
        end
        f1!(r, x) = sysOfEqs!(r, x, state, oldgrid, sgws)
        JA1!(Jx, x) = sysOfEqs_derivatives_update!(Jx, x, state, oldgrid, sgws)
        (status, results, info) = mcp_solve!(PathNLS(), f1!, JA1!, J, lb, ub, fx, silent=true, convergence_tolerance=1e-4)
        if status != 1
            @show status
            error("sparsegrids: solution update failed")
        end
        @views X[:, i] .= results
    end
    nothing
end

function NLsolve_solve!(X, lb, ub, fx, J, states, oldgrid, sgws, ftol, show_trace, chunk)
    for i in chunk
        @views begin 
            state = states[:, i]
            x = X[:, i]
        end
        f2!(r, x) = sysOfEqs!(r, x, state, oldgrid, sgws)
        JA2!(Jx, x) = sysOfEqs_derivatives_update!(Jx, x, state, oldgrid, sgws)
        df = OnceDifferentiable(f2!, JA2!, x, fx, J)
        let res
            try
                res = NLsolve.mcpsolve(df, lb, ub, x, ftol = ftol, show_trace = show_trace)
            catch e
                @show lb
                @show ub
                @show e
                @show x
                x .= evaluate(oldgrid, state)
                continue
            else
                if !res.f_converged
                    @show res.f_converged
                    error("sparsegrids: solution update failed")
                end
                x .= res.zero
            end 
        end
    end
    nothing
end

function NonlinearSolver_solve!(X, states, J, oldgrid, method, sgws, show_trace, chunk)
    f3!(r, x, state) = sysOfEqs!(r, x, state, oldgrid, sgws)
    JA3!(Jx, x, state) = sysOfEqs_derivatives_update!(Jx, x, state, oldgrid, sgws)
    nlf = NonlinearFunction(f3!, jac = JA3!, jac_prototype = J)
    for i in chunk
        @views begin 
            state = states[:, i]
            x = X[:, i]
        end
        nlp = NonlinearProblem(nlf, x, state)
        res = NonlinearSolve.solve(nlp, method, show_trace = Val(show_trace))
        @views x .= res.u
    end
    nothing
end
