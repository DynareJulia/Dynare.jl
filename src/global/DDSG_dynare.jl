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

include("DDSG.jl")

function AE_evaluate(grid,X)

    if hasproperty(grid, :is_ddsg) 
        return DDSG_evaluate(grid,X=X)
    else
        if length(X)==getNumDimensions(grid)
            return evaluate(grid,X)
        else
            return evaluateBatch(grid,X) 
        end
    end
end

function AE_evaluate!(Y,grid,X)
    if hasproperty(grid, :is_ddsg) 
        DDSG_evaluate!(grid,Y=Y,X=X)
    else
        evaluateBatch!(Y,grid,X) 
    end
end

function AE_differentiate!(J,grid,x)
    if hasproperty(grid, :is_ddsg) 
        DDSG_differentiate!(grid,J=J,x=x)
    else
        differentiate!(J,grid,x) 
    end
end

function AE_copyGrid(grid)
    if hasproperty(grid, :is_ddsg) 
        return DDSG(grid)
    else
        return copyGrid(grid) 
    end
end

function AE_getPoints(grid)
    if hasproperty(grid, :is_ddsg) 
        return DDSG_nodes(grid)
    else
        return getPoints(grid) 
    end
end

function AE_getNumPoints(grid)SG_NLsolve!
    if hasproperty(grid, :is_ddsg) 
        return sum(grid.grid_points)
    else
        return getPoints(grid) 
    end

end

function F_EQ!(; Y_buffer,X, oldgrid, sgws)

    map(s -> fill!(s.J, 0.0), sgws)

    @unpack system_variables, lb, ub, sgoptions, J, fx,
            sgmodel = sgws[1]
    @unpack mcp, method, solver, ftol, show_trace = sgoptions
    @unpack exogenous, parameters = sgmodel
    
    # Get the points that require function values

    nstates = size(X, 2)
    
    # Get the number of points that require function values
    # Array for intermediate update step
    n = length(system_variables)
    # TODO: check the case with non-autocorrelated shocks
    fill!(exogenous, 0.0)

    state = view(X, :, 1)

    # Default solver
    if isnothing(solver)
        solver = mcp ? NLsolver : NonlinearSolver 
    end

    # Solving nonlinear problem
    SG_NLsolve!(Y_buffer, lb, ub, fx, J, X, parameters, oldgrid, sgws, solver, method, ftol, show_trace )
    
    # Add the new function values to newgrid
    return view(Y_buffer, :, 1:nstates)
end

function F_EQ2(;X, oldgrid, sgws)

    @unpack system_variables, lb, ub, sgoptions, J, fx, sgmodel = sgws[1]
    @unpack mcp, method, solver, ftol, show_trace = sgoptions
    @unpack exogenous, parameters = sgmodel
    
    # Get the points that require function values

    #Y_buffer = Matrix{Float64}(undef, oldgrid.dof, size(X, 2))
    #Y_buffer  = DDSG_evaluate(oldgrid,X=X)


    nstates = size(X, 2)
    
    # Get the number of points that require function values
    # Array for intermediate update step
    n = length(system_variables)
    # TODO: check the case with non-autocorrelated shocks
    fill!(exogenous, 0.0)

    state = view(X, :, 1)

    # Default solver
    if isnothing(solver)
        solver = mcp ? NLsolver : NonlinearSolver 
    end

    #print("state size ",size(state))

    #DDSG_print(oldgrid)



    Y_buffer = zeros(n,nstates)

    # println(n)
    # println(nstates)
    # println(size(Y_buffer))
    # readline()

    # Solving nonlinear problem
    #println("-----> X ",X,size(X))
    #println("-----> Before ",norm(Y_buffer)/(prod(size(Y_buffer))),size(Y_buffer))
    # println("Debug - Input Types:")
    # println("Y_buffer: ", typeof(Y_buffer) , " ", size(Y_buffer))
    # println("lb: ", lb)
    # println("ub: ", ub)
    # println("fx: ", fx)
    # println("J: ", size(J))
    # println("X: ", size(X))
    # println("parameters: ", parameters)
    # println("oldgrid: ", typeof(oldgrid))
    # println("sgws: ", typeof(sgws))
    # println("solver: ", typeof(solver))
    # println("method: ", typeof(method))
    # println("ftol: ", typeof(ftol))
    # println("show_trace: ", typeof(show_trace))
    
    println("-----> X ",X)
    SG_NLsolve!(Y_buffer, lb, ub, fx, J, X, parameters, oldgrid, sgws, solver, method, ftol, show_trace)
    #println("-----> After ",norm(Y_buffer)/(prod(size(Y_buffer))),size(Y_buffer))
    println("-----> F(X) ",Y_buffer)
   #  readline()

    # Add the new function values to newgrid
    return Y_buffer
end

function DDSGapproximation(; context::Context=context,
    dimRef=-1,
    k_max = 1,
    ftol=1e-5,
    iterRefStart=25,
    gridDepth=2,
    gridOrder=1,
    gridRule="localp",
    maxiter=300,
    maxRef=1,
    maxRefLevel=gridDepth + maxRef,
    mcp=false,
    method=NewtonRaphson(),
    numstart=0,
    savefreq=10,
    scaleCorrInclude=[],
    scaleCorrExclude=[],
    show_trace=false,
    solver=mcp ? NLsolver : NonlinearSolver,
    surplThreshold=1e-3,
    tol_ti=1e-4,
    drawsnbr=10000,
    typeRefinement="classic",
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
        zeros(2 * endogenous_nbr),
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
            gridDomain[i, 1] = l[2].min
            gridDomain[i, 2] = l[2].max
        end
    end

    println("DDSG gridDim=$gridDim gridOut=$gridOut gridDepth=$gridDepth gridOrder=$gridOrder k_max=$k_max maxRefLevel=$maxRefLevel")

    ################################################################################
    #                        Adaptivity parameters                                 #
    ################################################################################

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
    forward_points = Matrix{Float64}(undef, gridDim, nodesnbr)
    forward_residuals = Matrix{Float64}(undef, length(forward_block.equations), nodesnbr)
    forward_variable_jacobian = Matrix{Float64}(undef, length(forward_block.equations), length(ids.forward_in_system))
    J = zeros(n, n)
    fx = zeros(n)
    evalPt = Vector{Float64}(undef, nstates)
    future_policy_jacobian = zeros(length(ids.state_in_system), length(ids.forward_in_system))
    policy_jacobian = Matrix{Float64}(undef, gridDim, n)
    dyn_endogenous_vector = [zeros(3 * endogenous_nbr) for i in 1:nodesnbr]
    shift_dyn_endogenous_vector = zeros(2 * endogenous_nbr, nodesnbr)
    forward_equations_nbr = length(forward_block.equations)
    M1 = Matrix{Float64}(undef, forward_equations_nbr, length(ids.state_in_system))
    policyguess = zeros(gridOut, nodesnbr)
    tmp_state_variables = Vector{Float64}(undef, length(dynamic_state_variables))
    sev = Sev(zeros(3 * endogenous_nbr), zeros(nstates))
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


    #sg_grid0, aPoints, aNum, nPols = set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain)
    #grid = AE_copyGrid(grid0)

    M, N, P = block_first_order_approximation(context, sgws[1])
    #sg_polGuess = guess_policy(context, aNum, nPols, aPoints, sgws[1], M, N, P)
    # polGuess = guess_policy(context, aNum, nPols, aPoints, dynamic_state_variables, system_variables, lb, ub, backward_block, forward_block, sgmodel)
    #loadNeededPoints!(sg_grid0, sg_polGuess)

    #################################
    # AE Edit - Start
    ################################
    l_max_min = maxRefLevel
    F_guess = X_in -> guess_policy(context, size(X_in, 2), size(X_in, 1), X_in, sgws[1], M, N, P)


    println("---- Building ... ")
    # Make copy of ddsg0
    grid = DDSG(dim=gridDim, dof=gridOut,l_min=l_max_min,l_max=l_max_min, order=gridOrder, rule=gridRule,k_max=k_max, domain=gridDomain)
    DDSG_init!(grid)
    DDSG_build!(grid,F=F_guess,X0=grid.centroid,refinement_type=typeRefinement,refinement_tol=surplThreshold,scale_corr_vec=Float64.(scaleCorr) )    

    println("---- Done Building ... ")

    my_points = DDSG_nodes(grid)
    FX_val = F_guess(my_points)
    FX_est = DDSG_evaluate(grid,X=my_points)
    ERR    = FX_val-FX_est
    ERR_MAX = maximum(abs.(ERR))

    # println("X=",my_points,size(my_points),"\n")
    # println("FX=",FX_val,size(FX_val),"\n")
    # println("FX_est=",FX_est,size(FX_est),"\n")
    # println("ERR=",ERR,size(ERR),"\n")
    # println("ERR_MAX=",ERR_MAX,"\n")



    # readline()



    #println("X0 =" , grid.X0 )
    #println("Y0 =" , grid.Y0 )

    #println("Eval  =" , DDSG_evaluate(grid,X=[0.5, 0.6, 0.1, 0.1]) )
    #readline()

    #################################
    # AE Edit - End
    ################################

    #scaleCorrMat = repeat(Float64.(scaleCorr), 1, Tasmanian.getNumLoaded(grid))

    #X_sample  = DDSG_nodes(grid)
    # Fill the matrix with random points within the domain
    num_points = 100
    Random.seed!(124)
    X_sample = Matrix{Float64}(undef, gridDim, num_points)
    for i in 1:size(gridDomain, 1)
        X_sample[i, :] .= (gridDomain[i, 1] .+ (gridDomain[i, 2] - gridDomain[i, 1]) .* rand(num_points))/2
    end


    F_old = F_guess
    total_time = 0
    context.timings["sparsegrids"] = Vector{Millisecond}(undef, 0)
    timings = context.timings["sparsegrids"]
    iter = 1
    while iter <= maxiter

        tin = now()
        map(s -> fill!(s.J, 0.0), sgws)
        F_solve = X -> F_EQ2(X=X, oldgrid=grid, sgws=sgws)

        t_build = -time()

        println("Building New Grid")
        newgrid = DDSG(dim=gridDim, dof=gridOut,l_min=l_max_min,l_max=l_max_min, order=gridOrder, rule=gridRule,k_max=k_max, domain=gridDomain)
        DDSG_init!(newgrid)
        DDSG_build!(newgrid,F=F_solve,X0=newgrid.centroid, refinement_type=typeRefinement,refinement_tol=surplThreshold,scale_corr_vec=Float64.(scaleCorr) )
        println("Finsihed Building New Grid")

        t_build = time()+t_build
        #println(" @ Time Build: $(t_build)")

        #pol      = DDSG_evaluate(grid,X=X_sample)
        pol       = DDSG_evaluate(newgrid,X=X_sample)
        pol_new   = F_EQ2(X=X_sample, oldgrid=newgrid, sgws=sgws)

        #println(" @ Norm:  $(norm(pol)) $(norm(pol_new)) $(norm(solved_v))  ")

        #user_input = readline()

        metric_sup = 0.0
        metric_L2  = 0.0
        for i in 1:grid.dof
            @views abs_diff = abs.(pol[i, :] - pol_new[i, :])
            @views metric_sup = max(metric_sup, maximum(abs_diff))
            @views metric_L2 += sum(abs_diff.^2)
        end
        metric_L2  = (metric_L2/(size(X_sample,2)*grid.dof))^0.5
        metric     = min(metric_L2, metric_sup)

        grid = AE_copyGrid(newgrid)

        iteration_walltime = now() - tin
        iter > 1 && (total_time += iteration_walltime.value)
        push!(timings, iteration_walltime)
        println(" + DDSG Iteration: $iter, Grid pts: $(AE_getNumPoints(grid)), Metric: $metric Computing time: $(now() -tin)")

        if (metric < tol_ti)
            break
        end
        iter += 1
        F_old = F_solve
    end
    average_time = total_time / (iter - 1)
    println("Last grid points: $(grid.grid_points)")
    println("Average iteration time (except first one): $average_time")

    results = context.results.model_results[1].sparsegrids
    results.average_iteration_time = average_time
    results.drawsnbr = drawsnbr
    results.ftol = ftol
    #results.grid = grid
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

