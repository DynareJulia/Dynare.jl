using PATHSolver, Dynare, Tasmanian, LinearAlgebra, SparseArrays
using FiniteDiff, Parameters
using Test

context = @dynare "models/global/irbc_small_irr_test.mod" "notmpterms";

dimRef = -1
ftol = 1e-5
iterRefStart = 25
gridDepth = 2
gridOrder = 1
gridRule = "localp"
maxiter = 300
maxRef = 1
maxRefLevel = gridDepth + maxRef
mcp = false
method = Dynare.NewtonRaphson()
numstart = 0
savefreq = 10
scaleCorrInclude = []
scaleCorrExclude = []
show_trace = false
solver = nothing
surplThreshold= 1e-3
tol_ti = 1e-4
TT = 10000
typeRefinement = "classic"

@testset "blocks" begin
    endogenous_variables = Dynare.get_endogenous(context.symboltable)
    model = context.models[1]
    endogenous_nbr = model.endogenous_nbr
    exogenous_nbr = model.exogenous_nbr
    work = context.work
    (dynamic_state_variables, predetermined_variables, system_variables,
     backward_block, forward_block, preamble_block) = Dynare.make_block_functions(context)
    @test dynamic_state_variables == [1, 4, 9, 12]
    @test system_variables == [1, 3, 4, 6, 7]

    lb = Vector{Float64}(undef, 0)
    ub = Vector{Float64}(undef, 0)
    bmcps = Vector{Vector{Int}}(undef, 0)
    # must be consistent with SysOfEqs()
    block_eqs = union(forward_block.equations, backward_block.equations)

    @test block_eqs == [2, 5, 1, 4, 7]
    # computes lb, ub, bmcps
    Dynare.block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables .+ endogenous_nbr)

    @test bmcps == [[1, 2], [1, 4]]
    @test lb == [-Inf, 0.0, -Inf, 0.0, -Inf]
    @test ub == [Inf, Inf, Inf, Inf, Inf]

    equation_xref_list, variable_xref_list = Dynare.xref_lists(context)
    params = context.work.params
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    nstates = length(dynamic_state_variables)

    sgmodel = Dynare.SGModel(repeat(steadystate, 3),
                             zeros(2*endogenous_nbr),
                             zeros(exogenous_nbr),
                             endogenous_nbr,
                             exogenous_nbr,
                             params,
                             steadystate,
                             []
                        )

    ids = Dynare.SGIndices(
        endogenous_nbr,
        model.i_fwrd_b,
        dynamic_state_variables,
        system_variables
    )
    forward_system_variables = model.i_fwrd_b
    @test forward_system_variables == [1, 2, 4, 5, 7]
    setdiff!(forward_system_variables, predetermined_variables .- endogenous_nbr)
    @test forward_system_variables == [1, 4, 7]

    gridDim = nstates
    gridOut = length(system_variables)
    limits = context.work.limits
    nl = length(limits)
    gridDomain = zeros(gridDim, 2)
    @test ids.state_variables == [1, 4, 2, 5]
    for l in limits
        i = findfirst(x -> x == context.symboltable[string(l[1])].orderintype, ids.state_variables)
        if !isnothing(i)
            gridDomain[i,1] = l[2].min
            gridDomain[i,2] = l[2].max
        end
    end

    grid0, aPoints, aNum, nPols =
        Dynare.set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain);
    @test aNum == 41
    @test nPols == 4

    ################################################################################
    #                        Adaptivity parameters                                 #
    ################################################################################

    grid = copyGrid(grid0)
    polGuess = Dynare.guess_policy(context, aNum, nPols, aPoints, dynamic_state_variables, system_variables, lb, ub, backward_block, forward_block, sgmodel)
    loadNeededPoints!(grid, polGuess)

    monomial = Dynare.MonomialPowerIntegration(exogenous_nbr)
    nodesnbr = length(monomial.nodes)
    
    # Scale correction in the refinement process:
    @views scaleCorr = Dynare.make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables])

    n = length(system_variables)
    backward_variable_jacobian = Matrix{Float64}(undef, length(backward_block.equations), length(system_variables)) 
    forward_equations_nbr = length(forward_block.equations)
    sgoptions = Dynare.SGOptions(mcp, method, solver, ftol, show_trace)
    residuals = Vector{Float64}(undef, length(system_variables))
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
    shift_dyn_endogenous_vector = zeros(2*endogenous_nbr)
    forward_equations_nbr = length(forward_block.equations)
    M1 = Matrix{Float64}(undef, forward_equations_nbr, length(ids.state_in_system))
    policyguess = zeros(gridOut, nodesnbr)
    tmp_state_variables = Vector{Float64}(undef, length(dynamic_state_variables))
    sev = Dynare.Sev(zeros(3*endogenous_nbr), zeros(nstates))
    sgws_ = Dynare.SGWS(monomial, dynamic_state_variables, system_variables, 
                bmcps, ids, lb, ub, backward_block, backward_variable_jacobian, 
                forward_block, preamble_block,
                residuals, forward_jacobian, forward_points, 
                forward_residuals, forward_variable_jacobian, J, fx,
                policy_jacobian, evalPt, 
                future_policy_jacobian, 
                dyn_endogenous_vector, M1, polGuess, tmp_state_variables,
                sev, sgmodel, sgoptions)
    if Threads.nthreads() > 1
        BLAS.set_num_threads(1)
        sgws = [deepcopy(sgws_) for i in 1:Threads.nthreads()]
    else
        sgws = [sgws]
    end

    @testset "forward looking equations residuals!" begin
        state = [0.95, 1.1, 0.1, -0.1]
        policy = [1, 0.1, 1, 0.1, 1.38]
        @unpack monomial, dynamic_state_variables, system_variables, bmcps, forward_block,
        preamble_block, sgmodel, forward_points, forward_residuals, dyn_endogenous_vector = sgws[1]
        @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel
        nodes = monomial.nodes
        numNodes = length(nodes)
        map(v -> fill!(v, 0.0), dyn_endogenous_vector)

        dyn_endogenous[dynamic_state_variables] .= state
        dyn_endogenous[system_variables .+ endogenous_nbr] .= policy 

        Dynare.set_dyn_endogenous_vector!(dyn_endogenous_vector, policy, state, grid, sgws[1])
        
        @test forward_points[:, 1] ≈ [1, 1, 0.095, -0.095 ]
        # 3) Determine relevant variables within the expectations operator
        X = evaluateBatch(grid, forward_points)
        fill!(exogenous, 0.0)
        fill!(forward_residuals, 0.0)
        for i in 1:numNodes
            # setting future variables
            @views begin
                dyn_endogenous_vector[i][2*endogenous_nbr .+ system_variables] .= X[:, i]
                forward_block.get_residuals!([], forward_residuals[:, i], dyn_endogenous_vector[i], exogenous, parameters, steadystate)
            end
        end
        xx = X[:, 1]
        @test dyn_endogenous_vector[1] ≈  [0.95, 0, 0, 1.1, 0, 0, 0,
                                              1, 0.1, 0.1, 1, -0.1, 0.1, 1.38,
                                              xx[1], 0.095, xx[2], xx[3], -0.095, xx[4], xx[5]]
        expected_residuals = forward_residuals*monomial.weights
        backward_residuals = zeros(length(system_variables) - length(expected_residuals))
        backward_block.get_residuals!([], backward_residuals, dyn_endogenous, exogenous, parameters, steadystate)
        
        residuals = vcat(expected_residuals, backward_residuals)
        Dynare.reorder!(residuals, bmcps) 

        residuals_1 = zeros(length(system_variables))
        Dynare.sysOfEqs!(residuals_1, policy, state, grid, sgws[1])

        @test  residuals ≈ residuals_1 
    end

    @testset "SysOfEqs" begin
        state = [1.0, 1.0, 0.1,  -0.1]
        policy = [1.0, 0.1, 1.0, 0.1, 1.38]
        Dynare.sysOfEqs!(residuals, policy, state, grid, sgws[1])
    end

    
    @testset "Tasmananian derivatives" begin

        state = [1.0, 1.0, 0, 0]
        Ds1 = zeros(4, 5)
        Ds2 = zeros(5, 4)
        differentiate!(Ds1, grid, state)
        resid = zeros(5)
        f(s) = evaluateBatch(grid, s)
        Ds2 = transpose(FiniteDiff.finite_difference_jacobian(f, state))
        @test Ds1 ≈ Ds2
    end

    @testset "forward_looking_equation_derivatives!" begin
#        forward_looking_equation_derivatives!(J, x, state, grid, nodes, weight, sgws)
        @unpack dynamic_state_variables, system_variables, ids, dyn_endogenous_vector,
            backward_block, forward_block, preamble_block, sgmodel, policy_jacobian, future_policy_jacobian, evalPt, M1 = sgws[1]
        @unpack dyn_endogenous, exogenous, endogenous_nbr, exogenous_nbr, parameters, steadystate = sgmodel

        node_id = 1
        @views dyn_endogenous = dyn_endogenous_vector[node_id]
        # policy function Jacobian
        Tasmanian.differentiate!(policy_jacobian, grid, dyn_endogenous[endogenous_nbr .+ dynamic_state_variables])
        # future policy Jacobian = policy Jacobian x derivatives of predetermined variables w.r. state
        Dynare.copy_from_submatrix!(future_policy_jacobian, policy_jacobian, ids.state_in_system, ids.forward_in_system)
        @test future_policy_jacobian == policy_jacobian[[1, 2], [1, 3, 5]]
        forward_block.update_jacobian!([], forward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
        @views forward_j = Matrix(forward_block.jacobian[:, 2*endogenous_nbr .+ system_variables[ids.forward_in_system]])
        mul!(M1, forward_j, transpose(future_policy_jacobian))
        @test M1 == forward_block.jacobian[:, [15, 18, 21]]*transpose(future_policy_jacobian)
                
        @views jacobian = Matrix(forward_block.jacobian[:, endogenous_nbr .+ system_variables])
        old_jacobian = copy(jacobian)        
        @views jacobian[:, ids.system_in_state] .+= M1
        @test jacobian[:, ids.system_in_state] ≈ old_jacobian[:, [1, 3]] .+ M1

        sgws1 = deepcopy(sgws[1])
        policy = dyn_endogenous[system_variables .+ endogenous_nbr]
        state = dyn_endogenous[dynamic_state_variables .+ endogenous_nbr]
        Dynare.sysOfEqs!(residuals, policy, state, grid, sgws1)
        J = zeros(length(forward_block.equations), length(system_variables))
        Dynare.forward_looking_equation_derivatives!(J, grid, dyn_endogenous_vector[node_id], 1.0, sgws1)
        @test J ≈ jacobian
    end

    @testset "sysOfEqs_derivatives_update!" begin
        x = [1, 0, 1, 0, 1.3]
        state = [1, 1, 0, 0]
        n = length(system_variables) 
        J = zeros(n,n)
        sgws1 = deepcopy(sgws[1])
        Dynare.sysOfEqs_derivatives_update!(J, x, state, grid, sgws[1])
        
        Jtarget = zeros(n,n)
        f!(r, y) = Dynare.sysOfEqs!(r, y, state, grid, sgws1)
        FiniteDiff.finite_difference_jacobian!(Jtarget, f!, x)
        @test J ≈ Jtarget
    end

    @testset "sysOfEqs_with_jacobian!" begin
        x = [1, 0, 1, 0, 1.3]
        state = [1, 1, 0, 0]
        n = length(system_variables) 
        resid = zeros(n)
        J = zeros(n,n)
        Dynare.sysOfEqs_with_jacobian!(resid, J, x, state, grid, sgws[1])

        rtarget = zeros(n)
        Jtarget = zeros(n,n)
        f!(r, y) = Dynare.sysOfEqs!(r, y, state, grid, sgws[1])
        Dynare.sysOfEqs!(rtarget, x, state, grid, sgws[1])
        FiniteDiff.finite_difference_jacobian!(Jtarget, f!, x)
        @test J ≈ Jtarget
    end

 #    global SGWS = sgws
#    global PREAMBLE_BLOCK = preamble_block
end

function test_diff(J, x, y, state, grid,
                   nPols, exogenous, nodes, weights, params, steadystate,
                   forward_equations_nbr, endogenous_nbr, exogenous_nbr,
                   dynamic_state_variables, system_variables, bmcps, dynamicws, model,
                   predetermined_variables, forward_system_variables,
                   forward_expressions_eqs, other_expressions_eqs, ids)
    y = repeat(steadystate, 3)
    T = dynamicws.temporary_values

    @views begin
        y[dynamic_state_variables] .= state
        y[system_variables] .= x
    end

    jacobian = get_dynamic_jacobian!(
        dynamicws,
        params,
        y,
        exogenous,
        steadystate,
        model,
        2,
    )
    y0 = copy(y)
    function f1(x)
        y[system_variables] .= x
        resid = zeros(forward_equations_nbr)
        forward_block(resid, T, y0, zeros(exogenous_nbr), params, steadystate)
        return resid
    end
    fjd = FiniteDiff.finite_difference_jacobian(f1, x)
    display(fjd)
    display(jacobian[forward_expressions_eqs, system_variables])

    function f2(y)
        resid = zeros(forward_equations_nbr)
        forward_block(resid, T, y, zeros(exogenous_nbr), params, steadystate)
        return resid
    end
    @show y
    fjd = FiniteDiff.finite_difference_jacobian(f2, y)
    display(fjd)
    display(Matrix(jacobian[forward_expressions_eqs, 1:21]))

    @show "Policy jacobian"
    evalPtr = y[endogenous_nbr .+ dynamic_state_variables]
    Jpolicy = zeros(4,5)
    differentiate!(Jpolicy, grid, evalPtr)
    display(transpose(Jpolicy))
    f3(x) = evaluateBatch(grid, x)
    Jpolicy_fd = FiniteDiff.finite_difference_jacobian(f3, evalPtr)
    display(Jpolicy_fd)

    @show "forward_looking residuals jacobian"
    @show "analytic"
    @show dynamic_state_variables[ids.state_in_system] .+ endogenous_nbr
    forward_equations_jacobian = copy(jacobian)
    @show weights[1]
    forward_equations_jacobian[forward_expressions_eqs, dynamic_state_variables[ids.state_in_system]] .+= jacobian[forward_expressions_eqs, forward_system_variables .+ 2*endogenous_nbr]*Jpolicy_fd[[1, 3, 5], 1:2]
    display(weights[1]*Matrix(forward_equations_jacobian[forward_expressions_eqs, system_variables .+ endogenous_nbr]))
    J = zeros(2, 5)
    forward_looking_equation_derivatives!(J, copy(y), x, exogenous, state, params, dynamicws, grid, nodes[1], weights[1], steadystate,
                                               model, forward_equations_nbr, state_variables, endogenous_nbr,
                                               exogenous_nbr, system_variables, forward_system_variables, forward_expressions_eqs,
                                               predetermined_variables, ids)
    display(J)
    println("finite difference")
    f4(x) = weights[1]*forward_looking_equation_residuals(x, state, params, grid, nodes,
                                                          steadystate, forward_equations_nbr,
                                                          state_variables, endogenous_nbr,
                                                          exogenous_nbr, system_variables)[1, :]
    jfd = FiniteDiff.finite_difference_jacobian(f4, x)
    display(jfd)

    @show "expected forward looking residuals jacobian "
    @show "analytic"
    display(forward_equations_derivatives(x, copy(y), nodes, weights, dynamicws, params, exogenous,
                                          steadystate, model,
                                          forward_equations_nbr, state, grid, state_variables, endogenous_nbr,
                                          exogenous_nbr, system_variables, forward_system_variables, forward_expressions_eqs,
                                          predetermined_variables, ids))
    @show "finite difference"
    function f5(x)
        expected_residuals = zeros(forward_equations_nbr)
        expectation_equations!(expected_residuals, T, x, state, params, grid, nodes, weights, steadystate, forward_equations_nbr, state_variables, endogenous_nbr, exogenous_nbr, system_variables)
        return expected_residuals
    end
    fjd = FiniteDiff.finite_difference_jacobian(f5, x)
    display(fjd)

    @show "sysOfEqs jacobian"
    @show "analytic"
    J = zeros(length(system_variables), length(system_variables))
    sysOfEqs_derivatives!(J, x, copy(y), state, grid,
    nPols, exogenous, nodes, weights, params, steadystate,
    forward_equations_nbr, endogenous_nbr, exogenous_nbr,
    state_variables, system_variables, bmcps, dynamicws, model,
    predetermined_variables, forward_system_variables,
    forward_expressions_eqs, other_expressions_eqs, ids)
    display(J)

    @show "finite differences"
    function f10(x)
        return sysOfEqs(x, T, copy(y), state, grid, getNumDimensions(grid), nodes, weights,
                      params, steadystate, forward_equations_nbr, endogenous_nbr,
                      exogenous_nbr, state_variables, system_variables, bmcps)
    end
    fjd = FiniteDiff.finite_difference_jacobian(f10, x)
    display(fjd)
end


