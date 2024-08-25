using Chairmarks, Dynare, Parameters, Tasmanian


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

endogenous_variables = Dynare.get_endogenous(context.symboltable)
model = context.models[1]
endogenous_nbr = model.endogenous_nbr
exogenous_nbr = model.exogenous_nbr
work = context.work

@show "make_block_functions"
s = @b (dynamic_state_variables, predetermined_variables, system_variables,
 backward_block, forward_block, preamble_block) = Dynare.make_block_functions(context)
@show s

(dynamic_state_variables, predetermined_variables, system_variables,
 backward_block, forward_block, preamble_block) = Dynare.make_block_functions(context)
lb = Vector{Float64}(undef, 0)
ub = Vector{Float64}(undef, 0)
bmcps = Vector{Vector{Int}}(undef, 0)
# must be consistent with SysOfEqs()
block_eqs = union(forward_block.equations, backward_block.equations)

# computes lb, ub, bmcps
@show "block_mcp!"
s = @b Dynare.block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables .+ endogenous_nbr)
@show s

Dynare.block_mcp!(lb, ub, bmcps, context, block_eqs, system_variables .+ endogenous_nbr)
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
setdiff!(forward_system_variables, predetermined_variables .- endogenous_nbr)

gridDim = nstates
gridOut = length(system_variables)
limits = context.work.limits
nl = length(limits)
gridDomain = zeros(gridDim, 2)
for l in limits
    i = findfirst(x -> x == context.symboltable[string(l[1])].orderintype, ids.state_variables) 
    if !isnothing(i)
        gridDomain[i,1] = l[2].min
        gridDomain[i,2] = l[2].max
    end
end

@show "set_initial_grid"
s = @b grid0, aPoints, aNum, nPols =
    Dynare.set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain);
@show s

grid0, aPoints, aNum, nPols =
    Dynare.set_initial_grid(gridDim, gridOut, gridDepth, gridOrder, gridRule, gridDomain);

################################################################################
#                        Adaptivity parameters                                 #
################################################################################

grid = Dynare.copyGrid(grid0)
polGuess = Dynare.guess_policy(context, aNum, nPols, aPoints, dynamic_state_variables, system_variables, lb, ub, backward_block, forward_block, sgmodel)
loadNeededPoints!(grid, polGuess)

monomial = Dynare.MonomialPowerIntegration(exogenous_nbr)
nodesnbr = length(monomial.nodes)

# Scale correction in the refinement process:
@views scaleCorr = Dynare.make_scaleCorr(scaleCorrInclude, scaleCorrExclude, endogenous_variables[system_variables])

n = length(system_variables)
forward_equations_nbr = length(forward_block.equations)
sgoptions = Dynare.SGOptions(mcp, method, solver, ftol, show_trace)
residuals = Vector{Float64}(undef, n)
backward_variable_jacobian = Matrix{Float64}(undef, length(backward_block.equations), length(system_variables)) 
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
forward_equations_nbr = length(forward_block.equations)
M1 = Matrix{Float64}(undef, forward_equations_nbr, length(ids.state_in_system))
tmp_state_variables = Vector{Float64}(undef, length(dynamic_state_variables))
sev = Dynare.Sev(zeros(3*endogenous_nbr), zeros(nstates))
sgws = Dynare.SGWS(monomial, dynamic_state_variables, system_variables, 
            bmcps, ids, lb, ub, backward_block, backward_variable_jacobian,
            forward_block, preamble_block,
            residuals, forward_jacobian, forward_points, 
            forward_residuals, forward_variable_jacobian, J, fx,
            policy_jacobian, evalPt, 
            future_policy_jacobian, 
            dyn_endogenous_vector, M1, polGuess, tmp_state_variables, sev, sgmodel, sgoptions)
scaleCorrMat = repeat(Float64.(scaleCorr), 1, getNumLoaded(grid))
            
iter = 1
# new policy guess to be computed
polGuess1 = copy(polGuess)
newgrid = Dynare.copyGrid(grid0) 
# Index of current grid level to control the number of refinements
ilev = gridDepth

fill!(sgws.J, 0.0)

grid0 = Dynare.copyGrid(grid)
polGuess0 = copy(polGuess1)

@show "ti_step"
s = @b Dynare.ti_step!(Dynare.copyGrid(newgrid), grid, polGuess1, sgws)
@show s

newgrid1 = Dynare.copyGrid(newgrid)
Dynare.ti_step!(newgrid1, grid, polGuess1, sgws)

@show "refine!"
s = @b Dynare.refine!(Dynare.copyGrid(newgrid1), polGuess1, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)
@show s

polGuess1 = Dynare.refine!(newgrid1, polGuess1, scaleCorr, scaleCorrMat, surplThreshold, dimRef, typeRefinement)

# Calculate (approximate) errors on tomorrow's policy grid
@show "policy_update"
s = @b Dynare.policy_update(grid, newgrid1, polGuess, polGuess1, length(system_variables))
@show s

metric, polGuess, grid = Dynare.policy_update(grid, newgrid1, polGuess, polGuess1, length(system_variables))

@unpack parameters = sgws.sgmodel
# Get the points that require function values
@show "getNeedPoints"
s = @b getNeededPoints(newgrid)
@show s

states = getNeededPoints(newgrid)
nstates = size(states, 2)

# Get the number of points that require function values
# Array for intermediate update step
n = length(system_variables)
# TODO: check the case with non-autocorrelated shocks
exogenous = zeros(exogenous_nbr)
fill!(exogenous, 0.0)

state = view(states, :, 1)

@show "closure definitions"
solver = Dynare.NonlinearSolver

@show solver

if isnothing(solver)
    solver = mcp ? Dynare.NLsolver : Dynare.NonlinearSolver 
end

# Time Iteration step
@show "Solution"
X = copy(polGuess0)
@show size(X)
s = @b    if solver == Dynare.NonlinearSolver
        Dynare.NonLinearSolver_solve!(copy(X), J, states, parameters, grid, sgws, method, ftol, show_trace )
    elseif solver == Dynare.NLsolver
        Dynare.NLsolve_solve!(copy(X), lb, ub, fx, J, states, grid, sgws, ftol, show_trace)
    elseif solver == Dynare.PATHSolver
        Dynare.PATHsolver_solve!(copy(X), lb, ub, states, fx, J, grid, sgws, ftol, show_trace)
    else
        error("Sparsegrids: unknown solver")
    end
@show s

    # Add the new function values to newgrid
#@show "loadNeededPoints"
#s = @b loadNeededPoints!(Dynare.copyGrid(newgrid), view(X, :, 1:nstates))
#@show s

@show "sysOfEqs"
residuals = zeros(length(system_variables))
policy = X[:, end]
state = states[:, end]
s = @b Dynare.sysOfEqs!(residuals, policy, state, grid, sgws)
@show s

nfwrd = length(forward_block.variables)
@show "set_dyn_endogenous_vector!"
dyn_endogenous_vector_1 = deepcopy(dyn_endogenous_vector)
s = @b Dynare.set_dyn_endogenous_vector!(dyn_endogenous_vector_1, policy, state, grid, sgws) 
@show s

@show "expectation_equations!"
s = @b (@views Dynare.expectation_equations!(residuals[1:nfwrd], policy, state, grid, sgws))
@show s

@show "backward_block.get_residuals!"
@unpack dyn_endogenous, exogenous, parameters, steadystate = sgws.sgmodel 
tempterms = []
s = @b (@views backward_block.get_residuals!(tempterms, residuals[nfwrd  + 1:end], dyn_endogenous, exogenous, parameters, steadystate))
@show s

@show "reorder!" 
s = @b Dynare.reorder!(residuals, bmcps)
@show s

@show "set_endogenous_variables!"
tempterms = []
shift_dyn_endogenous = zeros(2*endogenous_nbr)
s = @b Dynare.set_endogenous_variables!(tempterms, shift_dyn_endogenous, monomial.nodes[1], parameters, steadystate) evals=1
@show s

@show "forward_looking_equation_residuals!"
s = @b Dynare.forward_looking_equation_residuals!(forward_residuals, policy, state, grid, sgws)
@show s

@unpack J = sgws

#=
@show "sysOfEqs_derivatives!"
s = @b Dynare.sysOfEqs_derivatives!(J, policy, state, grid, sgws)
@show s
=#

@show "sysOfEqs_derivatives_update!"
s = @b Dynare.sysOfEqs_derivatives_update!(J, policy, state, grid, sgws)
@show s

@show "expectation_equations_derivatives"
J = zeros(2, length(system_variables))
s = @b Dynare.expectation_equations_derivatives(J, policy, state, grid, sgws)
@show s

@show "forward_looking_equation_derivatives!"
J = zeros(2, length(system_variables))
s = @b Dynare.forward_looking_equation_derivatives!(J, grid, dyn_endogenous, monomial.weights[1], sgws) evals=1
@show s



nothing

