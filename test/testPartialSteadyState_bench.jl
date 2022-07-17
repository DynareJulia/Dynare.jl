using Dynare
using ForwardDiff
using LinearAlgebra
using BenchmarkTools
using Optim

context = @dynare "E:\\DynareJulia\\example1ss_partial.mod";

# function representing the static model
static_fn = @timev context.dynarefunctions.static!.static!

# function computing the (partial) steady state
steady_state_fn = @timev context.dynarefunctions.steady_state!.steady_state!

# endogenous variables number
endo_nbr = @timev context.models[1].endogenous_nbr

# exogenous variables number
exo_nbr = @timev context.models[1].exogenous_nbr

# number of temporary variables
tmp_nbr = @timev context.dynarefunctions.static!.tmp_nbr[1]

# temporary variables
tmp = Vector{Float64}(undef, tmp_nbr)

# residuals
residuals = Vector{Float64}(undef, endo_nbr)

# exogenous variables
x = zeros(exo_nbr)

# workspace
ys = zeros(endo_nbr)

# parameters
params = context.work.params

# static function
f(y) = begin
    @show y
    static_fn(tmp, residuals, y, x, params)
    return residuals
end

# partial steady state function
g!(ys, y) = begin
    copy!(ys, y)
    steady_state_fn(ys, x, params)
end

# objective function for the minimization problem
# When the minimum of h(y) is zero we have the solution 
function h(y)
    g!(ys, y)
    res = f(ys)
    return sum(res.*res)
end

@show "guess values are all one"
y = @timev ones(endo_nbr)

@show "residuals"
@show @timev f(y)

@show "compute partial solution for steady state"
g!(ys, y)
@benchmark @show ys

@show "residuals with partial solution"
@benchmark @show f(y)

@show "objective function for the solution"
@benchmark @show h(y)

@show "Jacobian of the objective function"
df = zeros(endo_nbr, endo_nbr)
dh(y) = begin
    dg = ForwardDiff.jacobian(g!, ys, y)
    static_fn(tmp, residuals, df, y, x, params)
    return df*dg
end
@show dh(y)

results = optimize(h, [2.31, 0.0025, 0.256, 0.2635, 3.265, 5.265])

println("minimum = $(results.minimum) with in "*
"$(results.iterations) iterations")