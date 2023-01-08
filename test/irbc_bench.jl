using Dynare
using ForwardDiff
using LinearAlgebra
using Optim
using BenchmarkTools

context = @dynare "test/models/irbc/irbc_ss.mod" "-DN=2";

# function representing the static model
static_fn = context.dynarefunctions.static!.static!

# function computing the (partial) steady state
steady_state_fn = context.dynarefunctions.steady_state!.steady_state!

# endogenous variables number
endo_nbr = context.models[1].endogenous_nbr

# exogenous variables number
exo_nbr = context.models[1].exogenous_nbr

# number of temporary variables
tmp_nbr = sum(context.dynarefunctions.static!.tmp_nbr[1:2])

# temporary variables
tmp = Vector{Float64}(undef, (tmp_nbr))

# residuals
residuals = Vector{Float64}(undef, (endo_nbr))

# exogenous variables
x = zeros((exo_nbr))

# workspace
ys = zeros((endo_nbr))

# parameters
params = context.work.params

# static function
f(y) = begin
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
    return sum(res .* res)
end

@show "guess values are all one"
y = ones(endo_nbr)

@show "residuals"
@show f(y)

@show "compute partial solution for steady state"
g!(ys, y)
@show (ys)

@show "residuals with partial solution"
@show f(y)

@show "objective function for the solution"
@show h(y)

@show "Jacobian of the objective function"
df = zeros(endo_nbr, endo_nbr)
dh!(gh, y) = begin
    dg = ForwardDiff.jacobian(g!, ys, y)
    static_fn(tmp, residuals, df, y, x, params)
    mul!(gh, transpose(df * dg), residuals, 2, 0)
    return nothing
end
# @show dh(y)

@btime results = optimize(
    h,
    [
        2.31,
        0.02,
        0.25,
        0.26,
        3.26,
        5.26,
        2.25,
        3.65,
        6.56,
        8.98,
        5.25,
        2.25,
        5.65,
        6.56,
        7.24,
    ],
    BFGS(),
)
@btime results = optimize(
    h,
    dh!,
    [
        2.31,
        0.02,
        0.25,
        0.26,
        3.26,
        5.26,
        2.25,
        3.65,
        6.56,
        8.98,
        5.25,
        2.25,
        5.65,
        6.56,
        7.24,
    ],
    BFGS(),
)

println("minimum = $(results.minimum) with in " * "$(results.iterations) iterations")
