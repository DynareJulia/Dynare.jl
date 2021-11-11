using Dynare
using IterativeSolvers

include("perfectforesight/perfectforesight_solver_1.jl")
include("perfectforesight/linesearch.jl")

get_initial_dynamic_endogenous_variables! = Dynare.get_initial_dynamic_endogenous_variables!
get_terminal_dynamic_endogenous_variables! =
    Dynare.get_terminal_dynamic_endogenous_variables!
get_dynamic_endogenous_variables! = Dynare.get_dynamic_endogenous_variables!

struct NewtonWs
    fvec::Vector{Float64}
    fold::Vector{Float64}
    p::Vector{Float64}
    g::Vector{Float64}
    function NewtonWs(n)
        fvec = Vector{Float64}(undef, n)
        fold = Vector{Float64}(undef, n)
        p = similar(fvec)
        g = similar(fvec)
        new(fvec, fold, p, g)
    end
end

function check_convergency(f, tolf, n)
    return norm(f) / n < tolf
end

function get_direction!(
    p::AbstractVector{Float64},
    fvec::AbstractVector{Float64},
    jacobian::AbstractMatrix{Float64},
    algo::String,
)
    if algo == "LU"
        F = lu(jacobian)
        copy!(p, fvec)
        ldiv!(F, p)
    elseif algo == "GMRES"
        τ = 0.00000
        my_lu = ilu(jacobian, τ = τ)
        @show any(isnan.(jacobian))
        @show any(isnan.(fvec))
        fill!(p, 0.0)
        gmres!(
            p,
            jacobian,
            fvec,
            Pr = my_lu,
            initially_zero = true,
            verbose = true,
            maxiter = 4,
        )
    end
    lmul!(-1.0, p)
    return p
end

function newton_solver!(x, x0, func!, jac!, maxiter, tolf, tolx, direction_algo, ws)
    n = length(x)
    fvec = ws.fvec
    fold = ws.fold
    g = ws.g
    p = ws.p
    fvec = func!(fvec, x0)
    if !isreal(fvec)
        @show fvec
    end
    normf = norm(fvec)
    if normf / n < 0.1 * tolf
        return 0
    end
    jacobian = jac!(x0)
    iter = 0
    while iter < maxiter
        p = get_direction!(p, fvec, jacobian, direction_algo)
        mul!(g, transpose(jacobian), fvec)
        check, f = linesearch!(x, fvec, x0, g, p, func!; stpmx = 100, tolx = 1e-5)
        @debug "new f: $f"
        if check == 1
            return 1
        end
        #        x .+= p
        copy!(fold, fvec)
        if check_convergency(f, tolf, n)
            return (x, 0)
        end
        copy!(x0, x)
        jacobian = jac!(x0)
    end
    return 1
end

function get_first_order_solution!(context::Context)
    results = context.results.model_results[1]
    model = context.models[1]
    work = context.work
    options = Dynare.StochSimulOptions(Dict{String,Any}())
    endogenous = results.trends.endogenous_steady_state
    exogenous = results.trends.exogenous_steady_state
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both
    tmp_nbr = model.dynamic!.tmp_nbr::Vector{Int64}
    ws = Dynare.DynamicJacobianWs(
        model.endogenous_nbr,
        model.exogenous_nbr,
        ncol,
        sum(tmp_nbr[1:2]),
    )
    Dynare.compute_first_order_solution!(
        results.linearrationalexpectations,
        endogenous,
        exogenous,
        endogenous,
        work.params,
        model,
        ws,
        options,
    )
end

#context = @dynare "test/models/irbc/irbc1" "-DN=2" 
context = @dynare "test/models/example1/example1"

get_first_order_solution!(context)

md = context.models[1]
results = context.results.model_results[1]
n = md.endogenous_nbr
p = md.exogenous_nbr
A = zeros(n, n)
A[:, md.i_bkwrd_b] .= results.linearrationalexpectations.g1_1
B = copy(results.linearrationalexpectations.g1_2)
periods = 120
exogenous = zeros(periods, p)
exogenous[2, md.exogenous_nbr] = 0.5
window = 2
y = zeros(n, periods)
tmp1 = similar(y)
tmp2 = similar(y)
c = context.results.model_results[1].trends.endogenous_steady_state
simul_first_order_1!(y, zeros(n), A, B, exogenous, window, periods, tmp1, tmp2)
steadystate = context.results.model_results[1].trends.endogenous_steady_state
y .+= steadystate
residuals = zeros(n * (periods - 2))
dynamic_variables = zeros(md.n_bkwrd + md.n_fwrd + md.n_current + 2 * md.n_both)
temp_vec = context.work.temporary_values

initialvalues = y[1:n]
terminalvalues = y[(periods-1)*n+1:periods*n]
periods = 3
dy = zeros(n * periods)
params = context.work.params
x0 = view(vec(y), n+1:(periods+1)*n)
x = similar(x0)

JJ = Jacobian(context, periods)
ws_threaded = [Dynare.DynamicJacobianWs(context) for i = 1:Threads.nthreads()]


func!(resid::AbstractVector{Float64}, x0::AbstractVector{Float64}) = get_residuals!(
    resid,
    x0,
    initialvalues,
    terminalvalues,
    exogenous,
    dynamic_variables,
    steadystate,
    params,
    md,
    periods,
    temp_vec,
)

jac!(x0::AbstractVector{Float64}) = makeJacobian!(
    JJ,
    x0,
    initialvalues,
    terminalvalues,
    exogenous,
    context,
    periods,
    ws_threaded,
)


M = jac!(x0);
@show length(JJ.ss.I)
@show length(JJ.ss.J)
@show length(JJ.ss.V)
@show length(JJ.ss.klasttouch)
@show length(JJ.ss.csrrowptr)
@show length(JJ.ss.csrcolval)
@show length(JJ.ss.csrnzval)
@show length(JJ.ss.csccolptr)
maxiter = 50
tolf = 1e-8
tolx = 1e-5
ws = NewtonWs(periods * n)
@time newton_solver!(x, x0, func!, jac!, maxiter, tolf, tolx, "GMRES", ws)
