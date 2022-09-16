using BenchmarkTools
using Dynare
using SparseArrays

context = @dynare "../test/models/irbc/irbc_pf" "-DN=2"

periods = 300

m = context.models[1]
ddf = context.dynarefunctions
ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
tmp_nbr = ddf.dynamic!.tmp_nbr::Vector{Int64}
dynamic_ws = Dynare.DynamicWs(m.endogenous_nbr, m.exogenous_nbr, ncol, sum(tmp_nbr[1:2]))
perfect_foresight_ws = Dynare.PerfectForesightWs(context, periods)
X = perfect_foresight_ws.shocks
initialvalues = Dynare.get_dynamic_initialvalues(context)
guess_values = Dynare.perfect_foresight_initialization!(
    context,
    periods,
    "",
    X,
    perfect_foresight_ws,
    Dynare.steadystate,
    dynamic_ws,
)

results = context.results.model_results[1]
work = context.work
residuals = zeros(periods * m.endogenous_nbr)
dynamic_variables = dynamic_ws.dynamic_variables
temp_vec = dynamic_ws.temporary_values
steadystate = results.trends.endogenous_steady_state
terminalvalues = view(guess_values, :, periods)
params = work.params
JJ = perfect_foresight_ws.J

exogenous = perfect_foresight_ws.x

ws_threaded = [
    Dynare.DynamicWs(
        m.endogenous_nbr,
        m.exogenous_nbr,
        length(dynamic_variables),
        length(temp_vec),
    ) for i = 1:Threads.nthreads()
]

function f!(residuals, y)
    @debug "$(now()): start f!"
    Dynare.get_residuals!(
        residuals,
        vec(y),
        initialvalues,
        terminalvalues,
        exogenous,
        dynamic_variables,
        steadystate,
        params,
        m,
        ddf,
        periods,
        temp_vec,
    )
    @debug "$(now()): end f!"
end

function J!(
    A::SparseArrays.SparseMatrixCSC{Float64,Int64},
    y::AbstractVecOrMat{Float64},
)
    @debug "$(now()): start J!"
    A = Dynare.makeJacobian!(
        JJ,
        vec(y),
        initialvalues,
        terminalvalues,
        exogenous,
        context,
        periods,
        ws_threaded,
    )
    @debug count(!iszero, A) / prod(size(A))
    @debug "$(now()): end J!"
end

function fj!(residuals, JJ, y)
    f!(residuals, vec(y))
    J!(JJ, vec(y))
end

@debug "$(now()): start makeJacobian"
A0 = Dynare.makeJacobian!(
    JJ,
    vec(guess_values),
    initialvalues,
    terminalvalues,
    exogenous,
    context,
    periods,
    ws_threaded,
)
@debug "$(now()): end makeJacobian"
@debug "$(now()): start f!"
f!(residuals, vec(guess_values))
@debug "$(now()): end f!"
@debug "$(now()): start J!"
J!(A0, guess_values)
@debug "$(now()): end J!"
df = Dynare.OnceDifferentiable(f!, J!, vec(guess_values), residuals, A0)
@debug "$(now()): start nlsolve"

rr = copy(residuals)
F = Dynare.lu(A0)
b = @benchmark res = Dynare.nlsolve(df, vec(guess_values), method = :robust_trust_region, show_trace = false, ftol=cbrt(eps()));


