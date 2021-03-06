include("../src/perfectforesight/gmres_solver.jl")
include("../src/perfectforesight/perfectforesight_solvers.jl")
context = @dynare "../test/models/example1/example1.mod"
#context = @dynare "../test/models/irbc/irbc1.mod" "savemacro"

md = context.models[1]
endo_steadystate = context.results.model_results[1].trends.endogenous_steady_state
exo_steadystate = context.results.model_results[1].trends.exogenous_steady_state

function bench1(periods, preconditioner_window)
    res = zeros(periods*md.endogenous_nbr)
    res[[4, 3]] .= 0.001
    rout = zeros(periods*md.endogenous_nbr)
    ws = GmresWs(periods, preconditioner_window, context, "CR")
    ws.endogenous .= repeat(endo_steadystate, periods + 2)
    ws.exogenous .= repeat(exo_steadystate', periods + 2)
    gmres_solver!(rout, res, periods, preconditioner_window,
                  md, context.work, ws, verbose=true)
end

@time bench1(200, 20)
