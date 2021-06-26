using Krylov
using LinearOperators

@btime begin
NN = periods*n
ws = Dynare.GmresWs(periods, 50, context, "GS")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= zeros(periods + 2, md.exogenous_nbr)

LP = LinearOperator(Float64, NN, NN, false, false,
                    v -> Dynare.preconditioner!(rout, v, ws.P.g, ws.P.hh,
                                                ws.P.preconditioner_window,
                                                ws.P.periods))
dqgmres(ws.LREMap, res, N=LP, itmax=30, verbose=2)
end
