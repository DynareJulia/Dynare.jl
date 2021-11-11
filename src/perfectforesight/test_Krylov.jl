using Krylov
using LinearOperators

n = periods * nvar
LP = LineraOperator(
    Float64,
    n,
    n,
    false,
    false,
    v -> preconditioner!(v, x, ws.P.g, ws.P.hh, ws.P.preconditioner_window, ws.P.periods),
)
dqgmres(ws.LREMap, res, N = LP)
