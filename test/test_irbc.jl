using Dynare
using Test
context = @dynare "test/models/irbc/irbc1.mod" "savemacro" "-DN=2"


md = context.models[1]
endo_steadystate = context.results.model_results[1].trends.endogenous_steady_state
exo_steadystate = context.results.model_results[1].trends.exogenous_steady_state

periods = 2;
preconditioner_window = 2
res = zeros(periods*md.endogenous_nbr)
res[[4, 3]] .= 0.001
rout = zeros(periods*md.endogenous_nbr)
ws = Dynare.GmresWs(periods, preconditioner_window, context, "CR")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= repeat(exo_steadystate', periods + 2)

work = context.work
Dynare.get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
jacobian = Matrix{Float64}(undef, size(work.jacobian))
copy!(jacobian, work.jacobian)
Dynare.get_abc!(ws.a, ws.b, ws.c, jacobian, md)
P = Dynare.LREprecond(periods, preconditioner_window, ws.a, ws.b, ws.c, ws.g)

n = md.endogenous_nbr
time_vec = zeros(md.endogenous_nbr)
Dynare.get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
LREMap = Dynare.LinearMap(periods*n) do C, B
    @inbounds copyto!(ws.residuals, n + 1, B, 1, periods*n) 
    Dynare.jacobian_time_vec!(C, ws.dynamic_variables, ws.residuals, ws.endogenous, ws.exogenous,
                       ws.steadystate, ws.presiduals, ws.g, time_vec, context.work, md, periods)
end
JA = Dynare.Jacobian(context, periods)
A = Dynare.makeJacobian!(JA, ws.endogenous, ws.exogenous, context, periods)
x = randn(n*periods)
y = P\x
@test LREMap*y ≈ x
@test A*x ≈ LREMap*x
@show size(P.hh)
@show size(x)
target = [P.hh*x; ws.g*P.hh*x + P.hh[:, 1:n]*x[16:30]]
@test P\x ≈ target
x = zeros(n*periods)
x[[4,6]] .= 0.1
y = P\x
@test LREMap*y ≈ x

periods = 10
preconditioner_windows = 2
res = zeros(periods*md.endogenous_nbr)
res[[4, 6]] .= 0.1
rout = zeros(periods*md.endogenous_nbr)
ws = Dynare.GmresWs(periods, preconditioner_window, context, "CR")

Dynare.get_abc!(ws.a, ws.b, ws.c, jacobian, md)
P1 = Dynare.LREprecond(periods, preconditioner_window, ws.a, ws.b, ws.c, ws.g)
Dynare.get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
LREMap = Dynare.LinearMap(periods*n) do C, B
    @inbounds copyto!(ws.residuals, n + 1, B, 1, periods*n) 
    Dynare.jacobian_time_vec!(C, ws.dynamic_variables, ws.residuals, ws.endogenous, ws.exogenous,
                       ws.steadystate, ws.presiduals, ws.g, time_vec, context.work, md, periods)
end
y1 = P1\res

ws = Dynare.GmresWs(periods, preconditioner_window, context, "CR")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= zeros(periods + 2, md.exogenous_nbr)
Dynare.gmres!(rout, ws.LREMap, res, log=false,
       verbose=false, Pr=ws.P)
y = P1\rout
@test ws.LREMap*y ≈ res

Dynare.get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
jacobian = Matrix{Float64}(undef, size(work.jacobian))
copy!(jacobian, work.jacobian)
Dynare.get_abc!(ws.a, ws.b, ws.c, jacobian, md)
JA = Dynare.Jacobian(context, periods)
A = Dynare.makeJacobian!(JA, ws.endogenous, ws.exogenous, context, periods)
display([LREMap*rout A*rout res])
display([y1 rout])

periods = 400;
preconditioner_window = 3
res = zeros(periods*md.endogenous_nbr)
res[[3, 4]] .= 0.1
rout = zeros(periods*md.endogenous_nbr)
ws = Dynare.GmresWs(periods, preconditioner_window, context, "GS")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= repeat(exo_steadystate', periods + 2)

n = md.endogenous_nbr
work = context.work
function f1(periods)
@time    JA = Dynare.Jacobian(context, periods)
@time    A = Dynare.makeJacobian!(JA, ws.endogenous, ws.exogenous, context, periods)
@time    y = A\res
    return 0
end


function f2(periods, ws; verbose = false)
@time    JA = Dynare.Jacobian(context, periods)
@time    A = Dynare.makeJacobian!(JA, ws.endogenous, ws.exogenous, context, periods)
@time    Dynare.gmres!(rout, A, res, log=false,
                       verbose=verbose, Pr=ws.P)
@time    Dynare.ldiv!(ws.P, rout)
    return 0
end

function f3(periods, ws; verbose = false)
    ws.endogenous .= repeat(endo_steadystate, periods + 2)
    ws.exogenous .= zeros(periods + 2, md.exogenous_nbr)
@time    Dynare.gmres!(rout, ws.LREMap, res,Pr=ws.P,
                  log=false, verbose=verbose)
@time    Dynare.ldiv!(ws.P, rout)
    return 0
end

f1(periods)
fill!(rout, 0.0)
f2(periods, ws)
fill!(rout, 0.0)
f3(periods, ws)




