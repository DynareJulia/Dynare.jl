using BenchmarkTools
using LinearAlgebra.BLAS
include("../src/perfectforesight/gmres_solver.jl")
include("../src/perfectforesight/perfectforesight_solvers.jl")
#context = @dynare "../test/models/example1/example1.mod"
context = @dynare "test/models/irbc/irbc1.mod" "savemacro" "-DN=200"


md = context.models[1]
endo_steadystate = context.results.model_results[1].trends.endogenous_steady_state
exo_steadystate = context.results.model_results[1].trends.exogenous_steady_state
work = context.work
dynamic_variables =  
dynamic_variables = work.dynamic_variables
lli = md.lead_lag_incidence
endogenous = repeat(endo_steadystate, 3)
exogenous = repeat(exo_steadystate', 3)
period = 2
@btime get_dynamic_endogenous_variables!(dynamic_variables, endogenous,
                                         lli, md, period)
dynamic! = md.dynamic!.dynamic!
temporary_values = work.temporary_values
residuals = work.residuals
jacobian = work.jacobian
params = work.params
@btime fill!(jacobian, 0.0)
#=
include("models/irbc/irbc1Dynamic.jl")
@btime irbc1Dynamic.dynamicResidTT!(temporary_values,
                         dynamic_variables, exogenous, params, endo_steadystate, period)

@time irbc1Dynamic.dynamicResid!(temporary_values, residuals,
                       dynamic_variables, exogenous, params, endo_steadystate, period, true)


@btime irbc1Dynamic.dynamicG1TT!(temporary_values,
                      dynamic_variables, exogenous, params, endo_steadystate, period)

@btime irbc1Dynamic.dynamicG1!(temporary_values, jacobian,
                    dynamic_variables, exogenous, params, endo_steadystate, period, true)

@btime irbc1Dynamic.dynamic!(temporary_values, residuals, jacobian,
                  dynamic_variables, exogenous, params, endo_steadystate, period)
=#

@btime Base.invokelatest(dynamic!,
                         temporary_values,
                         residuals,
                         jacobian,
                         dynamic_variables,
                         exogenous,
                         params,
                         endo_steadystate,
                         period)

@btime compute_jacobian(work, dynamic_variables, exogenous,
                        endo_steadystate, md, period)
#=
periods = 200;
preconditioner_window = 7
res = zeros(periods*md.endogenous_nbr)
res[[3, 4]] .= 0.1
rout = zeros(periods*md.endogenous_nbr)
ws = GmresWs(periods, preconditioner_window, context, "GS")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= repeat(exo_steadystate', periods + 2)

n = md.endogenous_nbr
work = context.work
function f1(periods)
    work = context.work
    get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
    get_abc!(ws.a, ws.b, ws.c, work.jacobian, md)
    A = makeA(hcat(ws.a, ws.b, ws.c), ws.g, periods)
    y = A\res
end


function f2(periods, ws)
    ws.endogenous .= repeat(endo_steadystate, periods + 2)
    ws.exogenous .= zeros(periods + 2, md.exogenous_nbr)
    gmres!(rout, ws.LREMap, res, log=false,
       verbose=false, Pr=ws.P)
end


BLAS.set_num_threads(3)

z = zeros(n*periods)
y = zeros(n*periods)
jacobian_time_vec!(z, ws.dynamic_variables, ws.residuals,
                   ws.endogenous, ws.exogenous,
                   ws.steadystate, ws.presiduals, ws.g,
                   ws.temp_vec, context.work, md, periods)
=#

#=
@btime get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
@btime get_abc!(ws.a, ws.b, ws.c, work.jacobian, md)
@btime makeA(hcat(ws.a, ws.b, ws.c), ws.g, periods)
A = makeA(hcat(ws.a, ws.b, ws.c), ws.g, periods)
@btime y = A\res
@btime f1(periods)
@btime jacobian_time_vec!(z, ws.dynamic_variables, ws.residuals,
                   ws.endogenous, ws.exogenous,
                   ws.steadystate, ws.presiduals, ws.g,
                   ws.temp_vec, context.work, md, periods)
@btime ldiv!(z, ws.P, res)
@btime y = ws.LREMap*z
@btime f2(periods, ws)
=#

#=
for i = 1:4
    BLAS.set_num_threads(i)
    @show i
    @btime f1(periods)

    @btime f2(periods, preconditioner_window)
end
=#
