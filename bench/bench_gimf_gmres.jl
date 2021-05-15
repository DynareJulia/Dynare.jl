using BenchmarkTools
using Dynare
using SparseArrays
#context = @dynare "../../models/GIMF/gimf_gmres.mod" "nostrict"

context.options["stoch_simul"] = Dict()
context.options["stoch_simul"]["generalized_schur"] = Dict()
context.options["stoch_simul"]["cyclic_reduction"] = Dict()

md = context.models[1]
endo_steadystate = context.results.model_results[1].trends.endogenous_steady_state
exo_steadystate = context.results.model_results[1].trends.exogenous_steady_state

function bench1(periods, preconditioner_window, verbose=false)
    res = zeros(periods*md.endogenous_nbr)
    res[[4, 6]] .= 0.01
    rout = zeros(periods*md.endogenous_nbr)
    ws = Dynare.GmresWs(periods, preconditioner_window, context, "GS")
    ws.endogenous .= repeat(endo_steadystate, periods + 2)
    ws.exogenous .= repeat(exo_steadystate', periods + 2)
    Dynare.gmres!(rout, ws.LREMap, res,Pr=ws.P,
                  log=false, verbose=verbose)
    Dynare.ldiv!(ws.P, rout)
    return 0
end

function bench2(periods, preconditioner_window; verbose=false)
    res = zeros(periods*md.endogenous_nbr)
    res[[4, 6]] .= 0.01
    rout = zeros(periods*md.endogenous_nbr)
    ws = Dynare.GmresWs(periods, preconditioner_window, context, "GS")
    ws.endogenous .= repeat(endo_steadystate, periods + 2)
    ws.exogenous .= repeat(exo_steadystate', periods + 2)
    JA = Dynare.Jacobian(context, periods)
    ws_threaded = [Dynare.JacTimesVec(context) for i=1:Threads.nthreads()]
    A = Dynare.makeJacobian!(JA, ws.endogenous, ws.exogenous, context, periods, ws_threaded)
    Dynare.gmres!(rout, A, res, log=false,
           verbose=verbose, Pr=ws.P, initially_zero=true)
    Dynare.ldiv!(ws.P, rout)
    return 0
end


function bench3(periods, preconditioner_window; verbose=false)
    res = zeros(periods*md.endogenous_nbr)
    res[[4, 6]] .= 0.01
    rout = zeros(periods*md.endogenous_nbr)
    ws = Dynare.GmresWs(periods, preconditioner_window, context, "GS")
    ws.endogenous .= repeat(endo_steadystate, periods + 2)
    ws.exogenous .= repeat(exo_steadystate', periods + 2)
    JA = Dynare.Jacobian(context, periods)
    ws_threaded = [Dynare.JacTimesVec(context) for i=1:Threads.nthreads()]
    A = Dynare.makeJacobian!(JA, ws.endogenous, ws.exogenous, context, periods, ws_threaded)
    rout = A\res
    return 0
end

periods = 150
preconditioner_window = 10
algo = "CR"
m = context.models[1]
        n = m.endogenous_nbr
        a = Matrix{Float64}(undef, n, n)
        b = Matrix{Float64}(undef, n, n)
        c = Matrix{Float64}(undef, n, n)
        g = zeros(n, n)
        steadystate = Vector{Float64}(undef, n)
        steadystate_exo = Vector{Float64}(undef, m.exogenous_nbr)
        steadystate .= context.results.model_results[1].trends.endogenous_steady_state
        steadystate_exo .= context.results.model_results[1].trends.exogenous_steady_state
@btime        endogenous = repeat(steadystate, periods + 2)
@btime        exogenous = repeat(steadystate_exo',periods + 2)
        endogenous = repeat(steadystate, periods + 2)
        exogenous = repeat(steadystate_exo',periods + 2)
        lli = m.lead_lag_incidence
        dynamic_variables = zeros(nnz(sparse(lli)))
        temp_vec=Vector{Float64}(undef, n)
        work = context.work
@btime        LREWs = Dynare.LinearRationalExpectationsWs(
            algo, n, m.exogenous_nbr, m.exogenous_deterministic_nbr,
            m.i_fwrd_b, m.i_current, m.i_bkwrd_b, m.i_both, m.i_static)
LREWs = Dynare.LinearRationalExpectationsWs(
    algo, n, m.exogenous_nbr, m.exogenous_deterministic_nbr,
    m.i_fwrd_b, m.i_current, m.i_bkwrd_b, m.i_both, m.i_static)
@btime        LREresults = Dynare.LinearRationalExpectationsResults(n,
                                                       m.exogenous_nbr,
                                                       LREWs.backward_nbr)
        LREresults = Dynare.LinearRationalExpectationsResults(n,
                                                              m.exogenous_nbr,
                                                              LREWs.backward_nbr)
        options = context.options["stoch_simul"]
        if algo == "GS"
            options["generalized_schur"]["criterium"] = 1 + 1e-6
        else
            options["cyclic_reduction"] = Dict(["tol" => 1e-8])
        end

#@btime        Dynare.first_order_solver!(LREresults,
#                            algo,
#                            work.jacobian,
#                            options,
#                                         LREWs)
@show "first order"
results = LREresults
jacobian = work.jacobian
ws = LREWs
@btime Dynare.LinearRationalExpectations.remove_static!(jacobian, ws)
    if algo == "CR"
        @btime Dynare.LinearRationalExpectations.get_abc!(ws, jacobian)
        @btime Dynare.LinearRationalExpectations.cyclic_reduction!(ws.x, ws.c, ws.b, ws.a, ws.solver_ws, options["cyclic_reduction"]["tol"], 300)
        @btime begin
            vg = view(results.gs1, :, 1:ws.backward_nbr)
            vx = view(ws.x, ws.backward_indices_d, ws.backward_indices_d)
            copy!(vg, vx)
        end
        @btime begin
            vg = view(results.g1,ws.dynamic_indices, 1:ws.backward_nbr)
            vx = view(ws.x, :, ws.backward_indices_d)
            copy!(vg, vx)
        end
    elseif algo == "GS"
        @btime Dynare.LinearRationalExpectations.get_de!(ws, jacobian)
        @btime Dynare.LinearRationalExpectations.gs_solver!(ws.solver_ws, ws.d, ws.e, ws.backward_nbr, options["generalized_schur"]["criterium"])
        @btime results.gs1 .= ws.solver_ws.g1
        @btime for i = 1:ws.backward_nbr
            for j = 1:ws.backward_nbr
                x = ws.solver_ws.g1[j,i]
                results.g1[ws.backward_indices[j],i] = x
            end
            for j = 1:(ws.forward_nbr - ws.both_nbr)
                results.g1[ws.purely_forward_indices[j], i] =
                    ws.solver_ws.g2[ws.icolsE[ws.backward_nbr + j] - ws.backward_nbr, i]
            end
        end
    else
        error("Algorithm $algo not recognized")
    end
    if ws.static_nbr > 0
        Dynare.LinearRationalExpectations.add_static!(results, jacobian, ws)
    end
    #    A = view(jacobian, :, ws.backward_nbr + ws.current_nbr .+ (1:ws.forward_nbr))
    #    B = view(jacobian, :, ws.backward_nbr .+ ws.current_indices)
@btime begin
    vt = view(ws.temp8, :, 1:ws.forward_nbr)
    vj = view(jacobian, :, ws.backward_nbr + ws.current_nbr .+ (1:ws.forward_nbr))
    copy!(vt, vj)
end
@btime begin
    vt = view(ws.temp9, :, 1:ws.current_nbr)
    vj = view(jacobian, :, ws.backward_nbr .+ (1:ws.current_nbr))
    copy!(vt, vj)
end
    @btime Dynare.LinearRationalExpectations.make_lu_AGplusB!(ws.AGplusB, ws.temp8, results.g1_1, ws.temp9, ws)        
    @btime Dynare.LinearRationalExpectations.solve_for_derivatives_with_respect_to_shocks!(results, jacobian, ws)
@show "end first"
LB = LREWs.backward_indices
g1_1 = LREresults.g1_1
@btime        @inbounds for i = 1:LREWs.backward_nbr
    vg = view(g, :, LB[i])
    vg1 = view(g1_1, :, i)
    copy!(vg, vg1)
end
    if any(isnan.(g))
            throw(ArgumentError("NaN in g"))
        end
        residuals = zeros((periods+2)*n)
        presiduals = zeros(m.n_bkwrd + m.n_current + m.n_fwrd + 2*m.n_both)
@btime        ws_threaded = [Dynare.JacTimesVec(context) for i=1:Threads.nthreads()]
        ws_threaded = [Dynare.JacTimesVec(context) for i=1:Threads.nthreads()]
@btime        Dynare.get_jacobian!(ws_threaded[1], endogenous, exogenous, steadystate, m, 2)
        Dynare.get_jacobian!(ws_threaded[1], endogenous, exogenous, steadystate, m, 2)
        n = size(g, 1)
        if any(isnan.(ws_threaded[1].jacobian))
            throw(ArgumentError("NaN in jacobian"))
        end
@btime        @inbounds Dynare.get_abc!(a, b, c, ws_threaded[1].jacobian, m)
@btime        P = Dynare.LREprecond(periods, preconditioner_window, a, b, c, g)
P = Dynare.LREprecond(periods, preconditioner_window, a, b, c, g)
@btime        LREMap = Dynare.LinearMap(periods*n) do C, B
    @inbounds copyto!(residuals, n + 1, B, 1, periods*n)
    Dynare.jacobian_time_vec!(C, residuals, endogenous, exogenous,
                              steadystate, g, m, periods,
                              ws_threaded)
end
LREMap = Dynare.LinearMap(periods*n) do C, B
    @inbounds copyto!(residuals, n + 1, B, 1, periods*n)
    Dynare.jacobian_time_vec!(C, residuals, endogenous, exogenous,
                              steadystate, g, m, periods,
                              ws_threaded)
end
@btime        J = Dynare.Jacobian(context, periods)
    J = Dynare.Jacobian(context, periods)

#=
@btime bench1(150, 3)
@btime bench2(150, 3)
@btime bench3(150)
=#
