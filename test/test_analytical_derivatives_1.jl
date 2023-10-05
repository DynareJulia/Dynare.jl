include("../src/derivatives.jl")

using Dynare #Use analytical_derivatives branch
using LinearAlgebra
using Dynare.FastLapackInterface
using Test
using SparseArrays
using SuiteSparse
using GeneralizedSylvesterSolver
using FiniteDifferences

context = @dynare "models/analytical_derivatives/fs2000_sa.mod" "params_derivs_order=1" "notmpterms";

model = context.models[1]
endo_nbr = model.endogenous_nbr
exo_nbr = model.exogenous_nbr
param_nbr = model.parameter_nbr
dss_dp = zeros(endo_nbr, param_nbr)

# 1 - derivatives of steady state with respect to parameter
wss = Dynare.StaticWs(context)
params = context.work.params
results = context.results.model_results[1]
steadystate = results.trends.endogenous_steady_state
exogenous = results.trends.exogenous_steady_state 

dss_dx = Dynare.get_static_jacobian!(wss, params, steadystate, exogenous, model) 
(dfstatic_dp, gp) = Dynare.DFunctions.SparseStaticParametersDerivatives!(steadystate, exogenous, params)

SSDerivatives!(dss_dp, dss_dx, dfstatic_dp)

function fss(p)
    Dynare.evaluate_steady_state!(steadystate, exogenous, p)
    return steadystate
end

fd = central_fdm(5, 1)
dss_dp_target = FiniteDifferences.jacobian(fd, fss, params)[1]
@test dss_dp ≈ dss_dp_target

# 2 - derivatives of linearized coefficients with respect to parameters

wsd = Dynare.DynamicWs(context, order=2)
endogenous3 = repeat(steadystate, 3)
jacobian = Dynare.get_dynamic_jacobian!(
    wsd,
    params,
    endogenous3,
    exogenous,
    steadystate,
    model,
    2,
)

dA_dp, dB_dp, dC_dp = ABCderivatives(jacobian, context, wsd)

function get_jacobian_(params, context)
    context1 = deepcopy(context)
    context1.work.params .= params
    Dynare.compute_steady_state!(context1)
    steadystate = copy(context1.results.model_results[1].trends.endogenous_steady_state)
    wsd = Dynare.DynamicWs(context1, order=2)
    endogenous3 = repeat(steadystate, 3)
    exogenous = context.results.model_results[1].trends.exogenous_steady_state
    jacobian = Dynare.get_dynamic_jacobian!(
        wsd,
        params,
        endogenous3,
        exogenous,
        steadystate,
        model,
        2,
    )
    return jacobian
end

get_jacobian(params) = get_jacobian_(params, context)

fd = central_fdm(5, 1)
dJ_dp_target = FiniteDifferences.jacobian(fd, get_jacobian, params)[1]

for i = 1:param_nbr
    target = reshape(dJ_dp_target[:, i], endo_nbr, 3*endo_nbr + exo_nbr)
    @test isapprox(dA_dp[:, :, i], target[:, 2*endo_nbr .+ (1:endo_nbr)], atol = 0.05, rtol = 0.05)
    @test isapprox(dB_dp[:, :, i], target[:, endo_nbr .+ (1:endo_nbr)], atol = 0.05, rtol = 0.05)
    @test isapprox(dC_dp[:, :, i], target[:, (1:endo_nbr)], atol = 0.05, rtol = 0.05)
end

localapproximation!(irf=0, display=false)
g1 = context.results.model_results[1].linearrationalexpectations.g1_1
n = model.endogenous_nbr
params = context.work.params
jacobian = get_jacobian_(params, context)
A = jacobian[:, 2*n .+ (1:n)]
B = jacobian[:, n .+ (1:n)]
C = jacobian[: , 1:n]
state_indices = model.i_bkwrd_b
param_nbr = model.parameter_nbr

dX = d1_lre_solution(dA_dp, dB_dp, dC_dp, g1, A, B, C, state_indices, n, param_nbr)

index=1
function funX(param)
    context.work.params[index] = param
    localapproximation!(irf=0, display=false);
    return copy(context.results.model_results[1].linearrationalexpectations.g1_1)
end

fd = central_fdm(5, 1)
for i = 1:param_nbr
    global index = i
    old_param = params[i]
    dXi_target = FiniteDifferences.jacobian(fd, funX, params[i])
    context.work.params[i] = old_param
    @test isapprox(dX[i][:, state_indices], reshape(dXi_target[1], n, 4), atol = 0.01, rtol = 0.01)
end
