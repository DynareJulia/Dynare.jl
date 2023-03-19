# test estimation code agains DGSE Estimation.zip (Herbst Schorfheide, 2014 and using Dynare
using Test

context = Dynare.parser("nk", Dynare.CommandLineOptions())

symboltable = context.symboltable
varobs = context.work.observed_variables
has_trends = context.modfileinfo.has_trends
varobs_ids =
    [symboltable[v].orderintype for v in varobs if Dynare.is_endogenous(v, symboltable)]
model = context.models[1]
results = context.results.model_results[1]
lre_results = results.linearrationalexpectations
estimation_results = results.estimation

options = Dynare.EstimationOptions()

Yorig =
    Dynare.get_data("dsge1_data.csv", varobs, start = options.first_obs, last = options.last_obs)

observations = copy(Yorig)
nobs = size(observations, 2)
ssws = Dynare.SSWs(context, nobs, varobs)
ssws.Q .= model.Sigma_e
ssws.H .= context.work.Sigma_m
estimated_parameters = context.work.estimated_parameters

initial_values = Dynare.get_initial_value_or_mean(estimated_parameters)

Dynare.set_estimated_parameters!(context, initial_values)

@show initial_values
@show context.work.params
@show context.models[1].Sigma_e
@show ssws.Q

ll = Dynare.loglikelihood(initial_values, context, observations, ssws)
@show ssws.kalman_ws.v

@test isapprox(ll, -301.0627, atol=0.01) 

ep = context.work.estimated_parameters
@show initial_values
lprior = Dynare.logpriordensity(initial_values, ep)
@test isapprox(lprior, -11.9243, atol=0.01)

problem = Dynare.DSGELogPosteriorDensity(context, observations, 1, nobs)
@show problem(initial_values)
    
transformation = Dynare.DSGETransformation(ep)
transformed_problem = Dynare.TransformedLogDensity(transformation, problem)
transformed_density(θ) = -problem.f(collect(Dynare.TransformVariables.transform(transformation, θ)))
transformed_density_gradient!(g, θ) = (g = finite_difference_gradient(transformed_density, θ))

pit = Dynare.TransformVariables.inverse(transformation, tuple(initial_values...))
@show pit

pt = Dynare.TransformVariables.transform(transformation, pit)
@show pt

#Dynare.posterior_mode(context, observations)
