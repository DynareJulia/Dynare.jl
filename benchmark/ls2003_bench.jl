using BenchmarkTools
using Dynare

context = @dynare "../test/models/ls2003/ls2003.mod"

problem = Dynare.DSGELogPosteriorDensity(context, "../test/models/ls2003/data_ca1.csv", 8, 86)
transformation = Dynare.DSGETransformation(context.work.estimated_parameters)
transformed_problem = Dynare.TransformedLogDensity(transformation, problem)
transformed_density(θ) = -problem.f(collect(Dynare.TransformVariables.transform(transformation, θ)))
(p0, v0) = Dynare.get_initial_values(context.work.estimated_parameters)
ip0 = collect(Dynare.TransformVariables.inverse(transformation, tuple(p0...)))

@btime transformed_density(ip0)

@btime res = Dynare.optimize(
    transformed_density,
    $ip0,
    Dynare.LBFGS(),
)
res = Dynare.optimize(
    transformed_density,
    ip0,
    Dynare.LBFGS(),
    Dynare.Optim.Options(f_tol=1e-5)
)

(res, mode, tstdh, inv_hess) = Dynare.posterior_mode(context, datafile="../test/models/ls2003/data_ca1.csv", first_obs=8, last_obs=86, iterations=1000, show_trace=false)

model = Dynare.DensityModel(transformed_density)

spl = Dynare.RWMH(Dynare.MvNormal(zeros(length(mode)), Matrix(inv_hess)))

# Sample from the posterior.
@btime chain = Dynare.sample(model, spl, 10000,
                             init_params = mode,
                             param_names = context.work.estimated_parameters.name,
                             chain_type = Dynare.Chains)


