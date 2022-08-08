using BenchmarkTools
using Dynare

#context = @dynare "../test/models/fs2000/fs2000_bench.jl"

problem = Dynare.DSGENegativeLogPosteriorDensity(context, "../test/models/fs2000/fsdata_simul.csv", 1, 0)
transformation = Dynare.DSGETransformation(context.work.estimated_parameters)
transformed_problem = Dynare.TransformedLogDensity(transformation, problem)
transformed_density(θ) = problem.f(collect(Dynare.TransformVariables.transform(transformation, θ)))
(p0, v0) = Dynare.get_initial_values(context.work.estimated_parameters)
ip0 = collect(Dynare.TransformVariables.inverse(transformation, tuple(p0...)))

@btime problem.f(p0)
@btime transformed_density(ip0)
