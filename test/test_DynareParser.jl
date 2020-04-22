using Dynare


context = Dynare.Context()

parser("models/example1/example1", context)
println(context)
steady_state!(context)
println(context.results)
