module PardisoSolver

using Dynare, MKL, Pardiso

ps = MKLPardisoSolver()

@show "PardisoSolver"

Dynare.linear_solver!(::Dynare.PardisoLS, x, A, b) = Pardiso.solve!(ps, x, A, b)

end
