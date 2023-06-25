module PardisoSolver

using SparseArrays
using Dynare, MKL, Pardiso

ps = MKLPardisoSolver()

@show "PardisoSolver"

Dynare.linear_solver!(::Dynare.PardisoLS, x::Vector{Float64}, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}) = Pardiso.solve!(ps, x, A, b)

end
