module PardisoSolver

using SparseArrays
using Dynare, Pardiso

ps = MKLPardisoSolver()

Dynare.linear_solver!(::Dynare.PardisoLS, x::Vector{Float64}, A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}) = Pardiso.solve!(ps, x, A, b)

end
