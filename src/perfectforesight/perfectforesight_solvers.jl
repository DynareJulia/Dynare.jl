using FastLapackInterface
using FastLapackInterface.LinSolveAlgo
using IterativeSolvers
using LinearAlgebra
using LinearAlgebra.BLAS
using LinearMaps
using LinearRationalExpectations
import Base.\
import LinearAlgebra.mul!
import LinearAlgebra.ldiv!
import LinearAlgebra.BLAS: @blasfunc, BlasInt, BlasFloat, libblas
using LinearRationalExpectations
using SparseArrays

include("makeA.jl")
include("gmres_solver.jl")
