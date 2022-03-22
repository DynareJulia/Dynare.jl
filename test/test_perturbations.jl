using Test
using Dynare
using LinearAlgebra

@testset "Perturbations" begin

    n = 10 
    Sigma_e = rand(n,n)
    Sigma_e = transpose(Sigma_e)*Sigma_e
    sd = sqrt.(diag(Sigma_e))
    C = similar(Sigma_e)
    Dynare.correlation!(C, Sigma_e, sd)
    @test C .* sd .* transpose(sd) ≈ Sigma_e

    C = Dynare.correlation(Sigma_e)
    @test C .* sd .* transpose(sd) ≈ Sigma_e

    order = 4
    A = rand(n,n)
    AA = [zeros(n,n) for i = 1:order]
    Dynare.autocovariance!(AA, A, Sigma_e, zeros(n,n), zeros(n,n), order)
    @test AA[1]  ≈ A*Sigma_e
    @test AA[2]  ≈ A*A*Sigma_e
    @test AA[3]  ≈ A*A*A*Sigma_e
    @test AA[4]  ≈ A*A*A*A*Sigma_e
   
    VV = [zeros(n) for i = 1:order]
    Dynare.autocovariance!(VV, A, Sigma_e, zeros(n,n), zeros(n,n), order)
    @test VV[1]  ≈ diag(A*Sigma_e)
    @test VV[2]  ≈ diag(A*A*Sigma_e)
    @test VV[3]  ≈ diag(A*A*A*Sigma_e)
    @test VV[4]  ≈ diag(A*A*A*A*Sigma_e)
   
end
