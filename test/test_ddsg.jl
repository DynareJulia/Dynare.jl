# %%
import Pkg; Pkg.activate("./Dynare.jl/")
include("../src/Dynare.jl")
using .Dynare
using Test
using Tasmanian
using Statistics
# %%
context=Dynare.dynare("rbc.mod", "stoponerror")
# %%
α     = 1/3;
β     = 0.99;
ρ     = 0.95;
σ     = 0.01;
ss_l  = 1/3;
ss_k  = ss_l*(β*α)^(1/(1-α));
ss_y  = ss_k^α*ss_l^(1-α);
ss_c  = ss_y-ss_k;
ss_rk = α*ss_y/ss_k;
θ = (1-α)*(1-ss_l)*ss_y/(ss_l*ss_c);
# %%
c_pol(K,Z) = (1-α*β)*exp(Z)*K^α*ss_l^(1-α)
k_pol(K,Z) = α*β*exp(Z)*K^α*ss_l^(1-α)
y_pol(K,Z) = c_pol(K,Z)+k_pol(K,Z)
# %%
# initialPolGuess = Dynare.UserPolicyGuess(
#     function (x)
#         Z = x[1]
#         K = x[2]
#         return [
#             0.15,
#             0.05,
#             0.20,
#             1.
#         ]
#     end,
#     ["z", "k"],
#     ["c","k","y","rk"]
# )
# %%
# # Test of the user-provided initial policy guess
# initialPolGuess = Dynare.UserPolicyGuess(
#     function (x)
#         Z = x[1]
#         K = x[2]
#         return [
#             c_pol(K,Z),
#             k_pol(K,Z),
#             y_pol(K,Z),
#             α*y_pol(K,Z)/K
#         ]
#     end,
#     ["z", "k"],
#     ["c","k","y","rk"]
# )
# %%
# RBC Model
(ddsg, sgws) = Dynare.DDSGapproximation(tol_ti=1e-6,gridDepth=3,ftol=1e-8, polUpdateWeight=1., maxRef=0, k_max=1, surplThreshold=0.0);
# (ddsg, sgws) = Dynare.DDSGapproximation(tol_ti=1e-6,gridDepth=3, polUpdateWeight=1., maxRef=0, k_max=1, surplThreshold=0.0, show_trace=false);
# %%
# (SG_grid, sgws) = Dynare.sparsegridapproximation(tol_ti=1e-6,gridDepth=3,maxRef=0, polUpdateWeight=1., surplThreshold=0.0);
# %%
Z_vals = range(-2*σ/sqrt(1-ρ^2), 2*σ/sqrt(1-ρ^2), length=100)
K_vals = range(ss_k * 0.1, ss_k * 1.9, length=100) 
# %%
using Plots
using LaTeXStrings
c_theoretical = [c_pol(K, Z) for Z in Z_vals, K in K_vals]
pol_ddsg_1 = x->Dynare.interpolate(ddsg, x)
i_c = context.symboltable["c"].orderintype
i_c = findall(i_c .== context.models[1].i_dyn)[1]
c_ddsg_1 = [pol_ddsg_1([K, Z])[i_c] for Z in Z_vals, K in K_vals]
plot(K_vals, c_theoretical[50, :], label="Theoretical", title="Consumption Policy")
plot!(K_vals, c_ddsg_1[50, :], label="DDSG (K = 1)", linestyle=:dash)
ylabel!(L"C_t \mid Z_t = Z")
# %%
# RBC Model
# (ddsg, sgws) = Dynare.DDSGapproximation(tol_ti=1e-6,gridDepth=3, polUpdateWeight=0., maxRef=0, k_max=1, surplThreshold=0.);
# %%
@testset "DDSG" begin
    @testset "Constructors" begin
        @testset "Basic" begin
            dim=10
            k_max=10
            dof = 10
            l = 2
            d = zeros(dim,2);
            d[:,2] .= 1.
            centroids = 0.5*ones(dim)
            ddsg = Dynare.DDSG(dim,dof,l,l,k_max)
            @test ddsg.dim == dim
            @test ddsg.dof == dof
            @test ddsg.l_min == l
            @test ddsg.l_max == l
            @test ddsg.order == 1
            @test ddsg.rule == "localp"
            @test ddsg.domain == d
            @test ddsg.centroid == centroids
            @test length(ddsg.lookup) == k_max
            @test isempty(ddsg.grid)
            @test isempty(ddsg.coeff)
            @test isempty(ddsg.grid_points)
            @test isempty(ddsg.X0)
            @test isempty(ddsg.Y0)
            @test ddsg.coeff0 == 0.0
            @test !ddsg.is_ddsg
            k_max = dim-1
            ddsg = Dynare.DDSG(dim,dof,l,l,k_max)
            @test ddsg.is_ddsg
        end
        @testset "Copy" begin
            dim = 10
            k_max = 10
            dof = 10
            l = 2
            d = zeros(dim,2)
            d[:,2] .= 1.0
            centroids = 0.5 * ones(dim)
            ddsg = Dynare.DDSG(dim,dof,l,l,k_max)
            ddsg_copy = Dynare.DDSG(ddsg)
            @testset "Copied values match original" begin
                @test ddsg_copy.dim == ddsg.dim
                @test ddsg_copy.dof == ddsg.dof
                @test ddsg_copy.l_min == ddsg.l_min
                @test ddsg_copy.l_max == ddsg.l_max
                @test ddsg_copy.k_max == ddsg.k_max
                @test ddsg_copy.order == ddsg.order
                @test ddsg_copy.rule == ddsg.rule
                @test ddsg_copy.domain == ddsg.domain
                @test ddsg_copy.centroid == ddsg.centroid
                @test ddsg_copy.lookup == ddsg.lookup
                @test ddsg_copy.coeff == ddsg.coeff
                @test ddsg_copy.grid_points == ddsg.grid_points
                @test ddsg_copy.X0 == ddsg.X0
                @test ddsg_copy.Y0 == ddsg.Y0
                @test ddsg_copy.coeff0 == ddsg.coeff0
                @test ddsg_copy.is_ddsg == ddsg.is_ddsg
            end
            @testset "Shallow and deep copied fields (should have same content but independent memory)" begin
                @test ddsg_copy.domain !== ddsg.domain
                @test ddsg_copy.centroid !== ddsg.centroid
                @test ddsg_copy.lookup !== ddsg.lookup
                @test ddsg_copy.grid !== ddsg.grid
                @test ddsg_copy.coeff !== ddsg.coeff
                @test ddsg_copy.grid_points !== ddsg.grid_points
                @test ddsg_copy.X0 !== ddsg.X0
                @test ddsg_copy.Y0 !== ddsg.Y0
            end
            @testset "Deep copied fields (should have independent memory for the inner elements too)" begin
                # Lookup dictionary should be deep copied, so modifying the copy shouldn't affect original
                push!(ddsg_copy.lookup[3], ([1, 2, 3] => 42))
                @test ([1, 2, 3] in keys(ddsg_copy.lookup[3]))
                @test ([1, 2, 3] ∉ keys(ddsg.lookup[3]))  # Original should not be affected
            end
        end
    end
    @testset "Initialization" begin
        @testset "Uncut DDSG" begin
            dim = 2
            k_max = 2
            dof = 2
            l = 2
            ddsg = Dynare.DDSG(dim,dof,l,l,k_max)
            Dynare.DDSG_init!(ddsg)
            # Vectors
            @test length(ddsg.coeff) == 1
            @test length(ddsg.grid_points) == 1
            @test length(ddsg.grid) == 1
            @test ddsg.coeff == ones(1)
            @test ddsg.grid_points == zero(ddsg.grid_points)
            # Lookup tables
            @test ([1,2] in keys(ddsg.lookup[dim]))
            @test isempty(ddsg.lookup[1])
        end
        @testset "Cut DDSG" begin
            dim = 3
            k_max = 2
            dof = 2
            l = 2
            ddsg = Dynare.DDSG(dim,dof,l,l,k_max)
            Dynare.DDSG_init!(ddsg)
            nb = reduce(+,[binomial(dim,k) for k in 1:k_max])
            # Vectors
            @test length(ddsg.grid) == nb
            @test length(ddsg.coeff) == nb
            @test length(ddsg.grid_points) == nb
            @test length(ddsg.lookup) == k_max
            # Lookup tables
            @test ddsg.coeff == zeros(nb)
            @test ddsg.grid_points == zero(ddsg.grid_points)
            @test all([[k] in keys(ddsg.lookup[1]) for k in 1:k_max])
            @test all(length(values(ddsg.lookup[k])) == binomial(dim,k) for k in 1:k_max)
        end
    end
    @testset "Build and Evaluate" begin
        # %%
        dim = 20
        k_max = 1
        c = 1
        dof = 1
        l = 5
        F_poly(X::AbstractVector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
        F_poly(X::AbstractMatrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
        F(X::Vector{Float64}) = F_poly(X, c=c)
        F(X::Matrix{Float64}) = F_poly(X, c=c)
        ddsg = Dynare.DDSG(dim, dof, l, l, k_max)
        Dynare.DDSG_init!(ddsg)
        Dynare.DDSG_build!(
            ddsg,
            F,
            ddsg.centroid
        )
        @testset "Evaluate" begin
            num_points = 1000
            X_sample = rand(dim, num_points)
            Y_ddsg = Dynare.DDSG_evaluate(ddsg, X_sample)
            Y_exact = F(X_sample)
            @test sum(ddsg.grid_points) == reduce(+,[Tasmanian.getNumPoints(g) for g in ddsg.grid])
            @test mean(@. abs((Y_ddsg - Y_exact)/Y_exact)) < 1e-4
        end
        @testset "Build Helper Functions" begin
            v = [3]
            X = zeros(1,5)
            X[1,:] .= [0.5,0.0,1.0,0.25,0.75]
            Y = Dynare.evaluate_at_̄x_xᵥ(F, ddsg.X0, X, v)
            Dynare.update_ddsg_coefficients!(ddsg)
        end
    end
end
# %%