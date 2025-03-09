# %%
using Dynare
using Test
using Tasmanian
using Statistics
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
# # %%
# # LHS Graph: Plotting 
# nb_ddsg = Vector{Int}()
# errors_ddsg = Vector{Float64}()
# nb_sg = Vector{Int}()
# errors_sg = Vector{Float64}()
# # %%
# c = 1
# num_points = 1000
# F_poly(X::AbstractVector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
# F_poly(X::AbstractMatrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
# F(X::Vector{Float64}) = F_poly(X, c=c)
# F(X::Matrix{Float64}) = F_poly(X, c=c)
# dim = 20
# k_max = 1
# dof = 1
# X_sample = Matrix{Float64}(undef, dim, num_points)
# for i in 1:num_points
#     X_sample[:, i] .= ddsg.domain[:, 1] .+ (ddsg.domain[:, 2] - ddsg.domain[:, 1]).*rand(dim)
# end
# Y_exact = F(X_sample)
# for l in 1:5
#     ddsg = Dynare.DDSG(dim, dof, l, l, k_max)
#     Dynare.DDSG_init!(ddsg)
#     Dynare.DDSG_build!(
#         ddsg,
#         F,
#         ddsg.centroid
#     )
#     # %%
#     Y_ddsg = Dynare.DDSG_evaluate(ddsg, X_sample)
#     push!(errors_ddsg, mean(@. abs((Y_ddsg - Y_exact)/Y_exact)))
#     push!(nb_ddsg, sum(ddsg.grid_points))
#     # %%
#     sg = Tasmanian.TasmanianSG(dim, dof, l)
#     Tasmanian.makeLocalPolynomialGrid!(sg)
#     Tasmanian.setDomainTransform!(sg, ddsg.domain)  # Match grid to domain dimensions
#     # Refinement based on surplus coefficients
#     X = Tasmanian.getNeededPoints(sg)
#     N = size(X,2)
#     Y_val = F(X)
#     Tasmanian.loadNeededPoints!(sg, Y_val)
#     Y_sg = Tasmanian.evaluateBatch(sg,X_sample)
#     # %%
#     push!(errors_sg, mean(@. abs((Y_sg - Y_exact)/Y_exact)))
#     push!(nb_sg, Tasmanian.getNumPoints(sg))
# end
# # %%
# using LaTeXStrings
# lhs = plot(nb_sg, errors_sg, color = :black, label=L"SG",
#      markershape=:star6, markersize=5, linewidth = 4)
# plot!(nb_ddsg, errors_ddsg, color = :blue, label=L"DDSG~(\mathcal{K} = 1)",
#       markershape=:circle, markersize=5, linewidth = 4)
# plot!(xscale=:log10, yscale=:log10, minorgrid=true)
# xlabel!("Number of Grid Points")
# ylabel!("Error")

# # %%
# # RHS Graph: Plotting 
# nb_ddsg = Matrix{Int}(undef, 5, 3)
# errors_ddsg = Matrix{Float64}(undef, 5, 3)
# nb_sg = Vector{Int}()
# errors_sg = Vector{Float64}()
# # %%
# c = 3
# num_points = 1000
# F_poly(X::AbstractVector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
# F_poly(X::AbstractMatrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
# F(X::Vector{Float64}) = F_poly(X, c=c)
# F(X::Matrix{Float64}) = F_poly(X, c=c)
# dim = 20
# dof = 1
# X_sample = Matrix{Float64}(undef, dim, num_points)
# for i in 1:num_points
#     X_sample[:, i] .= ddsg.domain[:, 1] .+ (ddsg.domain[:, 2] - ddsg.domain[:, 1]).*rand(dim)
# end
# Y_exact = F(X_sample)
# # %%
# for l in 1:5
#     for k_max in 1:3
#         ddsg = Dynare.DDSG(dim, dof, l, l, k_max)
#         Dynare.DDSG_init!(ddsg)
#         Dynare.DDSG_build!(
#             ddsg,
#             F,
#             ddsg.centroid
#         )
#         # %%
#         Y_ddsg = Dynare.DDSG_evaluate(ddsg, X_sample)
#         errors_ddsg[l,k_max] = mean(@. abs((Y_ddsg - Y_exact)/Y_exact))
#         nb_ddsg[l,k_max] =  sum(ddsg.grid_points)
#     end
#     # %%
#     sg = Tasmanian.TasmanianSG(dim, dof, l)
#     Tasmanian.makeLocalPolynomialGrid!(sg)
#     Tasmanian.setDomainTransform!(sg, ddsg.domain)  # Match grid to domain dimensions
#     # Refinement based on surplus coefficients
#     X = Tasmanian.getNeededPoints(sg)
#     N = size(X,2)
#     Y_val = F(X)
#     Tasmanian.loadNeededPoints!(sg, Y_val)
#     Y_sg = Tasmanian.evaluateBatch(sg,X_sample)
#     # %%
#     push!(errors_sg, mean(@. abs((Y_sg - Y_exact)/Y_exact)))
#     push!(nb_sg, Tasmanian.getNumPoints(sg))
# end
# # %%
# using LaTeXStrings
# rhs = plot(nb_sg, errors_sg, color = :black, label=L"SG",
#      markershape=:star6, markersize=5, linewidth = 4)
# plot!(nb_ddsg[:,1], errors_ddsg[:,1], color = :blue, label=L"DDSG~(\mathcal{K} = 1)",
#       markershape=:circle, markersize=5, linewidth = 4)
# plot!(nb_ddsg[:,2], errors_ddsg[:,2], color = :red, label=L"DDSG~(\mathcal{K} = 2)",
#       markershape=:xcross, markersize=5, linewidth = 4, linestyle = :dash)
# plot!(nb_ddsg[:,3], errors_ddsg[:,3], color = :green, label=L"DDSG~(\mathcal{K} = 3)",
#       markershape=:rect, markersize=5, linewidth = 4, linestyle = :dot)
# plot!(xscale=:log10, yscale=:log10, minorgrid=true)
# xlabel!("Number of Grid Points")
# # %%
# plot(lhs, rhs, layout=(1, 2), link=:y)
# savefig("analytical_DDSG.pdf")
# # %%
# using Profile, Profile.Allocs
# Profile.Allocs.clear()  # Clear old profiling data
# Profile.Allocs.@profile Dynare.DDSG_evaluate(ddsg, X_sample)
# # %%
# --track-allocation=user
# # %%
# # # F_poly(X::Vector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
# # # F_poly(X::Matrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
# # F_poly(X::Vector{Float64}; c)::Vector{Float64} = [sum(sin, X)^c]
# # F_poly(X::Matrix{Float64}; c)::Matrix{Float64} = reduce(hcat, [F_poly(col;c=c) for col in eachcol(X)])
# # # %%
# # dim = 3
# # k_max = 2
# # c = 1
# # dof = 1
# # l = 2
# # l_min = l
# # l_max = l
# # num_points = 2
# # F(X::Vector{Float64}) = F_poly(X, c=c)
# # F(X::Matrix{Float64}) = F_poly(X, c=c)
# # # %%
# # ddsg = Dynare.DDSG(dim, dof, l, l, k_max)
# # Dynare.DDSG_init!(ddsg)
# # # %%
# # Dynare.DDSG_build!(
# #     ddsg,
# #     F,
# #     ddsg.centroid
# # )
# # # %%
# # v = [3]
# # X = zeros(1,5)
# # X[1,:] .= [0.5,0.0,1.0,0.25,0.75]
# # Dynare.evaluate_at_̄x_xᵥ(F, ddsg.X0, X, v)
# # # %%
# # @code_warntype Dynare.evaluate_at_̄x_xᵥ(F, ddsg.X0, X, v)
# # # %%
# # @code_warntype Dynare.DDSG_build!(
# #     ddsg,
# #     F,
# #     ddsg.centroid
# # )
# # # %%
# # @code_warntype Dynare.DDSG_refine!(ddsg.grid[2], 0., "classic", [])
# # # %%

# # # %%
# # X_sample = Matrix{Float64}(undef, dim, num_points)
# # for i in 1:num_points
# #     X_sample[:, i] .= grid.domain[:, 1] .+ (grid.domain[:, 2] - grid.domain[:, 1]) .* rand(dim)
# # end

# # # Y = DDSG_evaluate(grid,X=X_sample)
# # # num_points = sum(grid.grid_points)
# # #

# # # @code_warntype Dynare.DDSG_build!(
# # #     ddsg,
# # #     F=F,
# # #     X0=ddsg.centroid,
# # #     refinement_tol = 0.
