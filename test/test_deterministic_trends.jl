using Dynare
using Test


n = 3
steady_state = [1, 1.5, 3]
linear_coefs = [0.05, 0, -0.02]
quadratic_coefs = [0.01, 0, 0.005]

x0 = ones(n, n)
y = zeros(n, n)
trend = collect(0:n-1)
@testset "add trend" begin
    @testset "steady state" begin
        Dynare.add_steady_state!(y, x0, steady_state, 1)
        @test y ≈ 1 .+ repeat(steady_state, 1, n)
        Dynare.add_steady_state!(y, x0, steady_state, 2)
        @test y ≈ 1 .+ repeat(steady_state', n)
        copy!(y, x0)
        Dynare.add_steady_state!(y, steady_state, 1)
        @test y ≈ 1 .+ repeat(steady_state, 1, n)
        copy!(y, x0)
        Dynare.add_steady_state!(y, steady_state, 2)
        @test y ≈ 1 .+ repeat(steady_state', n)
    end
    @testset "linear trend" begin
        Dynare.add_linear_trend!(y, x0, steady_state, linear_coefs, 1)
        @test y ≈ 1 .+ repeat(steady_state) .+ linear_coefs .* trend'
        Dynare.add_linear_trend!(y, x0, steady_state, linear_coefs, 2)
        @test y ≈ 1 .+ repeat(steady_state', n) .+ linear_coefs' .* trend
        copy!(y, x0)
        Dynare.add_linear_trend!(y, steady_state, linear_coefs, 1)
        @test y ≈ 1 .+ repeat(steady_state, 1, n) .+ linear_coefs .* trend'
        copy!(y, x0)
        Dynare.add_linear_trend!(y, steady_state, linear_coefs, 2)
        @test y ≈ 1 .+ repeat(steady_state', n) .+ linear_coefs' .* trend
    end
    @testset "quadratic trend" begin
        Dynare.add_quadratic_trend!(y, x0, steady_state, linear_coefs, quadratic_coefs, 1)
        @test y ≈ 1 .+ repeat(steady_state, 1, n) .+ linear_coefs .* trend' .+ quadratic_coefs .* (trend.^2)'
        Dynare.add_quadratic_trend!(y, x0, steady_state, linear_coefs, quadratic_coefs, 2)
        @test y ≈ 1 .+ repeat(steady_state', n) .+ linear_coefs' .* trend .+ quadratic_coefs' .* trend.^2
        copy!(y, x0)
        Dynare.add_quadratic_trend!(y, steady_state, linear_coefs, quadratic_coefs, 1)
        @test y ≈ 1 .+ repeat(steady_state, 1, n) .+ linear_coefs .* trend' .+ quadratic_coefs .* (trend.^2)'
        copy!(y, x0)
        Dynare.add_quadratic_trend!(y, steady_state, linear_coefs, quadratic_coefs, 2)
        @test y ≈ 1 .+ repeat(steady_state', n) .+ linear_coefs' .* trend .+ quadratic_coefs' .* trend.^2
    end
end

@testset "remove trend" begin
    @testset "steady state" begin
        Dynare.remove_steady_state!(y, x0, steady_state, 1)
        @test y ≈ 1 .- repeat(steady_state, 1, n)
        Dynare.remove_steady_state!(y, x0, steady_state, 2)
        @test y ≈ 1 .- repeat(steady_state', n)
        copy!(y, x0)
        Dynare.remove_steady_state!(y, steady_state, 1)
        @test y ≈ 1 .- repeat(steady_state, 1, n)
        copy!(y, x0)
        Dynare.remove_steady_state!(y, steady_state, 2)
        @test y ≈ 1 .- repeat(steady_state', n)
    end
    @testset "linear trend" begin
        Dynare.remove_linear_trend!(y, x0, steady_state, linear_coefs, 1)
        @test y ≈ 1 .- repeat(steady_state) .- linear_coefs .* trend'
        Dynare.remove_linear_trend!(y, x0, steady_state, linear_coefs, 2)
        @test y ≈ 1 .- repeat(steady_state', n) .- linear_coefs' .* trend
        copy!(y, x0)
        Dynare.remove_linear_trend!(y, steady_state, linear_coefs, 1)
        @test y ≈ 1 .- repeat(steady_state, 1, n) .- linear_coefs .* trend'
        copy!(y, x0)
        Dynare.remove_linear_trend!(y, steady_state, linear_coefs, 2)
        @test y ≈ 1 .- repeat(steady_state', n) .- linear_coefs' .* trend
    end
    @testset "quadratic trend" begin
        Dynare.remove_quadratic_trend!(y, x0, steady_state, linear_coefs, quadratic_coefs, 1)
        @test y ≈ 1 .- repeat(steady_state, 1, n) .- linear_coefs .* trend' .- quadratic_coefs .* (trend.^2)'
        Dynare.remove_quadratic_trend!(y, x0, steady_state, linear_coefs, quadratic_coefs, 2)
        @test y ≈ 1 .- repeat(steady_state', n) .- linear_coefs' .* trend .- quadratic_coefs' .* trend.^2
        copy!(y, x0)
        Dynare.remove_quadratic_trend!(y, steady_state, linear_coefs, quadratic_coefs, 1)
        @test y ≈ 1 .- repeat(steady_state, 1, n) .- linear_coefs .* trend' .- quadratic_coefs .* (trend.^2)'
        copy!(y, x0)
        Dynare.remove_quadratic_trend!(y, steady_state, linear_coefs, quadratic_coefs, 2)
        @test y ≈ 1 .- repeat(steady_state', n) .- linear_coefs' .* trend .- quadratic_coefs' .* trend.^2
    end
end


