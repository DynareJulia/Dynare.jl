using AxisArrayTables
using Dynare
using ExtendedDates
using JLD2
using LinearAlgebra
using Test

@dynare "models/example1pf/example1pf_conditional"
context = load("models/example1pf/example1pf_conditional/output/example1pf_conditional.jld2", "context")
function make_f_J(context, scenario, periods)
    datafile = ""
    
    FI = Dynare.FlipInformation(context, periods, 1)
    m = context.models[1]
    ncol = m.n_bkwrd + m.n_current + m.n_fwrd + 2 * m.n_both
    tmp_nbr =  m.dynamic_tmp_nbr::Vector{Int64}
    dynamic_ws = Dynare.DynamicWs(m.endogenous_nbr,
                                  m.exogenous_nbr,
                                  sum(tmp_nbr[1:2]),
                                  m.dynamic_g1_sparse_colptr,
                                  m.dynamic_g1_sparse_rowval)

    perfect_foresight_ws = Dynare.PerfectForesightWs(context, periods)
    X = perfect_foresight_ws.shocks
    initial_values = Dynare.get_dynamic_initialvalues(context)
    terminal_values = Dynare.get_dynamic_terminalvalues(context, periods)
    guess_values = Dynare.perfect_foresight_initialization!(
        context,
        terminal_values,
        periods,
        datafile,
        X,
        perfect_foresight_ws,
        Dynare.steadystate,
        dynamic_ws,
    )

    results = context.results.model_results[1]
    work = context.work
    dynamic_variables = dynamic_ws.dynamic_variables
    temp_vec = dynamic_ws.temporary_values
    steadystate = results.trends.endogenous_steady_state
    params = work.params

    residuals = zeros(periods * m.endogenous_nbr)
    JJ = perfect_foresight_ws.J
    
    exogenous = perfect_foresight_ws.x
    ws_threaded = [
        Dynare.DynamicWs(
            m.endogenous_nbr,
            m.exogenous_nbr,
            length(temp_vec),
            m.dynamic_g1_sparse_colptr,
            m.dynamic_g1_sparse_rowval
        ) for i = 1:Threads.nthreads()
    ]
    nzval = similar(m.dynamic_g1_sparse_rowval, Float64)
    nzval1 = similar(nzval)

    f! = Dynare.make_pf_residuals(initial_values,
                                  terminal_values,
                                  exogenous,
                                  dynamic_variables,
                                  steadystate,
                                  params,
                                  m,
                                  periods,
                                  temp_vec,
                                  perfect_foresight_ws.permutationsR,
                                  FI.ix_stack
                                  )
    permutations = []
    J! = Dynare.make_pf_jacobian(Dynare.DFunctions.dynamic_derivatives!,
                                 initial_values,
                                 terminal_values,
                                 dynamic_variables,
                                 exogenous,
                                 periods,
                                 temp_vec,
                                 params,
                                 steadystate,
                                 m.dynamic_g1_sparse_colptr,
                                 nzval,
                                 m.endogenous_nbr,
                                 m.exogenous_nbr,
                                 permutations,
                                 FI,
                                 nzval1
                                 )
    return f!, J!, guess_values, residuals, JJ, exogenous
end

@testset verbose = true "Scenario" begin
    @testset "shocks!" begin
        @test context.work.scenario[Undated(1)][Undated(1)][:e] == (0.01 => Symbol())
        @test context.work.scenario[Undated(1)][Undated(2)][:u] == (0.015 => Symbol())
        @test context.work.scenario[Undated(1)][Undated(1)][:y] == (1.0 => Symbol(:u))
        
        Dynare.scenario!(infoperiod=Undated(2), name=:y, value=0.2, period=2, context = context, exogenous = :u)
        @test length(context.work.scenario) == 3
        @test context.work.scenario[2][Undated(2)] == Dict(:y => (0.2 => :u))
        
        Dynare.scenario!(infoperiod=Undated(2), name=:y, value=0.3, period=3, context = context, exogenous = :u)
        @test length(context.work.scenario[2]) == 2
        @test context.work.scenario[2][3][:y] == (0.3 => :u)
    end

    @testset "set_future_information" begin
        y = zeros(60)
        x = zeros(20)
        Dynare.set_future_information!(y, x, context, 10, 1)
        @test x[1] ≈ 0.01
        @test x[4] ≈ 0.015
        @test y[1] ≈ 1.0
        @test y[13] ≈ 1.0
    end
    
    @testset "FlipInformation" begin
        FI = Dynare.FlipInformation(context, 5, 2)
        target = repeat([false], 5*6)
        target[[1, 7]] .= true
        @test FI.flipped_variables == target 
        @test FI.ix_period[[1, 7]] == [20, 20]
        @test FI.ix_stack[[1, 7]] == [2, 4]
    end

    @testset "flip!" begin
        FI = Dynare.FlipInformation(context, 5, 2)
        x = collect(1:10)
        x0 = copy(x)
        y = zeros(30)
        Dynare.flip!(y, x, FI.ix_stack)
        @test y[1] == x0[2]
        @test y[7] == x0[4]
    end
    
    @testset "residuals and Jacobian" begin
        periods = 4
        f!, J!, guess_values, residuals, JJ, exogenous = make_f_J(context, context.work.scenario, periods)

        f!(residuals, guess_values)
        
        @test guess_values[1] == context.results.model_results[1].trends.endogenous_steady_state[1]
        @test exogenous[2] == 0
        @test guess_values[13] == context.results.model_results[1].trends.endogenous_steady_state[1]
        @test exogenous[6] == 0
        @test guess_values[14] == context.results.model_results[1].trends.endogenous_steady_state[2]

        J!(JJ, guess_values)
        
        J1target = zeros(24)
        J1target[6] = -1
        @test JJ[:, 1] == J1target
        J13target = zeros(24)
        J13target[18] = -1
        @test JJ[:, 13] == J13target
    end

    @testset "Newton step" begin
        scenario = Dict(1 => Dict(1 => Dict(:y => (1.0 => :e))))
        context.work.scenario = scenario

        periods = 4
        FI = Dynare.FlipInformation(context, periods, 1)

        f!, J!, guess_values, residuals, JJ, exogenous = make_f_J(context, context.work.scenario, periods)

        guess_values[1] = 0
        exogenous[1] = 1.0

        f!(residuals, guess_values)
        
        y, c, k, a, h, b = context.results.model_results[1].trends.endogenous_steady_state
        beta, rho, alpha, delta, theta, psi, tau = context.work.params
        y_star = 1.0
        @test residuals[1] ≈ c*theta*h^(1+psi) - (1-alpha)*y_star;
        @test residuals[2] ≈ (k - beta*(((exp(b)*c)/(exp(b)*c))
                                        *(exp(b)*alpha*y + (1 - delta)*k)))
        @test residuals[3] ≈ y_star - exp(a)*(k^alpha)*(h^(1 - alpha))
        @test residuals[4] ≈ k - (exp(b)*(y_star - c) + (1 - delta)*k)
        @test all(abs.(residuals[5:24]) .< 2*eps())

        J!(JJ, guess_values)

        Dy =  -JJ\residuals
    end

    @testset "Newton solution" begin
        scenario = Dict(1 => Dict(1 => Dict(:y => (1.1 => :e))))
        context.work.scenario = scenario

        periods = 100
        FI = Dynare.FlipInformation(context, periods, 1)

        f!, J!, guess_values, residuals, JJ, exogenous = make_f_J(context, context.work.scenario, periods)

        Dynare.set_future_information!(guess_values, exogenous, context, periods, 1)
        Dynare.flip!(guess_values, exogenous, FI.ix_stack)
        @test guess_values[1] == 0
        @test exogenous[1] == 1.1

        f!(residuals, guess_values)
        
        y, c, k, a, h, b = context.results.model_results[1].trends.endogenous_steady_state
        beta, rho, alpha, delta, theta, psi, tau = context.work.params
        y_star = 1.1
        @test residuals[1] ≈ c*theta*h^(1+psi) - (1-alpha)*y_star;
        @test residuals[2] ≈ (k - beta*(((exp(b)*c)/(exp(b)*c))
                                        *(exp(b)*alpha*y + (1 - delta)*k)))
        @test residuals[3] ≈ y_star - exp(a)*(k^alpha)*(h^(1 - alpha))
        @test residuals[4] ≈ k - (exp(b)*(y_star - c) + (1 - delta)*k)
        @test all(abs.(residuals[5:600]) .< 2*eps())

        J!(JJ, guess_values)

        Dy =  -JJ\residuals

        y = copy(guess_values)

        for i=1:4
            y += Dy
            f!(residuals, y)
            J!(JJ, y)
            Dy = -JJ\residuals
        end
    end

    @testset "Complex scenario" begin
        
        periods = 100
        FI = Dynare.FlipInformation(context, periods, 1)
        
        f!, J!, guess_values, residuals, JJ, exogenous = make_f_J(context, context.work.scenario, periods)

        Dynare.set_future_information!(guess_values, exogenous, context, periods, 1)
        Dynare.flip!(guess_values, exogenous, FI.ix_stack)

        f!(residuals, guess_values)
        
        J!(JJ, guess_values)

        Dy =  -JJ\residuals

        y = copy(guess_values)
        for i=1:4
            y += Dy
            f!(residuals, y)
            J!(JJ, y)
            Dy = -JJ\residuals
        end
        @test norm(residuals) < 1e-12
        e = AxisArrayTables.data(simulation(:e))
        u = AxisArrayTables.data(simulation(:u))
        y = AxisArrayTables.data(simulation(:y))
        @test e[1] != 0
        @test u[1] != 0
        @test u[2] == 0.015
        @test u[3] != 0
        @test all(u[4:100] .== 0)
        @test all(e[2:100] .== 0)
        @test y[1] == 1
        @test y[3] == 1
    end

end

nothing







