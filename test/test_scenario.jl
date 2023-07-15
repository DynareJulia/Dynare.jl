using Dynare
using ExtendedDates
using JLD2
using Test

@dynare "models/example1pf/example1pf_conditional"
context = load("models/example1pf/example1pf_conditional/output/example1pf_conditional.jld2", "context")

@testset verbose = true "Scenario" begin
    @testset "shocks!" begin
        @test context.work.scenario[Undated(1)][Undated(1)][:e] == (0.01 => Symbol())
        @test context.work.scenario[Undated(1)][Undated(2)][:u] == (0.015 => Symbol())
        @test context.work.scenario[Undated(1)][Undated(1)][:y] == (1.0 => Symbol(:u))
        @test context.work.scenario[Undated(1)][Undated(2)][:h] == (0.35 => Symbol(:e))
        
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
        Dynare.set_future_information(y, x, context, 10, 1)
        @test x[1] ≈ 0.01
        @test x[4] ≈ 0.015
        @test y[1] ≈ 1.0
        @test y[11] ≈ 0.35
        @test y[13] ≈ 1.0
    end
    
    @testset "FlipInformation" begin
        FI = Dynare.FlipInformation(context, 5, 2)
        target = repeat([false], 5*6)
        target[[1, 7]] .= true
        @test FI.flipped_variables == target 
        @test FI.flips_period[[1, 7]] == [20, 20]
        @test FI.flips_stack[[1, 7]] == [2, 4]
    end

    @testset "flip!" begin
        FI = Dynare.FlipInformation(context, 5, 2)
        x = collect(1:10)
        x0 = copy(x)
        y = zeros(30)
        Dynare.flip!(y, x, FI.flips_stack)
        @test y[1] == x0[2]
        @test y[7] == x0[4]
    end
    
    @testset "residuals" begin
        datafile = ""
        periods = 4
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
        residuals = zeros(periods * m.endogenous_nbr)
        dynamic_variables = dynamic_ws.dynamic_variables
        temp_vec = dynamic_ws.temporary_values
        steadystate = results.trends.endogenous_steady_state
        params = work.params
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

        # infoperiod == 2
        
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
                                      FI.flips_stack
                                      )

        f!(residuals, guess_values)
        
        @test guess_values[1] == 0
        @test exogenous[2] == context.results.model_results[1].trends.endogenous_steady_state[1]
        @test guess_values[11] == 0
        @test exogenous[3] == context.results.model_results[1].trends.endogenous_steady_state[5]
        @test guess_values[14] == context.results.model_results[1].trends.endogenous_steady_state[2]
        
        permutations = []
        J, permutations = Dynare.makeJacobian(m.dynamic_g1_sparse_colptr, m.dynamic_g1_sparse_rowval,
                                m.endogenous_nbr, m.exogenous_nbr, periods, permutations, FI)

        Dynare.updateJacobian!(J,
                               Dynare.DFunctions.dynamic_derivatives!,
                               guess_values,
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
                               nzval1)
        
        J1target = zeros(24)
        J1target[6] = -1
        @test J[:, 1] == J1target
        J11target = zeros(24)
        J11target[11] = -1
        @test J[:, 11] == J11target
        J13target = zeros(24)
        J13target[18] = -1
        Jtarget, permutations = Dynare.makeJacobian(m.dynamic_g1_sparse_colptr, m.dynamic_g1_sparse_rowval,
                                m.endogenous_nbr, periods, permutations)

        Dynare.flip!(guess_values, exogenous, FI.flips_stack)
        Dynare.updateJacobian!(Jtarget,
                               Dynare.DFunctions.dynamic_derivatives!,
                               guess_values,
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
                               nzval1)
    end

end

nothing







