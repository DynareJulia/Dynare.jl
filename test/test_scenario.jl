using Dynare
using ExtendedDates
using JLD2
using Test

context = load("models/example1pf/example1pf/output/example1pf.jld2", "context")

@testset verbose = true "Scenario" begin
    @testset "shocks!" begin
        Dynare.shocks!(name=:e, value=0.1, period=2, context = context)
        @test context.work.scenario == Dict(Undated(1) => Dict(:e => Dict(Undated(2) => (0.1 => Symbol()))))
        
        Dynare.shocks!(infoperiod=Undated(2), name=:y, value=0.2, period=2, context = context, exogenous = :u)
        @test length(context.work.scenario) == 2
        @test context.work.scenario[2][:y] == Dict(Undated(2) => (0.2 => :u))
        
        Dynare.shocks!(infoperiod=Undated(2), name=:y, value=0.3, period=3, context = context, exogenous = :u)
        @test length(context.work.scenario[2][:y]) == 2
        @test context.work.scenario[2][:y][3] == (0.3 => :u)
    end

    @testset "make_flips" begin
        @test Dynare.make_flips(context, Undated(1))[2] == [(1 => 2), (7 => 4)]
    end

    @testset "residuals" begin
        flips = Dynare.make_flips(context, Undated(1))

        datafile = ""
        periods = 4
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
        flips1 = flips[2]
        
        f! = Dynare.make_pf_residuals(
            initial_values,
            terminal_values,
            exogenous,
            dynamic_variables,
            steadystate,
            params,
            m,
            periods,
            temp_vec,
            perfect_foresight_ws.permutationsR,
            flips1
        )

        f!(residuals, guess_values)

        @test guess_values[1] == 0
        @test exogenous[2] == 0.2
        @test guess_values[7] == 0
        @test exogenous[4] == 0.3
        @test guess_values[13] == context.results.model_results[1].trends.endogenous_steady_state[1]

    end
end

nothing







