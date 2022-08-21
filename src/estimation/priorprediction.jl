mutable struct PriorPredictiveCheckResults
    domain::Int64
    other::Int64
    undetermined::Int64
    unstable::Int64
    function PriorPredictiveCheckResults(n, m)
        domain = 0
        other = 0
        undetermined = 0
        unstable = 0
        new(domain, other, undetermined, unstable)
    end
end

function priorpredictivecheck(context, iterations)
    model = context.models[1]
    results = context.results.model_results[1]
    work = context.work
    estimated_parameters = work.estimated_parameters
    model_parameters = work.params
    D = eltype(context.work.params)
    tmp_nbr = context.dynarefunctions.dynamic!.tmp_nbr
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both
    dynamicws = Dynare.DynamicWs(
            model.endogenous_nbr,
            model.exogenous_nbr,
            ncol,
            sum(tmp_nbr[1:2]),
        )

    parameter_nbr = length(estimated_parameters.prior)
    draws = Matrix{Float64}(undef, iterations, parameter_nbr)
    for (i, p) in enumerate(estimated_parameters.prior)
        @views draws[:, i] = rand(p, iterations)
    end

    @show maximum(draws, dims=1)
    @show minimum(draws, dims=1)
    results = PriorPredictiveCheckResults(iterations, length)
    failure = zeros(iterations)
    
    for i = 1:iterations
        @views set_estimated_parameters!(context, draws[i, :])
        #compute steady state and first order solution
        try
            Dynare.compute_stoch_simul!(
                context,
                dynamicws,
                model_parameters,
                StochSimulOptions(Dict{String, Any}());
                variance_decomposition = false,
            )
        catch e
            if isa(e, DomainError)
                failure[i] = 1
                results.domain += 1
            elseif isa(e, LinearRationalExpectations.UndeterminateSystemException)
                failure[i] = 1
                results.undetermined += 1
            elseif isa(e, LinearRationalExpectations.UnstableSystemException)
                failure[i] = 1
                results.unstable += 1
            else
                failure[i] = 1
                results.other += 1
            end
        end
    end

    failure_nbr =sum(failure)

    println("Share of DomainError: $(results.domain/iterations)")
    println("Share of UndeterminedError: $(results.undetermined/iterations)")
    println("Share of UnstableError: $(results.unstable/iterations)")
    println("Share of other error: $(results.other/iterations)")

    (nbplt, nr, nc, lr, lc, nstar) = pltorg(parameter_nbr)

    sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]

    k = 1
    nfig = 1
    for (i, p) in enumerate(estimated_parameters.name)
        @views z = hcat(draws[:, i], failure)
        @views sz = z[sortperm(z[:, 1]), :]
        @views y = cumsum(sz[:, 2])
        
        @views sp[k] = Plots.plot(sz[:, 1], y/results.undetermined, title = p)
        if k == nr*nc
            pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900))
            display(pl)
            graph_display(pl)
            savefig("PriorChecks$(nfig).png")
            k = 1
            nfig += 1
            sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]        end
        k += 1
    end
    if k < nr*nc
        pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900))
        display(pl)
        graph_display(pl)
        savefig("PriorChecks$(nfig).png")
    end
end
