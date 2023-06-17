mutable struct PriorPredictionResults
    domain::Int64
    other::Int64
    undetermined::Int64
    unstable::Int64
    steady_state::Array{Union{Missing,Float64}, 2}
    stdev::Array{Union{Missing, Float64}, 2}
    irfs::Vector{Dict{Symbol, AxisArrayTable}}
    function PriorPredictionResults(n::Int64, m::Model)
        domain = 0
        other = 0
        undetermined = 0
        unstable = 0
        steady_state = Array{Union{Missing, Float64}, 2}(missing, m.endogenous_nbr, n)
        stdev = Array{Union{Missing,Float64}, 2}(missing, m.endogenous_nbr, n)
        irfs = Vector{Dict{Symbol, AxisArrayTable}}(undef, 0)
        new(domain, other, undetermined, unstable, steady_state, stdev, irfs)
    end
end

function priorprediction(;context::Context=context, iterations::Int64=1000)
    estimated_parameters = context.work.estimated_parameters
    parameter_nbr = length(estimated_parameters)
    model = context.models[1]
    modfilepath = context.modfileinfo.modfilepath
    
    draws = Matrix{Float64}(undef, iterations, parameter_nbr)
    for (i, p) in enumerate(estimated_parameters.prior)
        @views draws[:, i] = rand(p, iterations)
    end
    (failure, results) = run_simulations(context, draws)

    display_priorprediction_results(draws, results, failure, iterations, parameter_nbr, estimated_parameters.name, get_endogenous(context.symboltable), modfilepath)
    plot_priorprediction_irfs(results.irfs, model, context.symboltable, "$(modfilepath)/graphs/priorprediction_irfs")

    draws = Matrix{Float64}(undef, iterations, parameter_nbr)
    for (i, p) in enumerate(estimated_parameters.prior)
        lb, ub = quantile(p, [0.001, 0.999])
        @views draws[:, i] = rand(Uniform(lb, ub), iterations)
    end
    (failure, results) = run_simulations(context, draws)

    display_priorprediction_checks(draws, results, failure, iterations, parameter_nbr, estimated_parameters.name, get_endogenous(context.symboltable), modfilepath)
end


function run_simulations(context, draws)
    iterations = size(draws, 1)    
    modfilepath = context.modfileinfo.modfilepath
    mkpath("$(modfilepath)/graphs")
    model = context.models[1]
    model_results = context.results.model_results[1]
    lre_results = model_results.linearrationalexpectations
    work = context.work
    estimated_parameters = work.estimated_parameters
    model_parameters = work.params
    D = eltype(context.work.params)
    dynamicws = Dynare.DynamicWs(context)

    parameter_nbr = length(estimated_parameters.prior)
    results = PriorPredictionResults(iterations, model)
    failure = zeros(iterations)
    
    for i = 1:iterations
        @views set_estimated_parameters!(context, draws[i, :])
        fill!(model_results.exogenous_steady_state, 0.0)
        #compute steady state and first order solution
        try
            compute_stoch_simul!(
                context,
                dynamicws,
                model_parameters,
                StochSimulOptions(Dict{String, Any}());
                variance_decomposition = true,
            )
            irfs!(context, 40)
            @views begin
                results.steady_state[:, i] .= model_results.trends.endogenous_steady_state
                results.stdev[:, i] .= robustsqrt.(diag(lre_results.endogenous_variance))
                push!(results.irfs, copy(model_results.irfs))
            end 
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
    return (failure, results)
end

function display_priorprediction_results(draws, results, failure, iterations, parameter_nbr, parameter_names, endogenous_names, modfilepath)
    failure_nbr =sum(failure)

    println("\nDISTRIBUTION OF MOMENTS\n")
    println("Share of DomainError: $(results.domain/iterations)")
    println("Share of UndeterminedError: $(results.undetermined/iterations)")
    println("Share of UnstableError: $(results.unstable/iterations)")
    println("Share of other error: $(results.other/iterations)")

    display_priorprediction_moments("STEADY STATE", results.steady_state, endogenous_names)
    display_priorprediction_moments("STANDARD DEVIATION", results.stdev, endogenous_names)
end

function display_priorprediction_checks(draws, results, failure, iterations, parameter_nbr, parameter_names, endogenous_names, modfilepath)
    failure_nbr =sum(failure)

    println("\nDISTRIBUTION OF COMPUTATION FAILURES\n")
    println("Share of DomainError: $(results.domain/iterations)")
    println("Share of UndeterminedError: $(results.undetermined/iterations)")
    println("Share of UnstableError: $(results.unstable/iterations)")
    println("Share of other error: $(results.other/iterations)")

    (nbplt, nr, nc, lr, lc, nstar) = pltorg(parameter_nbr)

    sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]

    k = 1
    nfig = 1
    for (i, p) in enumerate(parameter_names)
        @views z = hcat(draws[:, i], failure)
        @views sz = z[sortperm(z[:, 1]), :]
        @views y = cumsum(sz[:, 2])
        
        @views sp[k] = Plots.plot(sz[:, 1], y/results.undetermined, title = p, labels = false)
        k += 1
        if k > nr*nc || i == length(parameter_names)
            pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title = "Parameter and computation failure ($nfig)")
            graph_display(pl)
            savefig("$(modfilepath)/graphs/PriorChecks$(nfig).png")
            k = 1
            nfig += 1
            sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]
        end
    end
end

function display_priorprediction_moments(title, x, names)
    println("")
    rss = x
    data = Matrix{Any}(undef, length(names) + 1, 3)
    data[1,1] = "Variable"
    data[1,2] = "Mean"
    data[1,3] = "80% interval"
    for i = axes(rss, 1)
        srss = skipmissing(rss[i,:])
        data[i+1, 1] = names[i]
        data[i+1, 2] = mean(srss)
        data[i+1, 3] = ClosedInterval(quantile(srss, [0.1, 0.9])...)
    end
    dynare_table(data, title, columnheader = true) 
end