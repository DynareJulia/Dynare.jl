using Revise
using Dynare
using Plots
using KernelDensity
using Distributions

context = @dynare "test/models/ls2003/ls2003.mod"

function find(names::Vector, s)
    for (i, n) in enumerate(names)
        if n == s
            return i
        end
    end
end

has_posterior_mode(context) = length(context.work.estimated_parameters.posterior_mode)>0

function mcmc_diagnostics(chains, context, names)
    indices = [find(context.work.estimated_parameters.name, e) for e in names]
    iterations, column_nbr = size(chains.value.data)
    prior_pdfs = []
    n_plots = length(indices)
    posterior_pdfs = []
    posterior_x_axes = []
    prior_x_axes = []
    transformation = Dynare.DSGETransformation(context.work.estimated_parameters)
    transformed_data = Matrix{Float64}(undef, iterations, column_nbr - 1)
    
    for it in 1:iterations
        @views transformed_data[it, :] .= TransformVariables.transform(transformation, chains.value.data[it, 1:end-1])
    end
    
    for i in indices
        U = kde(transformed_data[:, i])
        prior = context.work.estimated_parameters.prior[i]
        m, v = mean(prior), var(prior)
        prior_x_axis = LinRange(m-15*v, m+15*v, length(U.x))
        prior_pdf = [pdf(prior, e) for e in prior_x_axis]
        push!(posterior_pdfs, U.density*(U.x[2] - U.x[1]))
        push!(posterior_x_axes, U.x)
        push!(prior_pdfs, prior_pdf/sum(prior_pdf))
        push!(prior_x_axes, prior_x_axis)
    end

    posterior_pdfs = hcat(posterior_pdfs...)
    posterior_x_axes = hcat(posterior_x_axes...)
    prior_pdfs = hcat(prior_pdfs...)
    prior_x_axes = hcat(prior_x_axes...)
    names = hcat(names...)

    f = plot(posterior_x_axes, posterior_pdfs, layout=n_plots, linecolor=:darkblue, legend=false, linewidth=3)
    plot!(f, prior_x_axes, prior_pdfs, layout=n_plots, title=names, linecolor=:darkgrey, legend=false, linewidth=3)
    if has_posterior_mode(context)
        posterior_modes = [context.work.estimated_parameters.posterior_mode[i] for i in indices]
        plot!(f, hcat(posterior_modes...), layout=n_plots, seriestype=:vline, line=(:dot,3), linecolor=:green, legend=false, linewidth=3)
    end
    f
end


function plot_priors(context, names, n_points = 100)
    indices = [find(context.work.estimated_parameters.name, e) for e in names]
    prior_pdfs = []
    n_plots = length(indices)
    prior_x_axes = []
    for i in indices
        prior = context.work.estimated_parameters.prior[i]
        m, v = mean(prior), var(prior)
        prior_x_axis = LinRange(m-15*v, m+15*v, n_points)
        prior_pdf = [pdf(prior, e) for e in prior_x_axis]
        push!(prior_pdfs, prior_pdf/sum(prior_pdf))
        push!(prior_x_axes, prior_x_axis)
    end
    prior_pdfs = hcat(prior_pdfs...)
    prior_x_axes = hcat(prior_x_axes...)
    names = hcat(names...)
    f = plot(prior_x_axes, prior_pdfs, layout=n_plots, title=names, linecolor=:darkgrey, legend=false, linewidth=3)
    f
end


f = mcmc_diagnostics(Dynare.my_chains, context, ["psi1", "psi2", "psi3", "rho_R"])
display(f)
#g = plot_priors(context, ["psi1", "psi2", "psi3", "rho_R"])
#display(g)
