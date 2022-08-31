using Revise
using Dynare
using Plots
using KernelDensity
using Distributions

context = @dynare "test/models/ls2003/ls2003.mod"
chains = mh_estimation(context, datafile = "test/models/ls2003/data_ca1.csv", first_obs=8, last_obs=86, iterations = 10000)

function find(names::Vector, s)
    for (i, n) in enumerate(names)
        if n == s
            return i
        end
    end
end

function mcmc_diagnostics(chains, context, names)

    indices = [find(context.work.estimated_parameters.name, e) for e in names]
    posterior_pdfs = []
    prior_pdfs = []
    n_plots = length(indices)
    x_axis = []

    for i in indices
        U = kde(chains.value.data[:, i])
        push!(posterior_pdfs, U.density*(U.x[2] - U.x[1]))
        prior_pdf = [pdf(context.work.estimated_parameters.prior[i], e) for e in U.x]
        push!(x_axis, U.x)
        push!(prior_pdfs, prior_pdf/sum(prior_pdf))
    end

    posterior_pdfs = hcat(posterior_pdfs...)
    prior_pdfs = hcat(prior_pdfs...)
    x_axis = hcat(x_axis...)

    names = hcat([e for e in context.work.estimated_parameters.name[1:n_plots]]...)

    f = plot(x_axis, prior_pdfs, layout=n_plots, title=names, linecolor=:darkgrey, legend=false, linewidth=3)
    plot!(f, x_axis, posterior_pdfs, layout=n_plots, linecolor=:darkblue, legend=false, linewidth=3)
end

mcmc_diagnostics(chains, context, ["psi1", "psi2", "rr", "tau"])