function display_data(;context=context,
                      datafile = "",  
                      data = AxisArrayTable(AxisArrayTables.AxisArray(Matrix(undef, 0, 0))),
                      detrend = false,
                      first_obs = 1,
                      last_obs =0,
                      moments = true,
                      plot = true,
                      table = false)
    varobs = context.work.observed_variables
    Yorig = get_data(datafile, varobs, start = first_obs, last = last_obs)
    title = "Moments of observed variables"
    data = vcat("Name", varobs)
    data = hcat(data, vcat("Mean", mean(Yorig, dims=2)))
    data = hcat(data, vcat("Stdev", std(Yorig, dims=2)))
    data = hcat(data, vcat("AC", [autocor(Vector{Float64}(Yorig[i, :]), [1]) for i in axes(Yorig, 1)]...))
    dynare_table(data, title, columnheader = true) 
    title = "Observed variables correlation"
    data = reshape(vcat("", varobs), 1, length(varobs) + 1)
    data = vcat(data, hcat(varobs, cor(Yorig', Yorig')))
    dynare_table(data, title, columnheader = true)
    plot_data(Yorig, varobs, context.modfileinfo.modfilepath)
end

function plot_data(data, varobs, modfilepath)
    mkpath("$(modfilepath)/graphs")
    (nbplt, nr, nc, lr, lc, nstar) = pltorg(length(varobs))

    sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]

    k = 1
    nfig = 1
    for (i, p) in enumerate(varobs)
        @views sp[k] = Plots.plot(data[i, :], title = p, labels = false)
        k += 1
        if k > nr*nc || i == length(varobs)
            plottitle = "Observed variables"
            nbplt > 1 && (plotttitle = "$plottitle ($nfig)")
            pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title = plottitle)
            graph_display(pl)
            savefig("$(modfilepath)/graphs/ObservedVariables$(nfig).png")
            k = 1
            nfig += 1
            sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]
        end
    end
end
