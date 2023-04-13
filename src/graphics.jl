using KernelDensity
using Plots

function graph_display(g)
    if get(ENV, "TERM_PROGRAM", "") == "vscode"
        display(g)
    end
end

function plot(
    variables::Vector{Vector{Float64}};
    bar_variable::Vector{Float64} = Vector{Float64}([]),
    first_period::Int64 = 1,
    plot_title::String = "Smoothed value",
    plot_legend::Tuple = (),
    plot_filename::String = "",
    plot_legend_position::Symbol = :topright,
)
    local myplot
    # deal first with bar variable
    if length(bar_variable) > 0
        variables = pushfirst!(variables, bar_variable)
    end
    for (i, v) in enumerate(variables)
        if length(plot_legend) > 0
            thislabel = plot_legend[i]
        else
            thislabel = ""
            plot_legend_position = false
        end
        if i == 1
            if length(bar_variable) > 0
                len = length(bar_variable)
                xb = collect(range(first_period, length = len))
                myplot = Plots.bar(
                    xb,
                    bar_variable,
                    label = thislabel,
                    legend = plot_legend_position,
                    title = plot_title,
                )
                twinx()
            else
                x = collect(range(first_period, length = length(v)))
                myplot = Plots.plot(
                    x,
                    v,
                    label = thislabel,
                    legend = plot_legend_position,
                    title = plot_title,
                    linewidth = 3,
                )
            end
        else
            x = collect(range(first_period, length = length(v)))
            myplot = Plots.plot!(x, v, label = thislabel, linewidth = 3)
        end
    end
    graph_display(myplot)
    if length(plot_filename) > 0
        Plots.savefig(plot_filename)
    end
end

function plot(
    variable::Vector{Float64};
    bar_variable::Vector{Float64} = Vector{Float64}([]),
    first_period::Int64 = 1,
    plot_title::String = "Smoothed value",
    plot_legend::Tuple = (),
    plot_filename::String = "",
    plot_legend_position::Symbol = :topright,
)
    local myplot
    # deal first with bar variable
    if length(plot_legend) > 0
        thislabel = plot_legend
    else
        thislabel = ""
        plot_legend_position = false
    end
    if length(bar_variable) > 0
        len = length(bar_variable)
        xb = collect(range(first_period, length = len))
        myplot = Plots.bar(
            xb,
            bar_variable,
            label = thislabel[1],
            legend = plot_legend_position,
            title = plot_title,
        )
        x = collect(range(first_period, length = length(variable)))
        lims = Plots.ignorenan_extrema(myplot[1].attr[:yaxis])
        m, M = extrema(variable)
        tb = (lims[2] - lims[1]) / (M - m)
        variable = tb * (variable .- m) .+ lims[1]
        Plots.plot!(x, variable, label = thislabel[2], title = plot_title, linewidth = 3)
        myplot = twinx()
        myplot = Plots.plot!(myplot, ylims = (m, M))
    else
        x = collect(range(first_period, length = length(variable)))
        myplot = Plots.plot(
            x,
            variable,
            label = thislabel[1],
            legend = plot_legend_position,
            title = plot_title,
            linewidth = 3,
        )
    end
    graph_display(myplot)
    if length(plot_filename) > 0
        Plots.savefig(plot_filename)
    end
end

function plot_irfs(irfs, model, symboltable, filepath)
    x = axes(first(irfs)[2])[1]
    endogenous_names = get_endogenous_longname(symboltable)
    exogenous_names = get_exogenous_longname(symboltable)
    for i = 1:model.exogenous_nbr
        exogenous_name = exogenous_names[i]
        (nbplt, nr, nc, lr, lc, nstar) = pltorg(model.original_endogenous_nbr)
        ivars = 1:nr*nc
        for p = 1:nbplt-1
            filename = "$(filepath)_$(exogenous_name)_$(p).png"
            plot_panel(
                x,
                Matrix(irfs[Symbol(exogenous_name)][:,ivars]),
                "Orthogonal shock to $(exogenous_name)",
                endogenous_names[ivars],
                nr,
                nc,
                nr * nc,
                filename,
            )
            ivars += nr * nc
        end
        ivars = ivars[1:nstar]
        filename = "$(filepath)_$(exogenous_name)_$(nbplt).png"
        plot_panel(
            x,
            Matrix(irfs[Symbol(exogenous_name)][:, ivars]),
            "Orthogonal shock to $(exogenous_name)",
            endogenous_names[ivars],
            lr,
            lc,
            nstar,
            filename,
        )
    end
end

function plot_priors(;context::Context=context)
    ep = context.work.estimated_parameters
    @assert length(ep.prior) > 0 "There is no defined priors"
    path = "$(context.modfileinfo.modfilepath)/graphs/"
    mkpath(path)
    filepath = "$(path)/Priors"
    nprior = length(ep.prior)
    X = zeros(100, nprior)
    Y = zeros(100, nprior)
    for (i, p) in enumerate(ep.prior)
        X[:, i] = range(StatsPlots.yz_args(p)..., 100)
        for j in axes(X, 1)
            Y[j, i] = pdf(p, X[j, i])
        end 
    end     
    (nbplt, nr, nc, lr, lc, nstar) = pltorg(nprior)
    ivars = collect(1:nr*nc)
    for p = 1:nbplt-1
        filename = "$(filepath)_$(p).png"
        plot_panel(
            X[:, ivars],
            Y[:, ivars],
            "Priors ($p)",
            ep.name[ivars],
            nr,
            nc,
            nr * nc,
            filename,
        )
        ivars .+= nr * nc
    end
    ivars = ivars[1:nstar]
    filename = "$(filepath)_$(nbplt).png"
    plot_panel(
        X[:, ivars],
        Y[:, ivars],
        "Priors ($nbplt)",
        ep.name[ivars],
        lr,
        lc,
        nstar,
        filename,
    )
end

function plot_panel(
    x,
    y,
    title,
    vnames,
    nr,
    nc,
    nstar,
    filename,
    kwargs...
)
    sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nr*nc]
    for i = 1:nstar
        if ndims(x) == 1
            xx = x
        else
            xx = x[:, i]
        end 
        yy = y[:, i]
        if all(yy .> 0)
            lims = (0, Inf)
        elseif all(yy .< 0)
            lims = (-Inf, 0)
        else
            lims = (-Inf, Inf)
        end
        sp[i] =
            Plots.plot(xx, yy, title = vnames[i], ylims = lims, kwargs...)
    end

    pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title= title)
    graph_display(pl)
    savefig(filename)
end

function pltorg(number)
    nrstar = 3
    ncstar = 3
    nstar = nrstar * ncstar
    nbplt = 0
    nr = 0
    nc = 0
    lr = 0
    lc = 0
    nbplt = ceil(number / nstar)
    (d, r) = divrem(number, nstar)
    if r == 0
        return (d, nrstar, ncstar, nrstar, ncstar, nstar)
    else
        (lr, lc) = pltorg_0(r)
        return (d + 1, nrstar, ncstar, lr, lc, r)
    end
end

function pltorg_0(number)
    @assert number < 10
    if number == 1
        lr = 1
        lc = 1
    elseif number == 2
        lr = 2
        lc = 1
    elseif number == 3
        lr = 3
        lc = 1
    elseif number == 4
        lr = 2
        lc = 2
    elseif number < 7
        lr = 3
        lc = 2
    elseif number < 10
        lr = 3
        lc = 3
    end
    return (lr, lc)
end

function plot_priorprediction_irfs(irfs, model, symboltable, filepath)
    x = axes(first(first(irfs))[2])[1]
    endogenous_names = get_endogenous_longname(symboltable)
    exogenous_names = get_exogenous_longname(symboltable)
    endogenous_nbr = model.original_endogenous_nbr
    for i = 1:model.exogenous_nbr
        exogenous_name = exogenous_names[i]
        (nbplt, nr, nc, lr, lc, nstar) = pltorg(model.original_endogenous_nbr)
        k = 1
        p = 1
        filename = "$(filepath)_$(exogenous_name)_$(p).png"
        sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nstar]
        while p <= endogenous_nbr
            pl = plot_panel_priorprediction_irfs(
                x,
                irfs,
                exogenous_name,
                "Orthogonal shock to $(exogenous_name)",
                endogenous_names[p]
            )
            if pl.n > 0 
                sp[k] = pl
            else
                k -= 1
            end 
            if k == nr*nc || p == endogenous_nbr
                pl = Plots.plot(sp..., layout = (nr, nc), size = (900, 900), plot_title = "Priorprediction: Orthogonal shock to $(exogenous_name)")
                graph_display(pl)
                savefig(filename)
                sp = [Plots.plot(showaxis = false, ticks = false, grid = false) for i = 1:nstar]
                k = 1
            end
            k += 1
            p += 1
        end
    end
end

function plot_panel_priorprediction_irfs(
    x,
    irfs,
    exogenous_name,
    title,
    endogenous_name,
)
    M = zeros(40, length(irfs))
    for (k, ir) in enumerate(irfs)
        m = Matrix(ir[Symbol(exogenous_name)][:, Symbol(endogenous_name)])
        @views M[:, k] .= m
    end
    Q = zeros(40, 5)
    for k in 1:40
        @views Q[k, :] .= quantile(M[k,:], [0.1, 0.3, 0.5, 0.7, 0.9])
    end
    colors = palette(:Blues_3)
    Q1 = view(Q,:, 1)
    Q2 = view(Q,:, 2)
    Q3 = view(Q,:, 3)
    Q4 = view(Q,:, 4)
    Q5 = view(Q,:, 5)
    ex = extrema(Q5)
    if abs(ex[1]-ex[2]) < 1e-8
        return Plots.plot()
    end
    p = Plots.plot(x, Q3, color = :Blue, labels="Median")
    plot!(p, x, Q1, fillrange=Q5, fillalpha = 0.3, color=:Blue,labels="80%", linealpha=0)
    plot!(p, x, Q2, fillrange=Q3, fillalpha = 0.5, color=:Blue,labels="40%", linealpha=0)
    plot!(p, title=endogenous_name)
    return p
end

function plot(aat::AxisArrayTable; label = (), title = "", filename = "")
    pl = Plots.plot(aat, label = label, title = title)
    graph_display(pl)
    savefig(filename)
end
