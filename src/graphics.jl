using Plots
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
    Plots.display(myplot)
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
        @show (m, M)
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
    Plots.display(myplot)
    if length(plot_filename) > 0
        Plots.savefig(plot_filename)
    end
end

function plot_irfs(y, model, symboltable, filename)
    x = 1:size(y[1], 2)
    endogenous_names = get_endogenous_longname(symboltable)
    exogenous_names = get_exogenous_longname(symboltable)
    for i = 1:model.exogenous_nbr
        (nbplt, nr, nc, lr, lc, nstar) = pltorg(length(endogenous_names))
        firstvar = 1
        for p in 1:nbplt - 1
            plot_irf_panel(x, y[i], endogenous_names, exogenous_names[i],
                           firstvar, nr, nc, "$(filename)_$(exogenous_names[i])", i)
            firstvar += nr*nc
        end
        @show nr nc lr lc
        plot_irf_panel(x, y[i], endogenous_names, exogenous_names[i], firstvar, lr, lc,
                       "$(filename)_$(exogenous_names[i])", nbplt)
    end
end

function plot_irf_panel(x, y, endogenous_names, exogenous_name, firstvar, nr, nc,
                        filename, panel_nbr)
    sp = []
    for i = 1:nr*nc
        ivar = firstvar + i - 1
        title = (i == 1) ? "Orthogonal shock to $(exogenous_name)" : ""
        push!(sp, Plots.plot(x, view(y, ivar, :), 
                             title = title,
                             label = endogenous_names[ivar]))
    end
    Plots.plot(sp..., layout = (nr, nc))
    savefig("$(filename)_$(panel_nbr).png")
end

function pltorg(number)
    nrstar = 3
    ncstar = 3
    nstar  = nrstar*ncstar
    nbplt  = 0
    nr     = 0
    nc     = 0
    lr     = 0
    lc     = 0
    if number == 1
        nbplt = 1
        lr    = 1
        lc    = 1
    elseif number == 2
        nbplt = 1
        lr    = 2
        lc    = 1
    elseif number == 3
        nbplt = 1
        lr    = 3
        lc    = 1
    elseif number == 4
        nbplt = 1
        lr    = 2
        lc    = 2
    elseif number == 5
        nbplt = 1
        lr    = 3
        lc    = 2
    elseif number == 6
        nbplt = 1
        lr    = 3
        lc    = 2
    elseif number == 7
        nbplt = 1
        lr    = 3
        lc    = 3
    elseif number == 8
        nbplt = 1
        lr    = 3
        lc    = 3
    elseif number == 9
        nbplt = 1
        lr    = 3
        lc    = 3
    else
        if number/nstar == round(number/nstar)
            nbplt = number/nstar
            nr    = nrstar
            nc    = ncstar
            lr    = nr
            lc    = nc
        else
            nbplt = ceil(number/nstar)
            nr    = nrstar
            nc    = ncstar
            remaining = number-(nbplt-1)*nstar
            if remaining == 1
                lr    = 1
                lc    = 1
            elseif remaining == 2
                lr    = 2
                lc    = 1
            elseif remaining == 3
                lr    = 3
                lc    = 1
            elseif remaining == 4
                lr    = 2
                lc    = 2
            elseif remaining == 5
                lr    = 3
                lc    = 2
            elseif remaining == 6
                lr    = 3
                lc    = 2
            elseif remaining == 7
                lr    = 3
                lc    = 3
            elseif remaining == 8
                lr    = 3
                lc    = 3
            end
        end
    end
    return (nbplt, nr, nc, lr, lc, nstar) 
end
