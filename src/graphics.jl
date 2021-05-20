using Plots
function plot(variables::Vector{Vector{Float64}};
              bar_variable::Vector{Float64} = [],
              first_period::Int64 = 1,
              plot_title::String = "Smoothed value",
              plot_legend::Tuple = (),
              plot_filename::String = "",
              plot_legend_position::Symbol = :topright)
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
                xb = collect(range(first_period, length=len)) 
                myplot = Plots.bar(xb, bar_variable, label = thislabel,
                                   legend=plot_legend_position,
                                   title = plot_title)
                twinx()
            else
                x = collect(range(first_period, length=length(v)))
                myplot = Plots.plot(x, v, label = thislabel,
                                    legend=plot_legend_position,
                                    title = plot_title,
                                    linewidth = 3)
            end
        else
            x = collect(range(first_period, length=length(v)))
            myplot = Plots.plot!(x, v,
                                 label = thislabel,
                                 linewidth = 3)
        end
    end
    Plots.display(myplot)
    if length(plot_filename) > 0
        Plots.savefig(plot_filename)
    end
end

function plot(variable::Vector{Float64};
              bar_variable::Vector{Float64} = [],
              first_period::Int64 = 1,
              plot_title::String = "Smoothed value",
              plot_legend::Tuple = (),
              plot_filename::String = "",
              plot_legend_position::Symbol = :topright)
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
        xb = collect(range(first_period, length=len))
        myplot = Plots.bar(xb, bar_variable, label = thislabel[1],
                           legend=plot_legend_position,
                           title = plot_title)
        x = collect(range(first_period, length=length(variable)))
        lims = Plots.ignorenan_extrema(myplot[1].attr[:yaxis])
        m, M = extrema(variable)
        tb = (lims[2] - lims[1])/(M-m)
        variable = tb*(variable .- m) .+ lims[1]
        Plots.plot!(x, variable, label = thislabel[2],
                             title = plot_title,
                    linewidth = 3)
        myplot = twinx()
        @show (m, M)
        myplot = Plots.plot!(myplot, ylims = (m, M))
    else
        x = collect(range(first_period, length=length(variable)))
        Plots.plot(x, variable, label = thislabel[1],
                   legend=plot_legend_position,
                   title = plot_title,
                   linewidth = 3)
    end
    Plots.display(myplot)
    if length(plot_filename) > 0
        Plots.savefig(plot_filename)
    end
end
