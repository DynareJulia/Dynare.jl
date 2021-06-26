using Plots
using GR
import Plots: twinx, Subplot, Plot, current, px, link_axes!


function twinx(sp::Subplot)
    sp[:right_margin] = max(sp[:right_margin], 30px)
    plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:yaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
    link_axes!(sp[:xaxis], twinsp[:xaxis])
    twinsp
end

twinx(plt::Plot = current()) = twinx(plt[1])

a = randn(20)
b = 3*randn(20) .+ 10

plt1 = Plots.bar(a, label="First label")
lims = Plots.ignorenan_extrema(plt1[1].attr[:yaxis])
m, M = extrema(b)
tb = (lims[2] - lims[1])/(M-m)
bb = tb*(b .- m) .+ lims[1]
plt2 = Plots.plot!(bb, label="Second label")
plt3 = twinx()
Plots.plot!(plt3, ylims=(m, M))
