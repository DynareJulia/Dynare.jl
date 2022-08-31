using JLD
using PrettyTables
using Dynare

function f(context)
    x = context.results.model_results[1].linearrationalexpectations.g1
    vx = view(x, :, 1:size(x, 2)-1)
    Dynare.display_stoch_simul(
        vx',
        "Coefficients of approximate solution function",
        context,
    )
end
