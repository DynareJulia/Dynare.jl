module AdvancedMHForwardDiffExt

if isdefined(Base, :get_extension)
    import AdvancedMH
    import DiffResults
    import ForwardDiff
else
    import ..AdvancedMH
    import ..DiffResults
    import ..ForwardDiff
end

function AdvancedMH.logdensity_and_gradient(model::AdvancedMH.DensityModel, params)
    res = DiffResults.GradientResult(params)
    ForwardDiff.gradient!(res, model.logdensity, params)
    return DiffResults.value(res), DiffResults.gradient(res)
end


end # module
