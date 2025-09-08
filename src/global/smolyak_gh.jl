using SparseGrids

"""
    SmolyakGHIntegration(D::Int, O::Int)

Smolyak(Gauss–Hermite) grid via SparseGrids.jl, mapped to the 𝓝(0,1)^d measure.
Returns nodes as Vector{Vector{Float64}} and weights as ProbabilityWeights.
"""
function SmolyakGHIntegration(D::Int, O::Int)
    x, w = SparseGrids.sparsegrid(D, O)  # nodes for ∫ g(x) e^{-‖x‖²} dx
    nodes_std   = [ (√2) .* xi for xi in x ]        # map to 𝓝(0,1): z = √2 x
    weights_std = (1 / (π^(D/2))) .* w              # and ω = w / π^(D/2)
    s = sum(weights_std); s ≠ 0 && (weights_std ./= s)
    return SmolyakGHIntegration(nodes_std, collect(weights_std))
end