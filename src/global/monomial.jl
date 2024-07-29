"""
    monomial nodes and weights
"""
struct MonomialPowerIntegration
    nodes::Matrix{Float64}
    weights::Vector{Float64}
end

function MonomialPowerIntegration(nShocks::Int)
    # Number of integration nodes
    numNodes = 2*nShocks^2 + 1

    z0 = zeros(numNodes, nShocks)

    # Deviations in one dimension (note that the origin is row one)
    for i1 in 1:nShocks
        z0[(i1 - 1)*2 + 2,i1] = 1.0
        z0[(i1 - 1)*2 + 3,i1] = -1.0
    end
    
    i0 = 0
    # Deviations in two dimensions
    for i1 in 1:nShocks
        for i2 in i1+1:nShocks
            z0[2*nShocks+2+i0*4, i1] = 1.0
            z0[2*nShocks+3+i0*4, i1] = 1.0
            z0[2*nShocks+4+i0*4, i1] = -1.0
            z0[2*nShocks+5+i0*4, i1] = -1.0
            z0[2*nShocks+2+i0*4, i2] = 1.0
            z0[2*nShocks+3+i0*4, i2] = -1.0
            z0[2*nShocks+4+i0*4, i2] = 1.0
            z0[2*nShocks+5+i0*4, i2] = -1.0
            i0 += 1
        end
    end
    # Nodes
    IntNodes = Matrix{Float64}(undef, numNodes, nShocks)
    @views begin
        IntNodes[1, :] .= zeros(nShocks)
        IntNodes[2:nShocks*2+1, :] .= z0[2:nShocks*2+1, :]*sqrt(2.0 + nShocks)
        IntNodes[nShocks*2+2:end, :] .= z0[nShocks*2+2:end, :]*sqrt((2.0+nShocks)/2.0)
    end 

    # Weights
    IntWeights = Vector{Float64}(undef, numNodes)

    IntWeights[1] = 2.0/(2.0+nShocks)
    IntWeights[2:nShocks*2+1] .= (4-nShocks)/(2*(2+nShocks)^2)
    IntWeights[nShocks*2+2:end] .= 1.0/(nShocks+2)^2

    return MonomialPowerIntegration(IntNodes, IntWeights)
end

function (monomial::MonomialPowerIntegration)(points::AbstractVector{Float64})
    return dot(points, monomial.weights)
end 
