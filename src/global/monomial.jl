"""
    monomial nodes and weights
"""
function MonomialPowerIntegration(nShocks::Int)
    # Number of integration nodes
    numNodes = 2*nShocks^2 + 1

    z0 = zeros(nShocks, numNodes)

    # Deviations in one dimension (note that the origin is row one)
    for i1 in 1:nShocks
        z0[i1, (i1 - 1)*2 + 2] = 1.0
        z0[i1, (i1 - 1)*2 + 3] = -1.0
    end
    
    i0 = 0
    # Deviations in two dimensions
    for i1 in 1:nShocks
        for i2 in (i1 + 1):nShocks
            z0[i1, 2*nShocks + 2 + i0*4, ] = 1.0
            z0[i1, 2*nShocks + 3 + i0*4, ] = 1.0
            z0[i1, 2*nShocks + 4 + i0*4, ] = -1.0
            z0[i1, 2*nShocks + 5 + i0*4, ] = -1.0
            z0[i2, 2*nShocks + 2 + i0*4, ] = 1.0
            z0[i2, 2*nShocks + 3 + i0*4, ] = -1.0
            z0[i2, 2*nShocks + 4 + i0*4, ] = 1.0
            z0[i2, 2*nShocks + 5 + i0*4, ] = -1.0
            i0 += 1
        end
    end
    # Nodes
    IntNodes = [Vector{Float64}(undef, nShocks) for i = 1:numNodes]
    fill!(IntNodes[1], 0.0)
    for i = 2:nShocks*2+1
        IntNodes[i] .= z0[:, i]*sqrt(2.0 + nShocks)
    end
    for i = nShocks*2 + 2:numNodes
        IntNodes[i] .= z0[:, i]*sqrt((2.0+nShocks)/2.0)
    end 

    # Weights
    IntWeights = Vector{Float64}(undef, numNodes)

    IntWeights[1] = 2.0/(2.0+nShocks)
    IntWeights[2:nShocks*2+1] .= (4-nShocks)/(2*(2+nShocks)^2)
    IntWeights[nShocks*2+2:end] .= 1.0/(nShocks+2)^2

    return MonomialPowerIntegration(IntNodes, IntWeights)
end

