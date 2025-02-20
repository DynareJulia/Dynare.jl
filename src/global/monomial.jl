"""
    MonomialPowerIntegration(nShocks::Int) -> MonomialPowerIntegration

Constructs a monomial quadrature rule for numerical integration in high-dimensional stochastic models.

# Arguments
- `nShocks::Int`  
  Number of independent stochastic shocks in the model.

# Returns
- `MonomialPowerIntegration`  
  A structure containing:
  - `nodes::Vector{Vector{Float64}}` → List of integration nodes (state space points).
  - `weights::Vector{Float64}` → Corresponding integration weights for numerical expectation computation.
"""
function MonomialPowerIntegration(nShocks::Int)
    # Number of integration nodes
    numNodes = 2*nShocks^2 + 1

    # Initialize node matrix
    z0 = zeros(nShocks, numNodes)

    # Deviations in one dimension (note that the origin is row one)
    for i1 in 1:nShocks
        z0[i1, (i1 - 1)*2 + 2] = 1.0
        z0[i1, (i1 - 1)*2 + 3] = -1.0
    end
    
    idx = 0
    # Deviations in two dimensions
    for i1 in 1:nShocks
        for i2 in (i1 + 1):nShocks
            baseIdx = 2 * nShocks + 2 + idx * 4
            z0[i1, baseIdx, ] = 1.0
            z0[i1, baseIdx+1, ] = 1.0
            z0[i1, baseIdx+2, ] = -1.0
            z0[i1, baseIdx+3, ] = -1.0
            z0[i2, baseIdx, ] = 1.0
            z0[i2, baseIdx+1, ] = -1.0
            z0[i2, baseIdx+2, ] = 1.0
            z0[i2, baseIdx+3, ] = -1.0
            idx += 1
        end
    end

    # Compute integration nodes
    IntNodes = [Vector{Float64}(undef, nShocks) for _ = 1:numNodes]
    fill!(IntNodes[1], 0.0)
    for i = 2:nShocks*2+1
        IntNodes[i] .= z0[:, i]*sqrt(2.0 + nShocks)
    end
    for i = nShocks*2 + 2:numNodes
        IntNodes[i] .= z0[:, i]*sqrt((2.0+nShocks)/2.0)
    end 

    # Compute integration weights
    IntWeights = Vector{Float64}(undef, numNodes)
    IntWeights[1] = 2.0/(2.0+nShocks)
    IntWeights[2:nShocks*2+1] .= (4-nShocks)/(2*(2+nShocks)^2)
    IntWeights[nShocks*2+2:end] .= 1.0/(nShocks+2)^2

    return MonomialPowerIntegration(IntNodes, IntWeights)
end

