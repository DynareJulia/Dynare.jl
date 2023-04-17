struct Ensemble{D} <: MHSampler
    n_walkers::Int
    proposal::D
end

function transition(sampler::Ensemble, model::DensityModelOrLogDensityModel, params)
    return [Transition(model, x) for x in params]
end

# Define the other sampling steps.
# Return a 2-tuple consisting of the next sample and the the next state.
# In this case they are identical, and for each walker they are either a new proposal
# (if accepted) or the previous proposal (if not accepted).
function AbstractMCMC.step(
    rng::Random.AbstractRNG,
    model::DensityModelOrLogDensityModel,
    spl::Ensemble,
    params_prev::Vector{<:Transition};
    kwargs...,
)
    # Generate a new proposal. Accept/reject happens at proposal level.
    transitions = propose(rng, spl, model, params_prev)
    return transitions, transitions
end

#
# Initial proposal
#
function propose(rng::Random.AbstractRNG, spl::Ensemble, model::DensityModelOrLogDensityModel)
    # Make the first proposal with a static draw from the prior.
    static_prop = StaticProposal(spl.proposal.proposal)
    mh_spl = MetropolisHastings(static_prop)
    return [propose(rng, mh_spl, model) for _ in 1:spl.n_walkers]
end

#
# Every other proposal
#
function propose(
    rng::Random.AbstractRNG,
    spl::Ensemble,
    model::DensityModelOrLogDensityModel,
    walkers::Vector{<:Transition},
)
    new_walkers = similar(walkers)

    n_walkers = spl.n_walkers
    uniform_sampler = Random.Sampler(rng, 1:(n_walkers - 1))

    for i in 1:n_walkers
        walker = walkers[i]
        idx = mod1(i + rand(rng, uniform_sampler), n_walkers)
        other_walker = idx < i ? new_walkers[idx] : walkers[idx]
        new_walkers[i] = move(rng, spl, model, walker, other_walker)
    end

    return new_walkers
end

#####################################
# Basic stretch move implementation #
#####################################
struct StretchProposal{P, F<:AbstractFloat} <: Proposal{P}
    proposal :: P
    stretch_length::F
end

StretchProposal(p) = StretchProposal(p, 2.0)

function move(
    rng::Random.AbstractRNG, 
    spl::Ensemble{<:StretchProposal},
    model::DensityModelOrLogDensityModel,
    walker::Transition,
    other_walker::Transition,
)
    # Calculate intermediate values
    proposal = spl.proposal
    n = length(walker.params)
    a = proposal.stretch_length
    z = ((a - 1) * rand(rng) + 1)^2 / a
    alphamult = (n - 1) * log(z)

    # Make new parameters
    y = @. other_walker.params + z * (walker.params - other_walker.params)

    # Construct a new walker
    new_walker = Transition(model, y)

    # Calculate accept/reject value.
    alpha = alphamult + new_walker.lp - walker.lp

    if -Random.randexp(rng) <= alpha
        return new_walker
    else
        return walker
    end
end

#########################
# Elliptical slice step #
# #########################

# struct EllipticalSlice{E} <: ProposalStyle
#     ellipse::E
# end

# function move(
#     # spl::Ensemble,
#     spl::Ensemble{Proposal{T,P}},
#     model::DensityModel,
#     walker::Transition,
#     other_walker::Transition,
# ) where {T<:EllipticalSlice,P}
#     # Calculate intermediate values
#     proposal = spl.proposal
#     n = length(walker.params)
#     nu = rand(proposal.type.ellipse)

#     u = rand()
#     y = walker.lp - Random.randexp()

#     theta = 2 * π * rand()

#     theta_min = theta - 2.0*π
#     theta_max = theta
    
#     f = walker.params
#     while true
#         stheta, ctheta = sincos(theta)

#         f_prime = f .* ctheta + nu .* stheta

#         new_walker = Transition(model, f_prime)

#         if new_walker.lp > y
#             return new_walker
#         else
#             if theta < 0 
#                 theta_min = theta
#             else
#                 theta_max = theta
#             end

#             theta = theta_min + (theta_max - theta_min) * rand()
#         end
#     end 
# end

#####################
# Slice and stretch #
#####################
# struct EllipticalSliceStretch{E, S<:Stretch} <: ProposalStyle
#     ellipse::E
#     stretch::S
# end

# EllipticalSliceStretch(e) = EllipticalSliceStretch(e, Stretch(2.0))

# function move(
#     # spl::Ensemble,
#     spl::Ensemble{Proposal{T,P}},
#     model::DensityModel,
#     walker::Transition,
#     other_walker::Transition,
# ) where {T<:EllipticalSliceStretch,P}
#     # Calculate intermediate values
#     proposal = spl.proposal
#     n = length(walker.params)
#     nu = rand(proposal.type.ellipse)

#     # Calculate stretch step first
#     subspl = Ensemble(spl.n_walkers, Proposal(proposal.type.stretch, proposal.proposal))
#     walker = move(subspl, model, walker, other_walker)

#     u = rand()
#     y = walker.lp - Random.randexp()

#     theta = 2 * π * rand()

#     theta_min = theta - 2.0*π
#     theta_max = theta
    
#     f = walker.params

#     i = 0
#     while true
#         i += 1
        
#         stheta, ctheta = sincos(theta)
        
#         f_prime = f .* ctheta + nu .* stheta

#         new_walker = Transition(model, f_prime)

#         # @info "Slice step" i f f_prime y new_walker.lp theta theta_max theta_min

#         if new_walker.lp > y
#             return new_walker
#         else
#             if theta < 0 
#                 theta_min = theta
#             else
#                 theta_max = theta
#             end

#             theta = theta_min + (theta_max - theta_min) * rand()
#         end
#     end 
# end
