abstract type Proposal{P} end

struct StaticProposal{issymmetric,P} <: Proposal{P}
    proposal::P
end
const SymmetricStaticProposal{P} = StaticProposal{true,P}

StaticProposal(proposal) = StaticProposal{false}(proposal)
function StaticProposal{issymmetric}(proposal) where {issymmetric}
    return StaticProposal{issymmetric,typeof(proposal)}(proposal)
end

struct RandomWalkProposal{issymmetric,P} <: Proposal{P}
    proposal::P
end
const SymmetricRandomWalkProposal{P} = RandomWalkProposal{true,P}

RandomWalkProposal(proposal) = RandomWalkProposal{false}(proposal)
function RandomWalkProposal{issymmetric}(proposal) where {issymmetric}
    return RandomWalkProposal{issymmetric,typeof(proposal)}(proposal)
end

# Random draws
Base.rand(p::Proposal, args...) = rand(Random.default_rng(), p, args...)
Base.rand(rng::Random.AbstractRNG, p::Proposal{<:Distribution}) = rand(rng, p.proposal)
function Base.rand(rng::Random.AbstractRNG, p::Proposal{<:AbstractArray})
    return map(x -> rand(rng, x), p.proposal)
end

# Densities
Distributions.logpdf(p::Proposal{<:Distribution}, v) = logpdf(p.proposal, v)
function Distributions.logpdf(p::Proposal{<:AbstractArray}, v)
    # `mapreduce` with multiple iterators requires Julia 1.2 or later
    return mapreduce(((pi, vi),) -> logpdf(pi, vi), +, zip(p.proposal, v))
end

###############
# Random Walk #
###############

function propose(
    rng::Random.AbstractRNG,
    proposal::RandomWalkProposal{issymmetric,<:Union{Distribution,AbstractArray}},
    ::DensityModelOrLogDensityModel
) where {issymmetric}
    return rand(rng, proposal)
end

function propose(
    rng::Random.AbstractRNG,
    proposal::RandomWalkProposal{issymmetric,<:Union{Distribution,AbstractArray}},
    model::DensityModelOrLogDensityModel,
    t
) where {issymmetric}
    return t + rand(rng, proposal)
end

function q(
    proposal::RandomWalkProposal{issymmetric,<:Union{Distribution,AbstractArray}},
    t,
    t_cond
) where {issymmetric}
    return logpdf(proposal, t - t_cond)
end

##########
# Static #
##########

function propose(
    rng::Random.AbstractRNG,
    proposal::StaticProposal{issymmetric,<:Union{Distribution,AbstractArray}},
    model::DensityModelOrLogDensityModel,
    t=nothing
) where {issymmetric}
    return rand(rng, proposal)
end

function q(
    proposal::StaticProposal{issymmetric,<:Union{Distribution,AbstractArray}},
    t,
    t_cond
) where {issymmetric}
    return logpdf(proposal, t)
end

############
# Function #
############

# function definition with abstract types requires Julia 1.3 or later
for T in (:StaticProposal, :RandomWalkProposal)
    @eval begin
        function (p::$T{issymmetric,<:Function})() where {issymmetric}
            return $T{issymmetric}(p.proposal())
        end
        function (p::$T{issymmetric,<:Function})(t) where {issymmetric}
            return $T{issymmetric}(p.proposal(t))
        end
    end
end

function propose(
    rng::Random.AbstractRNG,
    proposal::Proposal{<:Function},
    model::DensityModelOrLogDensityModel
)
    return propose(rng, proposal(), model)
end

function propose(
    rng::Random.AbstractRNG,
    proposal::Proposal{<:Function}, 
    model::DensityModelOrLogDensityModel,
    t
)
    return propose(rng, proposal(t), model)
end

function q(
    proposal::Proposal{<:Function},
    t,
    t_cond
)
    return q(proposal(t_cond), t, t_cond)
end

####################
# Multiple proposals
####################

function propose(
    rng::Random.AbstractRNG,
    proposals::AbstractArray{<:Proposal},
    model::DensityModelOrLogDensityModel,
)
    return map(proposals) do proposal
        return propose(rng, proposal, model)
    end
end
function propose(
    rng::Random.AbstractRNG,
    proposals::AbstractArray{<:Proposal},
    model::DensityModelOrLogDensityModel,
    ts,
)
    return map(proposals, ts) do proposal, t
        return propose(rng, proposal, model, t)
    end
end

@generated function propose(
    rng::Random.AbstractRNG,
    proposals::NamedTuple{names},
    model::DensityModelOrLogDensityModel,
) where {names}
    isempty(names) && return :(NamedTuple())
    expr = Expr(:tuple)
    expr.args = Any[:($name = propose(rng, proposals.$name, model)) for name in names]
    return expr
end

@generated function propose(
    rng::Random.AbstractRNG,
    proposals::NamedTuple{names},
    model::DensityModelOrLogDensityModel,
    ts,
) where {names}
    isempty(names) && return :(NamedTuple())
    expr = Expr(:tuple)
    expr.args = Any[
        :($name = propose(rng, proposals.$name, model, ts.$name)) for name in names
    ]
    return expr
end

"""
    logratio_proposal_density(proposal, state, candidate)

Compute the log-ratio of the proposal densities in the Metropolis-Hastings algorithm.

The log-ratio of the proposal densities is defined as
```math
\\log \\frac{g(x | x')}{g(x' | x)},
```
where ``x`` is the current state, ``x'`` is the proposed candidate for the next state,
and ``g(y' | y)`` is the conditional probability of proposing state ``y'`` given state
``y`` (proposal density).
"""
function logratio_proposal_density(proposal::Proposal, state, candidate)
    return q(proposal, state, candidate) - q(proposal, candidate, state)
end

# ratio is always 0 for symmetric proposals
logratio_proposal_density(::RandomWalkProposal{true}, state, candidate) = 0
logratio_proposal_density(::StaticProposal{true}, state, candidate) = 0

# type stable implementation for `NamedTuple`s
function logratio_proposal_density(
    proposals::NamedTuple{names}, states::NamedTuple, candidates::NamedTuple
) where {names}
    if @generated
        args = map(names) do name
            :(logratio_proposal_density(
                proposals[$(QuoteNode(name))],
                states[$(QuoteNode(name))],
                candidates[$(QuoteNode(name))],
            ))
        end
        return :(+($(args...)))
    else
        return sum(names) do name
            return logratio_proposal_density(
                proposals[name], states[name], candidates[name]
            )
        end
    end
end

# use recursion for `Tuple`s to ensure type stability
logratio_proposal_density(proposals::Tuple{}, states::Tuple, candidates::Tuple) = 0
function logratio_proposal_density(
    proposals::Tuple{<:Proposal}, states::Tuple, candidates::Tuple
)
    return logratio_proposal_density(first(proposals), first(states), first(candidates))
end
function logratio_proposal_density(proposals::Tuple, states::Tuple, candidates::Tuple)
    valfirst = logratio_proposal_density(first(proposals), first(states), first(candidates))
    valtail = logratio_proposal_density(
        Base.tail(proposals), Base.tail(states), Base.tail(candidates)
    )
    return valfirst + valtail
end

# fallback for general iterators (arrays etc.) - possibly not type stable!
function logratio_proposal_density(proposals, states, candidates)
    return sum(zip(proposals, states, candidates)) do (proposal, state, candidate)
        return logratio_proposal_density(proposal, state, candidate)
    end
end
