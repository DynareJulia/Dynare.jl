# AdvancedMH.jl

AdvancedMH.jl currently provides a robust implementation of random walk Metropolis-Hastings samplers.

Further development aims to provide a suite of adaptive Metropolis-Hastings implementations.

AdvancedMH works by allowing users to define composable `Proposal` structs in different formats.

## Usage

First, construct a `DensityModel`, which is a wrapper around the log density function for your inference problem. The `DensityModel` is then used in a `sample` call.

```julia
# Import the package.
using AdvancedMH
using Distributions
using MCMCChains

using LinearAlgebra

# Generate a set of data from the posterior we want to estimate.
data = rand(Normal(0, 1), 30)

# Define the components of a basic model.
insupport(θ) = θ[2] >= 0
dist(θ) = Normal(θ[1], θ[2])
density(θ) = insupport(θ) ? sum(logpdf.(dist(θ), data)) : -Inf

# Construct a DensityModel.
model = DensityModel(density)

# Set up our sampler with a joint multivariate Normal proposal.
spl = RWMH(MvNormal(zeros(2), I))

# Sample from the posterior.
chain = sample(model, spl, 100000; param_names=["μ", "σ"], chain_type=Chains)
```

Output:

```julia
Object of type Chains, with data of type 100000×3×1 Array{Float64,3}

Iterations        = 1:100000
Thinning interval = 1
Chains            = 1
Samples per chain = 100000
internals         = lp
parameters        = μ, σ

2-element Array{ChainDataFrame,1}

Summary Statistics

│ Row │ parameters │ mean     │ std      │ naive_se    │ mcse       │ ess     │ r_hat   │
│     │ Symbol     │ Float64  │ Float64  │ Float64     │ Float64    │ Any     │ Any     │
├─────┼────────────┼──────────┼──────────┼─────────────┼────────────┼─────────┼─────────┤
│ 1   │ μ          │ 0.156152 │ 0.19963  │ 0.000631285 │ 0.00323033 │ 3911.73 │ 1.00009 │
│ 2   │ σ          │ 1.07493  │ 0.150111 │ 0.000474693 │ 0.00240317 │ 3707.73 │ 1.00027 │

Quantiles

│ Row │ parameters │ 2.5%     │ 25.0%     │ 50.0%    │ 75.0%    │ 97.5%    │
│     │ Symbol     │ Float64  │ Float64   │ Float64  │ Float64  │ Float64  │
├─────┼────────────┼──────────┼───────────┼──────────┼──────────┼──────────┤
│ 1   │ μ          │ -0.23361 │ 0.0297006 │ 0.159139 │ 0.283493 │ 0.558694 │
│ 2   │ σ          │ 0.828288 │ 0.972682  │ 1.05804  │ 1.16155  │ 1.41349  │

```

### Usage with [`LogDensityProblems.jl`](https://github.com/tpapp/LogDensityProblems.jl)

Alternatively, you can define your model with the [`LogDensityProblems.jl`](https://github.com/tpapp/LogDensityProblems.jl) interface:

``` julia
using LogDensityProblems

# Use a struct instead of `typeof(density)` for sake of readability.
struct LogTargetDensity end

LogDensityProblems.logdensity(p::LogTargetDensity, θ) = density(θ)  # standard multivariate normal
LogDensityProblems.dimension(p::LogTargetDensity) = 2
LogDensityProblems.capabilities(::LogTargetDensity) = LogDensityProblems.LogDensityOrder{0}()

sample(LogTargetDensity(), spl, 100000; param_names=["μ", "σ"], chain_type=Chains)
```

## Proposals

AdvancedMH offers various methods of defining your inference problem. Behind the scenes, a `MetropolisHastings` sampler simply holds
some set of `Proposal` structs. AdvancedMH will return posterior samples in the "shape" of the proposal provided -- currently
supported methods are `Array{Proposal}`, `Proposal`, and `NamedTuple{Proposal}`. For example, proposals can be created as:

```julia
# Provide a univariate proposal.
m1 = DensityModel(x -> logpdf(Normal(x,1), 1.0))
p1 = StaticProposal(Normal(0,1))
c1 = sample(m1, MetropolisHastings(p1), 100; chain_type=Vector{NamedTuple})

# Draw from a vector of distributions.
m2 = DensityModel(x -> logpdf(Normal(x[1], x[2]), 1.0))
p2 = StaticProposal([Normal(0,1), InverseGamma(2,3)])
c2 = sample(m2, MetropolisHastings(p2), 100; chain_type=Vector{NamedTuple})

# Draw from a `NamedTuple` of distributions.
m3 = DensityModel(x -> logpdf(Normal(x.a, x.b), 1.0))
p3 = (a=StaticProposal(Normal(0,1)), b=StaticProposal(InverseGamma(2,3)))
c3 = sample(m3, MetropolisHastings(p3), 100; chain_type=Vector{NamedTuple})

# Draw from a functional proposal.
m4 = DensityModel(x -> logpdf(Normal(x,1), 1.0))
p4 = StaticProposal((x=1.0) -> Normal(x, 1))
c4 = sample(m4, MetropolisHastings(p4), 100; chain_type=Vector{NamedTuple})
```

## Static vs. Random Walk

Currently there are only two methods of inference available. Static MH simply draws from the prior, with no
conditioning on the previous sample. Random walk will add the proposal to the previously observed value.
If you are constructing a `Proposal` by hand, you can determine whether the proposal is a
`StaticProposal` or a `RandomWalkProposal` using

```julia
static_prop = StaticProposal(Normal(0,1))
rw_prop = RandomWalkProposal(Normal(0,1))
```

Different methods are easily composeable. One parameter can be static and another can be a random walk,
each of which may be drawn from separate distributions.

## Multiple chains

AdvancedMH.jl implements the interface of [AbstractMCMC](https://github.com/TuringLang/AbstractMCMC.jl/) which means sampling of multiple chains is supported for free:

```julia
# Sample 4 chains from the posterior serially, without thread or process parallelism.
chain = sample(model, RWMH(init_params), MCMCSerial(), 100000, 4; param_names=["μ","σ"], chain_type=Chains)

# Sample 4 chains from the posterior using multiple threads.
chain = sample(model, RWMH(init_params), MCMCThreads(), 100000, 4; param_names=["μ","σ"], chain_type=Chains)

# Sample 4 chains from the posterior using multiple processes.
chain = sample(model, RWMH(init_params), MCMCDistributed(), 100000, 4; param_names=["μ","σ"], chain_type=Chains)
```

## Metropolis-adjusted Langevin algorithm (MALA)

AdvancedMH.jl also offers an implementation of [MALA](https://en.wikipedia.org/wiki/Metropolis-adjusted_Langevin_algorithm) if the `ForwardDiff` and `DiffResults` packages are available. 

A `MALA` sampler can be constructed by `MALA(proposal)` where `proposal` is a function that
takes the gradient computed at the current sample. It is required to specify an initial sample `init_params` when calling `sample`.

```julia
# Import the package.
using AdvancedMH
using Distributions
using MCMCChains
using ForwardDiff
using StructArrays

using LinearAlgebra

# Generate a set of data from the posterior we want to estimate.
data = rand(Normal(0, 1), 30)

# Define the components of a basic model.
insupport(θ) = θ[2] >= 0
dist(θ) = Normal(θ[1], θ[2])
density(θ) = insupport(θ) ? sum(logpdf.(dist(θ), data)) : -Inf

# Construct a DensityModel.
model = DensityModel(density)

# Set up the sampler with a multivariate Gaussian proposal.
σ² = 0.01
spl = MALA(x -> MvNormal((σ² / 2) .* x, σ² * I))

# Sample from the posterior.
chain = sample(model, spl, 100000; init_params=ones(2), chain_type=StructArray, param_names=["μ", "σ"])
```

### Usage with [`LogDensityProblems.jl`](https://github.com/tpapp/LogDensityProblems.jl)

As above, we can define the model with the LogDensityProblems.jl interface.
We can implement the gradient of the log density function manually, or use [`LogDensityProblemsAD.jl`](https://github.com/tpapp/LogDensityProblemsAD.jl) to provide us with the gradient computation used in MALA.
Using our implementation of the `LogDensityProblems.jl` interface above:

```julia
using LogDensityProblemsAD
model_with_ad = LogDensityProblemsAD.ADgradient(Val(:ForwardDiff), LogTargetDensity())
sample(model_with_ad, spl, 100000; init_params=ones(2), chain_type=StructArray, param_names=["μ", "σ"])
```
