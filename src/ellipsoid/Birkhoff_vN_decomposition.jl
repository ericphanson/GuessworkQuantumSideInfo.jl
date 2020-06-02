"""
    PermutationIterator(D::AbstractMatrix; timer=nothing, rtol=1e-6)

Iterates over permutations and weights to form a decomposition
of a doubly stochastic matrix `D`, known as the Birkhoff–von Neumann
decomposition.

Optionally, specify the keyword arguments

* `timer`: a TimerOutputs.jl timer
* `rtol`: (relative tolerance) for comparing `α_i` to `1` to see if all the permutations have been found.

## Example

```julia
# Construct a `d` by `d` doubly stochastic matrix
d = 4
α_init = rand(5)
α_init = α_init / sum(α_init)
D = sum( α_init[i]*I(d)[randperm(d), :] for i = 1:5 )

# Reconstruct it from permutations
D_reconstruct = zero(D)
for (π_i, α_i) in PermutationIterator(D)
    P_i = I(d)[π_i, :]
    global D_reconstruct += α_i * P_i
end
@test D_reconstruct ≈ D
```
"""
@Base.kwdef struct PermutationIterator{T1, T2, TTimer, V}
    D::T1
    αs::V = eltype(D)[]
    rtol::T2 = 1e-6
    timer::TTimer = nothing
end

PermutationIterator(D::AbstractMatrix; kwargs...) = PermutationIterator(; D=D, kwargs...)

Base.IteratorSize(::PermutationIterator) = Base.SizeUnknown()
Base.IteratorEltype(::PermutationIterator) = Base.HasEltype()
Base.eltype(::PermutationIterator) = Tuple{Vector{Int}, Float64}

function Base.iterate(PI::PermutationIterator)
    state = PI.D
    empty!(PI.αs)
    iterate(PI, state)
end

function Base.iterate(PI::PermutationIterator, state)
    state === nothing && return nothing
    D = state
    @unpack rtol, timer, αs = PI

    # Find the maximum weight permutation
    if timer === nothing
        π = maximum_weight_perm(D)
    else
        @timeit timer "maximum_weight_perm" begin
        π = maximum_weight_perm(D)
        end
    end

    n = size(D,1)
    α0 = minimum( D[i,π[i]] for i = 1:n )
    c = isempty(αs) ? 1 : prod( 1 - α for α in αs)
    push!(αs, α0)

    if isapprox(α0, 1; rtol=rtol)
        # This was the last permutation
        return (π, c*α0), nothing
    else
        # Prepare the next doubly stochastic matrix
        # solve for `D1` so that
        # `D = α0 P + (1-α0) D1`
        D1 = copy(D)
        for i = 1:n
            D1[i, π[i]] -= α0
        end
        return (π, c*α0), D1/(1-α0)
    end
end

function maximum_weight_perm(D)
    assignment, cost = hungarian(-D)
    assignment
end
