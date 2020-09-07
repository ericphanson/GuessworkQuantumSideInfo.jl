## Notation
"""
    ket([T = Float64], i::Integer, d::Integer) -> SparseVector{Complex{T}}

Create a vector representing the `i`th computational basis vector in dimension `d`.

## Example

```jldoctest
julia> ket(1,2)
2-element SparseVector{Complex{Float64},Int64} with 1 stored entry:
  [1]  =  1.0+0.0im

julia> collect(ans)
2-element Array{Complex{Float64},1}:
 1.0 + 0.0im
 0.0 + 0.0im
```
"""
ket(::Type{T}, i::Integer, d::Integer) where {T} = sparsevec([i], [one(complex(T))], d)

"""
    bra([T = Float64], i::Integer, d::Integer) -> SparseVector{Complex{T}}'

Create a dual vector representing the bra associated to `i`th computational basis vector in dimension `d`.

## Example

```jldoctest
julia> bra(1,2)
1×2 LinearAlgebra.Adjoint{Complex{Float64},SparseVector{Complex{Float64},Int64}}:
 1.0-0.0im  0.0-0.0im

julia> collect(ans)
1×2 Array{Complex{Float64},2}:
 1.0-0.0im  0.0-0.0im
```
"""
bra(::Type{T}, i::Integer, d::Integer) where {T} = ket(T, i, d)'

ket(i::Integer, d::Integer) = ket(Float64, i, d)
bra(i::Integer, d::Integer) = bra(Float64, i, d)

"""
    dm(ψ::AbstractVector)

Creates the density matrix version of a pure state `ψ` via the outer product.
"""
dm(ψ::AbstractVector) = ψ' ⊗ ψ

if VERSION < v"1.2-pre"
    (I::UniformScaling)(n::Integer) = Diagonal(fill(I.λ, n))
end

const ⊗ = kron

"""
    BB84_states([T::Type = Float64])

Generates the BB84 states `|0⟩`, `|1⟩`, `|-⟩`, and `|+⟩`, for use in [`guesswork`](@ref) or other functions. The
numeric type can be optionally specified by the argument.
"""
function BB84_states end

BB84_states() = BB84_states(Float64)

function BB84_states(::Type{T}) where {T<:Number}
    ketplus = (ket(T, 1, 2) + ket(T, 2, 2)) / sqrt(T(2))
    ketminus = (ket(T, 1, 2) - ket(T, 2, 2)) / sqrt(T(2))
    ketzero = ket(T, 1, 2)
    ketone = ket(T, 2, 2)
    return dm.([ketzero, ketone, ketminus, ketplus])
end

"""
    iid_copies(ρBs::AbstractVector{<:AbstractMatrix}, n::Integer) -> Vector{Matrix}

Create a vector of all states of the form ``ρ_1 \\otimes \\dotsm \\otimes ρ_n``
where the ``ρ_i`` range over the set `ρBs`.
"""
function iid_copies(ρBs::AbstractVector{<:AbstractMatrix}, n::Integer)
    inds = cartesian_power(length(ρBs), n)
    [foldl(⊗, ρBs[collect(I)]) for I in inds]
end

## Random states

"""
    simplexpt(unif)

Takes a vector of length `d-1` of numbers between `0.0` and `1.0` and converts it a point on the standard `d` dimensional simplex.
"""
function simplexpt(unif)
    d = length(unif) + 1
    T = eltype(unif)
    w = zeros(T, d + 1)
    w[2:d] .= sort(unif)
    w[d+1] = one(T)
    diff(w)
end

"""
    randprobvec([rng, T=Float64], d)

Generates points of type `T`, uniformly at random on the standard `d-1` dimensional simplex using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).

## Example

```julia
julia> randprobvec(3)
3-element Array{Float64,1}:
 0.24815974900033688
 0.17199716455672287
 0.5798430864429402 

```
"""
randprobvec(rng::AbstractRNG, ::Type{T}, d::Integer) where T = simplexpt(rand(rng, T, d - 1))
randprobvec(::Type{T}, d::Integer) where T = randprobvec(Random.GLOBAL_RNG, T, d)
randprobvec(d::Integer) = randprobvec(Float64, d)

"""
    randunitary([rng, T=Float64], d)

Generates a unitary matrix of dimension `d` at random according to the Haar measure, using an algorithm described by Maris Ozols in ["How to generate a random unitary matrix"](http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20%5Bpaper%5D.pdf).
"""
function randunitary(rng::AbstractRNG, ::Type{T}, d::Integer) where T
    # `randn` for `BigFloat`'s isn't supported yet
    # (https://github.com/JuliaLang/julia/issues/17629)
    # so we convert afterwards instead (losing some randomness).
    rg1 = T.(randn(rng, d, d))
    rg2 = T.(randn(rng, d, d))
    RG = rg1 + im * rg2
    Q, R = qr(RG)
    r = diag(R)
    L = diagm(0 => r ./ abs.(r))
    return Q * L
end

randunitary(::Type{T}, d::Integer) where {T} = randunitary(Random.GLOBAL_RNG, T, d)
randunitary(d::Integer) = randunitary(Float64, d)

"""
    randdm([rng, T = Float64], d)

Generates a density matrix with numeric type `Complex{T}`, of dimension `d` at
random.

## Example

```julia
julia> randdm(2)
2×2 Array{Complex{Float64},2}:
 0.477118+0.0im        0.119848-0.0371569im
 0.119848+0.0371569im  0.522882+0.0im      
```
"""
function randdm(rng::AbstractRNG, ::Type{T}, d::Integer) where T
    eigs = diagm(0 => randprobvec(rng, T, d))
    U = randunitary(rng, T, d)
    ρ = U * eigs * (U')
    return Matrix((ρ + ρ') / 2)
end

randdm(::Type{T}, d::Integer) where T = randdm(Random.GLOBAL_RNG, T, d)
randdm(d) = randdm(Float64, d)
