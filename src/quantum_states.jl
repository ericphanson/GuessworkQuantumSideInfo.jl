## Notation
ketbra(::Type{T}, i::Integer, j::Integer, d::Integer) where {T} =
    sparse([i], [j], [one(complex(T))], d, d)
ket(::Type{T}, i::Integer, d::Integer) where {T} = sparsevec([i], [one(complex(T))], d)
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

Generates the BB84 states for use in [`guesswork`](@ref) or other functions. The
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

six_states() = six_states(Float64)

function six_states(::Type{T}) where {T<:Number}
    ketplus = (ket(T, 1, 2) + ket(T, 2, 2)) / sqrt(T(2))
    ketminus = (ket(T, 1, 2) - ket(T, 2, 2)) / sqrt(T(2))
    ketzero = ket(T, 1, 2)
    ketone = ket(T, 2, 2)
    ketyplus = ketzero + im * ketone
    ketyminus = ketzero - im * ketone
    return ρBs = dm.([ketzero, ketone, ketminus, ketplus, ketyplus, ketyminus])
end

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
    randprobvec([T=Float64], d)

Generates points of type `T`, uniformly at random on the standard `d-1` dimensional simplex using an algorithm by [Smith and Tromble](http://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf).
"""
randprobvec(T::Type, d::Integer) = simplexpt(rand(T, d - 1))
randprobvec(d::Integer) = randprobvec(Float64, d)

"""
    randunitary([T=Float64], d)

Generates a unitary matrix of dimension `d` at random according to the Haar measure, using an algorithm described by Maris Ozols in ["How to generate a random unitary matrix"](http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20%5Bpaper%5D.pdf).
"""
function randunitary(T::Type, d::Integer)
    rg1 = T.(randn(d, d))
    rg2 = T.(randn(d, d))
    RG = rg1 + im * rg2
    Q, R = qr(RG)
    r = diag(R)
    L = diagm(0 => r ./ abs.(r))
    return Q * L
end

randunitary(d::Integer) = randunitary(Float64, d)

"""
    randdm(T, d)

Generates a density matrix of dimension `d` at random.
"""
function randdm(T::Type, d::Integer)
    eigs = diagm(0 => randprobvec(T, d))
    U = randunitary(T, d)
    ρ = U * eigs * (U')
    return Matrix((ρ + ρ') / 2)
end
randdm(d) = randdm(Float64, d)

"""
    randpure(T, d)

Generates a pure state density matrix of dimension `d` at random.
"""
function randpure(T::Type, d::Integer)
    U = randunitary(T, d)
    x = T.(randn(d))
    normalize!(x)
    return dm(U * x)
end
randpure(d::Integer) = randpure(Float64, d)
