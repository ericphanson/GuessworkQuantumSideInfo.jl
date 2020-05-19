
function make_ellipsoid_problem(p::AbstractVector{T}, ρBs::AbstractVector{<:AbstractMatrix}, c::AbstractVector = T.(1:length(p));
    max_retries = 50,
    max_time = Inf,
    num_constraints = Inf,
    verbose::Bool = false,
    num_steps_per_SA_run::Integer = length(p)^2 * 500,
    mutate! = rand_rev!,
    debug = false,) where {T}

    length(p) == length(ρBs) ||
    throw(ArgumentError("Length of prior and vector of side information must match J"))
    J = length(p)
    dB = size(ρBs[1], 1)
    all(ρB -> size(ρB) == (dB, dB), ρBs) ||
    throw(ArgumentError("All side-information states must be square matrices of the same size"))

    if (length(c) < length(p))
        throw(ArgumentError("Need `length(c) >= length(p)`."))
    end

    EllipsoidProblem(p, ρBs, dB, c, max_retries, max_time, verbose, num_steps_per_SA_run, mutate!, debug)
end

struct EllipsoidProblem{T, TρBs, Tc, Tm}
    p::Vector{T}
    ρBs::TρBs
    dB::Int
    c::Tc
    max_retries::Int
    max_time::Float64 
    verbose::Bool
    num_steps_per_SA_run::Integer
    mutate!::Tm
    debug::Bool 
end


function find_violation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, dB, ρBs, max_retries, num_steps_per_SA_run, mutate! = prob
    J = length(p)

    R = g -> sum(c[N(invperm(g), x)] * p[x] * ρBs[x] for x in eachindex(p, ρBs))

    # We optimize over the inverse of the permutations we are interested in
    # that way we don't need to take the inverse until the end.
    Rg_inv! = let ρBs_subnormalized = [p[x] * ρBs[x] for x = 1:J], c = c
        (R, g_inv) -> begin
            for x = 1:J
                R .+= ρBs_subnormalized[x] .* c[g_inv[x]]
            end
            R
        end
    end

    # Initialize some variables
    π = collect(1:J)
    scratch_π = copy(π)
    R_scratch::Matrix{Complex{T}} = zeros(Complex{T}, dB, dB)
    total_time_so_far::Float64 = 0.0

    f = let Y::Matrix{Complex{T}} = Y
        π -> begin
            R_scratch .= -1 .* Y
            eigmin(Hermitian(Rg_inv!(R_scratch, π)))::T
        end
    end


    # Choose `fval` as something `> 0` to enter the loop.
    fval::T = 1.0
    n_tries::Int = 0
    SA_time = @elapsed while fval > 0 && n_tries < max_retries
        n_tries += 1
        # initialize `π` to hold a random permutation, then use the simulated annealing algorithm to try to minimize `f`, starting from this permutation. Repeat until we find a violated constraint, or we've tried enough times.
        randperm!(π)
        fval = SA!(π, scratch_π, f, mutate!, num_steps_per_SA_run)
    end

    if fval < 0
        new_constraint = copy(π)
        return R(new_constraint) - Y
    else
        return zero(Y)
    end
end



using LinearAlgebra
const tol = Ref(1e-4)



function ellipsoid_solve(f::EllipsoidProblem, x, P, ϵ)
    dB = f.dB
    n = length(x)
    stepsize = inv(n+1)
    iter = 0
    f_best = Inf
    l_best = -Inf

    while true
        iter += 1
        g, fval, feasible = subgradient(f, x)
        γ = sqrt(dot(g, P, g))
        @info "Iteration" iter γ fval f_best l_best
        fval = real(fval)
        γ = real(γ)
        if feasible
            f_best = min(f_best, fval)
            l_best = max(l_best, fval - γ)
        end
        feasible && γ ≤ ϵ && return (x=x, P=P, Y = herm(x))
        g̃ = (1/γ)*g

        Pg̃ = P * g̃
        x = x - stepsize*Pg̃
        P = (n^2/(n^2-1)) *(P - (2*stepsize)* Pg̃ * transpose(Pg̃))
    end
end

struct Orth end
struct PSD end

function normal_cone_element(::Orth, λ)
    z = zeros(length(λ))
    for i in eachindex(λ)
        if real(λ[i]) <= tol[]
            z[i] = -1
        end
    end
    return z
end

function normal_cone_element(::PSD, P)
    λ, Q = eigen(P)
    Q * Diagonal(normal_cone_element(Orth(), λ)) * Q'
end


function subgradient(f::EllipsoidProblem, x)
    dB = f.dB
    Y = herm(x)
    if !(Y ≈ Y')
        @show Y - Y'
        error()
    end
    V = find_violation(f, Y)
    if V != zero(V)
        elt = normal_cone_element(PSD(), V)
        @info "Found violation" V
        @info "" eigvals(V)
        @info "" elt
        @info "" eigvals(elt)
        return -invherm(elt), Inf, false
    else
        return invherm(-I(dB)), -tr(Y), true
    end
    
end



# helper for `herm`
function f(x, y, i, j)
    if i == j
        return x/2 + im*zero(y)
    elseif i > j
        return x + im * y
    elseif i < j
        return zero(x) + im*zero(y)
    end
end

# helper for invherm
function finv(z, i, j)
    if i <= j
        return real(z)
    elseif i > j
        return imag(z)
    end
end

# Create a complex `d` by `d` Hermitian matrix from a `d^2`-dimensional real vector `v`
function herm(v::AbstractVector)
    d = isqrt(length(v))
    @assert d^2 == length(v)
    m = reshape(v, d, d)
    m = [ f(m[i,j], m[j,i], i, j) for i = 1:d, j = 1:d]
    m = m + m'
    return m
end

# Recover the vector representation
function invherm(A::AbstractMatrix)
    d = size(A,1)
    m = [ finv(A[i,j], i, j) for i = 1:d for j = 1:d ]
    return vec(m)
end
