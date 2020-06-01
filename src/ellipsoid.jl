function make_ellipsoid_problem(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    c::AbstractVector = T.(1:length(p)),
    max_retries = 50,
    max_time = Inf,
    num_constraints = Inf,
    verbose::Bool = false,
    num_steps_per_SA_run::Integer = length(p)^2 * 500,
    mutate! = rand_rev!,
    debug = false,
) where {T}

    length(p) == length(ρBs) ||
        throw(ArgumentError("Length of prior and vector of side information must match J"))
    J = length(p)
    dB = size(ρBs[1], 1)
    all(ρB -> size(ρB) == (dB, dB), ρBs) ||
        throw(ArgumentError("All side-information states must be square matrices of the same size"))

    if (length(c) < length(p))
        throw(ArgumentError("Need `length(c) >= length(p)`."))
    end

    EllipsoidProblem(
        p,
        ρBs,
        dB,
        c,
        max_retries,
        max_time,
        verbose,
        num_steps_per_SA_run,
        mutate!,
        debug,
    )
end

struct EllipsoidProblem{T,TρBs,Tc,Tm}
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

@enum FeasibilityStates FEASIBLE INFEASIBLE UNKNOWN

function simplex_relaxation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, ρBs = prob
    J = length(p)
    c_tot = sum(c)
    c_max = c_tot - (J - 1) * c[1]
    for x = 1:J
        extremal_pt = p[x] * ρBs[x] * c_max + sum(p[y] * ρBs[y] * c[1] for y = 1:J if y ≠ x)
        if eigmin(Hermitian(extremal_pt - Y)) < 0
            # simplex relaxation is not feasible; discrete problem may or may not be
            return UNKNOWN
        end
    end
    # simplex relaxation is feasible, so discrete problem is too
    return FEASIBLE
end

function find_violation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, dB, ρBs, max_retries, num_steps_per_SA_run, mutate! = prob
    J = length(p)

    # Easy check for feasibility:
    # We want to know if `Y <= R(π)` for all
    # permutations `π`. But since `R(π) >= ρB*c[1]` for any `π`,
    # as a quick check we can see if `Y <= ρB*c[1]`.
    ρB = sum(p[x] * ρBs[x] for x in eachindex(p))
    if eigmin(Hermitian(c[1] * ρB - Y)) >= 0
        return FEASIBLE, zero(Y)
    end

    # Strictly stronger check, but takes more time to compute (`J` eigmin's rather than 1)
    if simplex_relaxation(prob, Y) == FEASIBLE
        return FEASIBLE, zero(Y)
    end

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
        # found a cut, definitely infeasible
        return INFEASIBLE, R(new_constraint) - Y
    else
        # could not find a cut but also did not prove feasibility
        return UNKNOWN, zero(Y)
    end
end



using LinearAlgebra
const tol = Ref(1e-4)

function default_init(f::EllipsoidProblem{T}) where T
    @unpack p, ρBs, c, dB = f
    ρB = sum(p[x] * ρBs[x] for x in eachindex(p))
    Y_0 = c[1] * ρB
    x_0 = invherm(Y_0)
    
    # In semidefinite order,
    # 0 <= Y - Y_0 <= c[end]*I(dB^2) - Y_0
    Y_bound = opnorm(c[end]*I(dB) - Y_0)

    # Hence `||Y-Y_0||_2 <= dB^2 * ||Y - Y_0||_∞ = dB^2*Y_bound`
    P_0 = Y_bound * dB^2 * Matrix{T}(I, dB^2, dB^2)

    return (x = x_0, P = P_0)
end

function guesswork_ellipsoid(p, ρBs; tol=1e-3, kwargs...)
    prob = make_ellipsoid_problem(p, ρBs; kwargs...)
    results = ellipsoid_algorithm(prob; ϵ=tol, default_init(prob)...)
    (optval = tr(results.Y), Y = results.Y, p = p, ρBs = ρBs, c = prob.c, J = length(p), K = length(p))
end

using SpecialFunctions
function unit_ball_volume(n)
    π^(n/2) / gamma(n/2 + 1)
end


function ellipsoid_algorithm(
    f::EllipsoidProblem{T};
    x,
    P,
    ϵ,
) where {T}
    dB = f.dB
    n = length(x)
    stepsize = inv(n + 1)
    iter = 0
    f_best = Inf
    l_best = -Inf
    vol0 = unit_ball_volume(n)
    while true
        iter += 1
        g, fval, feasible = subgradient(f, x)
        γ = sqrt(dot(g, P, g))
        if f.verbose
            vol = vol0 * sqrt(det(P))
            @info "Iteration" iter vol γ fval f_best l_best
        end
        fval = real(fval)
        γ = real(γ)
        if feasible
            f_best = min(f_best, fval)
            l_best = max(l_best, fval - γ)
        end
        feasible && γ ≤ ϵ && return (x = x, P = P, Y = herm(x))
        g̃ = (1 / γ) * g

        Pg̃ = P * g̃
        x = x - stepsize * Pg̃

        P = (n^2 / (n^2 - 1)) * (P - (2 * stepsize) * Pg̃ * transpose(Pg̃))
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
        f.verbose && @show Y - Y'
        error("`Y` not Hermitian")
    end
    state, V = find_violation(f, Y)
    if f.verbose
        @info state
    end
    if V != zero(V)
        elt = normal_cone_element(PSD(), V)
        # if f.verbose
        #     @info "Found violation" V
        #     @info "" eigvals(V)
        #     @info "" elt
        #     @info "" eigvals(elt)
        # end
        return -invherm(elt), Inf, false
    else
        return invherm(-I(dB)), -tr(Y), true
    end

end


# We use the following basis for the Hermitian matrices
# of dimension `d`, ℍ_d:
# { E_jj for j = 1:d }
#   ∪ { (E_{jk} + E_{kj})/sqrt(2) for 1 ≤ j < k ≤ d }
#   ∪ { (E_{jk} - E_{kj})/sqrt(2) for 1 ≤ j < k ≤ d }
# where `E_{jk}` is the matrix with a 1 in the (j,k) position
# and zero elsewhere.
# This basis has the advantage that the induced map from ℝ^{d^2} to ℍ_d
# is an isometry w/r/t/ the 2-norm on ℝ^{d^2} and the 2-norm (Frobenius norm) on ℍ_d.
# This implies in particular that if ||Y* - Y_0||_∞ <= x
# then `||x* - x_0||_2 = ||Y^* - Y_0||_2 <= d^2 ||Y* - Y_0||_∞ <= d^2 x`
# By choosing `Y_0 = 0`, using `0 <= Y^* <= c[end]*I`, we find `x=c`.
# Thus, we may take P_0 = d^2 c * I.
function invherm(M::AbstractMatrix)
    d = size(M, 1)
    @assert size(M) == (d, d)
    v = zeros(real(eltype(M)), d^2)
    sqrt2 = sqrt(eltype(v)(2))
    c = 0
    for i = 1:d
        c += 1
        v[c] = real(M[i, i])
    end
    for k = 1:d
        for j = 1:k-1
            c += 1
            v[c] = sqrt2 * real(M[j, k])
            c += 1
            v[c] = sqrt2 * imag(M[j, k])
        end
    end
    return v
end

function herm(v::AbstractVector)
    d = isqrt(length(v))
    @assert d^2 == length(v)
    M = zeros(complex(eltype(v)), d, d)
    sqrt2 = sqrt(eltype(M)(2))
    c = 0
    for i = 1:d
        c += 1
        M[i, i] = v[c]
    end
    for k = 1:d
        for j = 1:k-1
            c += 1
            M[j, k] = v[c] / sqrt2
            c += 1
            M[j, k] += im * v[c] / sqrt2
        end
    end
    return Hermitian(M)
end
