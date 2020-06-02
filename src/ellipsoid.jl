@enum FeasibilityStates FEASIBLE INFEASIBLE UNKNOWN

using Logging, JuMP, MathOptInterface
const MOI = MathOptInterface
using TimerOutputs

function default_init(p, ρBs, c, dB, init_noise) 
    T = eltype(p)
    ρB = sum(p[x] * ρBs[x] for x in eachindex(p))

    Y_0 = c[1] * ρB
    x_0 = invherm(Y_0)

    # In semidefinite order,
    # 0 <= Y - Y_0 <= c[end]*I(dB^2) - Y_0
    Y_bound = opnorm(c[end]*I(dB) - Y_0)

    noise = rand(T, dB^2) - rand(T, dB^2)
    noise = init_noise * noise ./ norm(noise)

    # Hence `||Y-Y_0||_2 <= dB^2 * ||Y - Y_0||_∞ = dB^2*Y_bound`
    # Hence `||Y-Y_0 - noise||_2 <= init_noise + ||Y-Y_0||_2 <= init_noise +  dB^2*Y_bound`
    P_0 = (init_noise + Y_bound * dB^2) * Matrix{T}(I, dB^2, dB^2)

    return (x = x_0 + noise, P = P_0)
end

function guesswork_ellipsoid(p, ρBs; kwargs...)
    prob = EllipsoidProblem(p, ρBs; kwargs...)
    return ellipsoid_algorithm!(prob)
end

"""
    struct EllipsoidProblem

Stores settings and all mutable state of the ellipsoid method. This allows, e.g.
```julia
results = guesswork_ellipsoid(p, ρBs; tol=1e-3, nl_solver = ...)
# inspect `out.optval`, etc
# Continue to solve with a tighter solution tolerance:
results2 = ellipsoid_algorithm!(results.prob; tol=1e-4)
```
"""
Base.@kwdef struct EllipsoidProblem{T1,T,TρBs,Tc,Tm,TT,L, TTimer, NLS, TTrace, I, F}
    x::Vector{T1}
    x_best::Vector{T1}
    P::Matrix{T1}
    p::Vector{T}
    ρBs::TρBs
    dB::Int
    c::Tc
    max_retries::Int
    max_time::Float64
    verbose::Bool
    num_steps_per_SA_run::Integer
    mutate!::Tm
    deepcut::Bool
    debug::Bool
    trace::Bool
    tracelog::TTrace
    tol::TT
    normal_cone_tol::TT
    perm_tol::TT
    logger::L
    timer::TTimer
    nl_solver::NLS
    timer_log_interval::Millisecond
    cuts::Vector{Vector{Int}}
    iter::I
    f_best::F
end

function EllipsoidProblem(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    nl_solver,
    c::AbstractVector = T.(1:length(p)),
    max_retries = 2,
    max_time = Inf,
    num_constraints = Inf,
    verbose::Bool = true,
    num_steps_per_SA_run::Integer = length(p)^2 * 500,
    mutate! = rand_rev!,
    debug = false,
    trace = true,
    deepcut = true,
    timer_log_interval::Millisecond = Millisecond(1)*1e4,
    x = nothing,
    P = nothing,
    tol = 1e-3,
    normal_cone_tol = 1e-4,
    perm_tol = 1e-4,
    init_noise = 1e-6,
    logger = verbose ? ConsoleLogger() : NullLogger(),
    timer = TimerOutput(),
    cuts = Vector{Int}[],
    iter = Ref(1),
    f_best = Ref(T(Inf))
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

    if x === nothing || P === nothing
        @unpack x, P = default_init(p, ρBs, c, dB, init_noise) 
    end

    return EllipsoidProblem(
        x=x,
        x_best=copy(x),
        P=P,
        p=p,
        ρBs=ρBs,
        dB=dB,
        c=c,
        max_retries=max_retries,
        max_time=max_time,
        verbose=verbose,
        num_steps_per_SA_run=num_steps_per_SA_run,
        mutate! = mutate!,
        debug = debug,
        trace = trace,
        tol = Ref(tol),
        normal_cone_tol = Ref(normal_cone_tol),
        perm_tol = Ref(perm_tol),
        logger = logger,
        timer=timer,
        nl_solver = nl_solver,
        tracelog = [],
        deepcut = deepcut,
        timer_log_interval=timer_log_interval,
        cuts=cuts,
        iter=iter,
        f_best=f_best,
    )
end

function simplex_relaxation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, ρBs = prob
    J = length(p)
    c_tot = sum(c)
    c_max = c_tot - (J - 1) * c[1]
    for x = 1:J
        extremal_pt = p[x] * ρBs[x] * c_max + sum(p[y] * ρBs[y] * c[1] for y = 1:J if y ≠ x)
        e = eigmin(Hermitian(extremal_pt - Y))
        if e < 0
            # simplex relaxation is not feasible; discrete problem may or may not be
            return UNKNOWN
        end
    end
    # simplex relaxation is feasible, so discrete problem is too
    return FEASIBLE
end

include("permutations.jl")

function true_feasiblity(prob, Y)
    # only real case for now
    @unpack p, ρBs, c, dB, timer, nl_solver, perm_tol, tracelog, trace = prob
    @timeit timer "formulation" begin
        J = length(p)
        A_re = real.(ρBs) .* p
        A_im = imag.(ρBs) .* p
        Y_re = real.(Y)
        Y_im = imag.(Y)
        m = Model(nl_solver)
        @variable(m, 0 <= D[i=1:J,j=1:J] <= 1)
        @constraint(m, sum(D, dims = 1) .== 1)
        @constraint(m, sum(D, dims = 2) .== 1)
        @variable(m, -1 <= ψ_re[i=1:dB] <= 1)
        @variable(m, -1 <= ψ_im[i=1:dB] <= 1)

        @NLexpression(m, R_re[l=1:dB,k=1:dB], sum(c[j]*A_re[i][l,k]*D[i,j] for j = 1:J, i = 1:J))
        @NLexpression(m, R_im[l=1:dB,k=1:dB], sum(c[j]*A_im[i][l,k]*D[i,j] for j = 1:J, i = 1:J))
        @NLexpression(m, t1, sum(ψ_re[l] * R_re[l,k] * ψ_re[k] - ψ_re[l] * R_im[l,k] * ψ_im[k] + ψ_im[l]*R_im[l,k] *ψ_re[k] + ψ_im[l] * R_re[l,k]*ψ_im[k] for l = 1:dB, k = 1:dB))
        @NLexpression(m, t2, sum(ψ_re[l] * Y_re[l,k] * ψ_re[k] - ψ_re[l] * Y_im[l,k] * ψ_im[k] + ψ_im[l]*Y_im[l,k] *ψ_re[k] + ψ_im[l] * Y_re[l,k]*ψ_im[k] for l = 1:dB, k = 1:dB))
        @NLobjective(m, Min, t1 - t2)
    end
    @timeit timer "optimize!" JuMP.optimize!(m)
    @assert termination_status(m) == MOI.OPTIMAL

    
    obj_value = objective_value(m)
    feasible = obj_value >= 0
    @info "true_feasiblity" obj_value feasible
    if objective_value(m) >= 0
        return (status=FEASIBLE, RY=nothing, π=nothing)
    else

        PI = PermutationIterator(D=value.(D), rtol=perm_tol[], timer=timer)
        perm_iter = 0
        for (π, α) in PI
            perm_iter += 1
            R = sum( c[π[x]] * p[x] * ρBs[x] for x in eachindex(p, ρBs))
            RY = Hermitian(R - Y)
            eigmin_RY = eigmin(RY)
            @info "true_feasiblity" perm_iter eigmin_RY
            if eigmin_RY < 0
                return (status=INFEASIBLE, RY=RY, π=π)
            end
        end
        error("This should never be reached! Could not find cut yet somehow infeasible")
    end
end    

function SA_find_cut(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, dB, ρBs, max_retries, num_steps_per_SA_run, mutate!, timer = prob
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

    f = let Y::Matrix{Complex{T}} = Y
        π -> begin
            R_scratch .= -1 .* Y
            eigmin(Hermitian(Rg_inv!(R_scratch, π)))::T
        end
    end


    # Choose `SA_violation` as something `> 0` to enter the loop.
    SA_violation::T = 1.0
    n_tries::Int = 0
    @timeit timer "SA" begin
        while SA_violation > 0 && n_tries < max_retries
            n_tries += 1
            # initialize `π` to hold a random permutation,
            # then use the simulated annealing algorithm to try to minimize `f`,
            # starting from this permutation. Repeat until we find a violated constraint,
            # or we've tried enough times.
            randperm!(π)
            SA_violation = SA!(π, scratch_π, f, mutate!, num_steps_per_SA_run)
        end
    end

    found_cut = SA_violation < 0
    @info "SA_find_cut" n_tries SA_violation found_cut

    if found_cut
        # found a cut, definitely infeasible
        return (status = INFEASIBLE, RY = Hermitian(R(π) - Y), π=π)
    else
        # Could not find a cut, inconclusive
        return (status=UNKNOWN, RY=nothing, π=nothing)
    end
end


function find_violation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, ρBs, timer = prob
    J = length(p)

    # Simple feasability check
    @timeit timer "simplex_relaxation" begin
        state = simplex_relaxation(prob, Y)
    end
    if state == FEASIBLE
        return (status=FEASIBLE, RY=nothing, π=nothing, method=:simplex_relaxation)
    end

    # Now try to find a cut via simulated annealing
    @timeit timer "SA_find_cut" begin
        @unpack status, RY, π = SA_find_cut(prob, Y)
    end
    if status == INFEASIBLE
        return (status=state, RY=RY, π=π, method=:simulated_annealing)
    end
    
    # last resort: globally solve a nonlinear problem
    @timeit timer "true_feasiblity" begin
        result = true_feasiblity(prob, Y)
    end

    return (; result..., method=:global_optimization)
end

using SpecialFunctions
function unit_ball_volume(n)
    π^(n/2) / gamma(n/2 + 1)
end

function ellipsoid_algorithm!(f::EllipsoidProblem, newtol::Union{Number, Nothing}=nothing)
    with_logger(f.logger) do
        _ellipsoid_algorithm!(f, newtol)
    end
end

function _ellipsoid_algorithm!(f::EllipsoidProblem{T}, newtol::Union{Number, Nothing}=nothing) where {T}
    @unpack dB, x, P, timer, normal_cone_tol, verbose, tracelog, trace, deepcut = f
    @unpack iter, f_best, x_best = f
    if newtol !== nothing
        f.tol[] = newtol
    end
    @unpack timer_log_interval = f
    t_log = now()
    ϵ = f.tol[]
    dB = f.dB
    n = length(x)
    vol0 = unit_ball_volume(n)
    while true
        if verbose
            @timeit timer "Display timer" begin
                if now() - t_log >= timer_log_interval
                    t_log = now()
                    @show timer
                end
            end
        end

        @timeit timer "yield" yield()
        @timeit timer "find_violation" begin
            @unpack status, RY, π, method = find_violation(f, herm(x))
        end
        feasible = status == FEASIBLE
        if !feasible
            push!(f.cuts, π)
            # Constraint cut
            g = -invherm(normal_cone_element(PSD(), RY; tol = normal_cone_tol[]))
            γ = sqrt(dot(g, P, g))
            constraint_violation = -eigmin(RY)
            α = deepcut ? min(constraint_violation / γ, .9) : 0
        else
            # objective cut
            g = invherm(-I(dB))
            γ = sqrt(dot(g, P, g))
            f_val = -tr(herm(x))
            if f_best[] > f_val
                f_best[] = f_val
                x_best .= x
            end
            α = deepcut ? (f_val - f_best[]) / γ : 0
        end
        
        vol = vol0 * sqrt(det(P))
        if !feasible
            f_val = 0
        end
        @info "Iteration" iter[] vol γ feasible f_val f_best[] method α
        @timeit timer "Trace logging" begin
            trace && push!(tracelog, (iter=iter[], method=method, status=status, γ=γ, f_best=f_best[], vol=vol, f_val = f_val, maxradii=eigmax(P), α=α, π=π, x = copy(x), P = copy(P)))
        end
        if feasible && γ ≤ ϵ
            Y = herm(x)
            return (optval = tr(Y), Y=Y, p = f.p, ρBs = f.ρBs, c = f.c,
                    J = length(f.p), K = length(f.p), prob=f)
        end

        # Update ellipsoid center
        @timeit timer "update center" begin
            g̃ = (1 / γ) * g
            Pg̃ = P * g̃
            x .= x - (1 + n*α)*inv(n + 1) * Pg̃
        end

        # Update ellipsoid shape
        @timeit timer "update shape" begin
            P .= (n^2 / (n^2 - 1)) * (1 - α^2)*(P - 2*(1 + n*α) / ((n+1)*(1+α)) * Pg̃ * transpose(Pg̃))
        end
        iter[] += 1
    end
end

struct Orthant end
struct PSD end

function normal_cone_element(::Orthant, λ::AbstractVector; tol::Number)
    z = zeros(length(λ))
    for i in eachindex(λ)
        if real(λ[i]) <= tol
            z[i] = -1
        end
    end
    return z
end

function normal_cone_element(::PSD, P::AbstractMatrix; tol::Number)
    λ, Q = eigen(P)
    Q * Diagonal(normal_cone_element(Orthant(), λ; tol=tol)) * Q'
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
