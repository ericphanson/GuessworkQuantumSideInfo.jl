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

    noise = randn(dB^2)
    noise = init_noise * noise ./ norm(noise)

    # Hence `||Y-Y_0||_2 <= dB^2 * ||Y - Y_0||_∞ = dB^2*Y_bound`
    # Hence `||Y-Y_0 - noise||_2 <= init_noise + ||Y-Y_0||_2 <= init_noise +  dB^2*Y_bound`
    P_0 = (init_noise + Y_bound * dB^2) * Matrix{T}(I, dB^2, dB^2)

    return (x = x_0 + noise, P = P_0)
end

function guesswork_ellipsoid(p, ρBs; kwargs...)
    prob = make_ellipsoid_problem(p, ρBs; kwargs...)
    results = ellipsoid_algorithm!(prob)
    (optval = tr(results.Y), p = p, ρBs = ρBs, c = prob.c, J = length(p), K = length(p), prob=prob, results...)
end


Base.@kwdef struct EllipsoidProblem{T1,T,TρBs,Tc,Tm,TT,L, LS, TTimer, NLS}
    x::Vector{T1}
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
    debug::Bool
    trace::Bool
    tol::TT
    logger::L
    linear_solver::LS
    timer::TTimer
    nl_solver::NLS
end

function make_ellipsoid_problem(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    linear_solver,
    nl_solver,
    c::AbstractVector = T.(1:length(p)),
    max_retries = 50,
    max_time = Inf,
    num_constraints = Inf,
    verbose::Bool = false,
    num_steps_per_SA_run::Integer = length(p)^2 * 500,
    mutate! = rand_rev!,
    debug = false,
    trace = true,
    x = nothing,
    P = nothing,
    tol = 1e-3,
    init_noise = 1e-4,
    logger = ConsoleLogger(),
    timer = TimerOutput(),
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
        logger = logger,
        linear_solver=linear_solver,
        timer=timer,
        nl_solver = nl_solver,
    )
end


@enum FeasibilityStates FEASIBLE INFEASIBLE UNKNOWN

function simplex_relaxation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, ρBs, logger = prob
    J = length(p)
    c_tot = sum(c)
    c_max = c_tot - (J - 1) * c[1]
    for x = 1:J
        extremal_pt = p[x] * ρBs[x] * c_max + sum(p[y] * ρBs[y] * c[1] for y = 1:J if y ≠ x)
        e = eigmin(Hermitian(extremal_pt - Y))
        if e < 0
            # simplex relaxation is not feasible; discrete problem may or may not be
            with_logger(logger) do
                state = UNKNOWN
                @info "simplex_relaxation" state
            end
            return UNKNOWN
        end
    end
    with_logger(logger) do
        state = FEASIBLE
        @info "simplex_relaxation" state
    end
    # simplex relaxation is feasible, so discrete problem is too
    return FEASIBLE
end

include("permutations.jl")

function true_feasiblity(prob, Y)
    # only real case for now
    @unpack p, ρBs, c, dB, logger, linear_solver, timer, nl_solver = prob
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
    status = obj_value >= 0 ? FEASIBLE : INFEASIBLE
    with_logger(logger) do
        @info "true_feasiblity" obj_value status
    end
    if objective_value(m) >= 0
        return FEASIBLE, zero(Y)
    else

        PI = PermutationIterator(D=value.(D), solver = linear_solver)
        perm_iter = 0
        for P in PI
            perm_iter += 1
            R = sum( (P*c)[x] * p[x] * ρBs[x] for x in eachindex(p, ρBs))
            RY = Hermitian(R - Y)
            eigmin_RY = eigmin(RY)
            with_logger(logger) do
                @info "true_feasiblity" perm_iter eigmin_RY
            end
            if eigmin_RY < 0
                return INFEASIBLE, RY
            end
        end
        error("Could not find cut yet somehow infeasible")
    end
end    

function SA_find_cut(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, dB, ρBs, max_retries, num_steps_per_SA_run, mutate!, logger, timer = prob
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
    @timeit timer "SA" begin
        while fval > 0 && n_tries < max_retries
            n_tries += 1
            # initialize `π` to hold a random permutation, then use the simulated annealing algorithm to try to minimize `f`, starting from this permutation. Repeat until we find a violated constraint, or we've tried enough times.
            randperm!(π)
            fval = SA!(π, scratch_π, f, mutate!, num_steps_per_SA_run)
        end
    end

    with_logger(logger) do
        @info "SA_find_cut" n_tries fval
    end

    if fval < 0
        # found a cut, definitely infeasible
        return INFEASIBLE, R(π) - Y
    else
        # Could not find a cut, inconclusive
        return UNKNOWN, zero(Y)
    end
end


function find_violation(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, ρBs, logger, timer = prob
    J = length(p)

    # Easy check for feasibility:
    # We want to know if `Y <= R(π)` for all
    # permutations `π`. But since `R(π) >= ρB*c[1]` for any `π`,
    # as a quick check we can see if `Y <= ρB*c[1]`.
    ρB = sum(p[x] * ρBs[x] for x in eachindex(p))
    @timeit timer "easycheck" begin
        eigmin_easycheck = eigmin(Hermitian(c[1] * ρB - Y))
    end
    with_logger(logger) do
        @info "easycheck" eigmin_easycheck
    end
    if eigmin_easycheck >= 0
        return FEASIBLE, zero(Y)
    end

    # Strictly stronger check, but takes more time to compute (`J` eigmin's rather than 1)
    @timeit timer "simplex_relaxation" begin
        state = simplex_relaxation(prob, Y)
    end
    if state == FEASIBLE
        return FEASIBLE, zero(Y)
    end

    # Now try to find a cut via simulated annealing
    @timeit timer "SA_find_cut" begin
        state, V = SA_find_cut(prob, Y)
    end
    if state == INFEASIBLE
        return state, V
    end
    
    # last resort: globally solve a nonlinear problem
    @timeit timer "true_feasiblity" begin
        result = true_feasiblity(prob, Y)
    end

    return result
end



using LinearAlgebra
const tol = Ref(1e-4)

using SpecialFunctions
function unit_ball_volume(n)
    π^(n/2) / gamma(n/2 + 1)
end

using Dates
function ellipsoid_algorithm!(
    f::EllipsoidProblem{T};

) where {T}
    @unpack dB, x, P, logger, timer = f
    t_log = now()
    ϵ = f.tol[]
    dB = f.dB
    n = length(x)
    iter = 0
    f_best = Inf
    l_best = -Inf
    vol0 = unit_ball_volume(n)
    x_best = copy(x)
    log_interval = Millisecond(1)*1e4
    while true
        @timeit timer "Display timer" begin
            if now() - t_log >= log_interval
                t_log = now()
                @show timer
            end
        end

        @timeit timer "yield" yield()
        iter += 1
        @timeit timer "find_violation" begin
            state, V = find_violation(f, herm(x))
        end
        feasible = state == FEASIBLE
        if !feasible
            # Constraint cut
            g = -invherm(normal_cone_element(PSD(), V))
            γ = sqrt(dot(g, P, g))
            constraint_violation = -eigmin(V)
            # α = min(constraint_violation / γ, .5)
            α = 0
            with_logger(logger) do
                @info "constraint cut" γ g constraint_violation α
            end
        else
            # objective cut
            g = invherm(-I(dB))
            γ = sqrt(dot(g, P, g))
            f_val = -tr(herm(x))
            if f_best > f_val
                f_best = f_val
                x_best .= x
            end
            α = (f_val - f_best) / γ
            α = 0
            with_logger(logger) do
                @info "objective cut" γ g f_val α
            end
        end
        
        with_logger(logger) do
            vol = vol0 * sqrt(det(P))
            if !feasible
                f_val = 0
            end
            @info "Iteration" iter vol γ feasible f_val f_best
        end

        feasible && γ ≤ ϵ && return (x = x, P = P, Y = herm(x_best), x_best = x_best)

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
