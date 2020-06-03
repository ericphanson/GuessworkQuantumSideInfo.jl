@enum FeasibilityStates FEASIBLE INFEASIBLE UNKNOWN

"""
    guesswork_ellipsoid(p, ρBs; nl_solver, kwargs...) -> NamedTuple
 
Solves the guesswork SDP via a custom implementation of the ellipsoid algorithm.
When provided a global optimizer for the `nl_solver` keyword (e.g. via EAGO.jl),
this computes the value guesswork SDP. Keyword arguments are:

* `nl_solver`: choice of nonlinear optimizer to use
* `c::AbstractVector = T.(1:length(p))`: choice of cost vector
* `deepcut = true`: whether to use "deep cuts" (can be faster)
* `tol = 1e-3`: solve tolerance (final value should be within `tol` of true solution)
* `max_time::TimePeriod`: provide an approximate upper limit on the duration of the solve (defaults to unlimited)
* `verbose::Bool = true`: whether to print timing information every so often to the terminal
* `timer_log_interval::Millisecond = Millisecond(1)*1e4`: how often to print timing information
* `logger = verbose ? ConsoleLogger() : NullLogger()`: provide a logger (such as via TensorBoardLogger.jl) to log
  information during a run.
* `normal_cone_tol = 0.0`: absolute tolerance to decide if a number is negative for the purpose of computing an element in the normal cone
* `perm_tol = 1e-4`: a relative tolerance to decide if a number is `1` for the purpose of finding permutations in the support of a doubly stochastic matrix
* `init_noise = 1e-6`: an amount of noise to add to the initial point to reduce symmetries
* `trace = true`: whether or not to store a trace of the center and shape of the ellipsoid, along with other parameters, at each step
* `max_SA_retries = 2`: a parameter to decide how many times to run a simulated annealing algorithm to find cuts before resorting to global optimization
* `num_steps_per_SA_run::Integer = length(p)^2 * 500`: how many steps to perform during each simulated annealing run
* `mutate! = rand_rev!`: mutation function for updating during the simulated annealing run
* `x = nothing`: initial center point of the ellipse (automatically chosen by [`default_init`](@ref) if left unspecified)
* `P = nothing`: initial shape of the ellipse (automatically chosen by [`default_init`](@ref) if left unspecified)

To aid in resuming a previous run (after hitting a time limit or tightening `tol`), the following parameters may be passed:

* `x_best = nothing`: 
* `timer = TimerOutput()`: 
* `cuts = Vector{Int}[]`: 
* `iter = Ref(1)`: 
* `f_best = Ref(T(Inf))`: 
* `tracelog = []`: 

If `results = guesswork_ellipsoid(p, ρBs; nl_solver = ...)`, then `results.prob` holds the previous solution
object. This run may be continued by `ellipsoid_algorithm!(results.prob; tol=...)` or fields may be passed
here, e.g. `guesswork_ellipsoid(p, ρBs; x_best = results.prob.x_best, ...)`, or to the
[`EllipsoidProblem`](@ref) constructor, e.g. `EllipsoidProblem(x_best = results.prob.x_best, ...)`.
"""
function guesswork_ellipsoid(p, ρBs; kwargs...)
    prob = EllipsoidProblem(p, ρBs; kwargs...)
    return ellipsoid_algorithm!(prob)
end

"""
    default_init(p, ρBs, c, dB, init_noise) -> NamedTuple

Choose an initial center `x` and shape `P` of the ellipse
such that the solution is guaranteed to be inside the ellipse.
The initial point is chosen to correspond to `c[1]*ρB + noise` where
`ρB` is the average state, and `noise` is a small perturbation
whose norm is governed by `init_noise`.
"""
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

"""
    struct EllipsoidProblem

Stores settings and all mutable state of the ellipsoid method. This allows, e.g.
```julia
results = guesswork_ellipsoid(p, ρBs; tol=1e-3, nl_solver = ...)
# inspect `out.optval`, etc
# Continue to solve with a tighter solution tolerance:
results2 = ellipsoid_algorithm!(results.prob; tol=1e-4)
```

See [`guesswork_ellipsoid`](@ref) for the possible keyword arguments
for constructing an `EllipsoidProblem`, which may be constructed by e.g.
`EllipsoidProblem(p, ρBs; nl_solver = ..., kwargs...)`.
"""
Base.@kwdef struct EllipsoidProblem{T1,T,TρBs,Tc,Tm,TT,L, TTimer, NLS, TTrace, I, F}
    x::Vector{T1}
    x_best::Vector{T1}
    P::Matrix{T1}
    p::Vector{T}
    ρBs::TρBs
    dB::Int
    c::Tc
    verbose::Bool
    max_SA_retries::Int
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
    max_time::Base.RefValue{Millisecond}
end

function EllipsoidProblem(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    nl_solver,
    c::AbstractVector = T.(1:length(p)),
    max_SA_retries = 2,
    verbose::Bool = true,
    num_steps_per_SA_run::Integer = length(p)^2 * 500,
    mutate! = rand_rev!,
    debug = false,
    trace = true,
    deepcut = true,
    timer_log_interval::Millisecond = Millisecond(1)*1e4,
    x = nothing,
    x_best = nothing,
    P = nothing,
    tol = 1e-3,
    normal_cone_tol = 0.0,
    perm_tol = 1e-4,
    init_noise = 1e-6,
    logger = verbose ? ConsoleLogger() : NullLogger(),
    timer = TimerOutput(),
    cuts = Vector{Int}[],
    iter = Ref(1),
    f_best = Ref(T(Inf)),
    tracelog = [],
    max_time::TimePeriod = Millisecond(typemax(Int)),
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

    if !issorted(c)
        throw(ArgumentError("`c` must be increasing."))
    end
    if x === nothing || P === nothing
        @unpack x, P = default_init(p, ρBs, c, dB, init_noise)
    end
    if x_best === nothing
        x_best = copy(x)
    end

    max_time = Ref(convert(Millisecond, max_time))

    return EllipsoidProblem(
        x=x,
        x_best=x_best,
        P=P,
        p=p,
        ρBs=ρBs,
        dB=dB,
        c=c,
        max_SA_retries=max_SA_retries,
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
        tracelog = tracelog,
        deepcut = deepcut,
        timer_log_interval=timer_log_interval,
        cuts=cuts,
        iter=iter,
        f_best=f_best,
        max_time = max_time,
    )
end

function unit_ball_volume(n)
    π^(n/2) / gamma(n/2 + 1)
end

"""
    ellipsoid_algorithm!(f::EllipsoidProblem; tol=nothing, max_time=nothing)

Run the ellipsoid algorithm on the [`EllipsoidProblem`](@ref) `f`. Optionally, pass
the keyword argument `tol` or `max_time` to update these parameters of `f` before
running.
"""
function ellipsoid_algorithm!(f::EllipsoidProblem; tol::Union{Number, Nothing}=nothing, max_time::Union{TimePeriod, Nothing} = nothing)
    if tol !== nothing
        f.tol[] = tol
    end

    if max_time !== nothing
        f.max_time[] = convert(Millisecond, max_time)
    end

    with_logger(f.logger) do
        _ellipsoid_algorithm!(f)
    end
end

function _ellipsoid_algorithm!(f::EllipsoidProblem{T}) where {T}
    @unpack dB, x, P, timer, normal_cone_tol, verbose, tracelog, trace, deepcut = f
    @unpack iter, f_best, x_best, max_time, timer_log_interval = f

    tol = f.tol[]
    t_log = t_init = now()
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

        # Might help interrupt if you are pressing ctrl+c
        @timeit timer "yield" yield()

        # All the hard work happens here
        @timeit timer "separation_oracle" begin
            @unpack status, RY, π, method = separation_oracle(f, herm(x))
        end

        feasible = status == FEASIBLE

        if !feasible
            # Constraint cut
            push!(f.cuts, π)
            g = -invherm(normal_cone_element(PSD(), RY; tol = normal_cone_tol[]))
            γ = sqrt(dot(g, P, g))
            constraint_violation = -eigmin(RY)
            α = deepcut ? min(constraint_violation / γ, .9) : 0
        else
            # objective cut
            g = invherm(-I(dB)) # derivative of -tr(Y)
            γ = sqrt(dot(g, P, g))
            f_val = -tr(herm(x))
            if f_best[] > f_val
                f_best[] = f_val
                x_best .= x
            end
            α = deepcut ? (f_val - f_best[]) / γ : 0
        end
        
        # How big is our ellipsoid?
        vol = vol0 * sqrt(det(Hermitian(P)))

        if !feasible
            f_val = NaN
        end

        @info "Iteration" iter[] vol γ feasible f_val f_best[] method α
        @timeit timer "Trace logging" begin
            maxradii = eigmax(Hermitian(P))
            trace && push!(tracelog, (iter=iter[], method=method, status=status, γ=γ, f_best=f_best[], vol=vol, f_val = f_val, maxradii=maxradii, α=α, π=π, x = copy(x), P = copy(P)))
        end

        # Exit condition
        if feasible && γ ≤ tol
            Y = herm(x)
            return (status = MOI.OPTIMAL, optval = tr(Y), Y=Y,
                    p = f.p, ρBs = f.ρBs, c = f.c,
                    J = length(f.p), K = length(f.p), prob=f)
        end

        # Early stopping if max time exceeded
        if now() - t_init >= max_time[]
            Y = herm(x)
            if verbose
                println("Warning: hit maximum time, stopping early.")
            end
            return (status =  MOI.TIME_LIMIT, optval = tr(Y), Y=Y,
            p = f.p, ρBs = f.ρBs, c = f.c,
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

"""
    invherm(M::AbstractMatrix) -> Vector

Creates a real vector `v` representation of a
complex Hermitian matrix `M`. These are related
by a linear isometry (with respect to the 2-norm).
See also [`herm`](@ref).

## Example

```julia
M = Hermitian(rand(4,4) + im*rand(4,4))
v = invherm(M)

norm(v, 2) ≈ norm(M, 2)
```
"""
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

"""
    herm(M::AbstractVector) -> Hermitian

Creates the complex Hermitian matrix `M` represented
by the real vector `v`. 
See also [`invherm`](@ref).

## Example

```julia
v = rand(16)
M = herm(M)

norm(v, 2) ≈ norm(M, 2)
```
"""
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
