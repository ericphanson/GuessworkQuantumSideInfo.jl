"""
    guesswork_upper_bound(
        p::AbstractVector{T},
        ρBs::AbstractVector{<:AbstractMatrix};
        make_solver,
        K::Integer = length(p),
        c::AbstractVector = T[1:K..., 10_000],
        max_retries = 50,
        max_time = Inf,
        num_constraints = Inf,
        verbose::Bool = false,
        num_steps_per_SA_run::Integer = length(p)^2 * 500,
    ) where {T<:Number} -> NamedTuple

Computes an upper bound to the guesswork problem associated to the c-q state
specified by `p` and `ρBs`, as in [`guesswork`](@ref). The number of allowed
guesses `K`, and a custom cost vector `c` may be optionally passed. If the
keyword argument `verbose` is set to true, information is printed about each
iteration of the algorithm.

The keyword argument `make_solver` is required, and must pass a *function that
creates a solver instances*. For example, instead of passing `SCSSolver()`, pass
`() -> SCSSolver()`. This is needed because the algorithm used in
`guesswork_upper_bound` solves a sequence of SDPs, not just one.

The algorithm has three termination criteria which are controlled by keyword
arguments. The algorithm stops when any of the following occur:

* `max_retries` simulated annealing attempts fail to find a violated constraint.
* `num_constraints` constraints have been added to the dual SDP
* The total runtime of the algorithm is projected to exceed `max_time` on the next iteration.

By default, `max_retries` is set to 50, while `num_constraints` and `max_time` are set to infinity.

Lastly, the keyword argument `num_steps_per_SA_run` controls the runtime of the
simulated annealing algorithm. Increase `num_steps_per_SA_run` to search longer
for a violated constraint within a given simulated annealing run.
"""
function guesswork_upper_bound(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    make_solver,
    K::Integer = length(p),
    c::AbstractVector = T[1:K..., 10_000],
    max_retries = 50,
    max_time = Inf,
    num_constraints = Inf,
    verbose::Bool = false,
    num_steps_per_SA_run::Integer = length(p)^2 * 500,
    mutate! = rand_rev!,
    debug = false,
) where {T<:Number}

    length(p) == length(ρBs) ||
    throw(ArgumentError("Length of prior and vector of side information must match J"))
    J = length(p)
    dB = size(ρBs[1], 1)
    all(ρB -> size(ρB) == (dB, dB), ρBs) ||
    throw(ArgumentError("All side-information states must be square matrices of the same size"))

    if max_retries == Inf && max_time == Inf && num_constraints == Inf
        throw(ArgumentError("All three termination criteria (`max_retries`, `max_time`, and `num_constraints`) are infinite; algorithm would never terminate."))
    end

    constraints::Vector{Vector{Int}} = []
    convex_problem, convex_Y, convex_R = make_problem(p, ρBs, c, dB, constraints, T)
    solve_time::Float64 = 0.0
    convex_Y.value = Matrix{Complex{T}}(I, dB, dB)
    upper_bound::T = 0.0

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

    # run the main loop
    # verbose is passed as a "value type" so logging statements can be compiled out when it is set to false.
    upper_bound = _loop(
        upper_bound,
        R_scratch,
        constraints,
        π,
        scratch_π,
        convex_R,
        convex_problem,
        convex_Y,
        Rg_inv!,
        num_steps_per_SA_run,
        total_time_so_far,
        max_retries,
        max_time,
        num_constraints,
        mutate!,
        make_solver,
        T,
        Val(verbose),
    )

    input_data = (J = J, K = K, c = c, p = p, ρBs = ρBs, dB = dB)

    return (
        optval = upper_bound,
        status = convex_problem.status,
        Y = evaluate(convex_Y),
        povm_outcomes = invperm.(constraints),
        input_data...,
    )
end

# We put the inner loop behind a function barrier
function _loop(
    upper_bound::T,
    R_scratch,
    constraints,
    π,
    scratch_π,
    convex_R,
    convex_problem,
    convex_Y,
    Rg_inv!,
    num_steps_per_SA_run,
    total_time_so_far,
    max_retries,
    max_time,
    num_constraints,
    mutate!,
    make_solver,
    ::Type{T},
    ::Val{verbose},
) where {T,verbose}

    while length(constraints) < num_constraints
        f = let Y::Matrix{Complex{T}} = evaluate(convex_Y)
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
        total_time_so_far += SA_time

        if fval < 0
            # add the permutation to the list of constraints and resolve the problem.
            new_constraint = copy(π)

            push!(constraints, new_constraint)
            update_problem!(convex_problem, convex_Y, new_constraint, convex_R)

            solve_time = @elapsed solve!(convex_problem, make_solver(); verbose = verbose)
            upper_bound::T = convex_problem.optval

            if verbose
                _log(
                    length(constraints),
                    fval,
                    upper_bound,
                    SA_time,
                    solve_time,
                    total_time_so_far,
                )
            end
        else
            # We've reached `max_retries` attempts to find a violated constraint.
            verbose && _log_retries(max_retries)
            break
        end

        # if we do another loop, it will take a round of SA plus a solve, so let's stop if
        # that will exceed the max time.
        projected_time = total_time_so_far + 1.1*(SA_time + solve_time)
        if projected_time > max_time
            if verbose
                _log_time(total_time_so_far, total_time_so_far + SA_time + solve_time, max_time)
            end
            break
        end

    end
    return upper_bound
end

@noinline function _log(
    lconstraints,
    violation,
    upper_bound,
    SA_time,
    solve_time,
    total_time_so_far,
)
    @info "Adding constraint $(lconstraints) " violation upper_bound SA_time solve_time total_time_so_far
end


@noinline function _log_retries(max_retries)
    @info "Unable to find violated constraint after repeating ($max_retries) simulated annealing runs. Change `max_retries` to adjust this behavior."
end

@noinline function _log_time(total_time_so_far, projected_time, max_time)
    if total_time_so_far < max_time
        @info "Next iteration is projected to exceed maximum time $(max_time); terminating here. Change `max_time` to adjust this behavior." total_time_so_far projected_time
    else
        @info "Reached maximum time $(max_time); terminating here. Change `max_time` to adjust this behavior." total_time_so_far
    end
end

function make_problem(p, ρBs, c, dB, constraints, ::Type{T}) where {T}
    R = g -> sum(c[N(invperm(g), x)] * p[x] * ρBs[x] for x in eachindex(p, ρBs))
    Y = ComplexVariable(dB, dB)
    povm_outcomes = constraints
    Y = ComplexVariable(dB, dB)
    constraints = Convex.Constraint[R(outcome) ⪰ Y for outcome in povm_outcomes]
    push!(constraints, Y' == Y)
    objective = real(tr(Y))
    return compatible_problem(Convex.maximize, objective, constraints, T), Y, R
end

function update_problem!(problem, Y, new_constraint, R)
    push!(problem.constraints, R(new_constraint) ⪰ Y)
end

function get_dual_and_val(
    problem,
    Y,
    make_solver,
    ::Val{dB},
    ::Val{verbose},
    ::Type{T},
) where {dB,T,verbose}
    solve!(problem, make_solver(); verbose = verbose)
    return problem.optval::T, SMatrix{dB,dB,Complex{T}}(evaluate(Y))
end


"""
SA!(π, scratch_π, f, mutate! = rand_rev!)

Try to minimize `f` over permutations of 1:`length(π)` via a simulated annealing algorithm, returning the value of `f` on the final choice of `π`. Modifies `π` to hold the optimal permutation at the end of the search, and uses `scratch_π` to hold proposal permutations. The optional argument `mutate!` should be a function

    mutate!(π_to_mutate, scratch_π, f, J) -> nothing

which modifies `π_to_mutate` in place, can use `scratch_π` to hold temporary permutations, and has access to `f` and `J`, the length of `π_to_mutate`. By default, `rand_rev!` is chosen, which is very fast and requires no evaluations of `f`. Other choices of `mutate!` include  `swap_all!` and `two_opt!` (uses ~ `n^2` evaluations of `f`), or a random choice at each step, via e.g.

    mutate! = (π_to_mutate, scratch_π, f, J) -> begin
        if rand() < .01
            two_opt!(π_to_mutate, scratch_π, f, J)
        else
            rand_rev!(π_to_mutate, scratch_π, f, J)
        end
    end
"""
function SA!(π, scratch_π, f, mutate!, num_steps)
    J = length(π)
    init_temp = exp(8)
    final_temp = exp(-6.5)
    cool_rate = (final_temp / init_temp)^(1 / (num_steps - 1))
    # divide by cool_rate so when we first multiply we get init_temp
    temp = init_temp / cool_rate
    proposed_π = scratch_π

    current_cost = f(π)
    for _ = 1:num_steps
        temp *= cool_rate
        mutate!(proposed_π, π, f, J)
        proposed_cost = f(proposed_π)
        @fastmath accept = proposed_cost < current_cost ? true :
            rand() < exp((current_cost - proposed_cost) / temp)
        if accept
            current_cost = proposed_cost
            π .= proposed_π
        else
            proposed_π .= π
        end
    end
    # try all possible swaps at the end just to try to really have a optimal value.
    fval = two_opt!(π, scratch_π, f, J)
    return fval
end


"""
    rand_rev!(π, scratch_π = nothing, f = nothing, n = length(π))

Modifies `π` by randomly reversing a section of the permutation. Has the arguments `scratch_π` and `f` so it can be used in `SA!`.
"""
@inline function rand_rev!(π, scratch_π = nothing, f = nothing, n = length(π))
    i, j = rand(1:n), rand(1:n)
    if i > j
        i, j = j, i
    end
    reverse!(π, i, j)
end

"""
    swap_all!(π, scratch_π, f, n = length(π))

Modifies `π` in place to find a lower value of `f` by swapping each of the `n*(n-1)/2` pairs of entries and choosing the optimal one.
"""
function swap_all!(π, scratch_π, f, n = length(π))
    @inline swap!(π, i, j) = π[i], π[j] = π[j], π[i]
    _mutate_two!(π, scratch_π, f, swap!, n)
end

"""
    two_opt!(π, scratch_π, f, n = length(π))

Modifies `π` in place to find a lower value of `f` by reversing each possible ordered section of  `π`  by trying all `n*(n-1)/2` possible pairs of endpoints and choosing the optimal one.
"""
function two_opt!(π, scratch_π, f, n = length(π))
    _mutate_two!(π, scratch_π, f, reverse!, n)
end

@inline function _mutate_two!(π, scratch_π, f, mut!, n)
    best_fval_so_far = Inf
    best_π_so_far = scratch_π
    for (i, j) in combinations(1:n, 2)
        mut!(π, i, j)
        fval = f(π)
        if fval < best_fval_so_far
            best_fval_so_far = fval
            best_π_so_far .= π
        end
        mut!(π, i, j)
    end
    π .= best_π_so_far
    return best_fval_so_far
end
