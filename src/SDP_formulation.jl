"""
    N(g, x) -> Int

Given a sequence of guesses `g`, returns the number of guesses needed to
correctly guess `x`. Returns `length(g)+1` if the correct answer is never
guesses.
"""
function N(g, x)
    i = findfirst(==(x), g)
    if i === nothing
        return length(g) + 1
    else
        return i
    end
end

"""
    cartesian_power(J, K) -> iterator

Returns an interator over the set of tuples of length `K` with elements in `1:J`.

# Example

```jldoctest
julia> collect(GuessworkQuantumSideInfo.cartesian_power(3,2))
9-element Array{Tuple{Int64,Int64},1}:
 (1, 1)
 (2, 1)
 (3, 1)
 (1, 2)
 (2, 2)
 (3, 2)
 (1, 3)
 (2, 3)
 (3, 3)
```
"""
function cartesian_power(J, K)
    sz = ntuple(i -> 1:J, K)
    inds = CartesianIndices(sz)
    return (Tuple(inds[y]) for y = 1:J^K)
end

"""
    make_povm_outcomes(J, K, remove_repetition)

Returns an iterator over tuples or vectors each representing a POVM outcome. If
`remove_repetition` is true, it returns ``J!/(J-K)!`` outcomes, vectors of
length `K` with elements in ``[1,...,J]`` such that each vector has no
repetition of elements. If `remove_repetition` is false, it returns all ``J^K``
elements of ``[1, ..., J]^K``.
"""
function make_povm_outcomes(J, K, remove_repetition::Bool)
    if remove_repetition
        return multiset_permutations(1:J, K)
    else
        return cartesian_power(J, K)
    end
end

"""
    compatible_problem(f, objective, constraints, numeric_type) -> Convex.Problem

Returns a `Convex.Problem` object corresponding to `f`, which should be either
`maximize` or `minimize`, the objective and constraints given as arguments, and
the `numeric_type`. Since only Float64 numeric types are supported on Convex
versions 0.12 and below, the argument `numeric_type` is only passed to `f` when
it is different from `Float64`.
"""
function compatible_problem(f, objective, constraints, numeric_type)
    if numeric_type != Float64
        f(objective, constraints; numeric_type = numeric_type)
    else
        f(objective, constraints)
    end
end

"""
    guesswork(
        p::AbstractVector{T},
        ρBs::AbstractVector{<:AbstractMatrix};
        solver,
        K::Integer = length(p),
        c = T[1:K..., 5_000],
        dual::Bool = false,
        remove_repetition::Bool = true,
        povm_outcomes = make_povm_outcomes(length(p), K, remove_repetition),
        verbose::Bool = true,
    )

Computes the guesswork for the c-q state specified by a probability vector `p`,
giving the distribution `X`, and `ρBs`, giving the associated quantum states.

The keyword arguments are as follows:

* `solver` is the only required keyword argument; an SDP solver such as SCS or
  MOSEK must be passed.
* `K` corresponds to the maximum number of allowed guesses. The number of
  variables in the primal SDP (and the number of constraints in the dual SDP)
  scales as `length(p)^K`.
* `c` may be given a custom cost vector. If `K < length(p)`, then `c` should be
  of length `K+1`. The last entry, `c[K+1]`, corresponds to the cost of not
  guessing the correct answer within `K` guesses.
* `dual` is a boolean variable indicating whether the primal or dual
  optimization problem should be solved.
* `remove_repetition` is a boolean variable defaulting to true, indicating
  whether repeated guesses of the same value should be removed; as long as `c`
  is increasing, this decreases the size of the SDP without affecting the
  optimal value.
* `povm_outcomes` should be an iterator (or vector) corresponding to the
  possible guessing orders. This defaults to all subsets of length `K` of
  `1:length(p)` without repetition.
* `verbose` is a boolean which indicates if warnings should be printed when the
  problem is not solved optimally.
"""
function guesswork(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    solver,
    K::Integer = length(p),
    c = T[1:K..., 5_000],
    dual::Bool = false,
    remove_repetition::Bool = true,
    povm_outcomes = make_povm_outcomes(length(p), K, remove_repetition),
    verbose::Bool = true,
    debug::Bool = false,
) where {T<:Number}
    length(p) == length(ρBs) ||
    throw(ArgumentError("Length of prior and vector of side information must match J"))

    J = length(p)
    dB = size(ρBs[1], 1)
    all(ρB -> size(ρB) == (dB, dB), ρBs) ||
    throw(ArgumentError("All side-information states must be square matrices of the same size"))

    R(g) = sum(c[N(g, x)] * p[x] * ρBs[x] for x in eachindex(p, ρBs))

    num_outcomes = length(povm_outcomes)

    if dual
        Y = ComplexVariable(dB, dB)

        constraints = Convex.Constraint[R(outcome) ⪰ Y for outcome in povm_outcomes]
        push!(constraints, Y' == Y)

        objective = real(tr(Y))

        problem = compatible_problem(Convex.maximize, objective, constraints, T)
        solve!(problem, solver; verbose = verbose)

        # introduced for debugging a problem with SCS with MathOptInterface
        if debug
            @show length(constraints)
            for outcome in povm_outcomes
                @show outcome
                @show eigmin(Hermitian(R(outcome) - evaluate(Y)))
            end
        end

        Y = evaluate(Y)

        # SDP dual values not implemented in Convex.jl yet
        # Es = [Matrix{ComplexF64}(problem.constraints[y].dual) for y = eachindex(povm_outcomes)]
        Es = missing

    else

        Es = [ComplexVariable(dB, dB) for _ = 1:num_outcomes]

        objective =
            real(sum([tr(Es[y] * R(outcome)) for (y, outcome) in enumerate(povm_outcomes)]))
        constraints = Convex.Constraint[E ⪰ 0 for E in Es]
        push!(constraints, sum(Es) == I(dB))

        problem = compatible_problem(Convex.minimize, objective, constraints, T)
        solve!(problem, solver; verbose = verbose)

        Es = evaluate.(Es)
        # SDP dual values not implemented in Convex.jl yet
        # Y = Matrix{ComplexF64}(problem.constraints[end].dual)
        Y = missing
    end

    input_data = (J = J, K = K, c = c, p = p, ρBs = ρBs, dB = dB)

    return (
        optval = problem.optval,
        status = problem.status,
        Y = Y,
        Es = Es,
        povm_outcomes = povm_outcomes,
        input_data...,
    )
end
