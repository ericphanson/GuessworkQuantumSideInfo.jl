
"""
    guesswork_lower_bound(
        p::AbstractVector{T},
        ρBs::AbstractVector{<:AbstractMatrix};
        solver,
        c = T[1:length(p)..., 10_000],
        verbose::Bool = false,
    )

See [`guesswork`](@ref) for the meaning of the arguments. Computes a lower bound
to the optimal expected number of guesses by solving a relaxed version of the
primal SDP. For `J` states, only needs `J^2` PSD variables subject to two linear
constraints.
"""
function guesswork_lower_bound(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix};
    solver,
    c = T[1:length(p)..., 10_000],
    verbose::Bool = false,
) where {T<:Number}
    length(p) == length(ρBs) ||
    throw(ArgumentError("Length of prior and vector of side information must match J"))
    J = length(p)
    dB = size(ρBs[1], 1)
    all(ρB -> size(ρB) == (dB, dB), ρBs) ||
    throw(ArgumentError("All side-information states must be square matrices of the same size"))

    ℰs = [ComplexVariable(dB, dB) for j = 1:J, k = 1:J]

    objective = real(sum(c[j] * tr(ℰs[j, k] * p[k] * ρBs[k]) for j = 1:J for k = 1:J))

    constraints = Convex.Constraint[ℰ ⪰ 0 for ℰ in ℰs] |> vec

    for j = 1:J
        push!(constraints, sum(ℰs[j, :]) == I(dB))
        push!(constraints, sum(ℰs[:, j]) == I(dB))
    end

    problem = compatible_problem(Convex.minimize, objective, constraints, T)

    solve!(problem, solver; verbose = verbose)

    ℰs = evaluate.(ℰs)

    data = (J = J, K = J, c = c, p = p, ρBs = ρBs, dB = dB)

    return (optval = problem.optval, status = problem.status, ℰs = ℰs, data...)
end
