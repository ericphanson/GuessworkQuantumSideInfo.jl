"""
    guesswork_MISDP(
        p::AbstractVector{T},
        ρBs::AbstractVector{<:AbstractMatrix},
        num_outcomes;
        solver,
        c = T.(1:length(p)),
        verbose::Bool = true,
    ) where {T<:Number} -> NamedTuple

Computes an approximation to the guesswork with `num_outcomes` possible guessing
orders for the c-q state specified by a probability vector `p`, giving the
distribution `X`, and `ρBs`, giving the associated quantum states. If
``num_outcomes > d_B^2``, where `d_B` is the dimension of the Hilbert space on
which the quantum states live, then this computes exactly the guesswork. Note
that one may not supply a parameter `K` in this case; in the mixed-integer SDP
formulation, the guesser is always allowed to make `length(p)` guesses. However,
a custom cost vector `c` may be supplied to choose how to penalize each possible
number of guesses required to get the correct answer.

The size of the mixed-integer SDP solved by this function grows polynomially in
`length(p)`, the dimension `d_B`, and `num_outcomes`, but mixed-integer SDPs are
not known to be efficiently computable in general.

The keyword argument `solver` must be supplied with a solver capable of solving
mixed-integer SDPs. Currently, this means
[Pajarito.jl](https://github.com/JuliaOpt/Pajarito.jl) which solves
mixed-integer SDPs by solving an alternating sequence of mixed-integer linear
programs and SDPs. See the documentation for an example. Note that the
performance of Pajarito depends on the performance of the underlying
mixed-integer linear solver and SDP solver it is given, and that commerical
solvers often have academic licenses and can be much more performant.

The keyword argument `verbose` prints the status, optimal value, optimal
guessing orders, and POVMs, in addition to warnings from Convex.jl if the
problem is not optimally solved.
"""
function guesswork_MISDP(
    p::AbstractVector{T},
    ρBs::AbstractVector{<:AbstractMatrix},
    num_outcomes;
    solver,
    c = T.(1:length(p)),
    verbose::Bool = true,
) where {T<:Number}
    length(p) == length(ρBs) ||
    throw(ArgumentError("Length of prior and vector of side information must match J"))

    J = length(p)
    dB = size(ρBs[1], 1)
    all(ρB -> size(ρB) == (dB, dB), ρBs) ||
    throw(ArgumentError("All side-information states must be square matrices of the same size"))

    constraints = Convex.Constraint[]

    Es = [Convex.ComplexVariable(dB, dB) for _ = 1:num_outcomes]
    for E in Es
        push!(constraints, E ⪰ 0)
    end
    πs = [reshape(Convex.Variable(J^2, :Bin), J, J) for _ = 1:num_outcomes]
    for j = 1:num_outcomes, i = 1:J
        push!(constraints, sum(πs[j][i, :]) == 1)
        push!(constraints, sum(πs[j][:, i]) == 1)
    end
    push!(constraints, sum(Es) == I(dB))

    Λs = []

    for k = 1:num_outcomes
        Λre = [Convex.Variable(J, J) for a = 1:dB, b = 1:dB]
        Λim = [Convex.Variable(J, J) for a = 1:dB, b = 1:dB]
        push!(Λs, Λre + im * Λim)
        π = πs[k]
        E = Es[k]
        E_U = dB / 2
        E_L = -1 * dB / 2
        for a = 1:dB, b = 1:dB
            Λab_re = Λre[a, b]
            Λab_im = Λim[a, b]
            push!(constraints, Λab_re >= π * E_L)
            push!(constraints, Λab_re <= π * E_U)
            push!(constraints, Λab_re >= real(E[a, b]) + π * E_U - E_U)
            push!(constraints, Λab_re <= real(E[a, b]) + π * E_L - E_L)

            push!(constraints, Λab_im >= π * E_L)
            push!(constraints, Λab_im <= π * E_U)
            push!(constraints, Λab_im >= imag(E[a, b]) + π * E_U - E_U)
            push!(constraints, Λab_im <= imag(E[a, b]) + π * E_L - E_L)
        end
    end

    ρBs_subnormalized = p .* ρBs
    vs = [
        [
            real(sum(
                ρBs_subnormalized[i][a, b] * Λ[b, a][j, i] for a = 1:dB for b = 1:dB
                for i = 1:J
            ))
            for j = 1:J
        ]
        for Λ in Λs
    ]

    objective = sum(v[j] * c[j] for j = 1:J for v in vs)

    for v in vs, j = 1:J-1
        push!(constraints, v[j] >= v[j+1])
    end

    prob = compatible_problem(Convex.minimize, objective, constraints, T)

    if verbose
        @info "Starting MISDP solve"
    end
    Convex.solve!(prob, solver; verbose = verbose)
    optval = prob.optval
    status = prob.status

    πs_perm = [Tuple(round.(Int, (evaluate(π) * collect(1:J)))) for π in πs]

    Es = evaluate.(Es)

    if verbose
        @info "MISDP solve"
        @show status
        @show optval
        @show πs_perm
        @show Es
    end

    data = (J = J, K = J, c = c, p = p, ρBs = ρBs, dB = dB)

    return (
        optval = optval,
        Es = Es,
        povm_outcomes = πs_perm,
        num_outcomes = num_outcomes,
        data...,
    )
end
