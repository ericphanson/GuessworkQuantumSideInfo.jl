# """

# Note `num_outcomes` is not `K`. (clarify...) 
# """
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
