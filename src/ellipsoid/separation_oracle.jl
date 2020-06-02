
function separation_oracle(prob::EllipsoidProblem{T}, Y) where {T}
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


# To check feasability, we wish to check that
# `S := [ sum( c[π[x]]*p[x]*ρBs[x] for x = 1:J ) - Y for π in permutations(1:J) ] ⊆ PSD`
# We relax this problem by considering a larger set which is simpler to construct:
# `S̃ := [ sum( c̃[i] *p[x]*ρBs[x] for x = 1:J ) - Y for c̃[i] >= 1, sum(c̃) == sum(c) ]`
# Since `S ⊆ S̃`, if `S̃ ⊆ PSD`, then `S` is too.
# Moreover, the latter is a polytope (as the composition of affine maps applied to 
# the standard simplex), so we may simply check that its extreme points lie in the PSD cone.
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

# Heuristically try to find a cut via simulated annealing
function SA_find_cut(prob::EllipsoidProblem{T}, Y) where {T}
    @unpack c, p, dB, ρBs, max_SA_retries, num_steps_per_SA_run, mutate!, timer = prob
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
        while SA_violation > 0 && n_tries < max_SA_retries
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


# We aim to solve the non-linear problem
# `min_π,ψ  <ψ, (sum( c[π[x]]*p[x]*ρBs[x] for x = 1:J ) - Y) ψ>`.
# We can relax the discrete constraints that
# π is a permutation to be instead a doubly stochastic matrix `D`
# `min_D,ψ  <ψ, (sum( (D*c)[x]*p[x]*ρBs[x] for x = 1:J ) - Y) ψ>`
# since if `sum( c[π[x]]*p[x]*ρBs[x] for x = 1:J ) >= Y` for all `π`,
# then sum( (D*c)[x]*p[x]*ρBs[x] for x = 1:J ) >= Y` for all `D`
# simply by averaging, using that `D` is the convex combination of permutations.
# If `nl_solver` is a global solver (i.e. solves the problem to global optimality),
# then we truly obtain whether or not the choice `Y` is feasible or not.
# Moreover, we may obtain `π` such that `c[π[x]]*p[x]*ρBs[x] for x = 1:J ) ≱ Y`
# by decomposing the resulting bistochastic matrix `D` into constituent permutations.
function true_feasiblity(prob, Y)
    # only real case for now
    @unpack p, ρBs, c, dB, timer, nl_solver, perm_tol, tracelog, trace = prob
    @timeit timer "formulation" begin
        J = length(p)
        # JuMP can only handle real numbers,
        # so we split into real and imaginary parts.
        A_re = real.(ρBs) .* p
        A_im = imag.(ρBs) .* p
        Y_re = real.(Y)
        Y_im = imag.(Y)

        m = Model(nl_solver)
        # `D` is doubly stochastic:
        @variable(m, 0 <= D[i=1:J,j=1:J] <= 1)
        @constraint(m, sum(D, dims = 1) .== 1)
        @constraint(m, sum(D, dims = 2) .== 1)

        # `ψ` is any vector; the magnitudes do not play a role for whether or not the objective is
        # non-negative, so we constraint them which makes the branch-and-bound solver much faster.
        @variable(m, -1 <= ψ_re[i=1:dB] <= 1)
        @variable(m, -1 <= ψ_im[i=1:dB] <= 1)

        # R = (D*c)[i]*p[i]*ρBs[i]
        @NLexpression(m, R_re[l=1:dB,k=1:dB], sum(c[j]*A_re[i][l,k]*D[i,j] for j = 1:J, i = 1:J))
        @NLexpression(m, R_im[l=1:dB,k=1:dB], sum(c[j]*A_im[i][l,k]*D[i,j] for j = 1:J, i = 1:J))

        # t1 = <ψ, R ψ>, exploiting that this quantity is real
        @NLexpression(m, t1, sum(ψ_re[l] * R_re[l,k] * ψ_re[k] - ψ_re[l] * R_im[l,k] * ψ_im[k] + ψ_im[l]*R_im[l,k] *ψ_re[k] + ψ_im[l] * R_re[l,k]*ψ_im[k] for l = 1:dB, k = 1:dB))
        # t2 = <ψ, Y ψ>, exploiting that this quantity is real
        @NLexpression(m, t2, sum(ψ_re[l] * Y_re[l,k] * ψ_re[k] - ψ_re[l] * Y_im[l,k] * ψ_im[k] + ψ_im[l]*Y_im[l,k] *ψ_re[k] + ψ_im[l] * Y_re[l,k]*ψ_im[k] for l = 1:dB, k = 1:dB))
        @NLobjective(m, Min, t1 - t2)
    end
    @timeit timer "optimize!" JuMP.optimize!(m)

    if termination_status(m) != MOI.OPTIMAL
        error("Problem was not solved optimally; status $(termination_status(m)).")
    end

    obj_value = objective_value(m)
    feasible = obj_value >= 0

    @info "true_feasiblity" obj_value feasible

    if objective_value(m) >= 0
        return (status=FEASIBLE, RY=nothing, π=nothing)
    else
        # Time to find a cut
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
