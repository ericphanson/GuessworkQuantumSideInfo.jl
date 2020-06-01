@Base.kwdef struct PermutationIterator{T1, T2, T3}
    D::T1
    solver::T2
    rtol::T3 = 1e-6
end

function Base.iterate(PI::PermutationIterator, state = (PI.D, 0))
    state === nothing && return nothing
    _it(PI, state)
end

function _it(PI, (D,c))
    @unpack solver, rtol = PI
    @unpack P,α0, D1, done = perms(D, solver, rtol)
    c += 1
    # done = done || c >= size(D,1)
    state = done ? nothing : (D1,c)
    return P, state
end


function perms(D::AbstractMatrix, solver, rtol)
    n = size(D,1)
    P = maximum_weight_perm(D, solver=solver)
    α0 = minimum( D[i,j] / P[i,j] for i = 1:n, j = 1:n if P[i,j] > 0 )
    done = isapprox(α0, 1; rtol=rtol)
    D1 = done ? D : (D - α0*P)/(1-α0)
    return (P=P, D1=D1, done=done, α0 =α0)
end

function maximum_weight_perm(D; solver)
    model = Model(with_optimizer(solver))
    n = size(D,1)
    @variable(model, P[1:n, 1:n] >= 0)
    @objective(model, Max, sum(P[i,j]*D[i,j] for i = 1:n, j=1:n))
    @constraint(model, sum(P, dims = 1) .<= 1)
    @constraint(model, sum(P, dims = 2) .<= 1)
    JuMP.optimize!(model)
    @assert termination_status(model) == MOI.OPTIMAL
    value.(P)
end
