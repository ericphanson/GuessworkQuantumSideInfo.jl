using GuessworkQuantumSideInfo
using SDPAFamily
using Test
using LinearAlgebra, GenericLinearAlgebra

# Fix a compatibility problem
LinearAlgebra.eigmin(A::Hermitian{Complex{BigFloat},Array{Complex{BigFloat},2}}) =
    minimum(real.(eigvals(A)))


function testPOVM(Es; tol = 1e-25)
    dB = size(first(Es), 1)
    @test sum(Es) ≈ complex(1.0) * I(dB) atol = tol

    for E in Es
        @test E ≈ E' atol = tol
    end
    @test all(isposdef, Hermitian(E) + tol * I(dB) for E in Es)
end

T = BigFloat
default_sdp_solver(T) = SDPAFamily.Optimizer{T}(presolve = true)

@testset "BB84" begin
    ρBs = GuessworkQuantumSideInfo.BB84_states(T)

    p = ones(T, 4) / 4
    output = guesswork(p, ρBs; solver = default_sdp_solver(T))
    testPOVM(output.Es)
    relaxed_output =
        guesswork_upper_bound(p, ρBs; num_constraints = 4, make_solver = () -> default_sdp_solver(T))

    true_val = (big(1) / big(4)) * (10 - sqrt(big(10)))
    @test output.optval ≈ true_val rtol = 1e-25

    @test output.optval ≈ relaxed_output.optval rtol = 1e-4

    lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver(T)).optval
    @test lb <= output.optval + 1e-4
end
