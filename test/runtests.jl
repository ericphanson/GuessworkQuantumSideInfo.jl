using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack
using GenericLinearAlgebra
using SCS
using GuessworkQuantumSideInfo: herm, invherm

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = false # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = true # incompatible with the current version of Pajarito
TEST_ELLIPSOID = false # does not pass yet since currently a heuristic

@info "Beginning tests with" TEST_MATLAB TEST_MISDP TEST_BB84_MISDP TEST_MOI TEST_ELLIPSOID

default_sdp_solver() = TEST_MOI ? SCS.Optimizer(verbose = 0, eps = 1e-6) : SCSSolver(verbose = 0, eps = 1e-6)

# bad solve error on 1.0 and nightly observed on CI (though not locally)
# possibly reflects problems with SCS, but not with this package,
# so we'll just use a relaxed tolerance here.
TOL = 1e-2

if TEST_MISDP
    using Pajarito, Cbc
    function misdp_solver(; verbose = false)
        sdp_solver = default_sdp_solver()
        mip_solver = CbcSolver(loglevel = 0)
        PajaritoSolver(
            cont_solver = sdp_solver,
            mip_solver = mip_solver,
            mip_solver_drives = false,
            use_mip_starts = true,
            solve_relax = false,
            log_level = verbose ? 3 : 0,
        )
    end
end

if TEST_MATLAB
    include("test_matlab.jl")
end

include("test_utilities.jl")

dB = 2
ketplus = (ket(1, dB) + ket(2, dB)) / sqrt(2)
ketminus = (ket(1, dB) - ket(2, dB)) / sqrt(2)
ketzero = ket(1, dB)
ketone = ket(2, dB)


Random.seed!(5)

@testset "GuessworkQuantumSideInfo.jl" begin

    @testset "Hermitian basis" begin
        for T in (Float64, BigFloat)
            for d in (2,3,5)
                v = rand(T, d^2) - rand(T, d^2)
                M = herm(v)
                # Isometry
                @test norm(v, 2) ≈ norm(M, 2)
                # Maps to Hermitian
                @test M isa Hermitian
                # Roundtrips
                @test v ≈ invherm(M)
                # Linear
                v2 = rand(T, d^2) - rand(T, d^2)
                λ = rand(T) - rand(T)
                @test herm(v + λ*v2) ≈ herm(v) + λ*herm(v2)

                # Test the other direction
                M = Hermitian(rand(T, d,d) - rand(T, d, d) + im*( rand(T, d,d) - rand(T, d, d)))
                v = invherm(M)
                @test norm(M, 2) ≈ norm(v,2)
                @test eltype(v) <: Real
                @test v isa AbstractVector
                @test M ≈ herm(v)
                M2 = Hermitian(rand(T, d,d) - rand(T, d, d) + im*( rand(T, d,d) - rand(T, d, d)))
                λ = rand(T) - rand(T)
                @test invherm(M + λ*M2) ≈ invherm(M) + λ*invherm(M2)
            end
        end
    end

    @testset "Uninformative side information" begin
        ρBs = dm.([ketzero, ketzero])
        p = [0.5, 0.5]
        output = test_optimize(p, ρBs, 1.5)

        pmf = pmfN(output)
        @test pmf ≈ [0.5, 0.5] rtol = TOL

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + TOL

        ub = guesswork_upper_bound(p, ρBs; num_constraints = 2^2 + 1, make_solver = default_sdp_solver).optval
        @test output.optval <= ub + TOL

    end

    @testset "Orthogonal side information" begin
        ρBs = dm.([ketzero, ketone])
        p = [0.5, 0.5]
        output = test_optimize(p, ρBs, 1.0)

        pmf = pmfN(output)

        @test pmf ≈ [1.0, 0.0] rtol = TOL

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + TOL

        ub = guesswork_upper_bound(p, ρBs; num_constraints =  2^2 + 1, make_solver = default_sdp_solver).optval
        @test output.optval <= ub + TOL
    end

    @testset "Plus zero" begin
        ρBs = dm.([ketplus, ketzero])
        p = [0.5, 0.5]

        output = test_optimize(p, ρBs, cos(π / 8)^2 + 2 * sin(π / 8)^2)

        @test pmfN(output) ≈ [cos(π / 8)^2, sin(π / 8)^2] rtol = TOL

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + TOL

        ub = guesswork_upper_bound(p, ρBs; num_constraints =  2^2 + 1, make_solver = default_sdp_solver).optval
        @test output.optval <= ub + TOL
    end

    @testset "Three random qubits" begin
        ρBs = [randdm(2) for _ = 1:3]
        p = ones(3) / 3

        @testset "Basic tests" begin
            output = test_optimize(p, ρBs; test_MISDP = false)
            
            lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
            @test lb <= output.optval + TOL

            ub = guesswork_upper_bound(p, ρBs; num_constraints =  2^2 + 1, make_solver = default_sdp_solver).optval
            @test output.optval <= ub + TOL
        end

        @testset "Custom cost vector" begin
            c = cumsum(10*rand(3)) # increasing vector
            output_c = test_optimize(p, ρBs; test_MISDP = false, c = c)
            lb_c = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver(), c = c).optval
            @test lb_c <= output_c.optval + TOL
            ub_c = guesswork_upper_bound(p, ρBs; num_constraints =  2^2 + 1, make_solver = default_sdp_solver, c = c).optval
            @test output_c.optval <= ub_c + TOL
        end
    end

    @testset "Three random qutrits" begin
        ρBs = [randdm(3) for _ = 1:3]
        p = ones(3) / 3

        output = test_optimize(p, ρBs; test_MISDP = false)
        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + TOL

        ub = guesswork_upper_bound(p, ρBs; num_constraints =  3^2 + 1, make_solver = default_sdp_solver).optval
        @test output.optval <= ub + TOL
    end


    @testset "Four random qubits" begin
        ρBs = [randdm(2) for _ = 1:4]
        p = [0.25, 0.25, 0.25, 0.25]

        @testset "Basic tests" begin
            output = test_optimize(p, ρBs; test_MISDP = false)

            lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
            @test lb <= output.optval + TOL

            ub = guesswork_upper_bound(p, ρBs; num_constraints =  2^2 + 1, make_solver = default_sdp_solver).optval
            @test output.optval <= ub + TOL
        end

        @testset "K=2 guesses" begin
            output_K = test_optimize(p, ρBs; test_MISDP = false, K = 2)
            pmf = pmfN(output_K)
            @test length(pmf) == 3
            @test sum(pmf) ≈ 1.0 atol = TOL
            @test all(x -> x >= -TOL, pmf)
        end
    end

    @testset "BB84" begin
        ρBs = BB84_states()
        p = ones(4) / 4
        output = guesswork(p, ρBs; solver = default_sdp_solver())
        testPOVM(output.Es)

        if TEST_BB84_MISDP && TEST_MISDP
            @info "Starting BB84 MISDP"
            output_MISDP, t, _ = @timed guesswork_MISDP(p, ρBs, 2; solver = misdp_solver())
            @info "Finished BB84 MISDP in $(round(t;digits=3)) seconds."
            testPOVM(output_MISDP.Es)
            @test output_MISDP.optval ≈ output.optval rtol = TOL
        end

        relaxed_output = guesswork_upper_bound(p, ρBs; num_constraints = 4, make_solver = default_sdp_solver)

        @test output.optval ≈ relaxed_output.optval rtol = TOL

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + TOL

    end

    @testset "BB84 variant" begin
        ρBs =
            dm.([
                ketzero,
                ketone,
                (ketzero + im * ketone) / sqrt(2),
                (ketzero - im * ketone) / sqrt(2),
            ])
        J = length(ρBs)
        p = ones(J) / J
        output = test_optimize(p, ρBs; test_MISDP = false)
        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + TOL

        ub = guesswork_upper_bound(p, ρBs; num_constraints =  2^2 + 1, make_solver = default_sdp_solver).optval
        @test output.optval <= ub + TOL
    end

    @testset "Errors and logging" begin
        ρBs = BB84_states()
        J = length(ρBs)
        p = ones(J) / J
        @test_throws ArgumentError guesswork_upper_bound(p, ρBs; num_constraints = Inf, max_time=Inf, max_retries = Inf, make_solver = default_sdp_solver)
        @test_logs (:info, r"Adding constraint") match_mode=:any guesswork_upper_bound(p, ρBs; verbose=true, make_solver = default_sdp_solver, max_time = 2)

        two_states = [randdm(2) for _ = 1:2]
        p_3 = randprobvec(3)

        if TEST_MISDP
            @test_logs (:info, r"MISDP solve") match_mode=:any guesswork_MISDP(randprobvec(2), [randdm(2) for _ = 1:2], 5; verbose=true, solver = misdp_solver())
            @test_throws ArgumentError guesswork_MISDP(p_3, two_states, 5; solver = misdp_solver())
        end
        
        @test_throws ArgumentError guesswork(p_3, two_states; solver = default_sdp_solver())
        @test_throws ArgumentError guesswork_lower_bound(p_3, two_states; solver = default_sdp_solver())
        @test_throws ArgumentError guesswork_upper_bound(p_3, two_states; make_solver = default_sdp_solver)
    end


    @testset "Concavity of the guesswork" begin
        J = 3
        dB = 2
        p_1 = randprobvec(J)
        ρBs_1 = [randdm(dB) for _ = 1:J]

        p_2 = randprobvec(J)
        ρBs_2 = [randdm(dB) for _ = 1:J]

        λ = rand()

        p = λ .* p_1 + (1-λ) .* p_2
        ρBs = λ .* ρBs_1 + (1-λ) .* ρBs_2

        g_avg = guesswork(p, ρBs; solver = default_sdp_solver()).optval
        g_1 = guesswork(p_1, ρBs_1; solver = default_sdp_solver()).optval
        g_2 = guesswork(p_2, ρBs_2; solver = default_sdp_solver()).optval
        @test TOL + g_avg >= λ*g_1 + (1-λ)*g_2
    end

    @testset "Quantum states" begin
        for T in (BigFloat, Float64)
            for d in (2,3)
                ρ = randdm(T, d)
                @test eltype(ρ) == Complex{T}
                @test tr(ρ) ≈ one(T) rtol=1e-8
                @test all( x-> x >= -1e-8, eigvals(ρ))

                U = GuessworkQuantumSideInfo.randunitary(T, d)
                @test eltype(U) == Complex{T}
                @test U' * U ≈ I(d) atol=1e-6
                @test U * U' ≈ I(d) atol=1e-6

                k = ket(T, 1, d)
                @test k' ≈ bra(T, 1, d) rtol=1e-8

                p = randprobvec(T, d)
                @test eltype(p) == T
                @test sum(p) ≈ one(T) rtol=1e-8
                @test all( x-> x >= -1e-8, p)
            end

            ρBs = BB84_states(T)
            @test length(ρBs) == 4
            @test Set(iid_copies(ρBs, 2)) == Set([ ρ ⊗ σ for ρ in ρBs for σ in ρBs])
        end
        
        # Defaults
        @test BB84_states() ≈ BB84_states(Float64)
        @test BB84_states()  ≈ dm.([ketzero, ketone, ketminus, ketplus]) atol=1e-6

        @test ket(1,2) ≈ ket(Float64, 1, 2)
        @test bra(1, 2) ≈ bra(Float64, 1, 2)
        @test eltype(randdm(2)) == Complex{Float64}
        @test eltype(GuessworkQuantumSideInfo.randunitary(2)) == Complex{Float64}
        @test eltype(randprobvec(2)) == Float64
    end
end
