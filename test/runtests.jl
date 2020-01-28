using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack
using Combinatorics: multiset_permutations
using SCS

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = true # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = false # incompatible with the current version of Pajarito

@info "Beginning tests with" TEST_MATLAB TEST_MISDP TEST_BB84_MISDP TEST_MOI

default_sdp_solver() = TEST_MOI ? SCS.Optimizer(verbose = 0) : SCSSolver(verbose = 0)

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

    @testset "randpure" begin
        for d in (2, 3, 5)
            for _ = 1:5
                ψ = randpure(d)
                @test tr(ψ) ≈ 1.0
                @test tr(ψ^2) ≈ 1.0 atol = 1e-6
                @test imag.(eigvals(ψ)) ≈ zeros(d) atol = 1e-8
                @test minimum(real.(eigvals(ψ))) ≈ 0.0 atol = 1e-8
            end
        end

    end

    @testset "Uninformative side information" begin
        ρBs = dm.([ketzero, ketzero])
        p = [0.5, 0.5]
        output = test_optimize(p, ρBs, 1.5)

        pmf = pmfN(output)
        @test pmf ≈ [0.5, 0.5] rtol = 1e-4

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 2^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4

    end

    @testset "Orthogonal side information" begin
        ρBs = dm.([ketzero, ketone])
        p = [0.5, 0.5]
        output = test_optimize(p, ρBs, 1.0)

        pmf = pmfN(output)

        @test pmf ≈ [1.0, 0.0] rtol = 1e-4

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 2^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4
    end

    @testset "Plus zero" begin
        ρBs = dm.([ketplus, ketzero])
        p = [0.5, 0.5]

        output = test_optimize(p, ρBs, cos(π / 8)^2 + 2 * sin(π / 8)^2)

        @test pmfN(output) ≈ [cos(π / 8)^2, sin(π / 8)^2] rtol = 1e-4

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 2^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4
    end

    @testset "Three random qubits" begin
        ρBs = [randdm(2) for _ = 1:3]
        p = ones(3) / 3

        output = test_optimize(p, ρBs; test_MISDP = false)

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 2^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4
    end

    @testset "Three random qutrits" begin
        ρBs = [randdm(3) for _ = 1:3]
        p = ones(3) / 3

        output = test_optimize(p, ρBs; test_MISDP = false)
        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 3^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4
    end


    @testset "Four random qubits" begin
        ρBs = [randdm(2) for _ = 1:4]
        p = [0.25, 0.25, 0.25, 0.25]


        output = test_optimize(p, ρBs; test_MISDP = false)

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 2^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4
    end

    @testset "BB84" begin
        ρBs = GuessworkQuantumSideInfo.BB84_states()

        # make sure we have the right states...
        @test Set(ρBs) == Set(dm.([ketzero, ketone, ketplus, ketminus]))
        @test length(ρBs) == 4

        p = ones(4) / 4
        output = guesswork(p, ρBs; solver = default_sdp_solver())
        testPOVM(output.Es)

        if TEST_BB84_MISDP && TEST_MISDP
            @info "Starting BB84 MISDP"
            output_MISDP, t, _ = @timed guesswork_MISDP(p, ρBs, 2; solver = misdp_solver())
            @info "Finished BB84 MISDP in $(round(t;digits=3)) seconds."
            testPOVM(output_MISDP.Es)
            @test output_MISDP.optval ≈ output.optval rtol = 1e-4
        end

        relaxed_output = guesswork_upper_bound(p, ρBs, 4; make_solver = default_sdp_solver)

        @test output.optval ≈ relaxed_output.optval rtol = 1e-4

        lb = guesswork_lower_bound(p, ρBs; solver = default_sdp_solver()).optval
        @test lb <= output.optval + 1e-4

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
        @test lb <= output.optval + 1e-4

        ub = guesswork_upper_bound(p, ρBs, 2^2 + 1; make_solver = default_sdp_solver).optval
        @test output.optval <= ub + 1e-4
    end

    @testset "Upper bound" begin
        ρBs = GuessworkQuantumSideInfo.BB84_states()
        J = length(ρBs)
        p = ones(J) / J
        @test_throws ArgumentError guesswork_upper_bound(p, ρBs, Inf; max_time=Inf, max_retries = Inf, make_solver = default_sdp_solver)
        @test_logs (:info, r"Adding constraint") match_mode=:any guesswork_upper_bound(p, ρBs; verbose=true, make_solver = default_sdp_solver, max_time = 2)
    end


    @testset "Concavity of the guesswork" begin
        J = 3
        dB = 2
        p_1 = GuessworkQuantumSideInfo.randprobvec(J)
        ρBs_1 = [randdm(dB) for _ = 1:J]

        p_2 = GuessworkQuantumSideInfo.randprobvec(J)
        ρBs_2 = [randdm(dB) for _ = 1:J]

        λ = rand()

        p = λ .* p_1 + (1-λ) .* p_2
        ρBs = λ .* ρBs_1 + (1-λ) .* ρBs_2

        g_avg = guesswork(p, ρBs; solver = default_sdp_solver()).optval
        g_1 = guesswork(p_1, ρBs_1; solver = default_sdp_solver()).optval
        g_2 = guesswork(p_2, ρBs_2; solver = default_sdp_solver()).optval
        @test 1e-4 + g_avg >= λ*g_1 + (1-λ)*g_2
    end
end
