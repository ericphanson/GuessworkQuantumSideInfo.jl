using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack
using GenericLinearAlgebra
using SCS
using GuessworkQuantumSideInfo: vec_to_herm, herm_to_vec, PermutationIterator

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = false # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = true # incompatible with the current version of Pajarito
TEST_ELLIPSOID = false

default_sdp_solver() = TEST_MOI ? SCS.Optimizer(verbose = 0, eps = 1e-6) : SCSSolver(verbose = 0, eps = 1e-6)

if TEST_MATLAB
    include("test_matlab.jl")
end

# Performs various generic checks
using Aqua
Aqua.test_unbound_args(GuessworkQuantumSideInfo)
Aqua.test_undefined_exports(GuessworkQuantumSideInfo)

# This hangs:
# Aqua.test_ambiguities(GuessworkQuantumSideInfo)
# Probably https://github.com/JuliaLang/julia/issues/35909
rand_vec(T, n) = rand(T, n) - rand(T, n)
rand_herm(T, d) = Hermitian(rand(T, d,d) - rand(T, d, d) + im*( rand(T, d,d) - rand(T, d, d)))
@testset "Utilities" begin
    @testset "Hermitian basis" begin
        for T in (Float64, BigFloat)
            for d in (2,3,5)
                v = rand_vec(T, d^2)
                M = vec_to_herm(v)
                # Isometry
                @test norm(v, 2) ≈ norm(M, 2)
                # Maps to Hermitian
                @test M isa Hermitian
                # Roundtrips
                @test v ≈ herm_to_vec(M)
                # Linear
                v2 = rand_vec(T, d^2)
                λ = rand(T) - rand(T)
                @test vec_to_herm(v + λ*v2) ≈ vec_to_herm(v) + λ*vec_to_herm(v2)

                # Test the other direction
                M = rand_herm(T, d)
                v = herm_to_vec(M)
                @test norm(M, 2) ≈ norm(v,2)
                @test eltype(v) <: Real
                @test v isa AbstractVector
                @test M ≈ vec_to_herm(v)
                M2 = rand_herm(T, d)
                λ = rand(T) - rand(T)
                @test herm_to_vec(M + λ*M2) ≈ herm_to_vec(M) + λ*herm_to_vec(M2)

                # Test inner product
                M1 = rand_herm(T, d)
                M2 = rand_herm(T, d)
                v1 = herm_to_vec(M1)
                v2 = herm_to_vec(M2)
                @test tr(M1' * M2) ≈ dot(v1, v2)
                v1 = rand_vec(T, d^2)
                v2 = rand_vec(T, d^2)
                M1 = vec_to_herm(v1)
                M2 = vec_to_herm(v2)
                @test tr(M1' * M2) ≈ dot(v1, v2)

            end
        end
    end

    @testset "Birkhoff's theorem" begin
        for d in (2,3,5), T in (Float64, BigFloat)
            n = 5
            α_init = rand(T, n)
            α_init = α_init / sum(α_init)
            D = sum( α_init[i]*I(d)[randperm(d), :] for i = 1:n )

            # Reconstruct it from permutations
            D_reconstruct = zero(D)
            for (π_i, α_i) in PermutationIterator(D)
                P_i = I(d)[π_i, :]
                D_reconstruct += α_i * P_i
            end
            @test D_reconstruct ≈ D

            # Make sure `collect` works (tests length etc)
            results = collect(PermutationIterator(D))
            D_reconstruct2 = sum(α*I(d)[π, :] for (π, α) in results)
            @test D_reconstruct2 ≈ D
        end
    end
end

include("test_problems.jl")
