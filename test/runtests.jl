using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack
using GenericLinearAlgebra
using SCS
using GuessworkQuantumSideInfo: herm, invherm, PermutationIterator

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = false # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = true # incompatible with the current version of Pajarito
TEST_ELLIPSOID = false

default_sdp_solver() = TEST_MOI ? SCS.Optimizer(verbose = 0, eps = 1e-6) : SCSSolver(verbose = 0, eps = 1e-6)

if TEST_MATLAB
    include("test_matlab.jl")
end

@testset "Utilities" begin
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
        end
    end
end

include("test_problems.jl")
