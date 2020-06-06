using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack
using GenericLinearAlgebra
using SCS

using EAGO, JuMP

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = false # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = true # incompatible with the current version of Pajarito
TEST_ELLIPSOID = true

default_sdp_solver() = SCS.Optimizer(verbose = 0, eps = 1e-6)

function nl_solver()
    with_optimizer(() -> EAGO.Optimizer(verbosity=0))
end

using GuessworkQuantumSideInfo: EllipsoidProblem
using Serialization

@testset "Serialization" begin
    p = ones(2)/2
    ρBs = [ randdm(2) for _ = 1:2]
    prob = EllipsoidProblem(p, ρBs, nl_solver = nl_solver())
    local sol, re_prob
    mktemp() do path, io
        serialize(path, prob)
        sol = GuessworkQuantumSideInfo.ellipsoid_algorithm!(prob)
        re_prob = deserialize(path)
    end
    re_sol = GuessworkQuantumSideInfo.ellipsoid_algorithm!(re_prob)
    @test sol.prob.p ≈ re_sol.prob.p
    @test sol.prob.ρBs ≈ re_sol.prob.ρBs
    @test sol.optval ≈ re_sol.optval rtol=1e-2
    @test sol.prob.x ≈ re_sol.prob.x rtol=1e-2
end

include("../test_problems.jl")
