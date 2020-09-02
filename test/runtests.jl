using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack, StableRNGs
using GenericLinearAlgebra
using SCS

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = false # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = true # incompatible with the current version of Pajarito

default_sdp_solver() = TEST_MOI ? SCS.Optimizer(verbose = 0, eps = 1e-6) : SCSSolver(verbose = 0, eps = 1e-6)

if TEST_MATLAB
    include("test_matlab.jl")
end

include("test_problems.jl")
