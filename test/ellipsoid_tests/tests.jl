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

include("../test_problems.jl")
