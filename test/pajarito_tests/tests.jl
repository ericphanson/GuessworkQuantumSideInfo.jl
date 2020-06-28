using Test, GuessworkQuantumSideInfo, LinearAlgebra, Statistics, Random, UnPack
using GenericLinearAlgebra
using SCS

using Pajarito, Cbc

TEST_MATLAB = false # requires MATLAB and MATLAB.jl installed
TEST_MISDP = true # requires Pajarito or another MISDP solver
TEST_BB84_MISDP = false # takes ~100 seconds; requires TEST_MISDP
TEST_MOI = false # incompatible with the current version of Pajarito
TEST_ELLIPSOID = false

default_sdp_solver() = SCSSolver(verbose = 0, eps = 1e-6)

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

include("../test_problems.jl")
