module GuessworkQuantumSideInfo

using Dates, LinearAlgebra, Logging, SparseArrays # standard libraries
using Random: randperm! # standard library
using Convex # Julia SDP solvers / interfaces
using Combinatorics: multiset_permutations, combinations # used for POVM outcomes
using UnPack: @unpack # helper

# Used in the ellipsoid method implementation:
using Hungarian: hungarian # solves the perfect matching problem
using TimerOutputs: TimerOutput, @timeit # timing functionality
using SpecialFunctions: gamma # only used in `unit_ball_volume`
using JuMP, MathOptInterface # an optimization problem modelling language 
const MOI = MathOptInterface # common shorthand

# Basic functions for quantum states and examples
export ket, bra, ⊗, I, dm, randdm, randprobvec, iid_copies, BB84_states
include("quantum_states.jl")

# The basic SDP routine
export guesswork
include("SDP_formulation.jl")

export guesswork_lower_bound
include("lower_bound.jl")

# Re-run the SDP adding constraints one at a time
export guesswork_upper_bound
include("upper_bound.jl")

# reformulate as a mixed-integer SDP
export guesswork_MISDP
include("MISDP_formulation.jl")

# Utilities for analyzing the POVMs resulting from the SDPs
export pmfN
include("analyze_measurements.jl")

# implementation of the ellipsoid algorithm
export guesswork_ellipsoid
include("ellipsoid/ellipsoid_algorithm.jl")
include("ellipsoid/Birkhoff_vN_decomposition.jl")
include("ellipsoid/separation_oracle.jl")

end
