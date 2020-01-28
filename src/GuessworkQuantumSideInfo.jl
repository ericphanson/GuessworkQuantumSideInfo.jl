module GuessworkQuantumSideInfo

using LinearAlgebra, SparseArrays # standard libraries
using Random: randperm! # standard library
using Convex # Julia SDP solvers / interfaces
using Combinatorics: multiset_permutations, combinations # used for POVM outcomes
using UnPack # helper

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

end
