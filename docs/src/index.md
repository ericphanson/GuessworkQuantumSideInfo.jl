```@meta
CurrentModule = GuessworkQuantumSideInfo
```

# GuessworkQuantumSideInfo

See the [Examples](@ref) for some quick examples, [Computing the guesswork to high precision](@ref)
for an example solving a problem with a high-precision SDP solver,
[Using a mixed-integer SDP to find extremal strategies](@ref) for an example
using a mixed-integer SDP, or below for the documentation of the functions
provided by this package.

## Guesswork functions
```@docs
GuessworkQuantumSideInfo.guesswork
GuessworkQuantumSideInfo.guesswork_lower_bound
GuessworkQuantumSideInfo.guesswork_upper_bound
```

## Quantum states

```@docs
GuessworkQuantumSideInfo.ket
GuessworkQuantumSideInfo.bra
GuessworkQuantumSideInfo.BB84_states
GuessworkQuantumSideInfo.iid_copies
GuessworkQuantumSideInfo.randdm
GuessworkQuantumSideInfo.randprobvec
```

## Utilities

```@docs
GuessworkQuantumSideInfo.pmfN
```
