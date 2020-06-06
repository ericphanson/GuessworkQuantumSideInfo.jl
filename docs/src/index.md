```@meta
CurrentModule = GuessworkQuantumSideInfo
```

# GuessworkQuantumSideInfo

This is a package accompanying the preprint [*Guesswork with Quantum Side Information*](http://arxiv.org/abs/2001.03598).

See the [Examples](@ref) for some examples (or just below for a quick example),
[Computing the guesswork to high precision](@ref) for an example solving a
problem with a high-precision SDP solver, [Using a mixed-integer SDP to find
extremal strategies](@ref) for an example using a mixed-integer SDP, or below
for the documentation of the functions provided by this package.

## Quick example

Consider one party Alice who draws a random number in the set `[1,2,3,4]`
uniformly at random. If she draws `1` she sends another party, Bob, the quantum
state `|0⟩`; if she draws `2`, she sends `|1⟩`, if she draws `3` she sends
`|-⟩`, and finally if she draws `4`, she sends `|+⟩`. Bob, knowing this general
procedure but not which number Alice drew, aims to guess the value Alice drew by
performing experiments on the quantum state he was given. The average number of
guesses Bob needs in order to get the right answer, minimized over all quantum
strategies, is the so-called *guesswork with quantum side information*. This
package provides a means to compute this.

```@repl
using GuessworkQuantumSideInfo, SCS
p = [0.25, 0.25, 0.25, 0.25];
ketzero = ket(1, 2);
ketone = ket(2, 2);
ketminus = (ket(1, 2) - ket(2,2))/sqrt(2);
ketplus = (ket(1, 2) + ket(2,2))/sqrt(2);
ρBs = dm.([ ketzero, ketone, ketminus, ketplus  ])
output = guesswork(p, ρBs; solver = SCSSolver(verbose=false));
output.optval
```

## Guesswork functions

```@docs
GuessworkQuantumSideInfo.guesswork
GuessworkQuantumSideInfo.guesswork_lower_bound
GuessworkQuantumSideInfo.guesswork_upper_bound
GuessworkQuantumSideInfo.guesswork_ellipsoid
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

## Ellipsoid algorithm functions

```@docs
GuessworkQuantumSideInfo.default_init
GuessworkQuantumSideInfo.EllipsoidProblem
GuessworkQuantumSideInfo.invherm
GuessworkQuantumSideInfo.herm
GuessworkQuantumSideInfo.ellipsoid_algorithm!
GuessworkQuantumSideInfo.PermutationIterator
```
