# GuessworkQuantumSideInfo

[![Build Status](https://github.com/ericphanson/GuessworkQuantumSideInfo.jl/workflows/CI/badge.svg)](https://github.com/ericphanson/GuessworkQuantumSideInfo.jl/actions)
[![Coverage](https://codecov.io/gh/ericphanson/GuessworkQuantumSideInfo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ericphanson/GuessworkQuantumSideInfo.jl)
[![Stable documentation](https://img.shields.io/badge/documentation-stable-blue.svg)](https://ericphanson.github.io/GuessworkQuantumSideInfo.jl/stable)
[![Dev documentation](https://img.shields.io/badge/documentation-dev-blue.svg)](https://ericphanson.github.io/GuessworkQuantumSideInfo.jl/dev)
[![arXiv article](https://img.shields.io/badge/article-arXiv%3A2001.03598-B31B1B)](https://arxiv.org/abs/2001.03598)

This is a package accompanying the preprint [*Guesswork with Quantum Side Information*](https://arxiv.org/abs/2001.03598).

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

```julia
julia> using GuessworkQuantumSideInfo, SCS

julia> p = [0.25, 0.25, 0.25, 0.25];

julia> ketzero = ket(1, 2);

julia> ketone = ket(2, 2);

julia> ketminus = (ket(1, 2) - ket(2,2))/sqrt(2);

julia> ketplus = (ket(1, 2) + ket(2,2))/sqrt(2);

julia> ρBs = dm.([ ketzero, ketone, ketminus, ketplus  ])
4-element Array{Array{Complex{Float64},2},1}:
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]                                                              
 [0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]                                                              
 [0.4999999999999999 + 0.0im -0.4999999999999999 - 0.0im; -0.4999999999999999 + 0.0im 0.4999999999999999 + 0.0im]
 [0.4999999999999999 + 0.0im 0.4999999999999999 + 0.0im; 0.4999999999999999 + 0.0im 0.4999999999999999 + 0.0im]  

julia> output = guesswork(p, ρBs; solver = SCSSolver(verbose=false));

julia> output.optval
1.709431078700102

```

It turns out it takes `(1/4)*(10 - sqrt(10)) ≈ 1.71` guesses on average.
