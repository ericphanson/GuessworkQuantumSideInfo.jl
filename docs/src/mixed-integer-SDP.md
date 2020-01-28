# Using a mixed-integer SDP to find extremal strategies

Let us revisit [Example 2: the BB84 states](@ref) using a mixed-integer SDP.
This allows us to specify the number of non-zero POVM elements in the
measurement-- or, equivalently, the number of possible guessing orders the final
strategy chooses among.

We will use the [Pajarito.jl](https://github.com/JuliaOpt/Pajarito.jl) solver,
which solves mixed-integer SDPs by solving an alternating sequence of
mixed-integer linear programs and SDPs; see the [Pajarito
paper](https://arxiv.org/abs/1808.05290) for information on how Pajarito works.
Pajarito requires both an SDP solver and a mixed-integer linear solver; we will
use the open source solvers SCS and Cbc, respectively.

```@repl misdp
using GuessworkQuantumSideInfo
using Pajarito, Cbc, SCS # solvers
function misdp_solver(; verbose = false)
    sdp_solver = SCSSolver(verbose=0)
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

p = ones(4)/4
ρBs = BB84_states()
output = guesswork_MISDP(p, ρBs, 2; verbose=false, solver = misdp_solver());
output.optval
output.Es
```

We see that with only two measurement outcomes we recover the same optimal value
as the case without a constraint on  the number of measurement outcomes (in
[Example 2: the BB84 states](@ref)). We've thus found an extremal strategy (in
that the POVM associated to this strategy cannot be written as the convex
combination of two other POVMs). Often, the solutions returned by
[`guesswork`](@ref) do not return extremal POVMs, although the details depend on
the SDP solver used.

We can see what the associated guessing orders are:

```@repl misdp
output.povm_outcomes
```

To summarize, one optimal strategy for the case of `p` uniform, and `ρBs` given
by the four BB84 states, is to perform a projective measurement whose operators
are given by `output.Es` above, and then make guesses in one of the two orders
given by `output.povm_outcomes` (depending on which measurement outcome was
obtained).
