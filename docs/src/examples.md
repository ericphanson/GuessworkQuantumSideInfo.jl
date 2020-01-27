Let's load the package, an SDP solver, and define a simple plotting routine.
```@example ex
using GuessworkQuantumSideInfo
using SCS, Plots
get_sdp_solver() = SCSSolver(verbose=false)
plot_pmfN(data) = bar(pmfN(data); xlabel="Guess number", ylabel="Probability of guessing correctly", legend = false)
```

Next, we define some basic qubit states.

```@example ex
dB = 2
ketplus = (ket(1, dB) + ket(2,dB))/sqrt(2)
ketminus = (ket(1, dB) - ket(2,dB))/sqrt(2)
ketzero = ket(1, dB)
ketone = ket(2, dB)
```


## Example 1: A warmup with trivial examples

Let's consider the case with $J=2$ and both states are the same. The side
information is therefore completely uninformative, and the guesswork is

```math
G(X|B) = \frac{1}{2}\cdot 1 + \frac{1}{2} \cdot 2 = 1.5.
```

We can check this:

```@repl ex
p = [0.5, 0.5]
ρBs = dm.([ ketzero, ketzero  ])
output = guesswork(p, ρBs; solver = get_sdp_solver());
output.optval
```

We see the result agrees with `1.5`, as we expected. Likewise, if we choose the
two states as $ |0\rangle, |1\rangle$, we can get it in one guess every time, of
course, since they are orthogonal:

```@repl ex
p = [0.5, 0.5]
ρBs = dm.([ ketzero, ketone  ])
output = guesswork(p, ρBs; solver = get_sdp_solver());
output.optval
```

We can inspect the POVMs:
```@repl ex
output.Es
```

As we would expect, we (approximately) obtain the projection onto $|0 \rangle$
and the projection onto $|1 \rangle$.


## Example 2: the BB84 states

Let's consider the four states $|+ \rangle, |-\rangle, |0\rangle, |1\rangle$.

```@repl ex
p = [0.25, 0.25, 0.25, 0.25]
ρBs = dm.([ ketplus, ketminus, ketzero, ketone  ])
output = guesswork(p, ρBs; solver = get_sdp_solver());
output.optval
```

```@example ex
plot_pmfN(output)
```


Let's try the same example, but imposing a steep cost for the fourth guess.

```@repl ex
c = [1.0, 2.0, 3.0, 5000.0]
output = guesswork(p, ρBs; c = c, solver = get_sdp_solver());
output.optval
```

We see that the average number of guesses to get a correct answer has gone up.
However, inspecting the probability mass function for the number of guesses
under the optimal strategy

```@example ex
plot_pmfN(output)
```

we see that the probability that the probability of guessing correctly on the
fourth guess goes to almost zero.

## Example 3: two copies of the BB84 states

Let us consider two tensor copies of the BB84 states:

```@repl ex
p = ones(16)/16
bb84 = GuessworkQuantumSideInfo.BB84_states()
ρBs = GuessworkQuantumSideInfo.iid_copies(bb84, 2)
```

In this case, there are $16! = 20922789888000$ possible guessing orders, and
hence $16!$ variables in the primal formulation of the SDP, or $16!+1$
constraints in the dual form of the SDP. In either case, we can't even fit them
all into our computer's memory. Instead, we resort to bounds:

```@repl ex
lb_output = guesswork_lower_bound(p, ρBs, solver = get_sdp_solver());
lb_output.optval

ub_output = guesswork_upper_bound(p, ρBs; max_time = 30, make_solver = get_sdp_solver);
ub_output.optval
```
