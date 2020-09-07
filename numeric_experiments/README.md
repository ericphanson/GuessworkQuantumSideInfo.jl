*The following is taken from my PhD thesis (translated from LaTeX to Markdown by Pandoc 2.10).*

- Eric Hanson

# Numeric experiments

Some of the algorithms provided by GuessworkQuantumSideInfo.jl were compared on
a set of 12 test problems. Each problem has *p* ≡ *u*, the uniform
distribution,
for simplicity. The states are chosen as

1.  Two random qubit density matrices

2.  Two random qutrit density matrices

3.  Three pure qubits chosen equidistant within one plane of the Bloch
    sphere (the “Y-states”), i.e.
    cos(j 2pi/3) |0> + sin(j 2pi/3)|1>,  j = 1,2,3

4.  Three random qubit density matrices

5.  Three random qutrit density matrices

6.  The four BB84 states, |0>, |1>, and
    |±> = (|0> ± |1>) / sqrt{2}

as well as the “tensor-2” case of
{*ρ* ⊗ *σ* : *ρ*, *σ* ∈ *S*}
for each of the six sets *S* listed above, corresponding to the
guesswork problem with quantum side information associated to
*ρ*<sub>*XB*</sub><sup> ⊗ 2</sup>, where *ρ*<sub>*XB*</sub> is the
state associated to the original guesswork problem with quantum side
information. The random states were chosen uniformly at random
(i.e. according to the Haar measure).

The exponentially-large SDP formulation (and its dual), the
mixed-integer SDP algorithm, and the active set method were compared,
with several choices of parameters and underlying solvers. The
mixed-integer SDP formulation was evaluated with *M* = *d*<sub>*B*</sub>
(yielding an upper bound), *M* = *d*<sub>*B*</sub><sup>2</sup> + 1
(yielding the optimal value[1]), with the Pajarito mixed-integer SDP
solver , using Convex.jl (version 0.12.7) to formulate the problem.
Pajarito proceeds by solving mixed-integer linear problems (MILP) and
SDPs as subproblems, and thus uses both a MILP solver and an SDP solver
as subcomponents. Pajarito provides two algorithms: an iterative
algorithm, which alternates between solving MILP and SDP subproblems,
and solving a single branch-and-cut problem in which SDP subproblems are
solved via so-called lazy callbacks to add cuts to the mixed-integer
problem. The latter is called “mixed-solver drives” (MSD) in the
Pajarito documentation. We tested three configurations of Pajarito
(version 0.7.0):

* (c1)
Gurobi (version 9.0.3) as the MILP solver and MOSEK (version 8.1.0.82)
as the SDP solver, with Pajarito’s MSD algorithm

* (c2)
Gurobi as the MILP solver and MOSEK as the SDP solver, with Pajarito’s
iterative algorithm, with a relative optimality gap tolerance of 0,

* \(o\)
Cbc (version 2.10.3) as the MILP solver, and SCS (version 2.1.1) as the
SDP solver, with Pajarito’s iterative algorithm

Here, ‘c’ stands for commercial, and ‘o‘ for open-source. In
configuration c1, Gurobi was set to have with a relative optimality gap
tolerance of 10<sup> − 5</sup> and in c2, a relative optimality gap
tolerance of 0. In both configurations, Gurobi was given an absolute
linear-constraint-wise feasibility tolerance of 10<sup> − 8</sup>, and
an integrality tolerance of 10<sup> − 9</sup>. These choices of
parameters match those made in . Cbc was given an integrality tolerance
of 10<sup> − 8</sup>, and SCS’s (normalized) primal, dual residual and
relative gap were set to 10<sup> − 6</sup> for each problem. The default
parameters were used otherwise. Note the MSD option was not used with
Cbc, since the solver does not support lazy callbacks.

For the (exponentially large) SDP primal and dual formulations, the
problems were solved with both MOSEK and SCS, and likewise with the
active-set upper bound.

The active set method uses simulated annealing to iteratively add
violated constraints to the problem to find an upper bound, as described
in , and uses a maximum-time parameter *t*<sub>max</sub> to stop
iterating when the estimated time of finding another constraint to add
would cause the running time of to exceed the maximum-time[2]. This
provides a way to compare the improvement (or lack thereof) of running
the algorithm for more iterations. The algorithm also terminates when a
violated constraint cannot be found after 50 runs of simulated annealing
(started each time with different random initial conditions). Here, the
problems were solved with three choices of *t*<sub>max</sub>, 20 s,
60 s, and 240 s.

The exact answer was not known analytically for most of these problems,
so the average relative error was calculated by comparing to the mean of
the solutions (excluding the active-set method and the MISDP with
*M* = *d*<sub>*B*</sub>, which only give an upper bound in general). In
the BB84 case, in which the solution is known exactly (see ), the
solutions obtained here match the analytic value to a relative tolerance
of at least 10<sup> − 7</sup>.

The problems were run sequentially on a 4-core desktop computer (Intel
i7-6700K 4.00GHz CPU, with 16 GB of RAM, on Ubuntu-20.04 via Windows
Subsystem for Linux, version 2), via the programming language Julia
(version 1.5.1), with a 5 minute time limit.


| Algorithm                                 | Parameters                | average rel. error | average time | number solved | number timed out | number errored out |
|:------------------------------------------|:--------------------------|:-------------------|:-------------|:--------------|:-----------------|:-------------------|
| MISDP (*d*<sub>*B*</sub>)                 | Pajarito (c1)             | 0 %                | 23.74 s      | 6             | 6                | 0                  |
|                                           | Pajarito (c2)             | 0 %                | 24.03 s      | 6             | 6                | 0                  |
|                                           | Pajarito (o)              | 0 %                | 45.05 s      | 6             | 6                | 0                  |
| MISDP (*d*<sub>*B*</sub><sup>2</sup> + 1) | Pajarito (c1)             | 0 %                | 27.99 s      | 4             | 8                | 0                  |
|                                           | Pajarito (c2)             | 0 %                | 29.98 s      | 4             | 8                | 0                  |
|                                           | Pajarito (o)              | 0 %                | 24.71 s      | 1             | 11               | 0                  |
| SDP                                       | MOSEK                     | 0 %                | 8.99 s       | 8             | 3                | 1                  |
|                                           | SCS                       | 0 %                | 9.08 s       | 8             | 3                | 1                  |
| SDP (dual)                                | MOSEK                     | 0 %                | 8.74 s       | 8             | 3                | 1                  |
|                                           | SCS                       | 0 %                | 8.59 s       | 8             | 3                | 1                  |
| Active set upper bound (MOSEK)            | *t*<sub>max</sub> = 20 s  | 6.80 %             | 16.08 s      | 10            | 2                | 0                  |
|                                           | *t*<sub>max</sub> = 60 s  | 6.79 %             | 19.03 s      | 10            | 2                | 0                  |
|                                           | *t*<sub>max</sub> = 240 s | 6.80 %             | 26.17 s      | 10            | 2                | 0                  |
| Active set upper bound (SCS)              | *t*<sub>max</sub> = 20 s  | 6.09 %             | 33.24 s      | 12            | 0                | 0                  |
|                                           | *t*<sub>max</sub> = 60 s  | 6.09 %             | 35.78 s      | 12            | 0                | 0                  |
|                                           | *t*<sub>max</sub> = 240 s | 6.09 %             | 34.30 s      | 11            | 1                | 0                  |


Caption: Comparison of average relative error and average solve time for the 12
problems discussed above. A problem is considered “timed out” if an
answer is not obtained in 5 minutes, and “errored out” if the solution
was not obtained due to errors (such as running out of RAM). The average
relative error, which was rounded to two decimal digits, and the time
taken are calculated only over the problems which were solved by the
given algorithm and choice of parameters. “MISDP (*d*<sub>*B*</sub>)”
refers to the choice *M* = *d*<sub>*B*</sub>, and likewise “MISDP
(*d*<sub>*B*</sub><sup>2</sup> + 1)” refers to the choice
*M* = *d*<sub>*B*</sub><sup>2</sup> + 1.


| Algorithm                                 | Parameters                | Two random qubits                 | Two random qutrits                | Y-states                    |
|:------------------------------------------|:--------------------------|:----------------------------------|:----------------------------------|:----------------------------|
| MISDP (*d*<sub>*B*</sub>)                 | Pajarito (c1)             | 23.63 s, timeout                  | 23.60 s, timeout                  | 23.56 s, timeout            |
|                                           | Pajarito (c2)             | 22.99 s, timeout                  | 23.31 s, timeout                  | 23.21 s, timeout            |
|                                           | Pajarito (o)              | 23.47 s, timeout                  | 24.77 s, timeout                  | 26.15 s, timeout            |
| MISDP (*d*<sub>*B*</sub><sup>2</sup> + 1) | Pajarito (c1)             | 23.35 s, timeout                  | 33.87 s, timeout                  | 27.68 s, timeout            |
|                                           | Pajarito (c2)             | 22.85 s, timeout                  | 30.37 s, timeout                  | 32.93 s, timeout            |
|                                           | Pajarito (o)              | 24.71 s, timeout                  | timeout, timeout                  | timeout, timeout            |
| SDP                                       | MOSEK                     | 8.69 s, 8.84 s                    | 8.78 s, 9.23 s                    | 9.33 s, timeout             |
|                                           | SCS                       | 9.00 s, 8.90 s                    | 8.44 s, 11.22 s                   | 8.98 s, timeout             |
| SDP (dual)                                | MOSEK                     | 8.46 s, 8.63 s                    | 8.54 s, 8.83 s                    | 9.11 s, timeout             |
|                                           | SCS                       | 8.76 s, 8.32 s                    | 8.33 s, 9.20 s                    | 8.74 s, timeout             |
| Active set upper bound (MOSEK)            | *t*<sub>max</sub> = 20 s  | 8.76 s (19.5 %), 10.41 s (1.5 %)  | 8.89 s (19.5 %), timeout          | 9.72 s (0 %), 34.25 s (? %) |
|                                           | *t*<sub>max</sub> = 60 s  | 10.91 s (19.5 %), 10.41 s (1.9 %) | 8.87 s (19.5 %), timeout          | 9.66 s (0 %), 31.00 s (? %) |
|                                           | *t*<sub>max</sub> = 240 s | 9.47 s (19.5 %), 10.40 s (1.5 %)  | 8.90 s (19.5 %), timeout          | 9.81 s (0 %), 30.26 s (? %) |
| Active set upper bound (SCS)              | *t*<sub>max</sub> = 20 s  | 9.04 s (19.5 %), 10.92 s (1.5 %)  | 8.70 s (19.5 %), 101.06 s (1.1 %) | 9.23 s (0 %), 82.84 s (? %) |
|                                           | *t*<sub>max</sub> = 60 s  | 9.07 s (19.5 %), 10.22 s (1.9 %)  | 8.66 s (19.5 %), 32.94 s (1.2 %)  | 9.18 s (0 %), 50.79 s (? %) |
|                                           | *t*<sub>max</sub> = 240 s | 9.04 s (19.5 %), 10.02 s (1.9 %)  | 8.79 s (19.5 %), 22.69 s (1.2 %)  | 9.31 s (0 %), 37.36 s (? %) |

Caption: The individual
timings for each algorithm and choice of settings on problems (1)–(3),
and the corresponding “tensor-2” problems discussed at
<a href="#eq:tensor-2" data-reference-type="eqref" data-reference="eq:tensor-2">[eq:tensor-2]</a>.
For each algorithm, the running time of the original problem is given
followed by the running time on the “tensor-2” problem, e.g. the SDP
formulation with MOSEK on the two random qubits problem was solved in
8.69 seconds, and in 8.84 seconds for the corresponding tensor-2
problem. “timeout” is written whenever the problem was not solved within
5 minutes. For the active set algorithms, the relative error is also
given for each problem in parenthesis. Note that the MISDP formulation
with *M* = *d*<sub>*B*</sub> is also only known to be an upper bound,
but a relative error of less than 10<sup> − 5</sup> in each instance, so
the relative errors are omitted. Lastly, the relative error is written
as ? % in the case that only an upper bound was obtained.


| Algorithm                                 | Parameters                | Three random qubits           | Three random qutrits           | BB84 states                  |
|:------------------------------------------|:--------------------------|:------------------------------|:-------------------------------|:-----------------------------|
| MISDP (*d*<sub>*B*</sub>)                 | Pajarito (c1)             | 23.63 s, timeout              | 24.49 s, timeout               | 23.51 s, timeout             |
|                                           | Pajarito (c2)             | 23.17 s, timeout              | 26.50 s, timeout               | 25.01 s, timeout             |
|                                           | Pajarito (o)              | 27.11 s, timeout              | 95.01 s, timeout               | 73.80 s, timeout             |
| MISDP (*d*<sub>*B*</sub><sup>2</sup> + 1) | Pajarito (c1)             | 27.08 s, timeout              | timeout, timeout               | timeout, timeout             |
|                                           | Pajarito (c2)             | 33.76 s, timeout              | timeout, timeout               | timeout, timeout             |
|                                           | Pajarito (o)              | timeout, timeout              | timeout, timeout               | timeout, timeout             |
| SDP                                       | MOSEK                     | 9.35 s, timeout               | 8.82 s, timeout                | 8.87 s, error                |
|                                           | SCS                       | 9.09 s, timeout               | 8.46 s, timeout                | 8.54 s, error                |
| SDP (dual)                                | MOSEK                     | 9.14 s, timeout               | 8.55 s, timeout                | 8.62 s, error                |
|                                           | SCS                       | 8.86 s, timeout               | 8.25 s, timeout                | 8.23 s, error                |
| Active set upper bound (MOSEK)            | *t*<sub>max</sub> = 20 s  | 9.73 s (1.3 %), 32.29 s (? %) | 9.50 s (5.8 %), timeout        | 9.51 s (0 %), 27.72 s (? %)  |
|                                           | *t*<sub>max</sub> = 60 s  | 9.69 s (1.3 %), 24.12 s (? %) | 9.45 s (5.8 %), timeout        | 9.77 s (0 %), 66.47 s (? %)  |
|                                           | *t*<sub>max</sub> = 240 s | 9.69 s (1.3 %), 30.42 s (? %) | 9.45 s (5.8 %), timeout        | 9.49 s (0 %), 133.81 s (? %) |
| Active set upper bound (SCS)              | *t*<sub>max</sub> = 20 s  | 9.44 s (1.3 %), 35.66 s (? %) | 9.67 s (5.8 %), 84.43 s (? %)  | 9.09 s (0 %), 28.76 s (? %)  |
|                                           | *t*<sub>max</sub> = 60 s  | 9.44 s (1.3 %), 54.60 s (? %) | 9.11 s (5.8 %), 155.51 s (? %) | 9.29 s (0 %), 70.57 s (? %)  |
|                                           | *t*<sub>max</sub> = 240 s | 8.89 s (1.3 %), 75.43 s (? %) | 9.01 s (5.8 %), timeout        | 9.15 s (0 %), 177.64 s (? %) |

Caption: The
individual timings for each algorithm and choice of settings on problems
(4)–(6). See the previous caption for a description of the quantities shown. Here, “error”
means the solution was not obtained due to an error (such as running out
of memory).


[1] although *M* = *d*<sub>*B*</sub><sup>2</sup> suffices

[2] The maximum time can still be exceeded, since at least one iteration
must be performed and the estimate can be wrong.
