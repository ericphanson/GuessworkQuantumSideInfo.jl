var documenterSearchIndex = {"docs":
[{"location":"examples.html#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"Let's load the package, an SDP solver, and define a simple plotting routine.","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"using GuessworkQuantumSideInfo\nusing SCS, Plots\nget_sdp_solver() = SCSSolver(verbose=false)\nplot_pmfN(data) = bar(pmfN(data); xlabel=\"Guess number\", ylabel=\"Probability of guessing correctly\", legend = false)","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"Next, we define some basic qubit states.","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"dB = 2\nketplus = (ket(1, dB) + ket(2,dB))/sqrt(2)\nketminus = (ket(1, dB) - ket(2,dB))/sqrt(2)\nketzero = ket(1, dB)\nketone = ket(2, dB)","category":"page"},{"location":"examples.html#Example-1:-A-warmup-with-trivial-examples-1","page":"Examples","title":"Example 1: A warmup with trivial examples","text":"","category":"section"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"Let's consider the case with J=2 and both states are the same. The side information is therefore completely uninformative, and the guesswork is","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"G(XB) = frac12cdot 1 + frac12 cdot 2 = 15","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"We can check this:","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"p = [0.5, 0.5]\nρBs = dm.([ ketzero, ketzero  ])\noutput = guesswork(p, ρBs; solver = get_sdp_solver());\noutput.optval","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"We see the result agrees with 1.5, as we expected. Likewise, if we choose the two states as $ |0\\rangle, |1\\rangle$, we can get it in one guess every time, of course, since they are orthogonal:","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"p = [0.5, 0.5]\nρBs = dm.([ ketzero, ketone  ])\noutput = guesswork(p, ρBs; solver = get_sdp_solver());\noutput.optval","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"We can inspect the POVMs:","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"output.Es","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"As we would expect, we (approximately) obtain the projection onto 0 rangle and the projection onto 1 rangle.","category":"page"},{"location":"examples.html#Example-2:-the-BB84-states-1","page":"Examples","title":"Example 2: the BB84 states","text":"","category":"section"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"Let's consider the four states + rangle -rangle 0rangle 1rangle.","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"p = [0.25, 0.25, 0.25, 0.25]\nρBs = dm.([ ketplus, ketminus, ketzero, ketone  ])\noutput = guesswork(p, ρBs; solver = get_sdp_solver());\noutput.optval","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"plot_pmfN(output)","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"Let's try the same example, but imposing a steep cost for the fourth guess.","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"c = [1.0, 2.0, 3.0, 5000.0]\noutput = guesswork(p, ρBs; c = c, solver = get_sdp_solver());\noutput.optval","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"We see that the average number of guesses to get a correct answer has gone up. However, inspecting the probability mass function for the number of guesses under the optimal strategy","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"plot_pmfN(output)","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"we see that the probability that the probability of guessing correctly on the fourth guess goes to almost zero.","category":"page"},{"location":"examples.html#Example-3:-two-copies-of-the-BB84-states-1","page":"Examples","title":"Example 3: two copies of the BB84 states","text":"","category":"section"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"Let us consider two tensor copies of the BB84 states:","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"p = ones(16)/16\nρBs = iid_copies(BB84_states(), 2)","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"In this case, there are 16 = 20922789888000 possible guessing orders, and hence 16 variables in the primal formulation of the SDP, or 16+1 constraints in the dual form of the SDP. In either case, we can't even fit them all into our computer's memory. Instead, we resort to bounds:","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"lb_output = guesswork_lower_bound(p, ρBs, solver = get_sdp_solver());\nlb_output.optval\n\nub_output = guesswork_upper_bound(p, ρBs; max_time = 30, make_solver = get_sdp_solver);\nub_output.optval","category":"page"},{"location":"examples.html#A-closer-look-at-guesswork_upper_bound-1","page":"Examples","title":"A closer look at guesswork_upper_bound","text":"","category":"section"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"We can understand guesswork_upper_bound better by setting verbose=true. The algorithm computes a sequence of upper bounds by relaxing the dual problem. First, it removes all constraints, and chooses the dual variable Y as the identity matrix. Then it uses a simulated annealing algorithm to heuristically minimize lambda_textmin(R_vec g - Y) over the possible guessing orders vec g (which are permutations), to find a constraint that is \"maximally\" violated by this choice of Y. Then we add the corresponding constraint Y leq R_vec g to the dual problem and solve it again. This is repeated until either a fixed number of constraints is added, some number of simulated annealing runs fails to find another violated constraint, or a time limit is reached. In the following, we set a time limit of 30 seconds.","category":"page"},{"location":"examples.html#","page":"Examples","title":"Examples","text":"ub_output = guesswork_upper_bound(p, ρBs; verbose=true, max_time = 30, make_solver = get_sdp_solver);","category":"page"},{"location":"high-precision-example.html#Computing-the-guesswork-to-high-precision-1","page":"High precision example","title":"Computing the guesswork to high precision","text":"","category":"section"},{"location":"high-precision-example.html#","page":"High precision example","title":"High precision example","text":"In order to solve the SDP describing the guesswork in high-precision, we need to start a new Julia session. That's because Pajarito, the library used for mixed-integer SDPs, and SDPAFamily, the one for high-precision SDPs, each require different versions of Convex.jl, the optimization problem modeling library. The following example was run locally and the output copied to the documentation. First, add the packages Convex#master (or just Convex once version 0.13 of that package is released) and SDPAFamily, and then the following code.","category":"page"},{"location":"high-precision-example.html#","page":"High precision example","title":"High precision example","text":"julia> setprecision(2000) # set Julia's global BigFloat precision to 2000\n2000\n\njulia> using SDPAFamily\n\njulia> opt = () -> SDPAFamily.Optimizer{BigFloat}(\n                  presolve = true,\n                  params = (  epsilonStar = 1e-200, # constraint tolerance\n                              epsilonDash = 1e-200, # normalized duality gap tolerance\n                              precision = 2000 # arithmetric precision used in sdpa_gmp\n                  ))\n#186 (generic function with 1 method)\n\njulia> using GuessworkQuantumSideInfo\n\njulia> T = BigFloat\nBigFloat\n\njulia> ρBs = BB84_states(T);\n\njulia> p = ones(T, 4) / 4;\n\njulia> @time output = guesswork(p, ρBs; solver = opt());\n 54.613852 seconds (274.09 M allocations: 45.610 GiB, 11.66% gc time)\n\njulia> output.optval\n1.709430584957905167000276613891820366570111215168695793285623786801851390340190444663937972905174926203163178961789986212862785992386529996327370100824961524163769585705185014835212461631471665128968986016501876699676856588609960582022565322653047114497843315997252226645378373262132182609166891615169945992530274280324399117232937277795982220635506452810752194823768763057910726875757516626180726385923719763995534231714003266054518160879579903803264241437877679215965923661443029759736849138449576021864074135403089512757915961340265964663906514782565168514016103734338487088415453174248635495108648696\n\njulia> true_val = (big(1) / big(4)) * (10 - sqrt(big(10)))\n1.709430584957905167000276613891820366570111215168695793285623786801851390340190444663937972905174926203163178961789986212862785992386529996327370100824961524163769585705185014835212461631471665128968978670767605727382105884738907685530234490547847514542907072112898417861750220350793858949168259598827171415967762552027433367614096830955304662484652953430354884079216261880492214057666332006545121703580399192189551149445422988531788846349936605605074292838992743898005525991364002709162694336846983051837080992237890825551560768461069874907526272632644537498160412788839029317222677722739412986441086645\n\njulia> abs(output.optval - true_val) <= 10.0^(-200)\ntrue\n","category":"page"},{"location":"high-precision-example.html#","page":"High precision example","title":"High precision example","text":"Note that the output of the optimization solver matches true_val up to an error of at most 10^-200.","category":"page"},{"location":"high-precision-example.html#","page":"High precision example","title":"High precision example","text":"See also https://github.com/ericphanson/GuessworkQuantumSideInfo.jl/tree/master/test/high_precision_tests for a folder with a reproducible environment for running this kind of high-precision code.","category":"page"},{"location":"mixed-integer-SDP.html#Using-a-mixed-integer-SDP-to-find-extremal-strategies-1","page":"Mixed-integer SDP example","title":"Using a mixed-integer SDP to find extremal strategies","text":"","category":"section"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"Let us revisit Example 2: the BB84 states using a mixed-integer SDP. This allows us to specify the number of non-zero POVM elements in the measurement– or, equivalently, the number of possible guessing orders the final strategy chooses among.","category":"page"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"We will use the Pajarito.jl solver, which solves mixed-integer SDPs by solving an alternating sequence of mixed-integer linear programs and SDPs; see the Pajarito paper for information on how Pajarito works. Pajarito requires both an SDP solver and a mixed-integer linear solver; we will use the open source solvers SCS and Cbc, respectively.","category":"page"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"using GuessworkQuantumSideInfo\nusing Pajarito, Cbc, SCS # solvers\nfunction misdp_solver(; verbose = false)\n    sdp_solver = SCSSolver(verbose=0)\n    mip_solver = CbcSolver(loglevel = 0)\n    PajaritoSolver(\n        cont_solver = sdp_solver,\n        mip_solver = mip_solver,\n        mip_solver_drives = false,\n        use_mip_starts = true,\n        solve_relax = false,\n        log_level = verbose ? 3 : 0,\n    )\nend\n\np = ones(4)/4\nρBs = BB84_states()\noutput = guesswork_MISDP(p, ρBs, 2; verbose=false, solver = misdp_solver());\noutput.optval\noutput.Es","category":"page"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"We see that with only two measurement outcomes we recover the same optimal value as the case without a constraint on  the number of measurement outcomes (in Example 2: the BB84 states). We've thus found an extremal strategy (in that the POVM associated to this strategy cannot be written as the convex combination of two other POVMs). Often, the solutions returned by guesswork do not return extremal POVMs, although the details depend on the SDP solver used.","category":"page"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"We can see what the associated guessing orders are:","category":"page"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"output.povm_outcomes","category":"page"},{"location":"mixed-integer-SDP.html#","page":"Mixed-integer SDP example","title":"Mixed-integer SDP example","text":"To summarize, one optimal strategy for the case of p uniform, and ρBs given by the four BB84 states, is to perform a projective measurement whose operators are given by output.Es above, and then make guesses in one of the two orders given by output.povm_outcomes (depending on which measurement outcome was obtained).","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"CurrentModule = GuessworkQuantumSideInfo","category":"page"},{"location":"index.html#GuessworkQuantumSideInfo-1","page":"Home","title":"GuessworkQuantumSideInfo","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"This is a package accompanying the preprint Guesswork with Quantum Side Information.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"See the Examples for some examples (or just below for a quick example), Computing the guesswork to high precision for an example solving a problem with a high-precision SDP solver, Using a mixed-integer SDP to find extremal strategies for an example using a mixed-integer SDP, or below for the documentation of the functions provided by this package.","category":"page"},{"location":"index.html#Quick-example-1","page":"Home","title":"Quick example","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"Consider one party Alice who draws a random number in the set [1,2,3,4] uniformly at random. If she draws 1 she sends another party, Bob, the quantum state |0⟩; if she draws 2, she sends |1⟩, if she draws 3 she sends |-⟩, and finally if she draws 4, she sends |+⟩. Bob, knowing this general procedure but not which number Alice drew, aims to guess the value Alice drew by performing experiments on the quantum state he was given. The average number of guesses Bob needs in order to get the right answer, minimized over all quantum strategies, is the so-called guesswork with quantum side information. This package provides a means to compute this.","category":"page"},{"location":"index.html#","page":"Home","title":"Home","text":"using GuessworkQuantumSideInfo, SCS\np = [0.25, 0.25, 0.25, 0.25];\nketzero = ket(1, 2);\nketone = ket(2, 2);\nketplus = (ket(1, 2) + ket(2,2))/sqrt(2);\nketminus = (ket(1, 2) - ket(2,2))/sqrt(2);\nρBs = dm.([ ketplus, ketminus, ketzero, ketone  ])\noutput = guesswork(p, ρBs; solver = SCSSolver(verbose=false));\noutput.optval","category":"page"},{"location":"index.html#Guesswork-functions-1","page":"Home","title":"Guesswork functions","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"GuessworkQuantumSideInfo.guesswork\nGuessworkQuantumSideInfo.guesswork_lower_bound\nGuessworkQuantumSideInfo.guesswork_upper_bound","category":"page"},{"location":"index.html#GuessworkQuantumSideInfo.guesswork","page":"Home","title":"GuessworkQuantumSideInfo.guesswork","text":"guesswork(\n    p::AbstractVector{T},\n    ρBs::AbstractVector{<:AbstractMatrix};\n    solver,\n    K::Integer = length(p),\n    c = T[1:K..., 5_000],\n    dual::Bool = false,\n    remove_repetition::Bool = true,\n    povm_outcomes = make_povm_outcomes(length(p), K, remove_repetition),\n    verbose::Bool = true,\n)\n\nComputes the guesswork for the c-q state specified by a probability vector p, giving the distribution X, and ρBs, giving the associated quantum states.\n\nThe keyword arguments are as follows:\n\nsolver is the only required keyword argument; an SDP solver such as SCS or MOSEK must be passed.\nK corresponds to the maximum number of allowed guesses. The number of variables in the primal SDP (and the number of constraints in the dual SDP) scales as length(p)^K.\nc may be given a custom cost vector. If K < length(p), then c should be of length K+1. The last entry, c[K+1], corresponds to the cost of not guessing the correct answer within K guesses.\ndual is a boolean variable indicating whether the primal or dual optimization problem should be solved.\nremove_repetition is a boolean variable defaulting to true, indicating whether repeated guesses of the same value should be removed; as long as c is increasing, this decreases the size of the SDP without affecting the optimal value.\npovm_outcomes should be an iterator (or vector) corresponding to the possible guessing orders. This defaults to all subsets of length K of 1:length(p) without repetition.\nverbose is a boolean which indicates if warnings should be printed when the problem is not solved optimally.\n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.guesswork_lower_bound","page":"Home","title":"GuessworkQuantumSideInfo.guesswork_lower_bound","text":"guesswork_lower_bound(\n    p::AbstractVector{T},\n    ρBs::AbstractVector{<:AbstractMatrix};\n    solver,\n    c = T[1:length(p)..., 10_000],\n    verbose::Bool = false,\n)\n\nSee guesswork for the meaning of the arguments. Computes a lower bound to the optimal expected number of guesses by solving a relaxed version of the primal SDP. For J states, only needs J^2 PSD variables subject to two linear constraints.\n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.guesswork_upper_bound","page":"Home","title":"GuessworkQuantumSideInfo.guesswork_upper_bound","text":"guesswork_upper_bound(\n    p::AbstractVector{T},\n    ρBs::AbstractVector{<:AbstractMatrix};\n    make_solver,\n    c::AbstractVector = T.(1:length(p)),\n    max_retries = 50,\n    max_time = Inf,\n    num_constraints = Inf,\n    verbose::Bool = false,\n    num_steps_per_SA_run::Integer = length(p)^2 * 500,\n) where {T<:Number} -> NamedTuple\n\nComputes an upper bound to the guesswork problem associated to the c-q state specified by p and ρBs, as in guesswork. A custom cost vector c may be optionally passed. If the keyword argument verbose is set to true, information is printed about each iteration of the algorithm.\n\nThe keyword argument make_solver is required, and must pass a function that creates a solver instances. For example, instead of passing SCSSolver(), pass () -> SCSSolver(). This is needed because the algorithm used in guesswork_upper_bound solves a sequence of SDPs, not just one.\n\nThe algorithm has three termination criteria which are controlled by keyword arguments. The algorithm stops when any of the following occur:\n\nmax_retries simulated annealing attempts fail to find a violated constraint.\nnum_constraints constraints have been added to the dual SDP\nThe total runtime of the algorithm is projected to exceed max_time on the next iteration.\n\nBy default, max_retries is set to 50, while num_constraints and max_time are set to infinity.\n\nLastly, the keyword argument num_steps_per_SA_run controls the runtime of the simulated annealing algorithm. Increase num_steps_per_SA_run to search longer for a violated constraint within a given simulated annealing run.\n\n\n\n\n\n","category":"function"},{"location":"index.html#Quantum-states-1","page":"Home","title":"Quantum states","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"GuessworkQuantumSideInfo.ket\nGuessworkQuantumSideInfo.bra\nGuessworkQuantumSideInfo.BB84_states\nGuessworkQuantumSideInfo.iid_copies\nGuessworkQuantumSideInfo.randdm\nGuessworkQuantumSideInfo.randprobvec","category":"page"},{"location":"index.html#GuessworkQuantumSideInfo.ket","page":"Home","title":"GuessworkQuantumSideInfo.ket","text":"ket([T = Float64], i::Integer, d::Integer) -> SparseVector{Complex{T}}\n\nCreate a vector representing the ith computational basis vector in dimension d.\n\nExample\n\njulia> ket(1,2)\n2-element SparseVector{Complex{Float64},Int64} with 1 stored entry:\n  [1]  =  1.0+0.0im\n\njulia> collect(ans)\n2-element Array{Complex{Float64},1}:\n 1.0 + 0.0im\n 0.0 + 0.0im\n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.bra","page":"Home","title":"GuessworkQuantumSideInfo.bra","text":"bra([T = Float64], i::Integer, d::Integer) -> SparseVector{Complex{T}}'\n\nCreate a dual vector representing the bra associated to ith computational basis vector in dimension d.\n\nExample\n\njulia> bra(1,2)\n1×2 LinearAlgebra.Adjoint{Complex{Float64},SparseVector{Complex{Float64},Int64}}:\n 1.0-0.0im  0.0-0.0im\n\njulia> collect(ans)\n1×2 Array{Complex{Float64},2}:\n 1.0-0.0im  0.0-0.0im\n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.BB84_states","page":"Home","title":"GuessworkQuantumSideInfo.BB84_states","text":"BB84_states([T::Type = Float64])\n\nGenerates the BB84 states |0⟩, |1⟩, |-⟩, and |+⟩, for use in guesswork or other functions. The numeric type can be optionally specified by the argument.\n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.iid_copies","page":"Home","title":"GuessworkQuantumSideInfo.iid_copies","text":"iid_copies(ρBs::AbstractVector{<:AbstractMatrix}, n::Integer) -> Vector{Matrix}\n\nCreate a vector of all states of the form ρ_1 otimes dotsm otimes ρ_n where the ρ_i range over the set ρBs.\n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.randdm","page":"Home","title":"GuessworkQuantumSideInfo.randdm","text":"randdm([T = Float64], d)\n\nGenerates a density matrix with numeric type Complex{T}, of dimension d at random.\n\nExample\n\njulia> randdm(2)\n2×2 Array{Complex{Float64},2}:\n 0.477118+0.0im        0.119848-0.0371569im\n 0.119848+0.0371569im  0.522882+0.0im      \n\n\n\n\n\n","category":"function"},{"location":"index.html#GuessworkQuantumSideInfo.randprobvec","page":"Home","title":"GuessworkQuantumSideInfo.randprobvec","text":"randprobvec([T=Float64], d)\n\nGenerates points of type T, uniformly at random on the standard d-1 dimensional simplex using an algorithm by Smith and Tromble.\n\nExample\n\njulia> randprobvec(3)\n3-element Array{Float64,1}:\n 0.24815974900033688\n 0.17199716455672287\n 0.5798430864429402 \n\n\n\n\n\n\n","category":"function"},{"location":"index.html#Utilities-1","page":"Home","title":"Utilities","text":"","category":"section"},{"location":"index.html#","page":"Home","title":"Home","text":"GuessworkQuantumSideInfo.pmfN","category":"page"},{"location":"index.html#GuessworkQuantumSideInfo.pmfN","page":"Home","title":"GuessworkQuantumSideInfo.pmfN","text":"pmfN(data; tol = 1e-5) -> Vector\n\nCompute the probability mass function for the number of guesses N, given a strategy. The nth entry of the output vector vec gives the probability for guessing the correct answer on the nth try, for n = 1 : K. The last entry gives the probability of never guessing the correct answer.\n\ndata should be a NamedTuple with entries for p, ρBs, Es, K, and povm_outcomes\ntol provides a tolerance above which to warn about imaginary or negative probabilities.\n\n\n\n\n\n","category":"function"}]
}
