# Computing the guesswork to high precision

In order to solve the SDP describing the guesswork in high-precision, we need to
start a new Julia session. That's because Pajarito, the library used for
mixed-integer SDPs, and SDPAFamily, the one for high-precision SDPs, each
require different versions of Convex.jl, the optimization problem modeling
library. The following example was run locally and the output copied to the
documentation. First, add the packages `Convex#master` (or just `Convex` once
version 0.13 of that package is released) and `SDPAFamily`, and then the
following code.

```julia
julia> setprecision(2000) # set Julia's global BigFloat precision to 2000
2000

julia> using SDPAFamily

julia> opt = () -> SDPAFamily.Optimizer{BigFloat}(
                  presolve = true,
                  params = (  epsilonStar = 1e-200, # constraint tolerance
                              epsilonDash = 1e-200, # normalized duality gap tolerance
                              precision = 2000 # arithmetric precision used in sdpa_gmp
                  ))
#186 (generic function with 1 method)

julia> using GuessworkQuantumSideInfo

julia> T = BigFloat
BigFloat

julia> ρBs = BB84_states(T);

julia> p = ones(T, 4) / 4;

julia> @time output = guesswork(p, ρBs; solver = opt());
 54.613852 seconds (274.09 M allocations: 45.610 GiB, 11.66% gc time)

julia> output.optval
1.709430584957905167000276613891820366570111215168695793285623786801851390340190444663937972905174926203163178961789986212862785992386529996327370100824961524163769585705185014835212461631471665128968986016501876699676856588609960582022565322653047114497843315997252226645378373262132182609166891615169945992530274280324399117232937277795982220635506452810752194823768763057910726875757516626180726385923719763995534231714003266054518160879579903803264241437877679215965923661443029759736849138449576021864074135403089512757915961340265964663906514782565168514016103734338487088415453174248635495108648696

julia> true_val = (big(1) / big(4)) * (10 - sqrt(big(10)))
1.709430584957905167000276613891820366570111215168695793285623786801851390340190444663937972905174926203163178961789986212862785992386529996327370100824961524163769585705185014835212461631471665128968978670767605727382105884738907685530234490547847514542907072112898417861750220350793858949168259598827171415967762552027433367614096830955304662484652953430354884079216261880492214057666332006545121703580399192189551149445422988531788846349936605605074292838992743898005525991364002709162694336846983051837080992237890825551560768461069874907526272632644537498160412788839029317222677722739412986441086645

julia> abs(output.optval - true_val) <= 10.0^(-200)
true

```

Note that the output of the optimization solver matches `true_val` up to an
error of at most $10^{-200}$.

See also
<https://github.com/ericphanson/GuessworkQuantumSideInfo.jl/tree/master/test/high_precision_tests>
for a folder with a reproducible environment for running this kind of
high-precision code.
