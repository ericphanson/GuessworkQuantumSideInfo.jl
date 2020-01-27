# High precision tests

These tests use the high-precision SDP solver SDPA-GMP, via
[SDPAFamily.jl](https://github.com/ericphanson/SDPAFamily.jl). This currently
requires the `master` branch of Convex.jl which is incompatible with the
Pajarito solver used for mixed-integer SDPs (until that solver is updated to use
MathOptInterface).

To run the high-precision tests, call
```
julia --project=. -e 'include("tests.jl")'
```
from the command line, with the working directory set to this folder, or enter a
Julia session with
```
juila --project=.
```
and then include the tests file via `include("tests.jl")`.

The `Manifest.toml` file is included in this folder so that a set of compatible
versions of the necessary packages can be reproduced here.
