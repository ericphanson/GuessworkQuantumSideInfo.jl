using GuessworkQuantumSideInfo
using Documenter

previous_GKSwstype = get(ENV, "GKSwstype", "")
ENV["GKSwstype"] = "100"

makedocs(;
    modules=[GuessworkQuantumSideInfo],
    authors="Eric P. Hanson",
    repo="https://github.com/ericphanson/GuessworkQuantumSideInfo.jl/blob/{commit}{path}#L{line}",
    sitename="GuessworkQuantumSideInfo.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ericphanson.github.io/GuessworkQuantumSideInfo.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ericphanson/GuessworkQuantumSideInfo.jl",
)

ENV["GKSwstype"] = previous_GKSwstype;
