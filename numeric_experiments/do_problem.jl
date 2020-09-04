# https://github.com/JuliaLang/julia/issues/28679
function _redirect_stdout(dofunc::Function)
    nullfile = @static Sys.iswindows() ? "nul" : "/dev/null"
    open(nullfile, "w") do io
        redirect_stdout(dofunc, io)
    end
end

_redirect_stdout() do
    include(joinpath(@__DIR__, "common.jl"))
end

problem_idx = parse(Int, ARGS[1])
algo_idx = parse(Int, ARGS[2])

algo = algos[algo_idx]
prob = problems[problem_idx]


optval, elapsed = algo.f(prob)
println(optval)
println(elapsed)
