using ProgressMeter, CSV
function write_results(algo, prob, optval = NaN, elapsed_seconds = NaN; warnings = false, errors = false)
    row = (;algo.algo, algo.settings, prob.numeric_type, prob.problem, optval, elapsed_seconds, warnings, errors)
    CSV.write(joinpath(@__DIR__, "results.csv"), [row], append = isfile(joinpath(@__DIR__, "results.csv")))
end

isfile(joinpath(@__DIR__, "results.csv")) && error("`results.csv` already exists")

include("common.jl")

n_problems = length(problems)
n_algos = length(algos)
timeout = 120
for problem_idx in 1:n_problems, algo_idx in 1:n_algos
    algo = algos[algo_idx]
    prob = problems[problem_idx]
    out = Pipe()
    err = Pipe()
    proc = run(pipeline(`julia --project=$(@__DIR__) do_problem.jl $(algo_idx) $(problem_idx)`, stdout=out, stderr=err), wait=false)
    time_elapsed = 0.0
    meter = Progress(timeout*2, desc = "Problem $problem_idx/$n_problems, algorithm $algo_idx/$n_algos ")
    while time_elapsed < timeout && Base.process_running(proc)
        sleep(0.5)
        next!(meter)
        time_elapsed += 0.5
    end
    finish!(meter)
    if Base.process_running(proc)
        println("Timed out!")
        Base.kill(proc)
       
        write_results(algo, prob)
    else
        println("Finished!")
        close(err.in)
        close(out.in)
        errors_text = String(read(err))
        if !isempty(errors_text)
            @error errors_text
        end
        results = split(chomp(String(read(out))), '\n')
        warnings = false
        if length(results) > 2
            @warn join(results[1:end-2], "\n")
            warnings = true
        end
        @info results
        optval, elapsed = parse.(Float64, results[end-1:end])
        write_results(algo, prob, optval, elapsed; errors = !isempty(errors_text), warnings)
    end

end

df = CSV.read(joinpath(@__DIR__, "results.csv"))
