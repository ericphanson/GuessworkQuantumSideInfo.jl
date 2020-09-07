using ProgressMeter, CSV
function write_results(algo, prob, optval = NaN, elapsed_seconds = NaN; warnings = false, errors = false)
    row = (;algo.algo, algo.settings, prob.numeric_type, prob.problem, optval, elapsed_seconds, warnings, errors)
    CSV.write(joinpath(@__DIR__, "results.csv"), [row], append = isfile(joinpath(@__DIR__, "results.csv")))
end

isfile(joinpath(@__DIR__, "results.csv")) && error("`results.csv` already exists")

include("common.jl")

# The idea is these can be run in parallel, and the thunks that write to the CSV returned,
# then those could be executed serially.
# It's probably more robust to have each process write a separate one-line CSV
# and then collect them at the end, though.
function run_problem(algo_idx, problem_idx; verbose = true)
    algo = algos[algo_idx]
    prob = problems[problem_idx]
    out = Pipe()
    err = Pipe()
    proc = run(pipeline(`julia --project=$(@__DIR__) do_problem.jl $(algo_idx) $(problem_idx)`, stdout=out, stderr=err), wait=false)
    time_elapsed = 0.0
    if verbose
        meter = Progress(timeout*2, desc = "Problem $problem_idx/$n_problems, algorithm $algo_idx/$n_algos ")
    end
    while time_elapsed < timeout && Base.process_running(proc)
        sleep(0.5)
        verbose && next!(meter)
        time_elapsed += 0.5
    end
    verbose && finish!(meter)
    if Base.process_running(proc)
        verbose && println("Timed out!")
        Base.kill(proc)
       
        f = () -> write_results(algo, prob)
    else
        verbose && println("Finished!")
        close(err.in)
        close(out.in)
        errors_text = String(read(err))
        if !isempty(errors_text) && verbose
            @error errors_text
        end
        results = split(chomp(String(read(out))), '\n')
        warnings = false
        if length(results) > 2
            verbose && @warn join(results[1:end-2], "\n")
            warnings = true
        end
        verbose && @info results
        if length(results) < 2
            verbose && @error "Not enough returns, something went wrong"
            errors = true
            optval, elapsed = NaN, NaN
        else
            optval, elapsed = parse.(Float64, results[end-1:end])
        end
        f = () -> write_results(algo, prob, optval, elapsed; errors = !isempty(errors_text), warnings)
    end
    return f
end

n_problems = length(problems)
n_algos = length(algos)
timeout = 60*5
for problem_idx in 1:n_problems, algo_idx in 1:n_algos
    # To run just a subset of the problems, a filter can be put such as
    # the following, which were used when I added algorithms after the first run.
    ## algos[algo_idx].algo == "MISDP_dB" || continue
    ## algos[algo_idx].algo == "MISDP (dB^2)" || continue
    f = run_problem(algo_idx, problem_idx; verbose = true)
    f()
end


df = CSV.read(joinpath(@__DIR__, "results.csv"))
