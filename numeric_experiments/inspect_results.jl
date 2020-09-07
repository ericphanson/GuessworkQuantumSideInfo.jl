using CSV, DataFrames, Statistics
using Printf
df = CSV.read(joinpath(@__DIR__, "results.csv"))

## This bit is a little annoying:
# To check for OOMs, we need to manually look through the log (see `log.md`)
# Which was numbered according to the algorithm and problem index.
# However, later two more algorithms were added (`MISDP_dB` and `MISDP (dB^2)`)
# which can mess up the numbering. So first we will filter out those to mimick
# the state of the first run, and add flags for OOMs.
include("common.jl") # need to load the list of algorithms (`algos`)

algos_pre_dB = filter(algos) do a
    a.algo != "MISDP_dB" && a.algo != "MISDP (dB^2)"
end

# Looking through `log.md`, we find two issues:
# out of memory (OOM) errors during 4 solves, and unknown errors in 2 solves.
OOM = [(problem=10, algo=4), (problem=10, algo=5), (problem=10, algo=9),
       (problem=10, algo=10)]
unknown_error = [(problem=12, algo=4), (problem=12, algo=9)]

# Let us record the OOM errors in the dataframe.
inds = falses(size(df, 1))
for (p_idx, a_idx) in OOM
    inds .= inds .| ((df.problem .== problems[p_idx].problem) .& (df.algo .== algos_pre_dB[a_idx].algo) .& (df.settings .== algos_pre_dB[a_idx].settings))
end
df.OOM = inds
## End of problem-index specific content.

# Now we will filter out the `MISDP` algorithm (which had `dB^2 + 1` entries, which
# is 1 more than needed, and is replaced by the `MISDP (dB^2)` version.
df = filter(df) do row
    row.algo != "MISDP"
end

df.timeout = isnan.(df.optval) .& (!).(df.errors)

df.error_not_solved_not_timeout = df.errors .& isnan.(df.optval) .& (!).(df.timeout)
# Now we will caculate discrepencies in `optval`
# First we will filter to rows that claim to have a correct answer
row_should_be_correct(row) = !startswith(row.algo, "guesswork_upper_bound") && !isnan(row.optval) && row.algo != "MISDP_dB"

df_right = filter(row_should_be_correct, df)
gdf_right = groupby(df_right, "problem")

# Check for BB84(1)
bb84_true_val = (big(1) / big(4)) * (10 - sqrt(big(10)))

@assert all(gdf_right[4].problem .== Ref("BB84(1)")) # check we've got the right group
@assert maximum((gdf_right[4].optval .- bb84_true_val ) ./ bb84_true_val) < 1e-7 # not too big of an error

mean_solutions_df = combine(gdf_right, :optval => mean)
mean_solutions = Dict(mean_solutions_df.problem .=> mean_solutions_df.optval_mean)

# Now that we've calculated the mean solutions, we calculate the discrepencies

df.discrepency_from_mean = [ isnan(row.optval) ? NaN : !haskey(mean_solutions, row.problem) ? missing : abs(row.optval - mean_solutions[row.problem]) for row in eachrow(df)]

df.relative_discrepency_from_mean = [ isnan(row.optval) ? NaN : !haskey(mean_solutions, row.problem) ? missing : abs(row.optval - mean_solutions[row.problem])/mean_solutions[row.problem] for row in eachrow(df)]

# Let us check the maximum discrepency is not too bad
df_right = filter(row_should_be_correct, df)
findmax(df_right.relative_discrepency_from_mean) # (1.7202621117336732e-6, 45)

# Now let us generate the final table.
## Rows in table: algorithm with settings
## Columns in table: # solved, # timeout, # errored, avg relative_discrepency_from_mean, average time

df.algo_settings = df.algo .* Ref(" ") .* df.settings

function pretty_algo(row)
    if startswith(row.algo, "guesswork_upper_bound")
        return "Upper bound"
    end
    if row.algo == "MISDP_dB"
        return "MISDP (\$d_B\$)"
    elseif row.algo == "MISDP"
        return "MISDP (\$d_B^2 + 1\$)"
    elseif row.algo == "MISDP (dB^2)"
        return "MISDP (\$d_B^2\$)"
    elseif row.algo == "dual_SDP"
        return "SDP (dual)"
    end
    return row.algo
end

function pretty_settings(row)
    n = length("guesswork_upper_bound")

    if startswith(row.algo, "guesswork_upper_bound")
        j = findfirst(==('='), row.algo)
        max_time = row.algo[j+1:end-1]
        return row.settings * ", \$t_\\text{max}=$(max_time)\$" 
    end
    if row.settings == "Pajarito(Mosek, Gurobi, MSD=true)"
        return "Pajarito (c1)"
    elseif row.settings == "Pajarito(Mosek, Gurobi, MSD=false)"
        return "Pajarito (c2)"
    elseif row.settings == "Pajarito(SCS, Cbc, MSD=false)"
        return "Pajarito (o)"
    end
    return row.settings
end

df.algo2 = [ pretty_algo(row) for row in eachrow(df) ]
df.settings2 =  [ pretty_settings(row) for row in eachrow(df) ]

gdf = groupby(df, [:algo2, :settings2])

function mean_no_nan(collection)
    c = [x for x in collection if  !ismissing(x) && !isnan(x)]
    isempty(c) ? missing : mean(c)
end
float_fmt(f) = @sprintf("%2.2f", f)
pct_fmt(f) =  float_fmt(f*100) * "\\,\\%"
time_fmt(f) =  float_fmt(f) * "\\,s"

table1 = combine(gdf, :relative_discrepency_from_mean => pct_fmt ∘ mean_no_nan => "average relative error", :elapsed_seconds => time_fmt ∘ mean_no_nan => "average time", :optval => (x -> sum((!isnan).(x))) => "solved", :timeout => sum => "timeouts", :error_not_solved_not_timeout => sum => "errors")

sort!(table1, [:algo2, :settings2])
table1 = rename(table1, :algo2 => :Algorithm, :settings2 => :Parameters)

df.time_status = ifelse.(df.timeout, Ref(:timeout), ifelse.(df.error_not_solved_not_timeout, Ref(:error), df.elapsed_seconds))


df2 = copy(df)


function tol_time(row)
    disc = ismissing(row.relative_discrepency_from_mean) ? "(?\\,\\%)" : isnan(row.relative_discrepency_from_mean) ? "" : "("* pct_fmt(row.relative_discrepency_from_mean)*")"
    
    time = row.time_status isa Symbol ? string(row.time_status) : time_fmt(row.elapsed_seconds)
    disc == "" ? time : "$time $disc"
end

df2.value = [tol_time(row) for row in eachrow(df2) ]
df_times = unstack(df2, [:algo2, :settings2], :problem, :value)
problem_stems =  unique([ str[1:end-3] for str in df2.problem])

nms = names(df_times, Not(["algo2", "settings2"]))

for stem in problem_stems
    df_times[:, stem] =  ["$(row[1]), $(row[2])" for row in eachrow(df_times[:, nms[findall(startswith(stem), nms)]])]
end
table2 = df_times[:, ["algo2", "settings2", problem_stems...]]
sort!(table2, [:algo2, :settings2])
table2 = rename(table2, :algo2 => :Algorithm, :settings2 => :Parameters)

using Latexify
latexify(table1, env=:table, latex=false)
latexify(table2[:,["Algorithm", "Parameters", "2qubits", "2qutrits", "Y"] ], env=:table, latex=false)
latexify(table2[:,["Algorithm", "Parameters", "3qubits", "3qutrits", "BB84"] ], env=:table, latex=false)
