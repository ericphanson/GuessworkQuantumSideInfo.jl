using GuessworkQuantumSideInfo: ket, dm, BB84_states, iid_copies, randdm, guesswork,
                                guesswork_upper_bound, guesswork_MISDP
using StableRNGs, UnPack
using Mosek, Gurobi, Pajarito, Cbc, SCS

if !@isdefined(GRB_ENV)
    const GRB_ENV = Gurobi.Env()
end

function Y_states(::Type{T}) where {T}
    return [begin
                θⱼ = T(j) * (T(2) * π) / T(3)
                dm(cos(θⱼ) * ket(T, 1, 2) + sin(θⱼ) * ket(T, 2, 2))
            end
            for j in 1:3]
end

function two_random_qubits()
    rng = StableRNG(746)
    return [randdm(rng, Float64, 2) for _ in 1:2]
end

function three_random_qubits()
    rng = StableRNG(514)
    return [randdm(rng, Float64, 2) for _ in 1:3]
end

function two_random_qutrits()
    rng = StableRNG(948)
    return [randdm(rng, Float64, 3) for _ in 1:2]
end

function three_random_qutrits()
    rng = StableRNG(219)
    return [randdm(rng, Float64, 3) for _ in 1:3]
end

function iid_problem(states, n)
    J = length(states)
    p = ones(T, J^n) / J^n
    ρBs = iid_copies(states, n)
    (p=p, ρBs = ρBs)
end

T = Float64
problems = []
for n in 1:2
    append!(problems, [
                (iid_problem(two_random_qubits(), n)..., numeric_type=T, problem="2qubits($n)"),
                (iid_problem(three_random_qubits(), n)..., numeric_type=T, problem="3qubits($n)"),
                (iid_problem(Y_states(T), n)..., numeric_type=T, problem="Y($n)"),
                (iid_problem(BB84_states(T), n)..., numeric_type=T, problem="BB84($n)"),
                (iid_problem(two_random_qutrits(), n)..., numeric_type=T, problem="2qutrits($n)"),
                (iid_problem(three_random_qutrits(), n)..., numeric_type=T, problem="3qutrits($n)"),
                ])
end

algos = []

for (misdp_solver, settings) in [(PajaritoSolver(cont_solver=MosekSolver(LOG=0),
                     mip_solver=GurobiSolver(GRB_ENV, OutputFlag=0, IntFeasTol=1e-9, FeasibilityTol=1e-8, MIPGap=0),
                     mip_solver_drives=false, rel_gap = 1e-5,
                     log_level=0), "Pajarito(Mosek, Gurobi, MSD=false)"),
                     (PajaritoSolver(cont_solver=MosekSolver(LOG=0),
                     mip_solver=GurobiSolver(GRB_ENV, OutputFlag=0, IntFeasTol=1e-9, FeasibilityTol=1e-8, MIPGap=1e-5),
                     mip_solver_drives=true, rel_gap = 1e-5,
                     log_level=0), "Pajarito(Mosek, Gurobi, MSD=true)"),
     (PajaritoSolver(cont_solver=SCSSolver(verbose=false, eps=1e-6), mip_solver=CbcSolver(loglevel=0, integerT=1e-8),
                     mip_solver_drives=false, rel_gap = 1e-5,
                     log_level=0), "Pajarito(SCS, Cbc, MSD=false)")]
    push!(algos,
          (; f = (prob -> @timed guesswork_MISDP(prob.p, prob.ρBs; solver=misdp_solver, verbose=false).optval), algo = "MISDP",
                                                                                                                   settings))
end

for (sdp_solver, settings) in
    [((() -> MosekSolver(LOG=0)), "Mosek"), ((() -> SCSSolver(verbose=0, eps=1e-6)), "SCS")]
    push!(algos,
          (; f = (prob -> @timed guesswork(prob.p, prob.ρBs; solver=sdp_solver(), verbose=false).optval), algo = "SDP", settings))

    push!(algos,
          (; f = (prob -> @timed guesswork(prob.p, prob.ρBs; dual=true, solver=sdp_solver(), verbose=false).optval), algo = "dual_SDP", settings))

    for max_time in (20, 60, 60*4)
        push!(algos,
          (; f = (prob -> @timed guesswork_upper_bound(prob.p, prob.ρBs; make_solver=sdp_solver, verbose=false, max_time).optval), algo = "guesswork_upper_bound(max_time=$max_time)", settings))
    end
end
