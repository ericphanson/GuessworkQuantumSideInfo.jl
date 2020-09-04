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

function two_random_qubits(::Type{T}) where {T}
    rng = StableRNG(746)
    return [randdm(rng, Float64, 2) for _ in 1:2]
end

function problem1(::Type{T}) where {T}
    p = ones(T, 4) / 4
    ρBs = BB84_states(T)
    return (p=p, ρBs=ρBs)
end

function problem2(::Type{T}) where {T}
    p = ones(T, 3) / 3
    ρBs = Y_states(T)
    return (p=p, ρBs=ρBs)
end

function problem_qubits(::Type{T}, n) where {T}
    ρBs_1 = two_random_qubits(T)
    p = ones(T, 2^n) / 2^n
    return (p=p, ρBs=iid_copies(ρBs_1, n))
end

problem3(T) = problem_qubits(T, 1)
problem4(T) = problem_qubits(T, 2)

T = Float64
problems = [(problem3(T)..., numeric_type=T, problem="problem 3"),
            (problem2(T)..., numeric_type=T, problem="problem 2"),
            (problem4(T)..., numeric_type=T, problem="problem 4"),
            (problem1(T)..., numeric_type=T, problem="problem 1")]

algos = []

for (misdp_solver, settings) in [(PajaritoSolver(cont_solver=MosekSolver(LOG=0),
                     mip_solver=GurobiSolver(GRB_ENV, OutputFlag=0),
                     mip_solver_drives=false, use_mip_starts=true, solve_relax=false,
                     log_level=0), "Pajarito(Mosek, Gurobi)"),
     (PajaritoSolver(cont_solver=SCSSolver(verbose=false), mip_solver=CbcSolver(loglevel=0),
                     mip_solver_drives=false, use_mip_starts=true, solve_relax=false,
                     log_level=0), "Pajarito(SCS, Cbc)")]
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

    push!(algos,
          (; f = (prob -> @timed guesswork_upper_bound(prob.p, prob.ρBs; make_solver=sdp_solver, verbose=false, max_time=80).optval), algo = "guesswork_upper_bound", settings))
end
