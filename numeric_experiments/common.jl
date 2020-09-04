using GuessworkQuantumSideInfo: ket, dm, BB84_states, iid_copies, randdm, guesswork, guesswork_upper_bound, guesswork_MISDP
using StableRNGs, UnPack
using Mosek, Gurobi, Pajarito, Cbc, SCS

if !@isdefined(GRB_ENV)
    const GRB_ENV = Gurobi.Env()
end

sdp_solver() = MosekSolver(LOG = 0)
# sdp_solver() = SCSSolver(verbose = 0, eps = 1e-6)
function misdp_solver(; verbose = false)
    PajaritoSolver(
        cont_solver = sdp_solver(),
        mip_solver = GurobiSolver(GRB_ENV, OutputFlag = 0),
        # mip_solver = CbcSolver(loglevel = 0),
        mip_solver_drives = false,
        use_mip_starts = true,
        solve_relax = false,
        log_level = verbose ? 3 : 0,
    )
end


function Y_states(::Type{T}) where {T}
    return [begin
                θⱼ = T(j) * (T(2) * π) / T(3)
                dm(cos(θⱼ) * ket(T, 1, 2) + sin(θⱼ) * ket(T, 2, 2))
            end
            for j in 1:3]
end

function two_random_qubits(::Type{T}) where T
    rng = StableRNG(746)
    return [randdm(rng, Float64, 2) for _ = 1:2]
end

function problem1(::Type{T}) where T
    p = ones(T, 4)/4
    ρBs = BB84_states(T)
    return (p=p, ρBs=ρBs)
end

function problem2(::Type{T}) where T
    p = ones(T, 3)/3
    ρBs = Y_states(T)
    return (p=p, ρBs=ρBs)
end

function problem_qubits(::Type{T}, n) where T
    ρBs_1 = two_random_qubits(T)
    p = ones(T, 2^n)/2^n
    return (p=p, ρBs=iid_copies(ρBs_1, n))
end

problem3(T) = problem_qubits(T, 1)
problem4(T) = problem_qubits(T, 2)

function algo1(prob)
    @unpack p, ρBs, numeric_type, problem = prob
    return @timed guesswork_MISDP(p, ρBs; solver = misdp_solver(verbose=false), verbose=false).optval
end

function algo2(prob)
    @unpack p, ρBs, numeric_type, problem = prob
    return @timed guesswork(p, ρBs; solver = sdp_solver(), verbose=false).optval
end

function algo3(prob)
    @unpack p, ρBs, numeric_type, problem = prob
    return @timed guesswork(p, ρBs; dual=true, solver = sdp_solver(), verbose=false).optval
end

function algo4(prob)
    @unpack p, ρBs, numeric_type, problem = prob
    return @timed guesswork_upper_bound(p, ρBs; make_solver = sdp_solver, verbose=false, max_time=80).optval
end


T = Float64
problems = [
    (problem3(T)..., numeric_type=T, problem = "problem 3"),
    (problem2(T)..., numeric_type=T, problem = "problem 2"),
    (problem4(T)..., numeric_type=T, problem = "problem 4"),
    (problem1(T)..., numeric_type=T, problem = "problem 1"),
]

algos = ((f=algo1, algo="MISDP", settings = "solver = Pajarito"),
        (f = algo2, algo="SDP", settings = "dual=false & solver = MosekSolver"),
        (f = algo3, algo="SDP", settings = "dual=true & solver = MosekSolver"),
        (f = algo4, algo="guesswork_upper_bound", settings = "solver = MosekSolver"))
