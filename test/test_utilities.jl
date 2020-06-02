Base.isapprox(x) = Base.Fix2(isapprox, x)

function testPOVM(Es; tol = TOL)
    dB = size(first(Es), 1)
    @test sum(Es) ≈ complex(1.0) * I(dB) atol = tol

    for E in Es
        @test E ≈ E' atol = tol
    end
    @test all(isposdef, Hermitian(E) + tol * I(dB) for E in Es)
end


function test_misdp_objective_value(data)
    @unpack povm_outcomes, ρBs, p, Es = data
    c = collect(1:length(povm_outcomes[1]))
    tot_cost = 0.0
    ρBs_tilde = p .* ρBs
    for y in eachindex(Es)
        E = Es[y]
        v = [tr(ρBs_tilde[x] * E) for x in eachindex(ρBs_tilde)]
        perm = povm_outcomes[y] |> collect
        cost = dot(c[invperm(perm)], v)
        tot_cost += cost
    end

    @test tot_cost ≈ data.optval rtol = TOL
end


get_SCS_solver() = TEST_MOI ? SCS.Optimizer(verbose = 0, eps = 1e-6) : SCSSolver(verbose = 0, eps = 1e-6)
function test_optimize(
    p,
    ρBs,
    true_opt_val = nothing;
    test_repetition = false,
    test_MISDP = true,
    kwargs...
)
    optvals = []
    test_data = []
    local current_output

    if TEST_ELLIPSOID
        current_output = guesswork_ellipsoid(p, ρBs; nl_solver = nl_solver(), tol=1e-3, verbose=false, kwargs...)
        push!(optvals, current_output.optval)
        push!(
            test_data,
            (
                opval = current_output.optval,
                solver = :ellipsoid,
                params = (tol=1e-3, kwargs...),
            ),
        )
    end
    
    for dual in (true, false)
        for remove_repetition in (test_repetition ? (false, true) : (true,))
            for (solver, solver_name) in ((get_SCS_solver(), :SCS),)
                current_output = guesswork(p, ρBs; solver = solver, dual = dual, remove_repetition = remove_repetition, kwargs...)
                push!(optvals, current_output.optval)
                push!(
                    test_data,
                    (
                        opval = current_output.optval,
                        solver = solver_name,
                        params = (dual = dual, remove_repetition = remove_repetition, kwargs...),
                    ),
                )

            end

            if TEST_MATLAB
                # `GuessworkQuantumSideInfo.guesswork_MATLAB` only exists if `include("test_matlab.jl")` is called.
                current_output =
                    GuessworkQuantumSideInfo.guesswork_MATLAB(p, ρBs; dual = dual, remove_repetition = remove_repetition, kwargs...)
                push!(optvals, current_output.optval)
                push!(
                    test_data,
                    (
                        opval = current_output.optval,
                        solver = :MATLAB,
                        params = (dual = dual, remove_repetition = remove_repetition, kwargs...),
                    ),
                )
            end
        end
    end

    if test_MISDP && TEST_MISDP
        dB = size(ρBs[1], 1)
        num_outcomes = min(factorial(length(p)), dB^2 + 1)
        if :c in keys(kwargs)
            current_output = guesswork_MISDP(
                p,
                ρBs,
                num_outcomes;
                c = kwargs.c,
                solver = misdp_solver(),
                verbose = false,
            )
        else
            current_output = guesswork_MISDP(
                p,
                ρBs,
                num_outcomes;
                solver = misdp_solver(),
                verbose = false,
            )
        end
        test_misdp_objective_value(current_output)
        # test_primal(current_output)
        push!(optvals, current_output.optval)
        push!(
            test_data,
            (
                opval = current_output.optval,
                solver = :MISDP,
                params = (num_outcomes = num_outcomes,),
            ),
        )

    end

    if true_opt_val === nothing
        true_opt_val = mean(optvals)
    end

    all_close_to_optimal = all(x -> isapprox(x, true_opt_val; rtol = TOL), optvals)

    # some helpful information in case of test failure
    if !all_close_to_optimal
        @error "`test_optimize` failure" true_opt_val
        for (i, test) in enumerate(test_data)
            @error "Solve $i results:" test
        end
    end
    @test all_close_to_optimal
    return current_output
end
