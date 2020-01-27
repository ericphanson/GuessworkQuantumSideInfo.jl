"""
    pmfN(data; tol = 1e-5) -> Vector

Compute the probability mass function for the number of guesses `N`, given a
strategy. The `n`th entry of the output vector `vec` gives the probability for
guessing the correct answer on the `n`th try, for `n = 1 : K`. The last entry
gives the probability of never guessing the correct answer.

* `data` should be a `NamedTuple` with entries for `p`, `ρBs`, `Es`, `K`, and
  `povm_outcomes`
* `tol` provides a tolerance above which to warn about imaginary
  or negative probabilities.

"""
function pmfN(data; tol = 1e-5)
    @unpack p, K = data
    T = eltype(p)

    probs = [_prob(n, data, tol, T) for n = 1:K]

    if K < length(p)
        prob_no_guess = 1 - sum(probs)

        if prob_no_guess < 0
            prob_no_guess < -tol && @warn "Probability negative " prob_no_guess
            prob_no_guess = zero(T)
        end

        push!(probs, prob_no_guess)
    end

    return probs
end

function _prob(n, data, tol, T)
    @unpack Es, p, ρBs, povm_outcomes = data
    prob = zero(T)
    for (y, outcome) in enumerate(povm_outcomes)
        for x in eachindex(p, ρBs)
            if findfirst(==(x), outcome) == n
                prob += p[x] * tr(Es[y] * ρBs[x])
            end
        end
    end
    if !isreal(prob)
        abs(imag(prob)) > tol && @warn "Imaginary probability" n prob
    end
    prob = real(prob)
    if prob < 0
        prob < -tol && @warn "Probability negative " n prob
        prob = zero(T)
    end
    return prob
end
