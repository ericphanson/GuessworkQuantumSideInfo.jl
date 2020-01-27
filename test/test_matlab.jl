using MATLAB

# Monkeypatch the MATLAB versions into the `GuessworkQuantumSideInfo` module
@eval GuessworkQuantumSideInfo begin

    @doc raw"""
        makeR(data)

    Compute a three-dimensional array `R`, such that the third index ranges from
    `y=1:J^k`, and `R[:,:,y]` is described by
    ``
    \sum_{x=1}^J p(x) \sum_{n\in L} c_n μ(n,x,\vec y)  ρ_B^x
    ``
    """
    function makeR(data; do_transpose::Bool)
        @unpack J, K, p, ρBs, c, dB, povm_outcomes = data

        d_out = length(povm_outcomes)

        R_function(g) = sum(c[N(g, x)] * p[x] * ρBs[x] for x in eachindex(p, ρBs))


        R = Array{ComplexF64}(undef, (dB, dB, d_out))
        if do_transpose
            for (y, outcome) in enumerate(povm_outcomes)
                R[:, :, y] .= transpose(R_function(outcome))
            end
        else
            for (y, outcome) in enumerate(povm_outcomes)
                R[:, :, y] .= R_function(outcome)
            end
        end
        R
    end


    function guesswork_MATLAB(
        p,
        ρBs;
        K = length(p),
        c = [1:K..., 10_000],
        dual::Bool = false,
        remove_repetition::Bool = true,
        povm_outcomes = make_povm_outcomes(length(p), K, remove_repetition),
    )
        length(p) == length(ρBs) ||
        throw(ArgumentError("Length of prior and vector of side information must match J"))

        J = length(p)
        dB = size(ρBs[1], 1)
        all(ρB -> size(ρB) == (dB, dB), ρBs) ||
        throw(ArgumentError("All side-information states must be square matrices of the same size"))

        d_out = length(povm_outcomes)

        data = (
            J = J,
            K = K,
            c = c,
            p = p,
            ρBs = ρBs,
            dB = dB,
            d_out = d_out,
            povm_outcomes = povm_outcomes,
        )


        if dual
            R = makeR(data; do_transpose = false)
            mat"""
            dB = double($(dB));
            d_out = double($(d_out));
            R = $(R);
            cvx_begin sdp quiet
            variable Y(dB, dB) hermitian
            dual variable E
            maximize(trace(Y))
            subject to
                E : R >= repmat(Y, [1, 1, d_out])
            cvx_end
            $(Y) = Y;
            $(E) = E;
            $(cvx_optval) = cvx_optval;
            $(cvx_status) = cvx_status;
            """
        else
            R = makeR(data; do_transpose = true)

            mat"""
            dB = double($(dB));
            d_out = double($(d_out));
            R = $(R);
            cvx_begin sdp quiet
            variable E(dB, dB, d_out) hermitian semidefinite
            dual variable Y
            minimize(real( sum( E(:) .* R(:) )))
            subject to
                Y : sum(E, 3) == eye(dB)
            cvx_end
            $(E) = E;
            $(Y) = Y;
            $(cvx_optval) = cvx_optval;
            $(cvx_status) = cvx_status;
            """
        end
        status = Symbol(cvx_status)

        if status != :Solved
            @error "CVX did not successfully solve the problem" status
        end
        Es = [E[:, :, j] for j = 1:d_out]
        return (
            optval = cvx_optval,
            status = status,
            Y = Y,
            Es = Es,
            povm_outcomes = Tuple.(povm_outcomes),
            data...,
        )
    end

end
