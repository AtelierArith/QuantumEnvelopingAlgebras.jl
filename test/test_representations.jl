function mat_identity(q, n)
    M = fill(zero(q), n, n)
    for i in 1:n
        M[i,i] = one(q)
    end
    return M
end

function matmul(A, B, q)
    n = size(A, 1)
    result = fill(zero(q), n, n)
    for i in 1:n, j in 1:n, k in 1:n
        result[i,j] += A[i,k]*B[k,j]
    end
    return result
end

function matpow(A, exponent, q)
    result = mat_identity(q, size(A, 1))
    for _ in 1:exponent
        result = matmul(result, A, q)
    end
    return result
end

function matadd(A, B)
    return [A[i,j] + B[i,j] for i in axes(A,1), j in axes(A,2)]
end

function matscale(s, A)
    return [s*A[i,j] for i in axes(A,1), j in axes(A,2)]
end

function representation_matrices(q, n)
    makezero() = fill(zero(q), n+1, n+1)
    E, F, K, Ki = makezero(), makezero(), makezero(), makezero()
    for j in 0:n
        K[j+1,j+1] = q^(n-2j)
        Ki[j+1,j+1] = q^(-n+2j)
        if j < n
            F[j+2,j+1] = one(q)
        end
        if j > 0
            E[j,j+1] = qinteger(q,j)*qinteger(q,n-j+1)
        end
    end
    return E, F, K, Ki
end

function representation(x, matrices)
    E, F, K, Ki = matrices
    n = size(E, 1)
    q = one(AA.base_ring(x))
    result = fill(zero(q), n, n)
    for (m, coeff) in pbw_terms(x)
        Kpower = m.b >= 0 ? matpow(K, m.b, q) : matpow(Ki, -m.b, q)
        term = matmul(matmul(matpow(F, m.a, q), Kpower, q), matpow(E, m.c, q), q)
        result = matadd(result, matscale(coeff, term))
    end
    return result
end

function check_representations(U, q)
    E, F, K, Ki = AA.gens(U)
    x = E*F + K^(-2) + 2*F
    y = E + Ki*F + K
    for n in 0:6
        matrices = representation_matrices(q, n)
        Em, Fm, Km, Kim = matrices
        rho = highest_weight_representation(U, n)
        for (generator, independent) in zip((E,F,K,Ki), matrices)
            image = rho(generator)
            for i in 1:n+1, j in 1:n+1
                @test image[i,j] == independent[i,j]
            end
        end
        @test rho(x*y) == rho(x)*rho(y)
        external_E = rho.E
        external_E[1,1] = one(q)
        @test rho(E)[1,1] == zero(q)
        @test_throws ParentMismatchError rho(AA.gens(uqsl2(AA.base_ring(U),q)[1])[1])
        identity = mat_identity(q, n+1)
        @test matmul(Km, Kim, q) == identity
        @test matmul(matmul(Km, Em, q), Kim, q) == matscale(q^2, Em)
        @test matmul(matmul(Km, Fm, q), Kim, q) == matscale(q^(-2), Fm)
        @test matadd(matmul(Em,Fm,q), matscale(-one(q),matmul(Fm,Em,q))) ==
            matscale(inv(q-inv(q)), matadd(Km, matscale(-one(q), Kim)))
        @test representation(x*y, matrices) == matmul(representation(x, matrices), representation(y, matrices), q)
        @test representation(casimir(U)*x, matrices) == representation(x*casimir(U), matrices)
    end
end
