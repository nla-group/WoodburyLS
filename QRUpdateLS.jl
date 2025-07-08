function QRUpdateLS(Q0, R0, U, V, b)
    m, n = size(Q0)
    r = size(U, 2)
    U0 = copy(U)

    # U0 .-= Q0 * (Q0' * U0)
    Q0tU = Q0'*U
    mul!(U0, Q0, Q0tU, -1.0, 1.0)
    
    Q1 = Matrix(qr!(U0).Q)
    
    # QTb = [Q0' * b; Q1' * b]
    QTb = Matrix{Float64}(undef, n + r, 1)
    mul!(view(QTb, 1:n, 1), Q0', b)
    mul!(view(QTb, n+1:n+r, 1), Q1', b)

    # W = [Q0' * U; Q1' * U]
    W = Matrix{Float64}(undef, n + r, r)
    copyto!(view(W, 1:n, 1:r), Q0tU)
    mul!(view(W, n+1:n+r, 1:r), Q1', U)

    # R = [R0; zeros(r, n)]
    R = Matrix{Float64}(undef, n + r, n)
    copyto!(view(R, 1:n, 1:n), R0)
    fill!(view(R, n+1:n+r, 1:n), 0.0)
    
    # Apply Givens rotations to R, W and QTb
    @inbounds for i ∈ 1:r
        @inbounds for k ∈ n-1+r:-1:i
            c, s, _ = LinearAlgebra.givensAlgorithm(W[k, i], W[k+1, i])
            BLAS.rot!(n, view(R, k, 1:n), n+r, view(R, k+1, 1:n), n+r, c, s)
            BLAS.rot!(r, view(W, k, 1:r), n+r, view(W, k+1, 1:r), n+r, c, s)
            BLAS.rot!(1, view(QTb, k, 1:1), 1, view(QTb, k+1, 1:1), 1, c, s)
        end
    end
    # R = R + WV'
    mul!(R, W, V', 1.0, 1.0)
    # Apply Householder reflections to R and QTb
    tmp = Vector{Float64}(undef, n)
    @inbounds for i ∈ 1:n
        μ = -sign(R[i,i])*norm(view(R, i:i+r, i))
        β = LAPACK.larfg!(view(R, i:i+r, i))
        LAPACK.larf!('L', view(R, i:i+r, i), β, view(R, i:i+r, i+1:n), view(tmp,1:n-i))
        LAPACK.larf!('L', view(R, i:i+r, i), β, view(QTb, i:i+r, 1:1), view(tmp,1:1))
        R[i,i] = μ
    end
    # x = R[1:n,1:n] \ QTb
    x = resize!(reshape(QTb, n+r), n)
    ldiv!(UpperTriangular(view(R, 1:n, 1:n)), x)
    return x
end
